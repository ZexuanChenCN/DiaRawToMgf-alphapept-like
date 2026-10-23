using System.Diagnostics;
using System.Globalization;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;

namespace DiaRawToMgf;

public static class Pipeline
{
    public const string Algorithm = "alphapept-0.5.3-port-v5";
    public static string Key(Options o)
    {
        var file=new FileInfo(o.Input);
        string text=JsonSerializer.Serialize(o with {Jobs=1,Worker=false,Resume=false})+
                    $"|{file.Length}|{file.LastWriteTimeUtc.Ticks}|{Algorithm}";
        if(o.FeatureCache.Length>0)
        {
            var paths=ReferenceCache.Paths(o.FeatureCache);
            foreach(var p in new[]{paths.Features,paths.Curves}){var f=new FileInfo(p);text+=$"|{p}|{f.Length}|{f.LastWriteTimeUtc.Ticks}";}
        }
        return Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(text)));
    }
    public static void Run(Options o)
    {
        Directory.CreateDirectory(o.Output);
        string done=Path.Combine(o.Output,"summary.json");
        string stem=Path.GetFileNameWithoutExtension(o.Input);
        string mgf=Path.Combine(o.Output,stem+".mgf");
        string key=Key(o);
        if(o.Resume && File.Exists(done))
        {
            using var s=JsonDocument.Parse(File.ReadAllText(done));
            if(s.RootElement.GetProperty("Fingerprint").GetString()==key && File.Exists(mgf) &&
               new FileInfo(mgf).Length==s.RootElement.GetProperty("MgfBytes").GetInt64() &&
               File.Exists(Path.Combine(o.Output,"mapping.tsv")))
            { Console.WriteLine($"[{stem}] Resume: existing completed output matches inputs/options");return; }
        }
        if(File.Exists(done) || File.Exists(mgf))
            throw new IOException($"Output already exists in {o.Output}; use a new output directory or --resume with identical options");
        using var logWriter=new StreamWriter(Path.Combine(o.Output,"run.log"),append:true,Encoding.UTF8);
        void Log(string text){lock(logWriter){string line=$"{DateTimeOffset.Now:O} [{stem}] {text}";logWriter.WriteLine(line);logWriter.Flush();Console.WriteLine(line);}}
        var total=Stopwatch.StartNew();
        File.WriteAllText(Path.Combine(o.Output,"options.json"),JsonSerializer.Serialize(o,JsonOptions));
        Log(o.FeatureCache.Length==0?$"Algorithm={Algorithm}; pure C# source port; no Python.":
            "Algorithm=cached-feature-association-v1; reuse supplied TSV feature/curve cache; no Python");
        var scans=new List<Ms1>();var ms2=new List<Ms2>();
        using(var raw=new RawReader(o.Input))
        {
            for(int scan=raw.First;scan<=raw.Last;scan++)
            {
                double rt=raw.Rt(scan);
                if(rt<o.RtStart)continue;if(rt>o.RtEnd)break;
                int level=raw.Level(scan);
                if(level==1)scans.Add(new Ms1(scan,rt,o.FeatureCache.Length==0?raw.Peaks(scan):[]));
                else if(level==2)ms2.Add(raw.Metadata(scan,rt));
                if(scan%20000==0)Log($"Read metadata/MS1: RAW scan {scan:N0}/{raw.Last:N0}");
            }
        }
        if(scans.Count<o.MinHillPoints || ms2.Count==0)
            throw new InvalidDataException("Insufficient MS1 or MS2 scans in the selected RT interval");
        double readSeconds=total.Elapsed.TotalSeconds;
        Log($"Read {scans.Count:N0} MS1, {ms2.Count:N0} MS2; {readSeconds:F1}s");
        var features=o.FeatureCache.Length==0?AlphaPeptFeatureFinder.Find(scans,o,Log):ReferenceCache.Load(o.FeatureCache,scans);
        double featureSeconds=total.Elapsed.TotalSeconds-readSeconds;
        Log($"Features: {features.Count:N0}; feature finding {featureSeconds:F1}s");
        WriteFeatures(o.Output,features,scans);
        var local=new Dictionary<int,double>[scans.Count];
        for(int i=0;i<local.Length;i++)local[i]=[];
        foreach(var f in features)foreach(var (scan,intensity) in f.Signal)local[scan][f.Id]=intensity;
        int pairs=0,matched=0,unmatched=0,capped=0;
        long availablePairs=0;
        var rt1=scans.Select(s=>s.Rt).ToArray();
        string partial=mgf+".partial";
        using(var raw=new RawReader(o.Input))
        using(var writer=TextWriter(partial))
        using(var map=TextWriter(Path.Combine(o.Output,"mapping.tsv")))
        using(var missing=TextWriter(Path.Combine(o.Output,"unmatched_ms2.tsv")))
        {
            map.WriteLine("record_id\traw_scan\tcandidate_rank\ttitle\tfeature_id\tprecursor_mz\tcharge\trt_min\twindow_lower\twindow_upper\tleft_ms1_scan\tright_ms1_scan\tleft_intensity\tright_intensity\tlocal_intensity\tall_candidate_count\tpeaks_before\tpeaks_after");
            missing.WriteLine("raw_scan\trt_min\twindow_lower\twindow_upper\treason");
            foreach(var m in ms2)
            {
                int right=ApNumerics.LowerBound(rt1,m.Rt),left=right-1;
                var candidates=Match(m,features,local,rt1,left,right);
                int allCandidateCount=candidates.Count;
                availablePairs+=candidates.Count;
                if(o.MaxCandidates>0 && candidates.Count>o.MaxCandidates){capped++;candidates.RemoveRange(o.MaxCandidates,candidates.Count-o.MaxCandidates);}
                if(candidates.Count==0)
                {unmatched++;missing.WriteLine(F($"{m.Scan}\t{m.Rt:R}\t{m.Lower:R}\t{m.Upper:R}\tno_locally_observed_feature"));continue;}
                Peak[] original=raw.Peaks(m.Scan),peaks=SelectPeaks(original,o.TopPeaks);
                if(peaks.Length==0)
                {unmatched++;missing.WriteLine(F($"{m.Scan}\t{m.Rt:R}\t{m.Lower:R}\t{m.Upper:R}\tno_valid_ms2_peaks"));continue;}
                matched++;
                string fragmentText=string.Concat(peaks.Select(p=>
                    (Math.Round(p.Mz,6)>0?p.Mz.ToString("F6",CultureInfo.InvariantCulture):p.Mz.ToString("R",CultureInfo.InvariantCulture))+" "+
                    (Math.Round(p.Intensity,3)>0?p.Intensity.ToString("F3",CultureInfo.InvariantCulture):p.Intensity.ToString("R",CultureInfo.InvariantCulture))+"\r\n"));
                for(int rank=0;rank<candidates.Count;rank++)
                {
                    var c=candidates[rank];var f=features[c.FeatureId];pairs++;
                    string title=$"{stem}.{m.Scan}.{m.Scan}.{f.Charge}.{rank}.dta";
                    WriteSpectrum(writer,title,m,f,fragmentText);
                    map.WriteLine(F($"{pairs}\t{m.Scan}\t{rank}\t{title}\t{f.Id}\t{f.Mz:R}\t{f.Charge}\t{m.Rt:R}\t{m.Lower:R}\t{m.Upper:R}\t{(left>=0?scans[left].Scan:-1)}\t{(right<scans.Count?scans[right].Scan:-1)}\t{c.LeftIntensity:R}\t{c.RightIntensity:R}\t{c.Intensity:R}\t{allCandidateCount}\t{original.Length}\t{peaks.Length}"));
                }
                if(matched%10000==0)Log($"Matched MS2 {matched:N0}; candidate spectra {pairs:N0}");
            }
        }
        Log("Validating MGF scan numbers, titles, counts, precursor fields and m/z order");
        ValidateMgf(partial,ms2,o.TopPeaks,pairs);
        File.Move(partial,mgf);
        var report=new {Algorithm=o.FeatureCache.Length==0?Algorithm:"cached-feature-association-v1",Fingerprint=key,Input=o.Input,Ms1Scans=scans.Count,Ms2Scans=ms2.Count,
            Features=features.Count,MatchedMs2=matched,UnmatchedMs2=unmatched,Spectra=pairs,
            CandidatePairsBeforeCap=availablePairs,Ms2WithCandidatesCapped=capped,
            ReadSeconds=readSeconds,FeatureSeconds=featureSeconds,
            ExportAndValidationSeconds=total.Elapsed.TotalSeconds-readSeconds-featureSeconds,
            TotalSeconds=total.Elapsed.TotalSeconds,MgfBytes=new FileInfo(mgf).Length,
            OriginalScanNumbers=true,PeakOrder="m/z ascending",FragmentDeconvolution=false,
            Validation="PASS",AlphaPeptPortTarget="0.5.3",FeatureCache=o.FeatureCache,
            PeakWorkingSetBytes=Process.GetCurrentProcess().PeakWorkingSet64};
        string temp=done+".partial";
        File.WriteAllText(temp,JsonSerializer.Serialize(report,JsonOptions));File.Move(temp,done,true);
        Log($"DONE {pairs:N0} spectra from {matched:N0} original MS2; {total.Elapsed.TotalMinutes:F2} min");
    }
    public static List<Candidate> Match(Ms2 m,List<Feature> features,Dictionary<int,double>[] local,double[] rt,int left,int right)
    {
        var ids=new HashSet<int>();
        if(left>=0)ids.UnionWith(local[left].Keys);
        if(right<local.Length)ids.UnionWith(local[right].Keys);
        var result=new List<Candidate>();
        foreach(int id in ids)
        {
            var f=features[id];
            if(m.Width<=0 || f.Mz<m.Lower || f.Mz>=m.Upper || m.Rt<f.RtStart || m.Rt>f.RtEnd)continue;
            double il=left>=0?local[left].GetValueOrDefault(id):0;
            double ir=right<local.Length?local[right].GetValueOrDefault(id):0;
            double strength;
            if(left>=0 && right<rt.Length)
            {double alpha=(m.Rt-rt[left])/(rt[right]-rt[left]);strength=il*(1-alpha)+ir*alpha;}
            else strength=il+ir;
            if(strength>0)result.Add(new Candidate(id,il,ir,strength));
        }
        return result.OrderByDescending(c=>c.Intensity).ThenBy(c=>c.FeatureId).ToList();
    }
    public static Peak[] SelectPeaks(Peak[] peaks,int top)=>peaks
        .Where(p=>double.IsFinite(p.Mz)&&p.Mz>0&&double.IsFinite(p.Intensity)&&p.Intensity>0)
        .OrderByDescending(p=>p.Intensity).ThenBy(p=>p.Mz).Take(top).OrderBy(p=>p.Mz).ToArray();
    public static void WriteSpectrum(TextWriter writer,string title,Ms2 m,Feature f,string fragments)
    {
        writer.WriteLine("BEGIN IONS");writer.WriteLine($"TITLE={title}");writer.WriteLine($"CHARGE={f.Charge}+");
        writer.WriteLine(F($"RTINSECONDS={m.Rt*60:F6}"));writer.WriteLine(F($"PEPMASS={f.Mz:F8}"));
        writer.WriteLine($"SCANS={m.Scan}");writer.Write(fragments);writer.WriteLine("END IONS");writer.WriteLine();
    }
    private static void WriteFeatures(string output,List<Feature> features,List<Ms1> scans)
    {
        using var f=TextWriter(Path.Combine(output,"ms1_features.tsv"));
        using var c=TextWriter(Path.Combine(output,"precursor_curves.tsv"));
        f.WriteLine("feature_id\tmz\tmz_std\tmz_most_abundant\tcharge\trt_start\trt_apex\trt_end\tfwhm\tn_isotopes\tmass\tms1_int_sum_apex\tms1_int_sum_area\tms1_int_max_apex\tms1_int_max_area");
        c.WriteLine("feature_id\tms1_scan\trt_min\tisotope_sum_intensity");
        foreach(var x in features)
        {
            double[] values=x.ReportValues.Length==14?x.ReportValues:[x.Mz,0,x.Mz,x.Charge,x.RtStart,x.RtApex,x.RtEnd,0,x.HillIds.Length,x.Mass,x.ApexIntensity,x.Area,0,0];
            f.WriteLine(x.Id+"\t"+string.Join('\t',values.Select(v=>v.ToString("R",CultureInfo.InvariantCulture))));
            foreach(var (si,it) in x.Signal)c.WriteLine(F($"{x.Id}\t{scans[si].Scan}\t{scans[si].Rt:R}\t{it:R}"));
        }
    }
    public static void ValidateMgf(string path,List<Ms2> ms2,int top,int expected)
    {
        var scans=ms2.ToDictionary(s=>s.Scan);var titles=new HashSet<string>();
        bool inside=false;int count=0,peaks=0,scan=-1,charge=0;double lastMz=0,precursor=0,rtSeconds=-1;string title="";
        foreach(string line in File.ReadLines(path))
        {
            if(line=="BEGIN IONS"){if(inside)throw new InvalidDataException("Nested spectrum");inside=true;peaks=0;scan=-1;lastMz=0;title="";charge=0;precursor=0;rtSeconds=-1;}
            else if(line.StartsWith("TITLE=")){title=line[6..];if(!titles.Add(title))throw new InvalidDataException("Duplicate TITLE");}
            else if(line.StartsWith("SCANS="))scan=int.Parse(line[6..],CultureInfo.InvariantCulture);
            else if(line.StartsWith("CHARGE="))charge=int.Parse(line[7..].TrimEnd('+'),CultureInfo.InvariantCulture);
            else if(line.StartsWith("PEPMASS="))precursor=double.Parse(line[8..],CultureInfo.InvariantCulture);
            else if(line.StartsWith("RTINSECONDS="))rtSeconds=double.Parse(line[12..],CultureInfo.InvariantCulture);
            else if(line=="END IONS")
            {
                if(!inside || peaks<1 || peaks>top || !scans.ContainsKey(scan))throw new InvalidDataException("Invalid scan/peak count");
                string suffix=title[..^4];var parts=suffix.Split('.');
                if(parts.Length<5 || int.Parse(parts[^4])!=scan || int.Parse(parts[^3])!=scan)
                    throw new InvalidDataException("TITLE does not preserve original scan");
                var metadata=scans[scan];
                if(charge<1 || int.Parse(parts[^2])!=charge || !double.IsFinite(precursor) ||
                   precursor<metadata.Lower-1e-7 || precursor>=metadata.Upper+1e-7 ||
                   !double.IsFinite(rtSeconds) || Math.Abs(rtSeconds-metadata.Rt*60)>1e-5)
                    throw new InvalidDataException("Invalid precursor, charge or retention time");
                count++;inside=false;
            }
            else if(inside && line.Length>0 && !line.Contains('='))
            {
                var p=line.Split(' ',StringSplitOptions.RemoveEmptyEntries);
                if(p.Length!=2)throw new InvalidDataException("Malformed peak");
                double mz=double.Parse(p[0],CultureInfo.InvariantCulture),it=double.Parse(p[1],CultureInfo.InvariantCulture);
                if(!double.IsFinite(mz)||!double.IsFinite(it)||mz<lastMz||mz<=0||it<=0)throw new InvalidDataException("Invalid or unsorted peak");
                lastMz=mz;peaks++;
            }
        }
        if(inside || count!=expected)throw new InvalidDataException("Truncated MGF / unexpected count");
    }
    public static StreamWriter TextWriter(string path)=>new(path,false,new UTF8Encoding(false),1024*1024){NewLine="\r\n"};
    public static string F(FormattableString value)=>value.ToString(CultureInfo.InvariantCulture);
    public static readonly JsonSerializerOptions JsonOptions=new(){WriteIndented=true};
}
