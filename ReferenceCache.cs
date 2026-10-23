using System.Globalization;

namespace DiaRawToMgf;

/// <summary>Optional language-neutral reference cache for isolated association comparisons.</summary>
public static class ReferenceCache
{
    public static (string Features,string Curves) Paths(string directory)
    {
        string features=Path.Combine(directory,"ms1_features.tsv"),curves=Path.Combine(directory,"precursor_curves.tsv");
        if(!File.Exists(features)||!File.Exists(curves))throw new FileNotFoundException("Feature cache needs ms1_features.tsv and precursor_curves.tsv (see README.txt)");
        return(features,curves);
    }
    public static List<Feature> Load(string directory,List<Ms1> scans)
    {
        var paths=Paths(directory);
        using var reader=new StreamReader(paths.Features);
        var header=reader.ReadLine()!.Split('\t');
        var cols=header.Select((name,index)=>(name,index)).ToDictionary(x=>x.name,x=>x.index);
        var result=new List<Feature>();
        string? line;
        while((line=reader.ReadLine())is not null)
        {
            var values=line.Split('\t');
            double Get(string key)=>double.Parse(values[cols[key]],CultureInfo.InvariantCulture);
            double Either(string a,string b)=>Get(cols.ContainsKey(a)?a:b);
            result.Add(new Feature{Id=(int)Get("feature_id"),Mz=Get("mz"),Charge=(int)Get("charge"),
                RtStart=Get("rt_start"),RtEnd=Get("rt_end"),RtApex=Get("rt_apex"),
                ApexIntensity=Either("ms1_int_sum_apex","intensity_apex"),Area=Either("ms1_int_sum_area","area_min"),
                HillIds=new int[(int)Get("n_isotopes")]});
        }
        result=result.OrderBy(f=>f.Id).ToList();
        for(int i=0;i<result.Count;i++)if(result[i].Id!=i)throw new InvalidDataException("Cache feature IDs must be contiguous");
        var lookup=scans.Select((s,index)=>(s.Scan,index)).ToDictionary(s=>s.Scan,s=>s.index);
        using var curves=new StreamReader(paths.Curves);
        var curveHeader=curves.ReadLine()!.Split('\t');
        var ci=curveHeader.Select((name,index)=>(name,index)).ToDictionary(x=>x.name,x=>x.index);
        while((line=curves.ReadLine())is not null)
        {
            var v=line.Split('\t');
            int scan=int.Parse(v[ci["ms1_scan"]],CultureInfo.InvariantCulture);
            if(!lookup.TryGetValue(scan,out int index))continue;
            int id=int.Parse(v[ci["feature_id"]],CultureInfo.InvariantCulture);
            double intensity=double.Parse(v[ci["isotope_sum_intensity"]],CultureInfo.InvariantCulture);
            if((uint)id>=(uint)result.Count||!double.IsFinite(intensity)||intensity<=0)throw new InvalidDataException("Invalid cached curve point");
            result[id].Signal.Add(index,intensity);
        }
        return result;
    }
}
