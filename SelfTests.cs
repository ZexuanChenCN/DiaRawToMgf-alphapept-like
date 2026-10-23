namespace DiaRawToMgf;

public static class SelfTests
{
    public static void Run()
    {
        static void Check(bool condition,string message){if(!condition)throw new Exception("SELF TEST: "+message);}
        Peak[] selected=Pipeline.SelectPeaks([new(200,2),new(100,10),new(150,10),new(90,1)],2);
        Check(selected.SequenceEqual([new Peak(100,10),new Peak(150,10)]),"Select intensity top N; write m/z ascending");
        var features=new List<Feature>{
            new(){Id=0,Mz=500,Charge=2,RtStart=0,RtEnd=2},
            new(){Id=1,Mz=501,Charge=2,RtStart=0,RtEnd=2},
            new(){Id=2,Mz=505,Charge=3,RtStart=0,RtEnd=2},
            new(){Id=3,Mz=500,Charge=2,RtStart=2,RtEnd=3}};
        Dictionary<int,double>[] local=[new(){[0]=10,[1]=2,[2]=5,[3]=8},new(){[0]=20,[1]=8}];
        var scan=new Ms2(77,1,500,4,0);
        var matched=Pipeline.Match(scan,features,local,[0,2],0,1);
        Check(matched.Count==2 && matched.Select(c=>c.FeatureId).Distinct().Count()==2,"Deduplicate left/right feature IDs; window/RT filter");
        Check(matched[0].FeatureId==0 && matched[0].Intensity==15,"Interpolate MS1 signal and rank candidates");
        using(var writer=new StringWriter())
        {
            Pipeline.WriteSpectrum(writer,"sample.77.77.2.0.dta",scan,features[0],"100.000000 10.000\r\n");
            Check(writer.ToString().Contains("SCANS=77") && writer.ToString().Contains(".77.77.2.0.dta"),"Preserve original RAW scan IDs");
        }
        var scans=new List<Ms1>();double[] shape=[1,4,10,4,1];
        for(int i=0;i<shape.Length;i++)scans.Add(new Ms1(100+i*9,i*0.01,[
            new(500,shape[i]*100),new(500+1.00335483507/2,shape[i]*50),
            new(500+1.00335483507,shape[i]*15)]));
        var q=ApInput.From(scans);
        var hills=new ApHills([0,5,10],[0,3,6,9,12,1,4,7,10,13],[],0,0);
        var stats=AlphaPeptStats.Calculate(q,hills,2);
        var report=AlphaPeptReport.Calculate(q,hills,stats,new ApPatterns([],[[0,1]],[2]),2);
        Check(report.Features.Count==1 && report.Features[0].Charge==2 && Math.Abs(report.Features[0].Mz-500)<1e-8,"Reported monoisotopic mz and charge");
        Check(report.Features[0].Signal.Count==5 && report.Features[0].Signal[2]==1500,"Isotope-sum chromatogram");
        Check(report.Lookup[6,0]==0 && report.Lookup[7,1]==1 && report.Lookup[8,0]==-1,"Centroid feature/isotope assignments");
        var rng=new ApRandom();Check(rng.Choice(5)==3&&rng.Choice(5)==4,"Pinned Numba random draw sequence");
        var random=new Random(42);
        foreach(int n in new[]{7,16,30,99,100,101,200})for(int trial=0;trial<30;trial++)
        {
            var values=new double[n,6];for(int i=0;i<n;i++){values[i,0]=random.Next(n/2+1);values[i,2]=trial%2==0?random.Next(n/2+1):random.NextDouble();}
            var active=Enumerable.Range(0,n).Where(_=>random.NextDouble()>0.2).ToHashSet();
            var remaining=Enumerable.Range(0,n).Where(active.Contains).ToArray();
            var pattern=ApNumerics.NumbaArgSort(remaining.Select(i=>values[i,0]).ToArray()).Select(i=>remaining[i]).ToArray();
            var areaSorted=Enumerable.Range(0,n).OrderByDescending(i=>values[i,2]).ToArray();
            var expected=ApNumerics.NumbaArgSort(pattern.Select(i=>values[i,2]).ToArray());Array.Reverse(expected);
            var actual=AlphaPeptIsotopes.SeedOrder(values,pattern,areaSorted,active,100);
            Check(actual.Take(100).SequenceEqual(expected.Take(100)),"Seed selection must preserve Numba ties and top100 boundary");
        }
        for(int trial=0;trial<2000;trial++)
        {
            int n=random.Next(2,250),seed=random.Next(n),z=random.Next(1,7);var values=new double[n,6];double mz=300;
            for(int i=0;i<n;i++){mz+=random.NextDouble()<0.1?3:random.Next(0,4)*0.1;values[i,0]=mz;values[i,1]=random.NextDouble()*0.2;}
            var pattern=Enumerable.Range(0,n).ToArray();var correlated=pattern.Where(i=>i!=seed&&random.NextDouble()<0.3).ToArray();var set=correlated.ToHashSet();
            var expected=new List<int>{seed};
            foreach(int direction in new[]{-1,1})
            {
                int index=1;
                for(int pos=seed+direction;pos>=0&&pos<n;pos+=direction)
                {
                    if(set.Contains(pos)&&AlphaPeptIsotopes.Spacing(values,seed,pos,z,-direction*index))
                    {if(direction>0)expected.Add(pos);else expected.Insert(0,pos);index++;}
                    if(Math.Abs(values[seed,0]-values[pos,0])>(AlphaPeptIsotopes.Delta+AlphaPeptIsotopes.DeltaStd)*index)break;
                }
            }
            Check(AlphaPeptIsotopes.GrowSparse(values,pattern,seed,z,correlated,5).SequenceEqual(expected),"Sparse growth must preserve the original stop condition, including wide SD and equal m/z");
        }
        Check(Options.Parse(["--help"]).Help,"Help parser");
        bool invalid=false;
        try{Options.Parse(["--input",".","--output","tmp","--top-peaks","0"]);}catch(ArgumentException){invalid=true;}
        Check(invalid,"Reject invalid CLI options");
        Console.WriteLine("PASS: top peaks; m/z order; original scan IDs; local signal; window/RT filters; feature-ID deduplication; synthetic feature; CLI validation");
    }
}
