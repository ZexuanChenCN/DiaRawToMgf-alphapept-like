using System.Runtime.InteropServices;
using System.Text.Json;

namespace DiaRawToMgf;

public static class AlphaPeptStageVerification
{
    public static ApInput Load(string directory)
    {
        var bounds=Read<long>(directory,"indices_ms1").Select(x=>checked((int)x)).ToArray();
        return new(bounds,Read<double>(directory,"mass_list_ms1"),Read<double>(directory,"int_list_ms1"),
            Read<double>(directory,"rt_list_ms1"),Read<int>(directory,"scan_numbers"));
    }
    public static T[] Read<T>(string directory,string name)where T:unmanaged
        =>MemoryMarshal.Cast<byte,T>(File.ReadAllBytes(Path.Combine(directory,name+".bin"))).ToArray();
    public static void Verify(string directory)
    {
        var q=Load(directory);var report=new List<object>();
        void Equal(string stage,string name,int[] actual)
        {
            var expected=Read<int>(directory,name);int limit=Math.Min(actual.Length,expected.Length);
            int mismatches=0,first=-1;
            for(int i=0;i<limit;i++)if(actual[i]!=expected[i]){mismatches++;if(first<0)first=i;}
            bool pass=actual.Length==expected.Length&&mismatches==0;
            report.Add(new{Stage=stage,Array=name,ExpectedLength=expected.Length,ActualLength=actual.Length,Mismatches=mismatches,FirstMismatch=first,Pass=pass});
            Console.WriteLine($"{name}: {(pass?"PASS":"FAIL")} expected={expected.Length}, actual={actual.Length}, mismatches={mismatches}, first={first}");
            if(!pass)throw new InvalidDataException($"Stage {stage}: {name} differs");
        }
        try
        {
            var first=AlphaPeptHills.Extract(q,2,8);
            Equal("first","first_ptr",first.Ptr);Equal("first","first_data",first.Data);Equal("first","first_counts",first.Counts);
            using var manifest=JsonDocument.Parse(File.ReadAllText(Path.Combine(directory,"manifest.json")));
            double expectedTol=manifest.RootElement.GetProperty("adaptive_tolerance").GetDouble();
            double tolerance=first.Median+first.Std*3;
            Console.WriteLine($"Adaptive tolerance: {tolerance:R}; reference: {expectedTol:R}");
            if(Math.Abs(tolerance-expectedTol)>1e-9)throw new InvalidDataException("Adaptive tolerance differs");
            var second=AlphaPeptHills.Extract(q,2,tolerance);
            Equal("second","second_ptr",second.Ptr);Equal("second","second_data",second.Data);Equal("second","second_counts",second.Counts);
            var dedup=AlphaPeptHills.Deduplicate(second);
            Equal("dedup","dedup_ptr",dedup.Ptr);Equal("dedup","dedup_data",dedup.Data);
            var split=dedup with{Ptr=AlphaPeptHills.Split(dedup,q.Intensity,1.3,1)};
            Equal("split","split_ptr",split.Ptr);
            var filtered=AlphaPeptHills.Filter(split,q.Intensity);
            Equal("filter","filtered_ptr",filtered.Ptr);Equal("filter","filtered_data",filtered.Data);
            var stats=AlphaPeptStats.Calculate(q,filtered,4);
            var expectedOrder=Read<long>(directory,"sortindex").Select(x=>(int)x).ToArray();
            if(!stats.Order.SequenceEqual(expectedOrder))throw new InvalidDataException("Stats sort order differs");
            var expectedUpper=Read<long>(directory,"upper").Select(x=>(int)x).ToArray();
            if(!stats.Upper.SequenceEqual(expectedUpper))throw new InvalidDataException("Stats upper bounds differ");
            var expectedStats=Read<double>(directory,"stats");double maxError=0;
            for(int i=0;i<stats.Order.Length;i++)for(int j=0;j<6;j++)maxError=Math.Max(maxError,Math.Abs(stats.Values[i,j]-expectedStats[i*6+j]));
            Console.WriteLine($"Bootstrap stats max absolute error: {maxError:R}");
            report.Add(new{Stage="stats",MaxAbsoluteError=maxError,SortOrderPass=true,UpperBoundsPass=true});
            if(maxError>1e-7)throw new InvalidDataException("Bootstrap stats differ");
            var groups=AlphaPeptIsotopes.PrePatterns(q,filtered,stats);
            var expectedPrePtr=Read<long>(directory,"pre_ptr");var expectedPreData=Read<long>(directory,"pre_data");
            if(groups.Length+1!=expectedPrePtr.Length)throw new InvalidDataException($"Pre pattern count differs: {groups.Length} vs {expectedPrePtr.Length-1}");
            for(int i=0;i<groups.Length;i++)if(!groups[i].Select(x=>(long)x).SequenceEqual(expectedPreData[(int)expectedPrePtr[i]..(int)expectedPrePtr[i+1]]))
                throw new InvalidDataException($"Pre pattern {i} differs");
            Console.WriteLine($"Pre isotope patterns: PASS ({groups.Length})");
            report.Add(new{Stage="pre_isotopes",Groups=groups.Length,Pass=true});
            var patterns=AlphaPeptIsotopes.Extract(q,filtered,stats,groups,4);
            var expectedPatterns=Read<long>(directory,"patterns");var expectedIsoIdx=Read<long>(directory,"isoidx");var expectedCharges=Read<long>(directory,"charges");
            int same=0;
            for(int i=0;i<Math.Min(patterns.Patterns.Length,expectedCharges.Length);i++)
                if(patterns.Charges[i]==expectedCharges[i]&&patterns.Patterns[i].Select(x=>(long)x)
                    .SequenceEqual(expectedPatterns[(int)expectedIsoIdx[i]..(int)expectedIsoIdx[i+1]]))same++;
            Console.WriteLine($"Isotope pattern agreement: {same}/{expectedCharges.Length}; C# count={patterns.Patterns.Length}");
            report.Add(new{Stage="isotopes",Same=same,Reference=expectedCharges.Length,Actual=patterns.Patterns.Length});
            if(same!=expectedCharges.Length||patterns.Patterns.Length!=expectedCharges.Length)throw new InvalidDataException("Isotope patterns differ");
            var final=AlphaPeptReport.Calculate(q,filtered,stats,patterns,4);
            var lines=File.ReadAllLines(Path.Combine(directory,"feature_table.tsv"));double reportError=0;
            for(int i=0;i<final.Features.Count;i++)
            {
                var values=lines[i+1].Split('\t').Select(x=>double.Parse(x,System.Globalization.CultureInfo.InvariantCulture)).ToArray();
                for(int j=0;j<14;j++)reportError=Math.Max(reportError,Math.Abs(values[j]-final.Rows[i,j]));
            }
            var expectedLookup=Read<int>(directory,"lookup");int lookupDifferences=0;
            for(int i=0;i<q.Mz.Length;i++)for(int j=0;j<2;j++)if(final.Lookup[i,j]!=expectedLookup[i*2+j])lookupDifferences++;
            Console.WriteLine($"Feature report max absolute error: {reportError:R}; lookup differences: {lookupDifferences}");
            report.Add(new{Stage="report",MaxAbsoluteError=reportError,LookupDifferences=lookupDifferences});
            if(reportError>1e-6||lookupDifferences!=0)throw new InvalidDataException("Feature report differs");
        }
        finally
        {
            File.WriteAllText(Path.Combine(directory,"csharp-stages.json"),JsonSerializer.Serialize(report,Pipeline.JsonOptions));
        }
    }
}
