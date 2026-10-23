// AlphaPept 0.5.3 feature_finding pipeline port. See THIRD_PARTY_NOTICES.txt.
namespace DiaRawToMgf;

public static class AlphaPeptFeatureFinder
{
    public static List<Feature> Find(List<Ms1> scans,Options o,Action<string> log)
    {
        var q=ApInput.From(scans);
        log($"AlphaPept 0.5.3 port: {q.Mz.Length:N0} MS1 centroids");
        var first=AlphaPeptHills.Extract(q,o.MaxGap,o.Ppm,o.MinHillPoints);
        double tolerance=first.Median+3*first.Std;
        log($"First hills {first.Ptr.Length-1:N0}; adaptive tolerance {tolerance:R} ppm");
        if(!double.IsFinite(tolerance))return [];
        first=null!;
        var hills=AlphaPeptHills.Extract(q,o.MaxGap,tolerance,o.MinHillPoints);
        log($"Second hills {hills.Ptr.Length-1:N0}");
        hills=AlphaPeptHills.Deduplicate(hills);
        hills=hills with{Ptr=AlphaPeptHills.Split(hills,q.Intensity,o.SplitLevel,1)};
        hills=AlphaPeptHills.Filter(hills,q.Intensity);
        log($"Filtered/split hills {hills.Ptr.Length-1:N0}");
        if(hills.Ptr.Length<2)return [];
        var stats=AlphaPeptStats.Calculate(q,hills,o.Threads);
        log("Bootstrap statistics complete");
        var groups=AlphaPeptIsotopes.PrePatterns(q,hills,stats,o.MinCharge,o.MaxCharge,5,o.Correlation);
        log($"Pre isotope groups {groups.Length:N0}");
        var patterns=AlphaPeptIsotopes.Extract(q,hills,stats,groups,o.Threads,o.MinCharge,o.MaxCharge,100,5,o.Correlation,o.IsotopeSplitLevel,o.EnvelopeScore,log);
        log($"Isotope patterns {patterns.Patterns.Length:N0}");
        return AlphaPeptReport.Calculate(q,hills,stats,patterns,o.Threads).Features;
    }
}
