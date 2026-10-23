using System.Globalization;

namespace DiaRawToMgf;

public sealed record Options
{
    public string Input { get; set; } = "";
    public string Output { get; set; } = "";
    public string FeatureCache { get; set; } = "";
    public bool Recursive { get; set; }
    public bool Worker { get; set; }
    public bool Resume { get; set; }
    public bool SelfTest { get; set; }
    public bool Help { get; set; }
    public int Jobs { get; set; } = 1;
    public int Threads { get; set; } = 4;
    public int TopPeaks { get; set; } = 250;
    public int MaxCandidates { get; set; } = 5;
    public int MinCharge { get; set; } = 1;
    public int MaxCharge { get; set; } = 6;
    public int MaxGap { get; set; } = 2;
    public int MinHillPoints { get; set; } = 3;
    public int MinIsotopes { get; set; } = 2;
    public int MaxIsotopes { get; set; } = 10;
    public double Ppm { get; set; } = 8;
    public double IsotopePpm { get; set; } = 8;
    public double Correlation { get; set; } = 0.6;
    public double EnvelopeScore { get; set; } = 0.6;
    public double SplitLevel { get; set; } = 1.3;
    public double IsotopeSplitLevel { get; set; } = 1.3;
    public double MinMs1Intensity { get; set; } = 0;
    public double RtStart { get; set; } = 0;
    public double RtEnd { get; set; } = double.MaxValue;

    public static Options Parse(string[] args)
    {
        var o = new Options();
        for (int i = 0; i < args.Length; i++)
        {
            string Value() => ++i < args.Length ? args[i] : throw new ArgumentException($"Missing value for {args[i-1]}");
            int Int() => int.Parse(Value(), CultureInfo.InvariantCulture);
            double Double() => double.Parse(Value(), CultureInfo.InvariantCulture);
            switch (args[i])
            {
                case "--input": o.Input = Value(); break;
                case "--output": o.Output = Value(); break;
                case "--feature-cache": o.FeatureCache = Value(); break;
                case "--recursive": o.Recursive = true; break;
                case "--worker": o.Worker = true; break;
                case "--resume": o.Resume = true; break;
                case "--self-test": o.SelfTest = true; break;
                case "--help": case "-h": o.Help = true; break;
                case "--jobs": o.Jobs = Int(); break;
                case "--threads": o.Threads = Int(); break;
                case "--top-peaks": o.TopPeaks = Int(); break;
                case "--max-candidates": o.MaxCandidates = Int(); break;
                case "--min-charge": o.MinCharge = Int(); break;
                case "--max-charge": o.MaxCharge = Int(); break;
                case "--max-gap": o.MaxGap = Int(); break;
                case "--min-hill-points": o.MinHillPoints = Int(); break;
                case "--min-isotopes": o.MinIsotopes = Int(); break;
                case "--max-isotopes": o.MaxIsotopes = Int(); break;
                case "--centroid-ppm": o.Ppm = Double(); break;
                case "--isotope-ppm": o.IsotopePpm = Double(); break;
                case "--isotope-correlation": o.Correlation = Double(); break;
                case "--envelope-score": o.EnvelopeScore = Double(); break;
                case "--hill-split-level": o.SplitLevel = Double(); break;
                case "--isotope-split-level": o.IsotopeSplitLevel = Double(); break;
                case "--min-ms1-intensity": o.MinMs1Intensity = Double(); break;
                case "--rt-start": o.RtStart = Double(); break;
                case "--rt-end": o.RtEnd = Double(); break;
                default: throw new ArgumentException($"Unknown option: {args[i]}");
            }
        }
        if (o.Help || o.SelfTest) return o;
        if(o.IsotopePpm!=8||o.MinIsotopes!=2||o.MaxIsotopes!=10||o.MinMs1Intensity!=0)
            throw new ArgumentException("--isotope-ppm, --min-isotopes, --max-isotopes and --min-ms1-intensity belong to the retired approximate engine; the AlphaPept port does not use these settings");
        if(!double.IsFinite(o.IsotopeSplitLevel)||o.IsotopeSplitLevel<=1)throw new ArgumentException("Invalid isotope split level");
        if (string.IsNullOrWhiteSpace(o.Input) || string.IsNullOrWhiteSpace(o.Output))
            throw new ArgumentException("--input and --output are required");
        if (o.Jobs < 1 || o.Threads < 1 || o.TopPeaks < 1 || o.MaxCandidates < 0 ||
            o.MinCharge < 1 || o.MaxCharge < o.MinCharge || o.MaxCharge > 20 || o.MaxGap < 0 ||
            o.MinHillPoints < 2 || o.MinIsotopes < 2 || o.MaxIsotopes < o.MinIsotopes || o.MaxIsotopes > 30 ||
            !double.IsFinite(o.Ppm) || o.Ppm <= 0 || !double.IsFinite(o.IsotopePpm) || o.IsotopePpm <= 0 ||
            !double.IsFinite(o.Correlation) || o.Correlation < 0 || o.Correlation > 1 ||
            !double.IsFinite(o.EnvelopeScore) || o.EnvelopeScore < 0 || o.EnvelopeScore > 1 ||
            !double.IsFinite(o.SplitLevel) || o.SplitLevel <= 1 || !double.IsFinite(o.MinMs1Intensity) ||
            o.MinMs1Intensity < 0 || !double.IsFinite(o.RtStart) || !double.IsFinite(o.RtEnd) ||
            o.RtStart < 0 || o.RtEnd <= o.RtStart)
            throw new ArgumentException("Invalid numeric option; see --help");
        o.Input = Path.GetFullPath(o.Input);
        o.Output = Path.GetFullPath(o.Output);
        if(o.FeatureCache.Length>0)o.FeatureCache=Path.GetFullPath(o.FeatureCache);
        return o;
    }
}

public readonly record struct Peak(double Mz, double Intensity);
public sealed record Ms1(int Scan, double Rt, Peak[] Peaks);
public sealed record Ms2(int Scan, double Rt, double Center, double Width, double Offset)
{
    // Offset is the displacement of the selected mass from the window center.
    public double Lower => Center - Width / 2 + Offset;
    public double Upper => Center + Width / 2 + Offset;
}
public readonly record struct TracePoint(int Ms1Index, double Mz, double Intensity);
public sealed class Hill
{
    public int Id;
    public List<TracePoint> Points = [];
    public double Mz;
    public double MzStd;
    public double Sum;
    public int Start => Points[0].Ms1Index;
    public int End => Points[^1].Ms1Index;
}
public sealed class Feature
{
    public double[] ReportValues = [];
    public int Id;
    public double Mz;
    public int Charge;
    public double RtStart;
    public double RtEnd;
    public double RtApex;
    public double Area;
    public double ApexIntensity;
    public double EnvelopeScore;
    public int[] HillIds = [];
    public SortedDictionary<int, double> Signal = [];
    public double Mass => Mz * Charge - Charge * 1.00727646687;
}
public readonly record struct Candidate(int FeatureId, double LeftIntensity, double RightIntensity, double Intensity);
