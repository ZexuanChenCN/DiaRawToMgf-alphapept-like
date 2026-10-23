using ThermoFisher.CommonCore.Data.Business;
using ThermoFisher.CommonCore.Data.FilterEnums;
using ThermoFisher.CommonCore.Data.Interfaces;
using ThermoFisher.CommonCore.RawFileReader;

namespace DiaRawToMgf;

public sealed class RawReader : IDisposable
{
    private readonly IRawDataPlus raw;
    public int First => raw.RunHeaderEx.FirstSpectrum;
    public int Last => raw.RunHeaderEx.LastSpectrum;
    public RawReader(string path)
    {
        raw = RawFileReaderAdapter.FileFactory(path);
        if (!raw.IsOpen) { raw.Dispose(); throw new IOException($"Cannot open RAW: {path}"); }
        raw.SelectInstrument(Device.MS, 1);
    }
    public int Level(int scan) => (int)raw.GetScanEventForScanNumber(scan).MSOrder;
    public double Rt(int scan) => raw.RetentionTimeFromScanNumber(scan);
    public Ms2 Metadata(int scan, double rt)
    {
        var e = raw.GetScanEventForScanNumber(scan);
        return new Ms2(scan, rt, e.GetMass(0), e.GetIsolationWidth(0), e.GetIsolationWidthOffset(0));
    }
    public Peak[] Peaks(int scan)
    {
        var stats = raw.GetScanStatsForScanNumber(scan);
        double[] mz, intensity;
        if (stats.IsCentroidScan)
        {
            var s = raw.GetSegmentedScanFromScanNumber(scan);
            mz = s.Positions; intensity = s.Intensities;
        }
        else
        {
            var s = raw.GetCentroidStream(scan, false);
            if (s.Masses is null || s.Masses.Length == 0)
                throw new InvalidDataException($"Scan {scan} is profile and has no centroid stream. Convert to centroid data first.");
            mz = s.Masses; intensity = s.Intensities;
        }
        if (mz.Length != intensity.Length) throw new InvalidDataException($"Scan {scan}: array length mismatch");
        var peaks = new List<Peak>(mz.Length);
        for (int i = 0; i < mz.Length; i++)
            if (double.IsFinite(mz[i]) && mz[i] > 0 && double.IsFinite(intensity[i]) && intensity[i] > 0)
                peaks.Add(new Peak(mz[i], intensity[i]));
        peaks.Sort((a, b) => a.Mz.CompareTo(b.Mz));
        return peaks.ToArray();
    }
    public void Dispose() => raw.Dispose();
}
