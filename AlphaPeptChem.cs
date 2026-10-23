// AlphaPept 0.5.3 chem.py/constants.py/check_averagine numerical adaptation (Apache-2.0).
namespace DiaRawToMgf;

public static class AlphaPeptChem
{
    public sealed record Envelope(Dictionary<int,double> Intensities,double SquaredSum);
    private sealed record Distribution(float Mass,int Count,double[] Intensities);
    private static readonly double[] Average=[4.9384,7.7583,1.3577,1.4773,0.0417]; // C H N O S insertion order
    private static readonly float[] Masses=[12f,(float)1.007940,(float)14.0030740052,(float)15.9949146221,(float)31.97207069];
    private static readonly float[][] Abundance=[
        [0.9893f,0.0107f,0f],[0.999885f,0.000115f,0f],[0.99636f,0.00364f],
        [0.99757f,0.00038f,0.00205f],[0.9499f,0.0075f,0.0425f,0.0001f]];
    private static Distribution Add(Distribution a,Distribution b)
    {
        var values=new double[a.Count+b.Count-1];
        for(int i=0;i<a.Count;i++)for(int j=0;j<b.Count;j++)values[i+j]+=a.Intensities[i]*b.Intensities[j];
        double maximum=values.Max();int count=values.Length;
        while(count>0&&values[count-1]<0.000001)count--;
        for(int i=0;i<values.Length;i++)values[i]/=maximum;
        return new((float)(a.Mass+b.Mass),count,values);
    }
    private static Distribution Multiply(Distribution a,int n)
    {
        if(n==1)return a;
        var result=new Distribution(0,1,[1.0]);var multiple=a;
        while(n>0){if((n&1)!=0)result=Add(result,multiple);multiple=Add(multiple,multiple);n>>=1;}
        return result;
    }
    public static (double[] Masses,double[] Intensity) MassToDist(double mass)
    {
        var counts=new int[5];double total=0;
        for(int i=0;i<5;i++){counts[i]=(int)Math.Round(mass/111.1254*Average[i]);total+=(double)counts[i]*Masses[i];}
        counts[1]+=(int)Math.Round((mass-total)/Masses[1]);
        var dist=new Distribution(0,1,[1.0]);
        for(int i=0;i<5;i++)
        {
            var isotope=new Distribution(Masses[i],Abundance[i].Length,Abundance[i].Select(x=>(double)x).ToArray());
            // x.add(isotopes[element]) starts from the identity, normalizing/pruning once.
            var x=Add(new Distribution(0,1,[1.0]),isotope);
            dist=Add(dist,Multiply(x,counts[i]));
        }
        // mass_to_dist uses len(intensities), including the pruned tail kept in storage.
        return (Enumerable.Range(0,dist.Intensities.Length).Select(i=>dist.Mass+i*AlphaPeptIsotopes.Delta).ToArray(),dist.Intensities);
    }
    public static double Check(double[,] stats,int[] pattern,int charge,Dictionary<double,Envelope>? cache=null)
    {
        var order=ApNumerics.NumbaArgSort(pattern.Select(i=>stats[i,0]*charge-charge*AlphaPeptIsotopes.Proton).ToArray());
        var observed=new Dictionary<int,double>();double minimum=double.PositiveInfinity;
        foreach(int i in order){int index=pattern[i];double mass=stats[index,0]*charge-charge*AlphaPeptIsotopes.Proton;minimum=Math.Min(minimum,mass);observed[(int)Math.Floor(mass)]=stats[index,2];}
        Envelope envelope;
        if(cache is null||!cache.TryGetValue(minimum,out envelope!))
        {
            var (masses,intensities)=MassToDist(minimum);var theoretical=new Dictionary<int,double>();
            for(int i=0;i<masses.Length;i++)theoretical[(int)Math.Floor(masses[i])]=intensities[i];
            double bb=0;foreach(double value in theoretical.Values)bb+=value*value;
            envelope=new(theoretical,bb);
            // Reusing a pure calculation never changes its operation order or output.
            if(cache is not null&&cache.Count<16384)cache[minimum]=envelope;
        }
        double ab=0,aa=0;
        foreach(var (key,value) in observed){aa+=value*value;ab+=value*envelope.Intensities.GetValueOrDefault(key);}
        return ab/Math.Sqrt(aa*envelope.SquaredSum);
    }
}
