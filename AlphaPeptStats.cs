// Port of AlphaPept 0.5.3 hill_stats/get_hill_data (Apache-2.0).
using System.Collections.Concurrent;
namespace DiaRawToMgf;

public sealed record ApStats(double[,] Values,int[] Order,int[] Upper,int[] ScanIndex,int[] FirstScan,int[] LastScan);
public static class AlphaPeptStats
{
    public static ApStats Calculate(ApInput q,ApHills hills,int threads,int boots=150,int maxBootSize=300)
    {
        int n=hills.Ptr.Length-1;var raw=new double[n,6];var scanIndex=new int[q.Mz.Length];
        for(int si=0;si<q.Rt.Length;si++)Array.Fill(scanIndex,si,q.Bounds[si],q.Bounds[si+1]-q.Bounds[si]);
        var choices=new ConcurrentDictionary<int,int[]>();
        Parallel.For(0,n,new ParallelOptions{MaxDegreeOfParallelism=threads},idx=>
        {
            int start=hills.Ptr[idx],end=hills.Ptr[idx+1],length=end-start,size=Math.Min(length,maxBootSize);
            var indices=choices.GetOrAdd(length,len=>{var rng=new ApRandom();var result=new int[boots*Math.Min(len,maxBootSize)];
                for(int i=0;i<result.Length;i++)result[i]=rng.Choice(len);return result;});
            var averages=new double[boots];double mean=0;
            for(int i=0;i<boots;i++)
            {
                double numerator=0,denominator=0;
                for(int j=0;j<size;j++){int peak=hills.Data[start+indices[i*size+j]];numerator+=q.Mz[peak]*q.Intensity[peak];denominator+=q.Intensity[peak];}
                averages[i]=numerator/denominator;mean+=averages[i];
            }
            mean/=boots;double delta=0;foreach(double value in averages)delta+=(mean-value)*(mean-value);
            double apex=0,area=0,rtMin=double.PositiveInfinity,rtMax=double.NegativeInfinity;
            for(int j=start;j<end;j++)
            {
                int p=hills.Data[j];double rt=q.Rt[scanIndex[p]],it=q.Intensity[p];
                apex=Math.Max(apex,it);rtMin=Math.Min(rtMin,rt);rtMax=Math.Max(rtMax,rt);
                if(j>start){int before=hills.Data[j-1];area+=(q.Rt[scanIndex[p]]-q.Rt[scanIndex[before]])*(q.Intensity[before]+it)/2;}
            }
            raw[idx,0]=mean;raw[idx,1]=Math.Sqrt(delta/(boots-1));raw[idx,2]=Math.Abs(area);raw[idx,3]=apex;raw[idx,4]=rtMin;raw[idx,5]=rtMax;
        });
        var starts=Enumerable.Range(0,n).Select(i=>raw[i,4]).ToArray();var order=ApNumerics.NumpyArgSort(starts);
        var values=new double[n,6];var sortedStarts=new double[n];
        for(int i=0;i<n;i++){for(int j=0;j<6;j++)values[i,j]=raw[order[i],j];sortedStarts[i]=values[i,4];}
        var upper=new int[n];for(int i=0;i<n;i++)upper[i]=ApNumerics.UpperBound(sortedStarts,values[i,5]);
        var firstScan=new int[n];var lastScan=new int[n];
        for(int i=0;i<n;i++){int h=order[i];firstScan[i]=scanIndex[hills.Data[hills.Ptr[h]]];lastScan[i]=scanIndex[hills.Data[hills.Ptr[h+1]-1]];}
        return new(values,order,upper,scanIndex,firstScan,lastScan);
    }
}
