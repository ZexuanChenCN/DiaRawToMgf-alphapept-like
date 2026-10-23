// C# adaptation of alphapept/feature_finding.py 0.5.3, Apache-2.0.
namespace DiaRawToMgf;

public sealed record ApInput(int[] Bounds,double[] Mz,double[] Intensity,double[] Rt,int[] ScanNumbers)
{
    public static ApInput From(List<Ms1> scans)
    {
        int n=checked(scans.Sum(s=>s.Peaks.Length));var bounds=new int[scans.Count+1];
        var mz=new double[n];var it=new double[n];int at=0;
        for(int i=0;i<scans.Count;i++)
        {foreach(var p in scans[i].Peaks){mz[at]=p.Mz;it[at++]=p.Intensity;}bounds[i+1]=at;}
        return new(bounds,mz,it,scans.Select(s=>s.Rt).ToArray(),scans.Select(s=>s.Scan).ToArray());
    }
}
public sealed record ApHills(int[] Ptr,int[] Data,int[] Counts,double Median,double Std);
public static class AlphaPeptHills
{
    public static ApHills Extract(ApInput q,int gap,double tolerance,int minLength=3)
    {
        var forward=new int[q.Mz.Length];Array.Fill(forward,-1);
        var backward=new bool[q.Mz.Length];var scores=new List<double>();
        for(int x=0;x<q.Rt.Length;x++)
        {
            int count=q.Bounds[x+1]-q.Bounds[x];var chosen=new int[count];Array.Fill(chosen,-1);
            for(int g=0;g<=gap;g++)
            {
                int y=x+g+1;if(y>=q.Rt.Length)break;
                var best=new double[count];Array.Fill(best,double.PositiveInfinity);
                var targets=new int[count];Array.Fill(targets,-1);
                int i=q.Bounds[x],j=q.Bounds[y];
                while(i<q.Bounds[x+1]&&j<q.Bounds[y+1])
                {
                    double diff=q.Mz[i]-q.Mz[j];
                    double delta=2e6*Math.Abs(diff)/(q.Mz[i]+q.Mz[j]);
                    int local=i-q.Bounds[x];
                    if(delta<tolerance&&best[local]>delta){best[local]=delta;targets[local]=j;}
                    if(diff>0)j++;else i++;
                }
                for(int k=0;k<count;k++)if(targets[k]>=0)
                {scores.Add(best[k]);if(chosen[k]<0)chosen[k]=targets[k];}
            }
            for(int k=0;k<count;k++)if(chosen[k]>=0)
            {forward[q.Bounds[x]+k]=chosen[k];backward[chosen[k]]=true;}
        }
        var ptr=new List<int>{0};var data=new List<int>();var lengths=new List<int>();
        for(int i=0;i<forward.Length;i++)if(forward[i]>=0&&!backward[i])
        {
            int length=1;for(int next=forward[i];forward[next]>=0;next=forward[next])length++;
            length++; // Above counts edges excluding terminal node.
            if(length<minLength)continue;
            for(int next=i;next>=0;next=forward[next])data.Add(next);
            lengths.Add(length);ptr.Add(data.Count);
        }
        double median=double.NaN,std=double.NaN;
        if(scores.Count>0)
        {
            var all=scores.ToArray();median=ApNumerics.Median(all,0,all.Length);
            double mean=all.Sum()/all.Length,sum=0;foreach(double score in all)sum+=(score-mean)*(score-mean);std=Math.Sqrt(sum/all.Length);
        }
        return new(ptr.ToArray(),data.ToArray(),lengths.ToArray(),median,std);
    }
    public static ApHills Deduplicate(ApHills h)
    {
        if(h.Data.Length==0)return h;
        var order=ApNumerics.NumbaArgSort(h.Counts.Select(x=>(double)x).ToArray());Array.Reverse(order);
        var taken=new bool[h.Data.Max()+1];var ptr=new List<int>{0};var data=new List<int>();
        foreach(int i in order)
        {
            bool free=true;for(int p=h.Ptr[i];p<h.Ptr[i+1];p++)if(taken[h.Data[p]])free=false;
            if(free){for(int p=h.Ptr[i];p<h.Ptr[i+1];p++)data.Add(h.Data[p]);ptr.Add(data.Count);}
            // Upstream marks points even when the candidate hill was rejected.
            for(int p=h.Ptr[i];p<h.Ptr[i+1];p++)taken[h.Data[p]]=true;
        }
        // Preserve upstream [:c] rather than [:c+1], including its dropped last hill.
        ptr.RemoveAt(ptr.Count-1);
        return h with{Ptr=ptr.ToArray(),Data=data.ToArray(),Counts=[]};
    }
    public static int[] Split(ApHills h,double[] intensity,double level,int window)
    {
        var ptr=h.Ptr.ToList();var check=Enumerable.Range(0,ptr.Count-1).ToList();
        while(check.Count>0)
        {
            var splits=new Dictionary<int,int>();
            foreach(int k in check)
            {
                int start=ptr[k],end=ptr[k+1];var trace=new double[end-start];
                for(int i=0;i<trace.Length;i++)trace[i]=intensity[h.Data[start+i]];
                ApNumerics.SmoothInPlace(trace,window);
                var minima=new List<int>();
                for(int i=2;i<trace.Length-2;i++)
                    if((trace[i-1]>trace[i]&&trace[i+1]>trace[i])||
                       (trace[i-1]>trace[i]&&trace[i+1]==trace[i]&&trace[i+2]>trace[i])||
                       (trace[i-2]>trace[i]&&trace[i-1]==trace[i]&&trace[i+1]>trace[i])||
                       (trace[i-2]>trace[i]&&trace[i-1]==trace[i]&&trace[i+1]==trace[i]&&trace[i+2]>trace[i]))minima.Add(i);
                foreach(int index in ApNumerics.NumbaArgSort(minima.Select(i=>trace[i]).ToArray()))
                {
                    int cut=minima[index];double left=double.NegativeInfinity,right=double.NegativeInfinity;
                    for(int i=0;i<cut;i++)left=Math.Max(left,trace[i]);for(int i=cut;i<trace.Length;i++)right=Math.Max(right,trace[i]);
                    if(trace[cut]==0||Math.Min(left,right)/trace[cut]>level){splits[k]=start+cut;break;}
                }
            }
            var nextPtr=new List<int>();var nextCheck=new List<int>();
            for(int k=0;k<ptr.Count-1;k++)
            {
                int n=nextPtr.Count;nextPtr.Add(ptr[k]);
                if(splits.TryGetValue(k,out int cut)){nextCheck.Add(n);nextPtr.Add(cut);nextCheck.Add(n+1);}
            }
            nextPtr.Add(ptr[^1]);ptr=nextPtr;check=nextCheck;
        }
        return ptr.ToArray();
    }
    public static ApHills Filter(ApHills h,double[] intensity,int large=40)
    {
        var largeIndices=Enumerable.Range(0,h.Ptr.Length-1).Where(i=>h.Ptr[i+1]-h.Ptr[i]>=large).ToArray();
        var remove=new HashSet<int>();
        for(int i=0;i<largeIndices.Length;i++)
        {
            int k=largeIndices[i],start=h.Ptr[k],end=h.Ptr[k+1];var trace=new double[end-start];
            for(int j=0;j<trace.Length;j++)trace[j]=intensity[h.Data[start+j]];
            double maximum=trace.Max();ApNumerics.SmoothInPlace(trace,1);
            // filter_hills indexes the flags, not large_peaks[flags]; preserve 0.5.3 behavior.
            if(!(maximum/trace[0]>2&&maximum/trace[^1]>2))remove.Add(i);
        }
        var ptr=new List<int>{0};var data=new List<int>();
        for(int i=0;i<h.Ptr.Length-1;i++)if(!remove.Contains(i))
        {for(int j=h.Ptr[i];j<h.Ptr[i+1];j++)data.Add(h.Data[j]);ptr.Add(data.Count);}
        // Upstream also retains data beyond the last pointer (deduplication's orphan tail).
        for(int j=h.Ptr[^1];j<h.Data.Length;j++)data.Add(h.Data[j]);
        return h with{Ptr=ptr.ToArray(),Data=data.ToArray(),Counts=[]};
    }
}
