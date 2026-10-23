// Port of AlphaPept 0.5.3 report_ (Apache-2.0).
namespace DiaRawToMgf;

public sealed record ApReport(double[,] Rows,int[,] Lookup,List<Feature> Features);
public static class AlphaPeptReport
{
    public static ApReport Calculate(ApInput q,ApHills hills,ApStats stats,ApPatterns patterns,int threads)
    {
        var rows=new double[patterns.Patterns.Length,14];var lookup=new int[q.Mz.Length,2];
        for(int p=0;p<q.Mz.Length;p++){lookup[p,0]=-1;lookup[p,1]=-1;}
        var features=new Feature[patterns.Patterns.Length];
        Parallel.For(0,features.Length,new ParallelOptions{MaxDegreeOfParallelism=threads},id=>
        {
            var pattern=patterns.Patterns[id];double mz=double.PositiveInfinity,sd=0,lo=double.PositiveInfinity,hi=double.NegativeInfinity;
            int abundant=pattern[0];
            foreach(int k in pattern){mz=Math.Min(mz,stats.Values[k,0]);sd+=stats.Values[k,1];lo=Math.Min(lo,stats.Values[k,4]);hi=Math.Max(hi,stats.Values[k,5]);if(stats.Values[k,2]>stats.Values[abundant,2])abundant=k;}
            var grid=new double[100];double step=(hi-lo)/99;
            for(int i=0;i<100;i++)grid[i]=lo+i*step;grid[99]=hi;
            var trace=new double[100];double most=double.NegativeInfinity,maxApex=0,maxArea=0;
            var feature=new Feature{Id=id,Mz=mz,Charge=patterns.Charges[id],HillIds=pattern.Select(k=>stats.Order[k]).ToArray()};
            for(int isotope=0;isotope<pattern.Length;isotope++)
            {
                int hill=stats.Order[pattern[isotope]],start=hills.Ptr[hill],end=hills.Ptr[hill+1];
                double sum=0,apex=0,area=0;
                for(int p=start;p<end;p++)
                {
                    int peak=hills.Data[p],si=stats.ScanIndex[peak];double it=q.Intensity[peak];
                    lookup[peak,0]=id;lookup[peak,1]=isotope;
                    feature.Signal[si]=feature.Signal.GetValueOrDefault(si)+it;sum+=it;apex=Math.Max(apex,it);
                    if(p>start){int before=hills.Data[p-1];area+=(q.Rt[si]-q.Rt[stats.ScanIndex[before]])*(q.Intensity[before]+it)/2;}
                }
                if(sum>most){most=sum;maxApex=apex;maxArea=area;}
                int cursor=start;
                for(int i=0;i<100;i++)
                {
                    double x=grid[i],first=q.Rt[stats.ScanIndex[hills.Data[start]]],last=q.Rt[stats.ScanIndex[hills.Data[end-1]]];
                    if(x<first||x>last)continue;
                    while(cursor+1<end&&q.Rt[stats.ScanIndex[hills.Data[cursor+1]]]<=x)cursor++;
                    int a=hills.Data[cursor];double value=q.Intensity[a];
                    if(cursor+1<end&&x!=q.Rt[stats.ScanIndex[a]])
                    {int b=hills.Data[cursor+1];double slope=(q.Intensity[b]-value)/(q.Rt[stats.ScanIndex[b]]-q.Rt[stats.ScanIndex[a]]);value=slope*(x-q.Rt[stats.ScanIndex[a]])+value;}
                    trace[i]+=value;
                }
            }
            int apexIndex=0;for(int i=1;i<100;i++)if(trace[i]>trace[apexIndex])apexIndex=i;
            int Nearest(int start,int end,double target){int best=start;for(int i=start+1;i<end;i++)if(Math.Abs(trace[i]-target)<Math.Abs(trace[best]-target))best=i;return best;}
            int left=apexIndex==0?0:Nearest(0,apexIndex,trace[apexIndex]/2),right=Nearest(apexIndex,100,trace[apexIndex]/2);
            int rtStart=apexIndex==0?0:Nearest(0,apexIndex,trace[apexIndex]*(1-0.95)),rtEnd=Nearest(apexIndex,100,trace[apexIndex]*(1-0.95));
            double sumArea=0;for(int i=rtStart+1;i<rtEnd;i++)sumArea+=(grid[i]-grid[i-1])*(trace[i-1]+trace[i])/2;
            feature.RtStart=grid[rtStart];feature.RtEnd=grid[rtEnd];feature.RtApex=grid[apexIndex];feature.ApexIntensity=trace[apexIndex];feature.Area=sumArea;
            double[] row=[mz,sd/pattern.Length,stats.Values[abundant,0],feature.Charge,feature.RtStart,feature.RtApex,feature.RtEnd,grid[right]-grid[left],pattern.Length,feature.Mass,feature.ApexIntensity,sumArea,maxApex,maxArea];
            feature.ReportValues=row;for(int j=0;j<14;j++)rows[id,j]=row[j];features[id]=feature;
        });
        return new(rows,lookup,features.ToList());
    }
}
