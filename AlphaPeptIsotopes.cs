// Literal numerical adaptation of AlphaPept 0.5.3 feature_finding.py (Apache-2.0).
namespace DiaRawToMgf;

public sealed record ApPatterns(int[][] Groups,int[][] Patterns,int[] Charges);
public static class AlphaPeptIsotopes
{
    internal static List<int> GrowSparse(double[,] values,int[] pattern,int seed,int charge,int[] correlatedPositions,double massRange)
    {
        int anchor=pattern[seed];double mz=values[anchor,0];var trail=new List<int>{anchor};
        foreach(int direction in new[]{-1,1})
        {
            int index=1,cursor=seed+direction;
            IEnumerable<int> candidates=direction>0?correlatedPositions.Where(p=>p>seed):correlatedPositions.Reverse().Where(p=>p<seed);
            foreach(int pos in candidates)
            {
                // The first non-correlating entry outside the original mass bound
                // would stop grow(). Keep that boundary even when skipping zeros.
                double limit=(Delta+DeltaStd)*index;int lo=direction>0?cursor:0,hi=direction>0?pattern.Length:cursor+1;
                while(lo<hi)
                {
                    int mid=(lo+hi)/2;bool beyond=Math.Abs(mz-values[pattern[mid],0])>limit;
                    if(direction>0){if(beyond)hi=mid;else lo=mid+1;}
                    else{if(beyond)lo=mid+1;else hi=mid;}
                }
                int boundary=direction>0?lo:lo-1;
                if(direction>0?pos>boundary:pos<boundary)break;
                int y=pattern[pos];
                if(Spacing(values,anchor,y,charge,-direction*index,massRange))
                {if(direction>0)trail.Add(y);else trail.Insert(0,y);index++;}
                if(Math.Abs(mz-values[y,0])>(Delta+DeltaStd)*index)break;
                cursor=pos+direction;
            }
        }
        return trail;
    }
    internal static int[] SeedOrder(double[,] values,int[] pattern,int[] areaSorted,HashSet<int> active,int seeds)
    {
        var top=areaSorted.Where(active.Contains).Take(seeds+1).ToArray();bool tied=false;
        for(int i=1;i<top.Length;i++)if(values[top[i-1],2]==values[top[i],2]){tied=true;break;}
        if(tied){var order=ApNumerics.NumbaArgSort(pattern.Select(x=>values[x,2]).ToArray());Array.Reverse(order);return order;}
        int Position(int hill)
        {
            double mz=values[hill,0];int lo=0,hi=pattern.Length;
            while(lo<hi){int mid=(lo+hi)/2;if(values[pattern[mid],0]<mz)lo=mid+1;else hi=mid;}
            while(pattern[lo]!=hill)lo++;return lo;
        }
        return top.Take(seeds).Select(Position).ToArray();
    }
    public const double Delta=1.00286864,DeltaStd=0.0109135,Proton=1.00727646687;
    public static double Correlate(ApInput q,ApHills hills,ApStats stats,int a,int b)
    {
        int h1=stats.Order[a],h2=stats.Order[b],s1=hills.Ptr[h1],e1=hills.Ptr[h1+1],s2=hills.Ptr[h2],e2=hills.Ptr[h2+1];
        var scan=stats.ScanIndex;
        if(scan[hills.Data[s1]]+3>scan[hills.Data[e2-1]]||scan[hills.Data[s2]]+3>scan[hills.Data[e1-1]])return 0;
        double aa=0,bb=0,ab=0;
        for(int p=s1;p<e1;p++){double v=q.Intensity[hills.Data[p]];aa+=v*v;}
        for(int p=s2;p<e2;p++){double v=q.Intensity[hills.Data[p]];bb+=v*v;}
        int i=s1,j=s2;
        while(i<e1&&j<e2)
        {
            int x=hills.Data[i],y=hills.Data[j];
            if(scan[x]==scan[y]){ab+=q.Intensity[x]*q.Intensity[y];i++;j++;}
            else if(scan[x]<scan[y])i++;else j++;
        }
        return ab/Math.Sqrt(aa*bb);
    }
    public static bool Spacing(double[,] stats,int a,int b,int charge,int isotopeIndex=0,double massRange=5)
    {
        double mass1=stats[a,0],mass2=stats[b,0],d1=stats[a,1]*massRange,d2=stats[b,1]*massRange;
        double error=isotopeIndex==0?Math.Abs(Math.Abs(mass2-mass1)-Delta/charge):Math.Abs(mass1-mass2-isotopeIndex*Delta/charge);
        return error<=Math.Sqrt(Math.Pow(DeltaStd/charge,2)+d1*d1+d2*d2);
    }
    private static int Lower(List<int> values,int x)
    {int lo=0,hi=values.Count;while(lo<hi){int m=(lo+hi)/2;if(values[m]<x)lo=m+1;else hi=m;}return lo;}
    public static int[][] PrePatterns(ApInput q,ApHills hills,ApStats stats,int minCharge=1,int maxCharge=6,double massRange=5,double cutoff=0.6)
    {
        int n=stats.Order.Length;double width=Delta+DeltaStd;
        var bins=new Dictionary<int,List<int>>();
        for(int i=0;i<n;i++){int b=(int)Math.Floor(stats.Values[i,0]/width);if(!bins.TryGetValue(b,out var list))bins[b]=list=[];list.Add(i);}
        var adjacency=new List<int>?[n];var insertion=new List<int>();
        for(int runner=0;runner<n;runner++)
        {
            int bin=(int)Math.Floor(stats.Values[runner,0]/width);var candidates=new List<int>();
            for(int b=bin-1;b<=bin+1;b++)if(bins.TryGetValue(b,out var values))
                for(int p=Lower(values,runner+1);p<values.Count&&values[p]<stats.Upper[runner];p++)candidates.Add(values[p]);
            candidates.Sort();
            foreach(int j in candidates)
            {
                if(Math.Abs(stats.Values[runner,0]-stats.Values[j,0])>width)continue;
                bool matches=false;for(int z=minCharge;z<=maxCharge;z++)if(Spacing(stats.Values,runner,j,z,0,massRange)){matches=true;break;}
                if(!matches||Correlate(q,hills,stats,runner,j)<=cutoff)continue;
                if(adjacency[runner]is null){adjacency[runner]=[];insertion.Add(runner);}
                if(adjacency[j]is null){adjacency[j]=[];insertion.Add(j);}
                adjacency[runner]!.Add(j);adjacency[j]!.Add(runner);
            }
        }
        var visited=new bool[n];var groups=new List<int[]>();
        foreach(int node in insertion)if(!visited[node])
        {
            var members=new List<int>();var queue=new Queue<int>();queue.Enqueue(node);visited[node]=true;
            while(queue.TryDequeue(out int current)){members.Add(current);foreach(int next in adjacency[current]!)if(!visited[next]){visited[next]=true;queue.Enqueue(next);}}
            members.Sort();groups.Add(members.ToArray());
        }
        return groups.OrderByDescending(g=>g.Length).ToArray();
    }
    public static ApPatterns Extract(ApInput q,ApHills hills,ApStats stats,int[][] groups,int threads,
        int minCharge=1,int maxCharge=6,int seeds=100,double massRange=5,double cutoff=0.6,double splitLevel=1.3,double envelopeScore=0.6,Action<string>? log=null)
    {
        var extracted=new List<(int[] Pattern,int Charge)>[groups.Length];int completed=0;
        log?.Invoke($"Largest isotope group: {(groups.Length>0?groups[0].Length:0):N0} hills");
        Parallel.For(0,groups.Length,new ParallelOptions{MaxDegreeOfParallelism=threads},groupIndex=>
        {
            var pre=groups[groupIndex].ToList();var found=new List<(int[],int)>();
            var neighbors=new Dictionary<int,int[]>();
            var envelopes=new Dictionary<double,AlphaPeptChem.Envelope>();
            var massSorted=ApNumerics.NumbaArgSort(pre.Select(x=>stats.Values[x,0]).ToArray()).Select(i=>pre[i]).ToArray();
            bool uniqueMass=true;for(int i=1;i<massSorted.Length;i++)if(stats.Values[massSorted[i-1],0]==stats.Values[massSorted[i],0]){uniqueMass=false;break;}
            var areaSorted=pre.OrderByDescending(i=>stats.Values[i,2]).ToArray();var active=pre.ToHashSet();
            int[] Correlated(int anchor)
            {
                if(neighbors.TryGetValue(anchor,out var cached))return cached;
                var positive=new List<int>();
                foreach(int other in groups[groupIndex])
                {
                    if(stats.FirstScan[other]+3>stats.LastScan[anchor])break;
                    if(other==anchor||stats.FirstScan[anchor]+3>stats.LastScan[other])continue;
                    if(Correlate(q,hills,stats,anchor,other)>cutoff)positive.Add(other);
                }
                return neighbors[anchor]=positive.ToArray();
            }
            while(pre.Count>1)
            {
                // With distinct keys the numerical sorted order is independent of the
                // quicksort partition path. Preserve the original sort whenever ties
                // could affect either hill order or any of the selected seeds.
                var pattern=uniqueMass?massSorted.Where(active.Contains).ToArray():
                    ApNumerics.NumbaArgSort(pre.Select(x=>stats.Values[x,0]).ToArray()).Select(i=>pre[i]).ToArray();
                var seedOrder=SeedOrder(stats.Values,pattern,areaSorted,active,seeds);
                var positions=pattern.Select((hill,pos)=>(hill,pos)).ToDictionary(x=>x.hill,x=>x.pos);
                int[]? champion=null;int championCharge=0,longest=0;double championIntensity=0;
                foreach(int seed in seedOrder.Take(seeds))
                {
                    int anchor=pattern[seed];
                    var correlated=Correlated(anchor).Where(positions.ContainsKey).Select(h=>positions[h]).Order().ToArray();
                    for(int z=minCharge;z<=maxCharge;z++)
                    {
                        var trail=GrowSparse(stats.Values,pattern,seed,z,correlated,massRange);
                        if(trail.Count<longest)continue;
                        int seedPosition=trail.IndexOf(anchor),left=0,right=trail.Count;
                        for(int i=1;i<trail.Count-1;i++)
                        {
                            double valley=stats.Values[trail[i],2];
                            if(!(stats.Values[trail[i-1],2]>valley&&stats.Values[trail[i+1],2]>valley))continue;
                            double lm=0,rm=0;
                            for(int j=0;j<i;j++)lm=Math.Max(lm,stats.Values[trail[j],2]);
                            for(int j=i+1;j<trail.Count;j++)rm=Math.Max(rm,stats.Values[trail[j],2]);
                            if(Math.Min(lm,rm)/valley>=splitLevel)
                            {if(i<seedPosition)left=i;else if(i>seedPosition&&right==trail.Count)right=i;}
                        }
                        var arr=trail.Skip(left).Take(Math.Min(trail.Count,right+1)-left).ToArray();
                        if(z*stats.Values[anchor,0]<1000)
                        {int max=0;for(int i=1;i<arr.Length;i++)if(stats.Values[arr[i],2]>stats.Values[arr[max],2])max=i;arr=arr[max..];}
                        double intensity=0;foreach(int i in arr)intensity+=stats.Values[i,2];
                        if(arr.Length>longest||(arr.Length==longest&&intensity>championIntensity))
                            if(AlphaPeptChem.Check(stats.Values,arr,z,envelopes)>envelopeScore)
                            {champion=arr;championCharge=z;longest=arr.Length;championIntensity=intensity;}
                    }
                }
                if(champion is null||champion.Length<=1)break;
                found.Add((champion,championCharge));var used=champion.ToHashSet();pre.RemoveAll(x=>used.Contains(x));active.ExceptWith(used);
                if(groupIndex==0&&found.Count%100==0)log?.Invoke($"Largest isotope group: {found.Count:N0} features extracted, {pre.Count:N0} hills remaining; distinct m/z keys={uniqueMass}");
            }
            extracted[groupIndex]=found;
            int count=Interlocked.Increment(ref completed);if(count%10000==0)log?.Invoke($"Isotope groups completed {count:N0}/{groups.Length:N0}");
        });
        var flat=extracted.SelectMany(x=>x).ToArray();
        return new(groups,flat.Select(x=>x.Pattern).ToArray(),flat.Select(x=>x.Charge).ToArray());
    }
}
