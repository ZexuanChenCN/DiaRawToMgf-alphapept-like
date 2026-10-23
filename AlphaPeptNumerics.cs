// Adapted from AlphaPept 0.5.3 (Apache-2.0) and Numba quicksort (BSD).
// See THIRD_PARTY_NOTICES.txt. Compatibility code intentionally preserves upstream ordering.
namespace DiaRawToMgf;

public static class ApNumerics
{
    public static int LowerBound(double[] values,double x)
    {int lo=0,hi=values.Length;while(lo<hi){int mid=(lo+hi)/2;if(values[mid]<x)lo=mid+1;else hi=mid;}return lo;}
    // NumPy 1.26.4 Windows aquicksort; its tie ordering differs from Numba's argsort.
    public static int[] NumpyArgSort(double[] a)
    {
        var r=Enumerable.Range(0,a.Length).ToArray();if(r.Length<2)return r;
        var stack=new Stack<(int Lo,int Hi,int Depth)>();stack.Push((0,r.Length-1,2*(int)Math.Log2(r.Length)));
        void Swap(int x,int y)=>(r[x],r[y])=(r[y],r[x]);
        void Heap(int lo,int hi)
        {
            int n=hi-lo+1;
            for(int root=n/2;root>=1;root--){int key=r[lo+root-1],i=root,j=2*i;
                while(j<=n){if(j<n&&a[r[lo+j-1]]<a[r[lo+j]])j++;if(!(a[key]<a[r[lo+j-1]]))break;r[lo+i-1]=r[lo+j-1];i=j;j*=2;}r[lo+i-1]=key;}
            for(int end=n;end>1;end--){Swap(lo,lo+end-1);int key=r[lo],i=1,j=2;
                while(j<end){if(j+1<end&&a[r[lo+j-1]]<a[r[lo+j]])j++;if(!(a[key]<a[r[lo+j-1]]))break;r[lo+i-1]=r[lo+j-1];i=j;j*=2;}r[lo+i-1]=key;}
        }
        while(stack.TryPop(out var item))
        {
            var (lo,hi,depth)=item;if(depth<0){Heap(lo,hi);continue;}
            while(hi-lo>15)
            {
                int mid=lo+((hi-lo)>>1);
                if(a[r[mid]]<a[r[lo]])Swap(mid,lo);if(a[r[hi]]<a[r[mid]])Swap(hi,mid);if(a[r[mid]]<a[r[lo]])Swap(mid,lo);
                double pivot=a[r[mid]];int i=lo,j=hi-1;Swap(mid,j);
                while(true){do{i++;}while(a[r[i]]<pivot);do{j--;}while(pivot<a[r[j]]);if(i>=j)break;Swap(i,j);}
                Swap(i,hi-1);depth--;
                if(i-lo<hi-i){stack.Push((i+1,hi,depth));hi=i-1;}else{stack.Push((lo,i-1,depth));lo=i+1;}
            }
            for(int i=lo+1;i<=hi;i++){int key=r[i],j=i;while(j>lo&&a[key]<a[r[j-1]]){r[j]=r[j-1];j--;}r[j]=key;}
        }
        return r;
    }
    public static int[] NumbaArgSort(double[] a)
    {
        var r=Enumerable.Range(0,a.Length).ToArray();
        if(r.Length<2)return r;
        var stack=new Stack<(int Lo,int Hi)>();stack.Push((0,r.Length-1));
        void Swap(int x,int y)=>(r[x],r[y])=(r[y],r[x]);
        while(stack.TryPop(out var pair))
        {
            var (lo,hi)=pair;
            while(hi-lo>=15)
            {
                int mid=(lo+hi)>>1;
                if(a[r[mid]]<a[r[lo]])Swap(lo,mid);
                if(a[r[hi]]<a[r[mid]])Swap(hi,mid);
                if(a[r[mid]]<a[r[lo]])Swap(lo,mid);
                double pivot=a[r[mid]];Swap(hi,mid);
                int i=lo,j=hi-1;
                while(true)
                {
                    while(i<hi&&a[r[i]]<pivot)i++;
                    while(j>=lo&&pivot<a[r[j]])j--;
                    if(i>=j)break;
                    Swap(i++,j--);
                }
                Swap(i,hi);
                if(hi-i>i-lo){if(hi>i)stack.Push((i+1,hi));hi=i-1;}
                else{if(i>lo)stack.Push((lo,i-1));lo=i+1;}
            }
            for(int i=lo+1;i<=hi;i++)
            {
                int k=r[i],j=i;double value=a[k];
                while(j>lo&&value<a[r[j-1]]){r[j]=r[j-1];j--;}
                r[j]=k;
            }
        }
        return r;
    }
    public static double Median(double[] x,int start,int length)
    {
        var a=new double[length];Array.Copy(x,start,a,0,length);Array.Sort(a);
        return length%2==0?(a[length/2-1]+a[length/2])/2:a[length/2];
    }
    public static void SmoothInPlace(double[] x,int window=1)
    {
        for(int i=0;i<x.Length;i++)
        {int lo=Math.Max(0,i-window),hi=Math.Min(x.Length,i+window+1);x[i]=Median(x,lo,hi-lo);}
        for(int i=0;i<x.Length;i++)
        {int lo=Math.Max(0,i-window),hi=Math.Min(x.Length,i+window+1);double sum=0;for(int j=lo;j<hi;j++)sum+=x[j];x[i]=sum/(hi-lo);}
    }
    public static int UpperBound(double[] x,double value)
    {int lo=0,hi=x.Length;while(lo<hi){int m=(lo+hi)/2;if(x[m]<=value)lo=m+1;else hi=m;}return lo;}
}

/// <summary>Numba's numpy RNG: MT19937 seeded per hill; masked-rejection randint.</summary>
public sealed class ApRandom
{
    private readonly uint[] mt=new uint[624];private int index=624;
    public ApRandom(uint seed=42)
    {mt[0]=seed;for(uint i=1;i<624;i++)mt[i]=unchecked(1812433253U*(mt[i-1]^(mt[i-1]>>30))+i);}
    private uint Next()
    {
        if(index>=624)
        {
            for(int i=0;i<624;i++)
            {uint y=(mt[i]&0x80000000U)|(mt[(i+1)%624]&0x7fffffffU);mt[i]=mt[(i+397)%624]^(y>>1)^((y&1)!=0?0x9908b0dfU:0);}
            index=0;
        }
        uint z=mt[index++];z^=z>>11;z^=(z<<7)&0x9d2c5680U;z^=(z<<15)&0xefc60000U;z^=z>>18;return z;
    }
    public int Choice(int n)
    {
        if(n==1)return 0;
        uint mask=(uint)n-1;mask|=mask>>1;mask|=mask>>2;mask|=mask>>4;mask|=mask>>8;mask|=mask>>16;
        uint value;do{value=Next()&mask;}while(value>=n);return (int)value;
    }
}
