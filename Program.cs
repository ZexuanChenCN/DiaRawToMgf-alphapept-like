using System.Collections.Concurrent;
using System.Diagnostics;
using System.Reflection;
using System.Text.Json;
using DiaRawToMgf;

return await Entry.Run(args);

internal static class Entry
{
    public static async Task<int> Run(string[] args)
    {
        try
        {
            if(args.Length==2 && args[0]=="--verify-alphapept-stages")
            {AlphaPeptStageVerification.Verify(args[1]);return 0;}
            if(args.Length==2 && args[0]=="--worker-json")
            {
                var worker=JsonSerializer.Deserialize<Options>(File.ReadAllText(args[1]))??throw new ArgumentException("Invalid worker config");
                Pipeline.Run(worker);return 0;
            }
            var o=Options.Parse(args);
            if(o.Help){Console.WriteLine(Help);return 0;}
            if(o.SelfTest){SelfTests.Run();return 0;}
            string[] files;
            if(File.Exists(o.Input))files=Path.GetExtension(o.Input).Equals(".raw",StringComparison.OrdinalIgnoreCase)?[o.Input]:throw new ArgumentException("Input file must be .raw");
            else if(Directory.Exists(o.Input))files=Directory.EnumerateFiles(o.Input,"*",o.Recursive?SearchOption.AllDirectories:SearchOption.TopDirectoryOnly)
                .Where(f=>Path.GetExtension(f).Equals(".raw",StringComparison.OrdinalIgnoreCase)).Order(StringComparer.Ordinal).ToArray();
            else throw new DirectoryNotFoundException(o.Input);
            if(files.Length==0)throw new ArgumentException("No .raw files found");
            if(o.FeatureCache.Length>0 && files.Length!=1)throw new ArgumentException("--feature-cache is a single-RAW comparison option, not a shared cache for a batch");
            Directory.CreateDirectory(o.Output);
            var children=new ConcurrentDictionary<int,Process>();
            using var cancel=new CancellationTokenSource();
            ConsoleCancelEventHandler handler=(_,e)=>{e.Cancel=true;cancel.Cancel();foreach(var child in children.Values)try{child.Kill(true);}catch(InvalidOperationException){}};
            Console.CancelKeyPress+=handler;
            var results=new ConcurrentBag<object>();
            int failures=0;
            try
            {
                await Parallel.ForEachAsync(files,new ParallelOptions{MaxDegreeOfParallelism=o.Jobs,CancellationToken=cancel.Token},async(file,ct)=>
                {
                    string relative=File.Exists(o.Input)?Path.GetFileName(file):Path.GetRelativePath(o.Input,file);
                    string directory=Path.Combine(o.Output,Path.ChangeExtension(relative,null));
                    Directory.CreateDirectory(directory);
                    var childOptions=o with{Input=Path.GetFullPath(file),Output=Path.GetFullPath(directory),Worker=true,Jobs=1};
                    string config=Path.Combine(directory,"worker-options.json");
                    await File.WriteAllTextAsync(config,JsonSerializer.Serialize(childOptions,Pipeline.JsonOptions),ct);
                    var start=new ProcessStartInfo{FileName=Environment.ProcessPath!,UseShellExecute=false,CreateNoWindow=true,
                        RedirectStandardOutput=true,RedirectStandardError=true};
                    if(Path.GetFileNameWithoutExtension(start.FileName).Equals("dotnet",StringComparison.OrdinalIgnoreCase))start.ArgumentList.Add(Assembly.GetExecutingAssembly().Location);
                    start.ArgumentList.Add("--worker-json");start.ArgumentList.Add(config);
                    using var child=Process.Start(start)??throw new IOException("Cannot start worker");
                    children[child.Id]=child;
                    using var childLog=new StreamWriter(Path.Combine(directory,"worker.log"),append:true){AutoFlush=true};
                    async Task Drain(StreamReader reader)
                    {
                        while(await reader.ReadLineAsync(ct) is { } line)
                        {lock(childLog){childLog.WriteLine(line);}Console.WriteLine(line);}
                    }
                    try{await Task.WhenAll(child.WaitForExitAsync(ct),Drain(child.StandardOutput),Drain(child.StandardError));}
                    finally{children.TryRemove(child.Id,out _);if(!child.HasExited)child.Kill(true);}
                    if(child.ExitCode!=0)Interlocked.Increment(ref failures);
                    results.Add(new{Raw=file,Output=directory,ExitCode=child.ExitCode});
                });
            }
            finally{Console.CancelKeyPress-=handler;}
            await File.WriteAllTextAsync(Path.Combine(o.Output,"batch-summary.json"),JsonSerializer.Serialize(results,Pipeline.JsonOptions));
            Console.WriteLine($"Batch completed: {files.Length-failures}/{files.Length} succeeded");
            return failures==0?0:1;
        }
        catch(OperationCanceledException){Console.Error.WriteLine("Cancelled");return 130;}
        catch(Exception e){Console.Error.WriteLine(e.ToString());return 1;}
    }
    private const string Help="""
    Pure C# DIA RAW -> MS1 features -> per-MS2 feature candidates -> MGF
    Usage: DiaRawToMgf --input RAW_FOLDER --output OUTPUT_FOLDER [options]
    Also accepts one .raw file. Original RAW scan numbers are always preserved.
      --recursive                  Include nested folders (default: false)
      --jobs N                     Concurrent RAW worker processes (default: 1)
      --threads N                  Feature computation threads per RAW (default: 4)
      --top-peaks N                Highest-intensity MS2 peaks, then m/z sorted (250)
      --max-candidates N           Per-MS2 local-intensity ranked candidates (5; 0=all)
      --centroid-ppm X             Hill centroid linking tolerance (8)
      --max-gap N                  Missing MS1 scans permitted in a hill (2)
      --min-hill-points N          Minimum MS1 trace points (3)
      --hill-split-level X         Peak/valley split ratio (1.3)
      --isotope-split-level X      Isotope abundance valley split ratio (1.3)
      --min-charge N               Minimum candidate charge (1)
      --max-charge N               Maximum candidate charge (6)
      --isotope-correlation X      Minimum hill trace cosine (0.6)
      --envelope-score X           Minimum CHNOS averagine cosine (0.6)
      --rt-start MIN --rt-end MIN  Optional RT interval; default entire acquisition
      --resume                     Skip completed files only when fingerprint matches
      --feature-cache DIR          Single-RAW comparison: reuse feature/curve TSV files, no Python
      --self-test                  Run synthetic pipeline correctness checks
    Default feature finding ports AlphaPept 0.5.3, including bootstrap and ordering.
    No DIA fragment deconvolution. No Python / Mono required.
    Requires .NET 10 AND ASP.NET Core 10 runtimes (both included in the .NET 10 SDK).
    """;
}
