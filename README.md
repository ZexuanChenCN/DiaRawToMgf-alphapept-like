# DiaRawToMgf — AlphaPept feature finding 的 C# 移植

面向 **Windows / Windows 服务器**的 Thermo DIA `.raw` 批处理工具：读取原始 MS1 → 提取 feature → 为每张原始 MS2 匹配前体候选 → 导出可用于 pFind 等搜索引擎的 MGF。

- **纯 C#，目标框架 .NET 10**；正常构建和运行无需 Python、AlphaPept、NumPy、Numba 或 Mono。
- feature finding 数值行为对齐 **AlphaPept 0.5.3**，包括 bootstrap、同位素簇选择及排序细节。
- 多个 RAW 用独立进程并行；单个 RAW 的部分计算使用多线程。
- 保留原始 MS2 `scan_num`；MS1/MS2 按 m/z 升序；MS2 默认取强度最高的 250 根峰，再按 m/z 升序输出。
- 完整 Astral 测试中，C# 与独立 Python 流程导出的 MGF **大小、SHA-256 和逐字节内容全部相同**。这是该测试文件及默认参数下的验证结果，不代表对任意数据和 AlphaPept 版本的保证。

## 1. 安装要求

在 Windows 电脑或服务器安装：

1. [Git for Windows](https://git-scm.com/downloads/win)。
2. [.NET 10 SDK x64](https://dotnet.microsoft.com/en-us/download/dotnet/10.0)；选择 **SDK**，它包含所需的 .NET 和 ASP.NET Core 运行时。

安装后重新打开 PowerShell，检查：

```powershell
git --version
dotnet --version
dotnet --list-runtimes
```

SDK 版本应为 `10.0.x`；运行时列表应包含 `Microsoft.NETCore.App 10.0.x` 和 `Microsoft.AspNetCore.App 10.0.x`。ASP.NET Core 是 Thermo SDK 的传递依赖，本程序是命令行工具，不启动 Web 服务。

首次构建需要访问 GitHub 的 `raw.githubusercontent.com` 以及 `api.nuget.org`。**处理 RAW 时无需联网。** 系统支持情况参见 [Microsoft 的 Windows 安装说明](https://learn.microsoft.com/en-us/dotnet/core/install/windows)。

## 2. 克隆、构建并自检

```powershell
git clone https://github.com/ZexuanChenCN/DiaRawToMgf-alphapept-like.git
cd DiaRawToMgf-alphapept-like
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\build.ps1
```

`build.ps1` 会依次：

1. 从 Thermo 官方仓库的固定历史提交下载 **RawFileReader / Data 8.0.6**，并验证 SHA-256。
2. 使用仓库中的 `NuGet.Config` 和 `packages.lock.json` 恢复锁定的依赖。
3. 发布 Release 版本到 `artifacts\publish`。
4. 运行 `--self-test`，失败时立即停止并报错。

这里的 Thermo NuGet 包来自官方仓库，不要求它们出现在公共 NuGet 源。不要直接替换为官方仓库当前最新版，否则不再是本次结果验证时的 SDK 版本。

也可以分步执行：

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\prepare-sdk.ps1
dotnet restore .\DiaRawToMgf.csproj --configfile .\NuGet.Config --locked-mode
dotnet publish .\DiaRawToMgf.csproj -c Release --no-restore --no-self-contained -o .\artifacts\publish
.\artifacts\publish\DiaRawToMgf.exe --self-test
```

## 3. 批量处理一个文件夹

### 构建时遇到 CET 错误

如果发布步骤显示 `Your Windows doesn't fully support CET`，说明 SDK 编译器的 EXE 启动遇到 Windows CET 兼容问题；此时发布未完成，因此 `artifacts\publish\DiaRawToMgf.exe` 可能不存在。

先更新仓库，再选择兼容构建：

```powershell
git pull --ff-only
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\build.ps1 -CetCompatibility
```

该选项使用同一 SDK 的 `csc.dll`，通过已能启动的 `dotnet.exe` 执行，关闭共享编译器，并为本次生成的应用设置 `CETCompat=false`。这只取消本程序输出 EXE 的 CET 标记，不修改已安装的 SDK、系统注册表或其他程序的保护设置；feature 算法、.NET 10 目标框架及数值参数不变。取消 CET 标记意味着该输出不启用这项硬件影子栈保护；正常支持 CET 的机器继续用默认构建即可。[Microsoft CET 兼容性说明](https://learn.microsoft.com/en-us/dotnet/core/compatibility/interop/9.0/cet-support)。

兼容构建也会运行自检。完成后按下面的普通命令启动 EXE 即可。若 `dotnet --version` 本身已经无法启动，应先解决 Windows/.NET 安装问题，而不是反复运行构建脚本。

### 开始转换

将下面的输入、输出路径替换为服务器上的实际路径：

```powershell
.\artifacts\publish\DiaRawToMgf.exe --input "D:\DIA\raw" --output "D:\DIA\mgf" --jobs 1 --threads 4 --resume
```

默认仅处理输入目录直接包含的 `.raw` 文件；包含子目录时加 `--recursive`：

```powershell
.\artifacts\publish\DiaRawToMgf.exe --input "D:\DIA\raw" --output "D:\DIA\mgf" --recursive --jobs 2 --threads 4 --resume
```

输入也可以是单个文件：

```powershell
.\artifacts\publish\DiaRawToMgf.exe --input "D:\DIA\raw\sample.raw" --output "D:\DIA\mgf" --max-candidates 5 --top-peaks 250
```

等价的 DLL 启动方式：

```powershell
dotnet .\artifacts\publish\DiaRawToMgf.dll --input "D:\DIA\raw" --output "D:\DIA\mgf" --jobs 1 --threads 4 --resume
```

程序为每个 RAW 建立独立子目录，递归处理时保留相对目录结构。PowerShell 中可通过 `$LASTEXITCODE` 查看退出码；全部成功返回 `0`，存在失败文件返回非零。

## 4. 并行和内存

| 参数 | 含义 | 建议起点 |
|---|---|---:|
| `--jobs` | 同时处理的 RAW 数量，每个文件一个独立进程 | 1 |
| `--threads` | 每个 RAW 的 feature 统计、同位素提取等计算线程上限 | 4 |

完整测试 RAW 约 8.1 GB，单进程峰值工作集约 **13.5 GiB**。相近数据建议先按每个并发文件约 **20 GiB RAM** 规划，并给系统及其他任务留余量。例如两文件并发先预留约 40 GiB。线程数不是所有阶段的固定 CPU 使用量，也不是进程数。

优先使用较快的本地磁盘，避免所有进程同时争用较慢的共享盘。建议先用一个文件测量服务器速度和内存，再增加 `--jobs`。输入 RAW 只读；运行日志、feature 表、曲线和 MGF 均写入指定输出目录，不会主动在 C 盘暂存 RAW 或这些中间数据。

## 5. 命令行参数

| 参数 | 默认值 | 说明 |
|---|---:|---|
| `--input PATH` | 必填 | RAW 文件夹或单个 Thermo RAW |
| `--output PATH` | 必填 | 输出根目录 |
| `--recursive` | 关闭 | 递归处理子目录 |
| `--jobs N` | 1 | 并行 RAW 进程数 |
| `--threads N` | 4 | 单文件计算线程上限 |
| `--max-candidates N` | 5 | 每张 MS2 按局部 MS1 强度保留的前体数；`0` 表示不限 |
| `--top-peaks N` | 250 | MS2 按强度选前 N 根峰，然后按 m/z 升序写出 |
| `--centroid-ppm X` | 8 | 第一轮 MS1 hill 连接容差；第二轮使用自适应容差 |
| `--max-gap N` | 2 | hill 连接允许缺失的 MS1 扫描数 |
| `--min-hill-points N` | 3 | 初始 hill 的最少点数 |
| `--hill-split-level X` | 1.3 | 色谱 hill 峰谷拆分阈值 |
| `--isotope-split-level X` | 1.3 | 同位素丰度轮廓的峰谷拆分阈值 |
| `--min-charge N` | 1 | 最低前体电荷 |
| `--max-charge N` | 6 | 最高前体电荷 |
| `--isotope-correlation X` | 0.6 | 同位素 hill 的色谱余弦相关性阈值 |
| `--envelope-score X` | 0.6 | CHNOS averagine 余弦相似度阈值 |
| `--rt-start MIN` | 不限制 | 起始保留时间，单位分钟 |
| `--rt-end MIN` | 不限制 | 结束保留时间，单位分钟 |
| `--resume` | 关闭 | 跳过输入、参数及完成记录匹配的文件 |
| `--self-test` | — | 运行合成数据和算法边界检查 |
| `--help` | — | 显示帮助 |

完整结果对齐测试使用 `--max-candidates 5 --top-peaks 250`，其余使用默认参数。先保留该配置；增大候选数会为同一张 MS2 增加前体假设，并不等于增加碎片解混。

bootstrap 固定使用 seed 42、150 次迭代和至多 300 个样本；同位素 seed 上限为 100。旧近似版本的 `--isotope-ppm`、`--min-isotopes`、`--max-isotopes`、`--min-ms1-intensity` 已退役，非默认值会报错，不应作为当前算法的调参选项。

## 6. 输出文件与 scan 格式

```text
mgf/
  batch-summary.json
  sample/
    sample.mgf
    ms1_features.tsv
    precursor_curves.tsv
    mapping.tsv
    unmatched_ms2.tsv
    summary.json
    options.json
    worker-options.json
    run.log
    worker.log
```

| 文件 | 内容 |
|---|---|
| `sample.mgf` | 最终通过格式检查的 MGF |
| `ms1_features.tsv` | feature ID 和 AlphaPept 的 14 列汇总值 |
| `precursor_curves.tsv` | 每个 feature 在各原始 MS1 scan 上的同位素强度和 |
| `mapping.tsv` | MGF 记录与原始 MS2、feature、候选排名、局部强度的对应关系 |
| `unmatched_ms2.tsv` | 没有可导出候选的原始 MS2 及原因 |
| `summary.json` | 成功标记、输入指纹、参数相关指纹、数量、耗时、峰值内存 |
| `run.log` / `worker.log` | 算法阶段及进程日志 |
| `batch-summary.json` | 各输入文件的输出路径和进程退出码 |

MGF 格式示例：

```text
TITLE=sample.77079.77079.2.0.dta
CHARGE=2+
RTINSECONDS=1234.567890
PEPMASS=500.12345678
SCANS=77079
```

`TITLE` 的两个 scan 字段以及 `SCANS` 始终为**原始 RAW MS2 scan**。最后的候选排名从 `0` 开始，仅用于区分同一张 MS2 的不同前体，不替换 scan。MGF 记录 ID 也不会替换 scan。

## 7. 前体匹配逻辑与限制

1. 从完整 MS1 数据提取色谱 hills，并组装同位素 feature。
2. 对每张原始 MS2，取其紧邻的前后两张 MS1 中**实际观测到**的 feature 的并集，按 feature ID 去重。
3. 要求 MS2 的 RT 落在 feature 起止范围内，feature 的单同位素 m/z 落在真实隔离窗口内。
4. 对两侧 MS1 的同位素强度和做插值；某侧不存在该 feature 时按零处理。
5. 按插值强度降序排名，强度相同时按 feature ID 排序，默认取前 5 个候选。
6. 每个候选使用相同原始 MS2 的 top250 碎片。

**不做 MS2 碎片色谱解混。** 这些 MGF 是带有 MS1 前体指派的 DIA 混合谱。没有 feature 的 MS2 不伪造母离子或电荷，其 scan 保存在 `unmatched_ms2.tsv`。输入必须是 Thermo `.raw` 文件；不支持 Waters `.raw` 文件夹。没有可用质心数据/厂商 centroid stream 的 profile scan 会明确报错。

## 8. 完整 RAW 验证结果

测试文件：`Astral_IPP-IP_HumanBCell-Lulinrong-20260723_DIA.raw`；采集约 67.64 分钟。C# 从 RAW 独立运行，未加载 Python feature 缓存；Python 参考文件在 C# 运行前已生成。

| 检查项目 | 结果 |
|---|---|
| MS1 / MS2 数量 | 7,859 / 157,168 |
| feature 数量 | 223,027；ID 和全部 14 列数值与参考表一致 |
| 输出 MGF 候选数 | 215,483；全部候选对、电荷、排名和 TITLE 一致 |
| 匹配 / 未匹配原始 MS2 | 96,183 / 60,985 |
| MS2 top250 峰数组 | 全部一致 |
| MS1 曲线点 | 5,660,503 个点全部对应；强度仅有浮点尾差 |
| 最终 MGF 文件大小 | 两份均为 1,129,259,943 字节 |
| 最终 MGF 字节比较 | 完全相同 |
| C# RAW 读取 | 24.07 秒 |
| C# feature finding | 92.71 秒 |
| C# 导出及格式检查 | 130.94 秒 |
| C# 总耗时 | 247.72 秒，即约 4.13 分钟 |
| 单进程峰值工作集 | 14,500,270,080 字节，约 13.5 GiB |

两份完整 MGF 的 SHA-256 均为：

```text
72693bf835c22e0446211b5fe441c6f641b9e3622d021bdb71c918fc29af45cb
```

映射表的前体 m/z 最大差约 `5.7e-14`、插值强度最大差约 `2.5e-6`；MS1 曲线强度最大差约 `1.2e-7`。这些文本/浮点尾差未改变最终 MGF。

详细记录：[feature 表比较](validation/alphapept_port_full_features.json)、[MGF 与曲线比较](validation/alphapept_port_full_mgf.json)、[完整运行统计](validation/alphapept_port_full_summary.json)、[最终文件哈希和来源记录](validation/mgf_byte_identity.json)。耗时是一次 Windows 实测，不是服务器吞吐保证；MGF 对齐也不代表已验证鉴定灵敏度或 FDR。

## 9. 更新、续跑与故障排查

更新源代码后重新构建：

```powershell
git pull --ff-only
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\build.ps1
```

- `--resume` 只跳过输入、参数和成功指纹匹配的文件，并检查 MGF 大小及映射表是否存在。更改输入、算法版本或参数后使用新的输出根目录。
- MGF 先写入 `.partial`，格式检查通过后才成为 `.mgf`。`summary.json` 表示该文件完成；失败日志不能当成成功输出。
- RAW 进程失败时先看该文件的 `worker.log`、`run.log`，再看根目录的 `batch-summary.json`。
- 提示缺少运行时：安装 .NET 10 SDK，或同时安装 .NET 10 和 ASP.NET Core 10 运行时。
- Thermo 包恢复失败：重新运行 `prepare-sdk.ps1`，检查官方 GitHub 下载地址和网络代理；使用这里的 `NuGet.Config`，不要跳过 SDK 准备步骤。
- 仅需要部署已构建程序时，复制**整个** `artifacts\publish` 文件夹到目标 Windows 服务器；不要只复制 `.exe` 或 `.dll`。目标机仍需对应运行时。
- 输入路径可含空格，但需要用引号包裹。运行账号必须能读取 RAW 并写入输出目录；计划任务账号应确认可访问相应磁盘/共享目录。

## 10. 可选开发验证

`qa/` 中的 Python 脚本只用于重新生成 AlphaPept 参考快照和比较结果，正常使用不需要安装这些依赖。参考环境固定为 AlphaPept 0.5.3、NumPy 1.26.4、Numba 0.63.1；该移植保留 NumPy Windows 排序行为。

不安装 Python，也可以比较两份 MGF 的大小及 SHA-256：

```powershell
.\scripts\compare-mgf.ps1 -PythonMgf "D:\reference\sample.mgf" -CSharpMgf "D:\output\sample.mgf"
```

中间阶段核验入口为 `--verify-alphapept-stages ORACLE_DIRECTORY`。详见 [QA 说明](qa/README.md) 和 [移植说明](docs/IMPLEMENTATION.md)。`--feature-cache DIR` 仅用于单 RAW 开发对照，不是常规批量运行参数。

## 许可和第三方来源

源码按 [Apache-2.0](LICENSE) 发布。AlphaPept、Numba 和 NumPy 相关移植及其原始许可证见 [THIRD_PARTY_NOTICES.txt](THIRD_PARTY_NOTICES.txt) 与 [licenses/](licenses/)。

Thermo RawFileReader/Data 是独立厂商依赖，仍适用 [Thermo 官方 SDK 许可](https://github.com/thermofisherlsms/RawFileReader/blob/345aafcb46c944b843ca49880b5644bff6285f61/License.doc)。本仓库不提交 RAW、MGF 数据或 Thermo 二进制包；准备脚本从官方固定版本下载并校验。此项目不是 AlphaPept 或 Thermo 的官方发行版。
