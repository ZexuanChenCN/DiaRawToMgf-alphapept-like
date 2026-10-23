# 开发对照工具（正常运行无需 Python）

- `alphapept_stage_oracle.py`：从 AlphaPept MS1 HDF 提取指定 RT 片段，使用安装的 AlphaPept 0.5.3 导出各阶段二进制快照。
- `compare_feature_tables.py`：比较 feature ID 和所有共享数值列。
- `compare_csharp_dia.py`：比较候选对、TITLE、电荷、排名、完整 MS2 峰数组及 MS1 曲线；`--mode alphapept-port` 在完整核验失败时返回非零。
- `clean_pfind_mgf.py`：比较脚本所使用的 MGF 解析器，同时保留早期格式清洗命令。

oracle 环境需 AlphaPept 0.5.3、NumPy 1.26.4、Numba 0.63.1、h5py、pandas 及相应依赖。这些脚本不被 C# 程序调用。

各脚本使用 `python .\qa\脚本名.py --help` 查看参数。oracle 生成后运行：

```powershell
.\artifacts\publish\DiaRawToMgf.exe --verify-alphapept-stages "D:\validation\oracle"
python .\qa\compare_feature_tables.py "D:\reference\ms1_features.tsv" "D:\output\sample\ms1_features.tsv" "D:\validation\features.json"
```

完整 MGF 对照中，`--reference` 指向 Python `ms2_centric_top5` 输出目录，`--candidate` 指向单 RAW 的 C# 输出目录。`--reference-features` 接受原始 feature TSV，`--reference-curves` 接受 Python 生成的 `precursor_ms1_curves.h5`。

仓库 `validation/` 仅包含小型结果报告，不包含实验 RAW、MGF、HDF 或大体积 oracle。原始测试摘要中的绝对本机输入路径已归一为 RAW 文件名；数值统计未修改。
