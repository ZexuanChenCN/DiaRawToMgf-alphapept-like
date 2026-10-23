# 移植说明

兼容目标是 AlphaPept 0.5.3、Numba 0.63.1 和 NumPy 1.26.4 的 Windows 数值/排序行为。MS1 feature finding 与本项目的 DIA MS2 匹配/导出是两个步骤；后者采用 README 中说明的相邻 MS1、局部强度 top5 策略。

| 源文件 | 对应计算 |
|---|---|
| `AlphaPeptHills.cs` | 两轮自适应 centroid 连接、hill 去重、拆分和过滤 |
| `AlphaPeptStats.cs` | seed 42 bootstrap、面积、峰值、RT 排序 |
| `AlphaPeptIsotopes.cs` | 连通分量、seed/电荷遍历、同位素簇选择 |
| `AlphaPeptChem.cs` | CHNOS averagine、卷积及相似度 |
| `AlphaPeptNumerics.cs` | Numba/NumPy 排序、MT19937、平滑和数值辅助函数 |
| `AlphaPeptReport.cs` | 100 点插值、feature 的 14 列统计及 centroid 归属 |
| `Pipeline.cs` | MS2 匹配、MGF/TSV 输出和格式检查 |

bootstrap 为每条 hill 重新使用 seed 42，迭代 150 次，每次至多 300 个样本。同位素电荷默认 1–6，seed 上限 100；间距使用原始 `1.00286864 / charge`、`0.0109135 / charge` 和五倍 bootstrap 标准差。

为复现固定版本结果，保留了上游在 hill 去重指针尾部和 large-hill filter 索引上的边界行为，不在兼容模式中顺手修正。Numba 与 NumPy 对并列元素的排序不同，分别实现对应算法。

优化只复用不变的纯计算：理论分布、相关性、可证明唯一的排序。存在影响结果的并列键时仍执行原排序。稀疏 grow 遍历必须保留原扫描序列中第一个超出质量范围的元素，包括该元素仍可能匹配并扩大范围的情况；合成测试覆盖该边界。

完整 Windows RAW 测试的 feature 表 14 列、最终 MGF 已与独立 Python 结果对齐。跨其他 AlphaPept 版本或不同数据，不据此声称普遍逐字节一致；修改数值或顺序相关代码后应重新进行阶段和完整文件对照。
