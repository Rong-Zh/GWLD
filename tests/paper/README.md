# 2023 论文演示与数值验证

`paper_2023.R` 按 Zhang et al. (2023) 的 Box 1 和方法部分编写，适配本地
GWLD 2.0。它使用包内 `duck$SNP`（500 个样本、染色体 26 的 60 个位点），
逐对验证 1,770 个位点对，并输出计算结果和图片。

## 来源与原文

- 论文：[GWLD: an R package for genome-wide linkage disequilibrium analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC10468308/)，
  DOI：[10.1093/g3journal/jkad154](https://doi.org/10.1093/g3journal/jkad154)。
- `reference/box1_original.txt`：从论文 XML 提取的 Box 1 原文，保留原始引号和
  输入路径占位符，供追溯，不能直接作为 R 脚本执行。
- `reference/jkad154_supplementary_data.zip`：论文原始补充材料，包含 C++ 使用手册
  `File_S1_G3-2023-404081.docx` 及结果表/图片
  `Supplemental_Tables_and_Figures_G3-2023-404081.docx`；其中没有独立的验证脚本。
  SHA-256：`988f8aee2d823d1f1a886ccd663354207373527d59d4e473ca9a758db11e3339`。
- 获取接口：[Europe PMC 正文 XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10468308/fullTextXML)、
  [补充材料](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10468308/supplementaryFiles)。
  获取日期：2026-11-05。论文原文依据 CC BY 4.0 保留并注明来源。

本目录的可运行脚本是新编写的复现版本，不是从补充材料提取出的原始验证代码。
它取消自动安装远程 GWLD 的步骤，并将 Box 1 的 `SnpPosition = Info$POS`
改为当前接口要求的染色体和坐标两列。

## 运行

先安装本地 `GWLD-R` 和论文中使用的交叉验证包 `genetics`、`infotheo`。
验证脚本会打印实际加载的 GWLD 版本与路径，不会自动安装或更新任何包。

```r
install.packages(c("genetics", "infotheo"))
```

在仓库根目录运行，默认使用 2 个线程，输出到 `build/paper-validation`：

```bash
Rscript tests/paper/paper_2023.R
```

如果使用本次 macOS 测试安装的项目本地 R 库：

```bash
R_LIBS="$PWD/build/R-library" Rscript tests/paper/paper_2023.R
```

可指定输出目录和线程数：

```bash
R_LIBS="$PWD/build/R-library" Rscript tests/paper/paper_2023.R build/paper-validation-4 4
```

## 验证内容

| 内容 | 参考或检查方式 | 绝对误差阈值 |
| --- | --- | --- |
| MI | `infotheo::mutinformation`，由自然对数单位换算为 bit | `1e-10` |
| RMI | 独立 R 表达式计算 GWLD 使用的 RMI 近似式 | `1e-10` |
| r²、D、D′ | `genetics::LD`，按主等位基因方向对齐 D | `1e-3` |
| 并行计算 | 五种方法的单线程与指定线程结果逐项完全一致 | 精确一致 |
| R 接口 | `GWLD` 与 `MI`、`RMI`、`LD` 返回的数值一致 | R 默认数值比较 |
| decay、circos | 直接计算与复用 RMI 矩阵两种路径一致 | 精确一致 |

LD 的 `1e-3` 阈值用于容纳 `genetics` 默认优化器的收敛精度与 GWLD 2.0
显式比较可行区间端点产生的差异；这不表示逐位完全相同。
MI 每对位点仅使用完整样本。RMI 保留 GWLD 的整个位点类别数约定。
参考包或断言失败会产生非零退出码，不会跳过交叉验证。

脚本还输出 MI、RMI、r² 的 Pearson 相关矩阵、三种热图，以及 Box 1 使用的
包内历史全基因组 circos 图。`circos-stored-duck.png` 来自保存的数据，
不是用染色体 26 的 60 个位点重新计算的全基因组结果；该子集的跨染色体结果为空。

## 输出与范围

- `checks.tsv`：参考方法、配对数、最大误差、阈值和通过状态。
- `pairs-1.tsv` 至 `pairs-5.tsv`：按 MI、RMI、r²、D、D′ 顺序保存每对参考值和误差。
- `correlations.tsv`、`decay.tsv`、`circos.tsv`：相关性及 RMI 配对结果。
- `results.rds`：基因型、元数据、五种矩阵和参考值。
- `heatmap-*.png`、`circos-stored-duck.png`：图片。
- `sessionInfo.txt`：R 环境和已加载包的版本信息。

论文完整鸭数据为 542 个样本、40,601 个 SNP；仓库中的包内数据是子集。
因此不要求精确重现补充表 S1 的相关系数，也不将此脚本视为全部论文结果的复现。
此脚本没有运行 PLINK、TASSEL 或人、猪、玉米全数据集的验证和性能测试。

2026-11-05 在 macOS arm64、R 4.6.1、本地 GWLD 2.0.0 上运行：五项数值检查、
并行和两种 decay/circos 路径均通过。MI 最大绝对误差约 `4.15e-15`，
RMI 约 `7.88e-15`；r²、D、D′ 分别约 `3.70e-4`、`7.64e-5`、`6.48e-4`。
