# GWLD 2.0

用于全基因组连锁不平衡分析，支持 RMI、MI、r²、D 和 D'。
C++ 命令行版与 R 包共用 `GWLD` 命名空间下的计算内核。
VCF、`.vcf.gz` 和 BCF 的读取也共用 `GWLD-R/inst/include/gwld/vcf.hpp`。

## C++ 命令行版

依赖：CMake ≥ 3.14、C++17 编译器、Armadillo 头文件、OpenMP。
VCF 读取需要 HTSlib 及其压缩库依赖，详见 C++ 使用说明。
ArgumentsParser 使用你自己的仓库，以 Git 子模块放在
`GWLD-C++/deps/ArgumentsParser`。CLI 不再链接 Armadillo、BLAS 或 LAPACK 动态库。

```bash
git submodule update --init --recursive
cmake -S GWLD-C++ -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure

# calculate defaults to a full RMI matrix; --method accepts rmi / mi / r2 / d / dprime
./build/GWLD calculate --bed-prefix example/example --method mi --threads 4 -o result

# Calculate within-chromosome pairs up to 300 kb apart without allocating a full matrix
./build/GWLD calculate --vcf example/example-3.vcf --decay --max-distance-kb 300 -t 4 -o decay

# Add --heatmap to also write the full matrix; decay/circos reuse its values
./build/GWLD calculate --ped-prefix example/example-2 --decay --heatmap -o result

# Use separate subcommands for format conversion and frequency summaries
./build/GWLD convert --bed-prefix example/example --format vcf -o converted
./build/GWLD convert --vcf example/example-3.vcf --format dosage -o genotypes
./build/GWLD summarize --vcf example/example-3.vcf -o frequencies
```

CLI 统一使用 `GWLD <子命令> [选项]`；长选项使用小写和连字符，线程为
`-t/--threads`，输出前缀为 `-o/--output-prefix`。旧的扁平调用方式和参数已移除。
数值基因型输出为 `.dosage.tsv`（缺失值为 `NA`），频率输出为 `.freq.tsv`。

非标准安装位置可通过 `-DGWLD_ARMADILLO_INCLUDE_DIR=/path/to/include` 指定。
没有 OpenMP 时可用 `-DGWLD_OPENMP=OFF` 构建，运行时使用单线程。
完整参数及迁移说明见 [C++ 2.0 使用说明](GWLD-C++/README.md)。

### macOS（Apple Silicon）

使用 Homebrew 的原生依赖和 LLVM 编译器：

```bash
brew install cmake pkg-config armadillo htslib llvm libomp
git submodule update --init --recursive
cmake -S GWLD-C++ -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER="$(brew --prefix llvm)/bin/clang++" \
  -DGWLD_USE_BUNDLED_HTSLIB=OFF \
  -DCMAKE_INSTALL_PREFIX="$PWD/build/install"
cmake --build build -j
ctest --test-dir build --output-on-failure
cmake --install build
./build/install/bin/GWLD --version
```

macOS 构建会启用 `_REENTRANT`，以使用系统提供的线程安全 `lgamma_r`。

## R 包

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Rhtslib")
devtools::install_github("Rong-Zh/GWLD/GWLD-R")
library(GWLD)
# Rows are samples and columns are variants, coded as 0 / 1 / 2 / NA
g <- matrix(c(0L, 1L, 2L, 0L, 1L, 2L), nrow = 3)
GWLD(g, method = "MI", cores = 2)
```

R 保留现有函数接口，并优化矩阵计算、decay 和 circos。
`vcf2matrix("input.vcf.gz", threads = 4)` 直接返回 `Info` 和 `0/1/2/NA` 整数型
`Genotype` 矩阵，保留原始样本名；这替代了旧版包含 GT 文本的矩阵。
需要字符表时继续使用 `read.vcf()`，两者走同一个 HTSlib 读取器。
R 包通过 Rhtslib 构建，不链接 CLI 的 Linux 静态库。
回归测试：安装包后运行 `Rscript GWLD-R/tests/core.R` 和 `Rscript GWLD-R/tests/vcf.R`。

`ped2matrix()` 与 VCF/BED 读取器一样，默认直接返回 `0/1/2/NA` 整数基因型，
无需再调用 `codegeno()`；需要原始等位基因字符串时指定 `genotype = "char"`。
`codegeno()` 默认识别 `/` 和 `|`，对已编码数值矩阵保持等位基因方向。
`HeatMap(input)` 可直接使用读取器返回的 `Genotype` 和 `Info`；已有 LD 矩阵时用
`HeatMap(input, LDmatrix = r2)` 直接绘图，不重复计算。
circos 绘图可直接使用 `Info` 和 `circos()` / `calc_circos()` 的字符配对结果，
无需手动修改列类型。完整流程见 [example/example.R](example/example.R)，
安装后可从仓库根目录执行 `source("example/example.R", chdir = TRUE)`。
新增流程回归测试：`Rscript GWLD-R/tests/workflow.R`。
测试数据位于 `GWLD-R/tests/fixtures`，打包后也可独立执行全部 R 测试：
`R CMD build GWLD-R`，然后 `R CMD check --no-manual GWLD_2.0.0.tar.gz`。

圈图使用内置 `color_legend()` 组件，不依赖 ComplexHeatmap；支持横向/纵向连续色标、
标题和密集标签自动避让。默认将色标放在右上角，并按字号与画布大小调整高度。
推荐用等宽高画布导出圈图；长方形画布因保持圆形比例仍会留白。
图例回归测试：`Rscript GWLD-R/tests/circos_legend.R`。

R 接口在转换成整数前会拒绝小数、无穷值和非法基因型；整数型和数值型的
`0/1/2/NA`（以及计算接口中的 `-1` 缺失值）均可使用。
`GWLD()` 与 `LD()`、`MI()`、`RMI()` 一样保留位点名称，并统一接受 `r2`、
`Dprime` 和小写方法名。缓存热图的矩阵名称、顺序必须与输入位点一致；
无名称的矩阵按明确传入的位点顺序绘制。`color` 也可直接传入颜色向量。

大数据的 R 配对计算可直接写文件，避免将全部配对放入 R 矩阵：

```r
calc_decay(geno012, Info, method = "RMI", cores = 4, file = "decay.tsv")
calc_circos(geno012, Info, method = "RMI", cores = 4, file = "circos.tsv")
```

`file` 输出有标题行，按批次写入；不传 `file` 时仍返回原来的字符矩阵。
配对数值保留双精度，矩阵和配对计算都支持批次间中断；中断后的文件保留已写入行。
R 包可用 `GWLD_OPENMP=OFF R CMD INSTALL GWLD-R` 构建单线程版本。
回归测试还包括 `Rscript GWLD-R/tests/regressions.R` 和
`python3 GWLD-R/tests/interrupt_test.py`；测试单线程构建时设置 `GWLD_TEST_OPENMP=OFF`。

## 性能与验证

LD 求解器现在检查可行区间内的全部驻点和边界；可证明似然凹性的普通情形使用
单调导数搜索，可能存在多峰的情形使用多项式实根隔离。
这修复了部分双杂合子占比高的数据中选错似然峰、甚至选错 D 符号的问题。
LD 数值可能与旧求解器不同；下面的计时表保留为更新前的历史基准。

论文 Box 1 原文、原始补充材料及适配 2.0 的可运行验证脚本见
[2023 论文验证说明](tests/paper/README.md)。安装 `genetics`、`infotheo` 后运行
`Rscript tests/paper/paper_2023.R`，可检查 MI、RMI 和 LD 数值及并行一致性，
并生成热图与结果表。

本机 WSL、i5-12400、GCC 13.3 Release、单线程、512 样本 × 384 位点、
25% 缺失率、固定随机种子 2026，三轮交替取中位数：

| 方法 | 2.0 标量计数路径 | 位集内核 | 加速 |
| --- | ---: | ---: | ---: |
| MI | 0.130 s | 0.0080 s | 16.3 倍 |
| RMI | 0.155 s | 0.0097 s | 16.0 倍 |
| r² | 0.152 s | 0.0286 s | 5.3 倍 |
| D | 0.155 s | 0.0290 s | 5.4 倍 |
| D' | 0.153 s | 0.0288 s | 5.3 倍 |

标量参考保留提交 `b96745b` 使用的逐样本计数路径。上述测试的最大数值差异为 0。
计时包含预处理和矩阵创建，不包含文件读写；实际加速取决于数据规模、缺失率、
CPU 和线程数。4096 样本 × 128 位点、1% 缺失、4 线程的另一组测试加速约 3.7–5.4 倍。

```bash
cmake -S GWLD-C++ -B build -DGWLD_BUILD_BENCHMARKS=ON
cmake --build build -j
./build/gwld_engine_benchmark 512 384 1 25
./build/gwld_engine_benchmark 4096 128 4 1
```

至少 128 样本、4 位点时，新内核将每个位点编码为三条位集，通过 AND 和位计数
ARM64 使用 NEON；GCC/Clang 的 x86 构建按 CPU/系统能力和样本规模选择
AVX-512 VPOPCNTDQ、AVX2、POPCNT 或通用计数路径，不要求全局启用 AVX。
MSVC 保留 POPCNT/通用回退。
位集额外占每样本每位点 3 bit，按 64 样本补齐。小数据使用标量路径。
RMI 的整数阶乘对数缓存按计算量选择启用，上限 8 MiB，使用原数学函数预计算。
支持的平台使用可重入 lgamma_r 进行并行 RMI 计算。
decay 按染色体建立排序索引，
仅计算距离窗口内的配对，支持输入坐标无序以及超过 32 位的坐标。
R 的 OpenMP 工作线程只访问 C++ 内存，结果按原始位点顺序输出。
纯格式转换及频率统计按需执行；按染色体输出逐组建立计算缓存。

字符基因型转换按列缓存不同的 R 字符串，直接访问列内存，不再复制列或逐单元格
构造字符串；多字符等位基因按完整 token 匹配。200 万个基因型的本机三轮交替
基准中，`recode`、`stats_freq`、矩阵转换约快 10 倍；普通有效输入与旧实现逐项核对。
复现：`Rscript GWLD-R/tests/benchmarks/codegeno.R`。边界回归测试见
`Rscript GWLD-R/tests/codegeno.R`。


字符编码的列内缓存现使用可复用的开放寻址表，减少重复基因型的哈希节点分配。
2026-11-07，在本机 1000 万个字符基因型、预热一次后六轮交替计时的中位数中，
相对上轮 2.0 的 `unordered_map` 实现，`recode`、`stats_freq`、矩阵转换分别再快
约 1.25、1.31、1.28 倍，输出逐项完全一致。
`tests/legacy` 中的编码和计算文件仅作为旧版比较基准；生产实现分别为
`GWLD-R/src/codeGenotype.cpp` 和 `GWLD-R/inst/include/gwld/`。
