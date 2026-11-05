# GWLD 2.0

用于全基因组连锁不平衡分析，支持 RMI、MI、r²、D 和 D'。
C++ 命令行版与 R 包共用 `GWLD` 命名空间下的计算内核。
VCF、`.vcf.gz` 和 BCF 的读取也共用 `GWLD-R/inst/include/gwld/vcf.hpp`。

## C++ 命令行版

依赖：CMake ≥ 3.14、C++11 编译器、Armadillo 头文件、OpenMP。
VCF 读取需要 HTSlib 及其压缩库依赖，详见 C++ 使用说明。
ArgumentsParser 使用你自己的仓库，以 Git 子模块放在
`GWLD-C++/deps/ArgumentsParser`。CLI 不再链接 Armadillo、BLAS 或 LAPACK 动态库。

```bash
git submodule update --init --recursive
cmake -S GWLD-C++ -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure

# 默认计算 RMI 全矩阵；--method 支持 RMI / MI / r2 / r^2 / D / D'
./build/GWLD --bfile example/example --method MI --threads 4 -o result

# 只计算同染色体、300 kb 内的配对，不分配全矩阵
./build/GWLD --vcf example/example-3.vcf --decay --max-dist 300 -p 4 -o decay

# 加 --heatmap 同时输出全矩阵，decay/circos 复用已经计算的结果
./build/GWLD --file example/example-2 --decay --heatmap -o result
```

非标准安装位置可通过 `-DGWLD_ARMADILLO_INCLUDE_DIR=/path/to/include` 指定。
没有 OpenMP 时可用 `-DGWLD_OPENMP=OFF` 构建，运行时使用单线程。
完整参数及迁移说明见 [C++ 2.0 使用说明](GWLD-C++/README.md)。

## R 包

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Rhtslib")
devtools::install_github("Rong-Zh/GWLD/GWLD-R")
library(GWLD)
# 行为样本、列为位点，编码为 0 / 1 / 2 / NA
g <- matrix(c(0L, 1L, 2L, 0L, 1L, 2L), nrow = 3)
GWLD(g, method = "MI", cores = 2)
```

R 保留现有函数接口，并优化矩阵计算、decay 和 circos。
`vcf2matrix("input.vcf.gz", threads = 4)` 直接返回 `Info` 和 `0/1/2/NA` 整数型
`Genotype` 矩阵，保留原始样本名；这替代了旧版包含 GT 文本的矩阵。
需要字符表时继续使用 `read.vcf()`，两者走同一个 HTSlib 读取器。
R 包通过 Rhtslib 构建，不链接 CLI 的 Linux 静态库。
回归测试：安装包后运行 `Rscript GWLD-R/tests/core.R` 和 `Rscript GWLD-R/tests/vcf.R`。

## 性能与验证

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
一次处理 64 个样本；x86 运行时检测 POPCNT，其他机器使用对应的兼容路径。
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

