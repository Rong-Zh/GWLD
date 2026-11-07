# GWLD C++ 2.0

构建命令见[项目 README](../README.md)。默认 Release 构建，程序输出在构建目录中。
公共计算与 VCF 读取代码位于 `../GWLD-R/inst/include/gwld`，命名空间为 `GWLD`。
`vcf.hpp` 是 CLI 和 R 唯一的 VCF 读取实现，使用 HTSlib 流式读取，不使用内存映射。

## 子命令与参数

统一使用 `GWLD <command> [options]`。任务使用子命令，任务内部的输出和配置使用
小写连字符选项。所有方法名也统一为小写：`rmi`、`mi`、`r2`、`d`、`dprime`。

| 子命令 | 功能 |
| --- | --- |
| `calculate` | LD 计算，可组合 heatmap、decay、circos 和按染色体输出 |
| `convert` | 格式转换，`--format vcf`（默认）或 `--format dosage` |
| `summarize` | 输出 REF/ALT 等位基因频率和基因型缺失率 |

`GWLD --help` 查看子命令，`GWLD calculate --help` 或 `GWLD help calculate` 查看具体选项。
`GWLD -V` / `GWLD --version` 显示版本。

### 公共选项

三个输入选项必须且只能指定一个：

| 参数 | 输入 |
| --- | --- |
| `-b, --bed-prefix PREFIX` | PLINK SNP-major `.bed`、`.bim`、`.fam` |
| `-p, --ped-prefix PREFIX` | PLINK `.ped`、`.map` |
| `-v, --vcf PATH` | 双等位、二倍体 VCF、`.vcf.gz`（gzip/BGZF）或 BCF；支持任意 GT 位置 |

三个子命令都支持 `-t, --threads N`（默认 1，必须为正数；也用于 BGZF 解压）、
`-o, --output-prefix PREFIX`（默认 `gwld`）和 `-h, --help`。

### calculate 选项

| 参数 | 作用与默认值 |
| --- | --- |
| `-m, --method NAME` | `rmi`（默认）、`mi`、`r2`、`d`、`dprime` |
| `-d, --decay` | 输出同染色体内的配对到 `.decay` |
| `--max-distance-kb KB` | decay 最大距离，默认 300 kb，允许 0，边界包含在内；显式指定时必须搭配 `--decay` |
| `-c, --circos` | 输出跨染色体且值大于阈值的配对到 `.circos` |
| `--circos-threshold VALUE` | circos 阈值，默认 0.1；显式指定时必须搭配 `--circos` |
| `--by-chromosome` | 每条染色体输出 `PREFIX.chrNAME.heatmap` |
| `--heatmap` | 与其他输出同时生成全矩阵 |

### 输出格式

`convert --format vcf` 写入标准 `.vcf`，保留样本名和 REF/ALT 编码方向。
`convert --format dosage` 写入 `.dosage.tsv`：首行为标题，前五列为
`CHROM POS ID REF ALT`，之后每列是一个样本；值为 ALT 拷贝数 `0/1/2`，
缺失值为 `NA`。该表不是 VCF。

`summarize` 写入 `.freq.tsv`，首行为标题，八列为
`CHROM POS ID REF ALT REF_FREQ ALT_FREQ MISSING_RATE`。
频率按非缺失样本计算，全部缺失的位点频率为 `nan`。

转换和频率统计各自只生成对应文件，不计算 LD。

```bash
./build/GWLD calculate --vcf example/example-3.vcf --method r2 --heatmap -o ld
./build/GWLD calculate --bed-prefix example/example --decay --circos --threads 4 -o pairs
./build/GWLD convert --ped-prefix example/example-2 --format vcf -o converted
./build/GWLD convert --vcf example/example-3.vcf --format dosage -o genotypes
./build/GWLD summarize --bed-prefix example/example -o frequency
```

上面的路径以仓库根目录为工作目录。

## 2.0 行为变化

- `calculate` 只指定输入和方法时，默认输出全矩阵；指定 decay、circos 或按染色体输出时，
  仅生成对应输出。需要全矩阵时添加 `--heatmap`。全矩阵仍需要 O(位点数²) 内存。
- decay 使用绝对物理距离，只计算同染色体、窗口内的配对；不要求输入按坐标排序。
  输出为 `ID_1 ID_2 Dist Value` 四列，无标题行，制表符分隔。
- circos 为 `CHR_1 POS_1 ID_1 CHR_2 POS_2 ID_2 Value` 七列，无标题行。
  空配对结果会生成空文件。线程数不会改变结果顺序。
- heatmap 保留上三角布局，下三角和对角线为 `nan`，制表符分隔。
- 无完整配对的位点对返回缺失值。单态位点的 D 为 0，r² 和 D' 为缺失值。
  LD 的极大似然估计比较全部驻点和端点；对可证明凹性的情形使用快速单调导数搜索，
  对可能多峰的情形隔离驻点多项式的全部实根，避免选错似然峰。
  RMI 的类别数和 LD 的边际频率保持旧算法约定；缺失样本按配对排除。
- VCF 保持 REF/ALT 的编码方向；BED 保持旧版 allele-1 dosage，REF 为 BIM 第六列，
  ALT 为第五列；PED 以高频等位基因为 REF，频率相同时按字典序确定。
  D 的符号依赖等位基因方向，跨格式比较 D 时应先统一方向。
- PLINK 样本名与 R 端统一为 `FID:IID`，保留家系身份。分隔符 `:` 和 `%`
  在各组件内分别转义为 `%3A` 和 `%25`。重复的家系/个体组合会报错；VCF 样本名原样保留。
- 数值文本输出使用可往返恢复双精度的有效位数。
- 旧的扁平调用及 `--2vcf`、`--code012`、`--allele-freq`、`--bfile`、`--file`、
  `--by-chrom`、`--max-dist`、`--threshold`、`--out-prefix` 已移除，不保留别名。
  线程短选项现在是 `-t`，`-p` 为 PED 输入前缀。CLI 方法值使用小写，R 接口不变。
- 格式错误、缺失输入、非法方法和线程数使用非零退出码。
  基因型计算不支持多等位或非缺失的单倍体 GT。BCF 的坐标范围受其格式限制。

## 开发与测试

`deps/ArgumentsParser` 是指向 `git@github.com:Rong-Zh/ArgumentsParser.git` 的子模块。
本次从远端拉取的 HEAD 为 `16d621e92e9adc416be1d0befa64b366116ff4b5`。
Git 子模块固定提交以保证构建可复现；更新时显式拉取并记录新的子模块提交。

Linux x86-64 默认链接 `deps/htslib` 中的 1.24 静态库，需要 pkg-config、zlib、
bzip2、liblzma、libdeflate 开发包。其他平台安装原生 HTSlib ≥ 1.15 后使用
`-DGWLD_USE_BUNDLED_HTSLIB=OFF`，由 pkg-config 提供编译参数。
R 包的同一个读取器通过 Rhtslib 构建，不依赖这份 Linux 静态库。

`ctest --test-dir build --output-on-failure` 覆盖独立数值参考、缺失值、
窗口边界、格式转换、输入校验和并行一致性。可用 `GWLD_BUILD_BENCHMARKS=ON`
构建 `gwld_benchmark`，同编译条件下对比 `tests/legacy` 中保留的旧计算实现。
`gwld_vcf_benchmark input.vcf` 单独比较旧流解析器与公共 HTSlib 读取器，交替执行
三轮取中位数，并核对所有基因型、样本、等位基因和位点信息。
`gwld_engine_benchmark [样本数=512] [位点数=384] [线程数=1] [缺失百分比=25]`
对比提交 `b96745b` 使用的标量计数路径与位集内核，包含预处理，三轮交替取中位数，
核对 MI、RMI、r²、D、D' 的全部结果。位集内核和 RMI 缓存的条件、内存开销及
本机结果见项目 README。`packed` 回归测试还会强制测试兼容计数路径，并覆盖
64 位边界、缺失值、空矩阵以及多线程结果。
本机 WSL、GCC 13.3 Release、单线程、512 样本 × 12,000 位点、41,836,454 字节
GT:DP 合成 VCF 的一次三轮测试中，旧流解析器中位数为 1.424 s，HTSlib 读取器
为 0.291 s（约 4.9 倍）。该基准只计读取和矩阵生成，包含打开文件，不计 LD 计算；
结果依赖存储、缓存、HTSlib 版本和字段内容。
旧版中英文手册作为 1.x 资料保留，2.0 以本说明和 `--help` 为准。


### 计算内核基准

启用 `GWLD_BUILD_BENCHMARKS=ON` 后，运行：

```bash
./build/gwld_kernel_benchmark 512 768 0 1
./build/gwld_kernel_benchmark 512 768 5 1
./build/gwld_kernel_benchmark 8192 384 5 1
```

参数依次为样本数、位点数、缺失百分比和线程数。使用固定随机种子，预热一次后
取三轮中位数，计时包括位点汇总、位集构建、缓存和全矩阵计算，不包括文件读写。
可传入第五个参数作为输出前缀，保存每种方法的列主序双精度矩阵，供版本间逐项比较。

内核使用受区间约束的 Newton 迭代加速可证明凹性的 LD 似然求解；困难情形
退回二分，多峰情形仍比较全部驻点和端点。MI 使用最多 8 MiB 的计数对数缓存；
RMI 缓存无缺失位点的边际统计量，含缺失的配对仍执行完整公式。
无缺失配对只需四次位集交集即可还原 3×3 列联表；ARM64 使用 NEON 累积计数，
保留其他架构的原生及兼容路径。打包后的对角计数同时提供位点统计量，避免再次扫描整数矩阵。

2026-11-07，macOS arm64、LLVM 22 Release，第一轮优化前后的本地 2.0 实现对比：

| 样本 × 位点 | 缺失 | 线程 | MI 加速 | RMI 加速 | r² 加速 |
| --- | --- | --- | --- | --- | --- |
| 512 × 768 | 0% | 1 | 3.16× | 3.96× | 5.41× |
| 512 × 768 | 5% | 1 | 2.62× | 1.06× | 5.42× |
| 8192 × 384 | 5% | 1 | 1.78× | 1.50× | 3.34× |
| 512 × 1536 | 0% | 4 | 3.09× | 3.79× | 5.36× |

上述合成数据的最大逐项差异低于 `1.3e-14`，缺失位置一致；加速幅度依赖样本规模、缺失率和硬件，并非所有数据都能取得相同收益。


### SIMD 分派与验证

共享内核在编译时保留多个专用函数，运行时按 CPU 和操作系统支持情况选择：

- ARM64：NEON，含无缺失配对的四交集路径。
- GCC/Clang x86：至少 8 个 64 位字时优先选择 AVX512F + VPOPCNTDQ；
  否则至少 16 个字且支持 AVX2 时使用 AVX2 字节查表计数；小配对使用 POPCNT。
- 不支持上述指令的 CPU：通用 SWAR；MSVC 当前保留 POPCNT/SWAR。

不使用全局 `-mavx2`、`-mavx512*` 或 `-march=native`，因此启动和分派代码
保持兼容。AVX2 每 31 个向量归并字节累加值，NEON 每 4095 个向量归并，
避免大样本计数溢出。向量尾部使用标量计数。

`gwld_simd_test` 会强制执行当前 CPU 支持的每条路径，并逐项比较计数表、
完整样本数和无缺失的四交集还原结果，覆盖 64/128/256/512 位边界及累加器归并边界。
可直接运行 `./build/gwld_simd_test`，或由 CTest 自动运行。
2026-11-07：ARM64 的 NEON/通用路径、macOS x86_64 Clang 的 POPCNT 回退、
Debian x86_64 GCC 14 编译及 QEMU 10 的 AVX2/POPCNT/无 POPCNT 回退测试通过。
AVX-512 专用函数已由 GCC 与 Clang 编译；当前测试机器和 QEMU 均不能执行
AVX-512，因此该路径的运行结果和真实 Intel/AMD 的 SIMD 加速幅度尚待硬件验证。


### 后续热点优化

在上述实现上，进一步复用 MI 的完整位点边际量，为有缺失数据的 RMI 缓存
整数及半整数对数，并以缩放计数计算边际平方和；查表值使用原 `std::log`
预计算，额外缓存最多 8 MiB，只在配对数量足以摊销且确有部分缺失列时建立。
无缺失位集内核返回四个计数，减少聚合对象的传递和无用边际计算；矩阵输出
按列填充缺失部分与有效部分，避免先全矩阵填充再覆盖上三角。
凹性已确认的 LD 似然仅比较可能具有非零似然的端点，复用内部解的评分；
多峰情形仍隔离全部驻点并比较端点。

以下比较的基线是本节优化开始前的 2.0，不是 1.x。相同编译条件下，
每次预热后取三轮中位数，前后版本再交替执行三次后取中位数；
计时包含打包、统计缓存、分配和矩阵计算，不含文件读写。
`hwe` 数据的位点间次要等位基因频率从 0.01 到 0.5 变化。

| 样本 × 位点 | 缺失 | 分布 | 线程 | MI 加速 | RMI 加速 | r² 加速 |
| --- | --- | --- | --- | --- | --- | --- |
| 512 × 1536 | 0% | uniform | 1 | 1.32× | 1.08× | 2.12× |
| 512 × 1536 | 5% | uniform | 1 | 0.96× | 1.36× | 1.96× |
| 512 × 1536 | 5% | hwe | 1 | 1.04× | 1.36× | 1.65× |
| 8192 × 768 | 5% | uniform | 1 | 0.99× | 1.08× | 1.26× |
| 512 × 1536 | 5% | uniform | 4 | 1.02× | 1.36× | 1.93× |

最大逐项差异低于 `2.7e-15`，缺失位置一致。部分 MI 输入没有收益，表中同时保留
略慢的测量值。900 组类别数及缺失边界检查、10 万组似然压力检查、OpenMP
开启/关闭、ARM64 实测和 x86 AVX2 仿真核心检查均通过。真实 x86 性能仍需在
相应硬件上测量。

复现混合等位基因频率数据（第五个参数 `-` 表示不保存二进制矩阵）：

```bash
./build/gwld_kernel_benchmark 512 1536 5 1 - hwe
```


### 工具函数优化

`tests/legacy/utils.cpp` 的兼容接口已优化，历史实现另存为
`utils_reference.cpp/.h`，计算基准继续链接历史版本。新版本避免输入字符串、
分隔符和矩阵列的复制，直接拆分文件字段，批量写出数值，并复用共享的整数
分类计数。生产 PLINK 解析与 `summarize` 命令也使用对应的快速拆分/计数路径。

本机 macOS arm64、LLVM 22 Release，预热后前后交替三轮取中位数：

| 操作 | 旧实现 | 优化后 | 加速 |
| --- | --- | --- | --- |
| split | 27.16 ms | 8.39 ms | 3.24× |
| production_fields | 14.01 ms | 4.54 ms | 3.08× |
| readfile | 20.85 ms | 6.57 ms | 3.17× |
| stats_freq | 6.32 ms | 1.96 ms | 3.23× |
| counts_allele | 6.46 ms | 2.87 ms | 2.25× |
| writearma | 131.54 ms | 31.21 ms | 4.21× |

这些是函数级基准，不代表整个命令行的加速倍数；数据规模在
`tests/utils_benchmark.cpp` 中固定。频率与计数逐项相同，正常矩形数据的输出
文件逐字节相同。数值写出保留这个旧接口原有的 6 位有效数字；生产 CLI
仍使用 17 位精度。C++14/旧标准库的格式化回退及 GCC/Clang 兼容性由测试覆盖。

```bash
cmake -S GWLD-C++ -B build -DGWLD_BUILD_BENCHMARKS=ON
cmake --build build --target gwld_utils_benchmark
./build/gwld_utils_benchmark build/utils-perf
ctest --test-dir build --output-on-failure
```
