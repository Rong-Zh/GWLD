# GWLD C++ 2.0

构建命令见[项目 README](../README.md)。默认 Release 构建，程序输出在构建目录中。
公共计算与 VCF 读取代码位于 `../GWLD-R/inst/include/gwld`，命名空间为 `GWLD`。
`vcf.hpp` 是 CLI 和 R 唯一的 VCF 读取实现，使用 HTSlib 流式读取，不使用内存映射。

## 输入与参数

三个输入选项必须且只能指定一个：

| 参数 | 输入 |
| --- | --- |
| `-b, --bfile PREFIX` | PLINK SNP-major `.bed`、`.bim`、`.fam` |
| `-f, --file PREFIX` | PLINK `.ped`、`.map` |
| `-v, --vcf PATH` | 双等位、二倍体 VCF、`.vcf.gz`（gzip/BGZF）或 BCF；支持任意 GT 位置 |

| 参数 | 作用与默认值 |
| --- | --- |
| `-m, --method NAME` | `RMI`（默认）、`MI`、`r2` / `r^2`、`D`、`D'` |
| `-p, --threads N` | 线程数，默认 1，必须为正数；也用于 BGZF 解压 |
| `-o, --out-prefix PREFIX` | 输出前缀，默认 `gwld` |
| `-d, --decay` | 输出同染色体内的配对到 `.decay` |
| `--max-dist KB` | decay 最大距离，默认 300 kb，允许 0，边界包含在内 |
| `-c, --circos` | 输出跨染色体且值大于阈值的配对到 `.circos` |
| `--threshold VALUE` | circos 阈值，默认 0.1 |
| `--by-chrom` | 每条染色体输出 `PREFIX.chrNAME.heatmap` |
| `--heatmap` | 与其他输出同时生成全矩阵 |
| `--2vcf` | 输出带样本名和标准头部的 `.vcf` |
| `--code012` | 输出旧版 VCF 形状的数值基因型表，保留 `.vcf` 后缀；它不是标准 VCF |
| `--allele-freq` | `.freq` 三列：REF 频率、ALT 频率、基因型缺失率 |
| `-h, --help` | 参数帮助 |
| `-V, --version` | 版本号 |

## 2.0 行为变化

- 只指定输入和方法时，默认输出全矩阵；指定 decay、circos、按染色体或转换任务时，
  仅生成对应输出。需要全矩阵时添加 `--heatmap`。全矩阵仍需要 O(位点数²) 内存。
- decay 使用绝对物理距离，只计算同染色体、窗口内的配对；不要求输入按坐标排序。
  输出为 `ID_1 ID_2 Dist Value` 四列，无标题行，制表符分隔。
- circos 为 `CHR_1 POS_1 ID_1 CHR_2 POS_2 ID_2 Value` 七列，无标题行。
  空配对结果会生成空文件。线程数不会改变结果顺序。
- heatmap 保留上三角布局，下三角和对角线为 `nan`，制表符分隔。
- 无完整配对的位点对返回缺失值。单态位点的 D 为 0，r² 和 D' 为缺失值。
  LD 的极大似然估计显式比较可行区间端点，修复完全连锁时的边界数值问题。
  RMI 的类别数和 LD 的边际频率保持旧算法约定；缺失样本按配对排除。
- VCF 保持 REF/ALT 的编码方向；BED 保持旧版 allele-1 dosage，REF 为 BIM 第六列，
  ALT 为第五列；PED 以高频等位基因为 REF，频率相同时按字典序确定。
  D 的符号依赖等位基因方向，跨格式比较 D 时应先统一方向。
- 使用 `--vcf`、`--decay` 等标准长参数写法，或 `-v`、`-d` 等短参数。
  旧手册中的 `-vcf` 等单横线长参数写法不再支持。
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
本机 WSL、GCC 13.3 Release、单线程、512 样本 × 12,000 位点、41,836,454 字节
GT:DP 合成 VCF 的一次三轮测试中，旧流解析器中位数为 1.424 s，HTSlib 读取器
为 0.291 s（约 4.9 倍）。该基准只计读取和矩阵生成，包含打开文件，不计 LD 计算；
结果依赖存储、缓存、HTSlib 版本和字段内容。
旧版中英文手册作为 1.x 资料保留，2.0 以本说明和 `--help` 为准。
