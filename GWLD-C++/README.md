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
  LD 的极大似然估计显式比较可行区间端点，修复完全连锁时的边界数值问题。
  RMI 的类别数和 LD 的边际频率保持旧算法约定；缺失样本按配对排除。
- VCF 保持 REF/ALT 的编码方向；BED 保持旧版 allele-1 dosage，REF 为 BIM 第六列，
  ALT 为第五列；PED 以高频等位基因为 REF，频率相同时按字典序确定。
  D 的符号依赖等位基因方向，跨格式比较 D 时应先统一方向。
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
