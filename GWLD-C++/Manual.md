# GWLD manual

## 安装

1. 依赖

* armadillo <http://arma.sourceforge.net/download.html>
* openmp
  
2. 使用Cmake安装

* cmake >= 3.9 <https://cmake.org/download/>

        $mkdir bulid && cd build
        $cmake ..
        $make
        $../bin/GWLD

3. 源代码安装

        $g++ -std=c++11 -o bin/GWLD src/main.cpp src/filedeal.cpp src/methods.cpp src/utils.cpp -I include -larmadillo -fopenmp

&nbsp;

## 1. 参数列表<!--Parameter list-->

| Short parameter | Long parameter | Parameter Type | Description | Default |
| :----: | :----: | :----: | :----: | :----: |
| -b | --bfile| string | Input plink format (\*.bed \*.bim and \*.fam) file's prefix | - |
| -v | --vcf | string | Input vcf format file name | - |
| -f | --file | string | Input plink format (\*.ped and \*.map) file's prefix | - |
| -m | --method | string | Calculation method between two sites | RMI |
| - | --by-chrom | - | Calculated by chromosome | -|
| - | --decay | - | Output decay according to the chosen method | - |
| - | --max-dist | int | Calculates the maximum distance between two sites while decay parameter selected  | - |
| - | --circos | - | Calculate between different chromosomes | - |
| - | --code012 | - | Recode genotype with -1, 0 ,1, 2. -1 for missing value | - |
| - | --allele-freq  | - | Calculate allele frequency per site | - |
| - | --2vcf | - | Convert plink format file to vcf format file | - |
| -p | --threads | int | threads number | 1 |
| -o | --out-prefix| string | output file name's prefix | gwld |
| -h | --help | - | print help message | - |

&nbsp;<!--添加空行-->

## 2. 主要功能

### 2.1 计算相同染色体上两个位点之间连锁不平衡的值, 可选择的计算方法有($D$, $D'$, $r^2$, $MI$, $RMI$)

        ./demo  -vcf test.vcf -m RMI --by-chrom -p 10 -o result

### 2.2 计算不同染色体上两个位点之间连锁不平衡的值, 可选择的计算方法有($MI$, $RMI$)

        ./demo  -vcf test.vcf -m RMI -p 10 -o result

### 2.3 选择不同方法计算整个基因组的LD decay, 默认计算两个位点之间的最大距离为300Kb, 可使用(--max-dist 参数设置)

        ./demo -vcf test.vcf -m RMI -decay -p 10 -o result

### 2.4 不同染色体间连锁不平衡可视化

        ./demo -vcf test.vcf -m RMI -circos -p 10 -o result

&nbsp;<!--添加空行-->

## 3. 其他功能

### 3.1 文件格式转换

1. Plink文件格式(\*.bed, \*.bim和 \*.fam)转换为VCF文件格式

        ./demo -bfile test -2vcf -o result

2. Plink文件格式(\*.ped and \*.map)转换为VCF文件格式

        ./demo -file test -m -2vcf -o result

### 3.2 统计染色体上的标记数量以及起始相对位置和终止相对位置

        ./demo -vcf test.vcf --chr-info -o result

### 3.3 统计染色体上双等位基因的等位基因频率

        ./demo -vcf test.vcf --allele-freq -o result

### 3.4 基因型重新编码为-1, 0, 1, 2(其中-1 表示为缺失值),并以VCF文件格式输出

1. Plink 格式(\*.bed, \*.bim和 \*.fam)

        ./demo -bfile test --code012 -o result
2. Plink 格式(\*.ped and \*.map)

        ./demo -file test --code012 -o result

3. VCF 格式

        ./demo -vcf test.vcf --code012 -o result
<!--
### 4. 项目结构和效果

* 项目结构
![本地路径](C:/Users/zerel/Desktop/截图_1.jpg "项目结构")
&nbsp;
* 运行效果
![本地路径](C:/Users/zerel/Desktop/截图_2.jpg "项目结构")
-->