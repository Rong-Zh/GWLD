# GWLD manual (an R package for genome-wide linkage disequilibrium analysis) 

## 1. Installation

1. Dependences

* armadillo <http://arma.sourceforge.net/download.html>
* openmp

2. Cmake

* cmake >= 3.9 <https://cmake.org/download/>

``` bash
    $ mkdir bulid && cd build
    $ cmake ..
    $ make
    $ ../bin/GWLD
```

1.3 compiling source codes
``` bash
    $ g++ -std=c++11 -o bin/GWLD src/main.cpp src/filedeal.cpp src/methods.cpp src/utils.cpp -I include -larmadillo -fopenmp
```
 
&nbsp;

## 1. 参数列表<!--Parameter list-->

| Short parameter | Long parameter | Parameter Type |                                   Description                                    | Default |
| :-------------: | :------------: | :------------: | :------------------------------------------------------------------------------: | :-----: |
|       -b        |    --bfile     |     string     |           Input plink format (\*.bed \*.bim and \*.fam) file's prefix            |    -    |
|       -v        |     --vcf      |     string     |                            Input vcf format file name                            |    -    |
|       -f        |     --file     |     string     |               Input plink format (\*.ped and \*.map) file's prefix               |    -    |
|       -m        |    --method    |     string     |                       Calculation method between two sites                       |   RMI   |
|        -        |   --by-chrom   |       -        |                             Calculated by chromosome                             |    -    |
|        -        |    --decay     |       -        |                   Output decay according to the chosen method                    |    -    |
|        -        |   --max-dist   |      int       | Calculates the maximum distance between two sites while decay parameter selected |    -    |
|        -        |    --circos    |       -        |                     Calculate between different chromosomes                      |    -    |
|        -        |   --code012    |       -        |              Recode genotype with -1, 0 ,1, 2. -1 for missing value              |    -    |
|        -        | --allele-freq  |       -        |                       Calculate allele frequency per site                        |    -    |
|        -        |     --2vcf     |       -        |                   Convert plink format file to vcf format file                   |    -    |
|       -p        |   --threads    |      int       |                                  threads number                                  |    1    |
|       -o        |  --out-prefix  |     string     |                            output file name's prefix                             |  gwld   |
|       -h        |     --help     |       -        |                                print help message                                |    -    |

&nbsp;<!--添加空行-->

## 2 Main function

To calculate the linkage disequilibrium (LD) measures, including D, D', $r^{2}$, mutual information (MI) and reduced MI (RMI), the following options can be chosen.

### 3.1 Loci on the same chromosome (D, D', $\mathbf{r}^{\mathbf{2}}$, MI and RMI) 

``` bash
    $ ./GWLD  -vcf test.vcf -m RMI --by-chrom -p 10 -o result
```
    
### 3.2 Loci on two different chromosomes (MI and RMI) 

``` bash
    $ ./GWLD  -vcf test.vcf -m RMI -p 10 -o result
```

### 3.3 LD decay at the global level (default, maximum distance 300kb between two loci, using the --max-dist parameter) 

``` bash
    $ ./GWLD -vcf test.vcf -m RMI -decay -p 10 -o result
```
### 3.4 Circos-like map visualization between different chromosomes 

``` bash
    $ ./GWLD -vcf test.vcf -m RMI -circos -p 10 -o result
```
 

## 4. Other functions 

### 4.1 Exchange between different input file formats 

Files in the plink format can be reformatted into variant calling format(VCF).

4.1.1 From Plink file formats (\*.bed, \*.bim and \*.fam) to VCF

``` bash
    $ ./GWLD -bfile test -2vcf -o result
```
4.1.2 Plink file formats (\*.ped and \*.map) to VCF

``` bash
    $ ./GWLD -file test -m -2vcf -o result
```
### 4.2 Numbers, start and end positions of markers on each chromosome 

``` bash
    $ ./GWLD -vcf test.vcf --chr-info -o result
```
### 4.3 Allele frequency {#allele-frequency .unnumbered}

``` bash
    $ ./GWLD -vcf test.vcf --allele-freq -o result
```
### 4.4 Genotypes recoded into 0, 1, 2, -1 (missing values), and output as VCF files 

4.4.1 Plink format (\*.bed, \*.bim和 \*.fam)

``` bash
    $ ./GWLD -bfile test --code012 -o result
```
4.4.2 Plink format (\*.ped and \*.map)

``` bash
    $ ./GWLD -file test --code012 -o result
```
4.4.3 VCF format

``` bash
    $ ./GWLD -vcf test.vcf --code012 -o result
```