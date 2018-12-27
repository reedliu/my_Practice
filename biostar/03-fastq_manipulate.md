# 对fasta/fastq进行一些小操作

> 刘小泽写于18.12.27
>
> 参考biostar handbook以及Wei Shen的SeqKit（处理fa/fq领域还比较出名）http://bioinf.shenwei.me/seqkit/，它可以处理压缩的或者解压的fa/fq

### 下载测试数据

```shell
wget http://data.biostarhandbook.com/reads/duplicated-reads.fq.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz
```

### 小操作们

#### 先看一下整体的数据统计

```shell
$ seqkit stat *.gz
file                      format  type     num_seqs    sum_len  min_len   avg_len  max_len
duplicated-reads.fq.gz    FASTQ   DNA        15,000  1,515,000      101       101      101
viral.1.1.genomic.fna.gz  FASTA   DNA             7    197,450    1,608  28,207.1  154,675
viral.1.protein.faa.gz    FASTA   Protein       113     56,244       38     497.7    3,122
```

#### 关于GC含量

使用`fx2tab` 将fq/fa的统计信息转为制表符分割的表格，(fasta默认包括ID、sequence；fastq默认包括ID、sequence、quality）；另外可以自己指定length、GC等信息

```shell
# 默认输出：ID、sequence
$ seqkit fx2tab viral.1.1.genomic.fna.gz
# 指定输出(输出name、GC、length)
$ seqkit fx2tab --name --only-id --gc --length 
viral.1.1.genomic.fna.gz

NC_024015.1			1608	44.22
NC_001798.2			154675	70.38
NC_030692.1			8908	50.10
NC_027892.1			8957	40.57
NC_029642.1			9006	39.88
NC_001607.1			8910	50.07
NC_001422.1			5386	44.76
# [--name可以用-n替代，--only-id可以用-i替代]
```

当然，除了GC含量，任何碱基的含量也可以得到

```shell
$ seqkit fx2tab -H -n -i -B a -B c -B ac viral.1.1.genomic.fna.gz
#name	seq	qual	a	c	ac
NC_024015.1			27.18	26.43	53.61
NC_001798.2			14.87	35.03	49.90
NC_030692.1			25.03	25.25	50.28
NC_027892.1			29.68	19.44	49.11
NC_029642.1			30.14	19.22	49.36
NC_001607.1			25.30	25.54	50.84
NC_001422.1			23.97	21.48	45.45
```

#### 如何根据ID得到fa/fq的子集？

这个还真是一个高频问题，许多时候就想要一个fastq/fasta中的某些ID的信息，利用`seqkit grep`可以方便解决

```shell
# 先做一个测试数据，包括一些ID信息(就是从duplicated-reads.fq.gz中随机取样)
$ seqkit sample --proportion 0.001 duplicated-reads.fq.gz | seqkit seq --name --only-id >id.txt
SRR1972739.2996
SRR1972739.3044
SRR1972739.3562
SRR1972739.4162
SRR1972739.5975
# 然后根据ID抽取其中的序列
$ seqkit grep --pattern-file id.txt duplicated-reads.fq.gz > dup_sub.fq.gz
# 用这个检查一下是不是这些ID
$ seqkit seq --name --only-id dup_sub.fq.gz
```

#### 找到质量差的碱基（非ATCG）并定位它们

先用conda安装好`csvtk` 

```shell
$ eqkit fx2tab -n -i -a xxx.fq.gz | csvtk -H -t -f 4 -r -i grep  -p "[^ATCG]"
# 因为这里的实例数据都是不错的，没有低质量(可能之前的序列结果又更新了)
# csvtk参数表示：
# -H: --no-header-row 
# -t: --tabs 
# -f: --fileds 
# -r: --use-regexp
# -i: --ignore-case
# grep -p: grep --pattern
```

假入找到了低质量碱基，我们想过滤掉

