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
$ seqkit fx2tab -n -i -a xxx.fq.gz | csvtk -H -t -f 4 -r -i grep  -p "[^ATCG]"
# 因为这里的实例数据都是不错的，没有低质量(可能之前的序列结果又更新了)。可以自己找一个包括N碱基的序列试试，替换掉xxx.fq.gz就好
# csvtk参数表示：
# -H: --no-header-row 
# -t: --tabs 
# -f: --fileds 
# -r: --use-regexp
# -i: --ignore-case
# grep -p: grep --pattern
```

假入找到了低质量碱基，我们想过滤掉

```shell
# 先找到低质量序列ID
$ seqkit fx2tab -n -i -a xxx.fq.gz | csvtk -H -t -f 4 -r -i grep  -p "[^ATCG]" | csvtk -H -t cut -f 1 > id2.txt
# 然后找到这些ID并去除
$ seqkit grep --pattern-file id2.txt --invert-match xxx.fq.gz >clean.fa
```

或者我们想看看这些低质量碱基所在序列的信息（包括基因ID、正负链、起始终止位点等）

```shell
$ seqkit grep --pattern-file id2.txt xxx.fq.gz | seqkit locate --ignore-case --only-positive-strand --pattern N+
# 其中--pattern可以多设置几个，就找出包含多种碱基的序列，比如可以再加一个 --pattern K+
# 这个+表示含有一个或者多个前面的碱基
```

#### 想要移除重复的序列

```shell
$ seqkit rmdup --by-seq --ignore-case duplicated-reads.fq.gz > rmdup.fq.gz
# 如果序列非常大的话，可以用 -m/--md5，利用md5信息代替原始序列信息，来减少内存占用
# 如果想根据序列名进行去重，可以使用参数 --by-name
```

#### 序列定位motif/subsequence/enzyme digest sites

假入现在有一个文件存放的是motifs / enzyme digest sites

```shell
$ cat enzymes.fa
>EcoRI
GAATTC
>MmeI
TCCRAC
>SacI
GAGCTC
>XcmI
CCANNNNNNNNNTGG #表示我们这里就像匹配到CCA和TGG，中间保证9个碱基就好
```

我们为了保证上面的`XcmI`这种情况也能成功进行匹配，需要定义一个`-d或者--degenerate`参数

```shell
seqkit locate --degenerate --ignore-case --pattern-file enzymes.fa viral.1.1.genomic.fna.gz
#或者
seqkit locate -d -i -f enzymes.fa viral.1.1.genomic.fna.gz
```

#### 将fasta序列排序

可以按照序列长度

```shell
$ seqkit sort --by-length viral.1.1.genomic.fna.gz > sorted.fa
# 对于大型的fasta文件，可以使用-2或者--two-pass减少内存使用
```

> 以上都是一些较为常用的小操作，对于更复杂的比如“拆分fasta”、“合并fasta“等，用到可以再去搜索