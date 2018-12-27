# 关于Entrez与SRA，需要了解的东西

> 刘小泽写于2018.12.27
>
> 本文原文来自Biostar handbook。NCBI应该是最常用到的数据库了，其中的Entrez数据库作为一个一级数据库，整合了PubMed和DNA、蛋白的序列、基因、结构、变异、基因表达等信息。因为数据量实在是太大了，因此我们必须了解怎么获取、下载数据。

### 简介

Entrez是一个法语词汇，我们一般称之为“en-trez”，为了方便获取其中的数据，NCBI提供了一个`Entrez E-utils`的web API，不过这个使用比较麻烦；另外，命令行的下载工具`Entrez Direct` 是经常使用的，简化了下载流程。

https://www.ncbi.nlm.nih.gov/books/NBK179288/

https://www.biostars.org/p/92671/

### 看一下Entrez的基础搜索

我们一般搜索序列是怎么搜索的呢？比如我们要搜索`Ebola virus 1976`，结果返回Ebola病毒的全基因组信息（https://www.ncbi.nlm.nih.gov/nuccore/AF086833.2），对于一条序列整个过程就是这么简单

![](https://upload-images.jianshu.io/upload_images/9376801-23407a581f1a2667.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

> 注意到：**accession number**是`AF086833`，**version number** 是`AF086833.2` 。第一个即便序列日后再做改动，它的编号还是保持不变；第二个随着修改版本而修改小数点后的版本号
>
> 另外之前还有一种索引号叫GI编号，但是NCBI自从2015年以后对于新加的序列都取消了这个编号信息

### 如何使用Entrez Direct？

##### 最常用到的命令就是`efetch`

```shell
# 得到Genbank序列
efetch -db=nuccore -format=gb -id=AF086833 >AF086833.gb
# 得到fasta序列
efetch -db=nuccore -format=fasta -id=AF086833 >AF086833.fa
# 还可以选取其中一段序列
efetch -db=nuccore -format=fasta -id=AF086833 -seq_start=1 -seq_stop=3 
```

另外还可以指定正负链

```shell
# 指定正链
efetch -db=nuccore -format=fasta -id=AF086833 -seq_start=1 -seq_stop=5 -strand=1
>AF086833.2:1-5 Ebola virus - Mayinga, Zaire, 1976, complete genome
CGGAC
# 指定负链
efetch -db=nuccore -format=fasta -id=AF086833 -seq_start=1 -seq_stop=5 -strand=2
>AF086833.2:c5-1 Ebola virus - Mayinga, Zaire, 1976, complete genome
GTCCG
# 可以看到，第二个命令得到的结果中 >AF086833.2:c5-1，与第一个的>AF086833.2:1-5不同，说明这个是反向互补片段
```

##### 利用`esearch`进行搜索

假如有一个项目序列号`PRJNA257197`，一般文章中都会提供这样的序列号，拿到它可以去下载一些信息，而不用登陆NCBI自己去找

```shell
esearch -help 查看帮助信息
# 检索的话需要提供-db 和 -query
 esearch -db nucleotide -query PRJNA257197
 esearch -db protein -query PRJNA257197
# 结果会返回像这样的“环境”信息，这个信息我们目前用不到，但是可以在其他的Entrez direct项目中可以用到(比如下面的序列获取)
<ENTREZ_DIRECT>
  <Db>nucleotide</Db>
<WebEnv>NCID_1_77384293_130.14.22.215_9001_1545880173_1681182876_0MetA0_S_MegaStore</WebEnv>
  <QueryKey>1</QueryKey>
  <Count>249</Count>
  <Step>1</Step>
</ENTREZ_DIRECT>
# 想要得到相关的序列
esearch -db nucleotide -query PRJNA257197 | efetch -format=fasta > genome.fa
esearch -db protein -query PRJNA257197 | efetch -format=fasta > protein.fa
```

我们知道了存在这个xml文件，就可以利用它提取更多的内容

```shell
efetch -db taxonomy -id 9606,7227,10090 -format xml | xtract -Pattern Taxon -first TaxId ScientificName GenbankCommonName Division
# 结果得到
9606	Homo sapiens	human	Primates
7227	Drosophila melanogaster	fruit fly	Invertebrates
10090	Mus musculus	house mouse	Rodents
```

### 如何得到多个SRA数据信息？

说到SRA（Short Read Archive），你应该立刻想到原始数据都存在里面，我们处理数据的第一步也就是下载SRA文件，一般来讲，都会用到几个甚至十几个SRA原始数据，而它们的存储也是有自己的结构的。

- **第一级：NCBI的BioProject**：`PRJN...` （例如：`PRJNA257197`），里面存储了研究的目的等信息，一个project一般与多个sample和多个dataset相关；
- **第二级：NCBI的BioSample**： `SAMN...`或者`SRS...`（例如：`SAMN03254300`） ，存储了生物原材料的信息，每一个样本都有自己的属性
- **第三级：SRA Experiment**：`SRX...`存储了特定样本的测序文库信息等
- **第四级：SRA Run**：`SRR...`或`ERR...` (例如：`SRR1553610`)，存储利用某种测序手段得到的原始数据

结构知道了，那么怎么得到一个BioProject下面的SRA信息呢？虽然不会全部用到，但真的要一个个去网站搜索吗？

利用上面的`esearch + efetch`，我们可以得到一些SRA的信息

```shell
# 有了SRP信息(就是)
esearch -db sra -query PRJNA257197 | efetch -format runinfo > info.csv
# 存储了几十项信息，包括了run，ProjectID,Sample,BioSample等等，用到什么提取什么
```

另外这里我们保存成csv文件，就是为了方便后续命令行处理（一般不用excel处理）

> 好像在生信界，Excel的名声并不是很好。这篇文章可以称的上是**公开吐槽**了

![](https://upload-images.jianshu.io/upload_images/9376801-756a0d918e1d11e6.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

然后我们可以提取第一列看看这个project中的SRR编号

```shell
# 先统计一下有多少SRR数据
cat info.csv | cut -d ',' -f1 | grep SRR | wc -l # 结果有891个
# 还可以根据其他条件过滤一些SRR (例如根据日前)
cat info.csv | grep '2014-08-19' | cut -d ',' -f1 | grep SRR | head -10 > ids.txt

```

然后就可以根据这些SRR号去下载

- 方法一：使用`fastq-dump` ，直接下载并转fastq，还可以指定下载前多少行

  `cat ids.txt | xargs -n1 fastq-dump --split-3 -X 10000 $1 `

- 方法二：如果SRA太大，使用`ascp`高速下载器，下载sra，然后用`fastq-dump`拆分成fastq

  ```shell
  cat ids.txt | while read i;do ascp -v -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T -l200m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR155/$i/$i.sra ./;echo "Downloaded $i" done
  ```
