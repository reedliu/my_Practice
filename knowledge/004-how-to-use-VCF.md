# 变异信息那些事

> 刘小泽写于18.12.31
>
> 上篇：主要了解VCF的背景知识；
> 中篇：VCF怎么来？怎么进行一些小操作来获取内部信息？
> 下篇：得到变异位点需要做些什么
>
> 一般我们会从WES的上游得到SNP、InDel等信息，这些重要的信息都保存在VCF中，那么怎么对这些变异进行提取、评估与解释呢？一起来学习一下

### VCF（Variant Call Format）是什么？

之前也写过一篇相关的，这次想要更深层次去了解它https://www.jianshu.com/p/957efb50108f

> 我们知道，variant calling(找变异)的过程发生在alignment(比对)之后，那么肯定流程更加复杂，因此variant calling得到的结果也要更加精炼、内容更加丰富。于是VCF文件接手了这个棘手的工作。了解VCF，对于想要另辟蹊径发现新研究内容的人来说，真的是一块宝藏，就看你怎么挖掘了。

#### 主要内容

得到一个VCF文件，首先看到的就是它的**Header（表头）**，如下：（其实有非常非常多的头信息…这里只写几行）

```shell
##fileformat=VCFv4.1
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##GATKCommandLine.HaplotypeCaller=<ID=HaplotypeCaller,Version=3.4-3-gd1ac142,Date="Mon May 18 17:36:4
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##contig=<ID=chr1,length=249250621,assembly=b37>
##reference=file:human_genome_b37.fasta
```

第一行是VCF的版本信息，看似没用，但实际上我们在分析其他数据时，并不能保证一直使用最新的VCF格式，因此检查下VCF版本确保后续提取正确

`FILTER`行是说过滤了什么内容；

`FORMAT`和`INFO` 相当于变异位点的注释信息；

`CommandLine` 是说使用的call variant工具信息；

`contig`和`reference` 是当重复别人数据时，恰好没告诉你数据来源，这时就可以参考这个

接下来才是重点：**Records信息**

包含至少8列tab分割的常规信息（用来描述变异位点），第9列及以后表示各个样本的变异信息（可以包括成百上千个样本）

前9列信息包括：

- CHROM：变异发生的chromosome或者contig

- POS：变异发生的基因组坐标（对缺失来讲，显示的是缺失开始的位置）

- ID：一般是dbSNP的ID（可有可无）

- REF：Forward strand（正链）上的参考等位基因
  补充一下基础知识：
  ![image.png](https://upload-images.jianshu.io/upload_images/9376801-6e1b30707190400c.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

  等位基因（allele）：一对同源染色体相同位置上控制某一性状的不同形态的基因，可以简单想成控制同一性状的不同基因

- ALT：正链上对应的改变的等位基因（可以有多个）

- QUAL：`REF/ALT` 变异位点存在的概率（和FASTQ的质量值、SAM的MAPQ一样，都是Phred值 `-10 * log(1-p)`）

- FILTER：数据一般都要经过适当的过滤后才能继续使用variant callset（变异位点数据集），关于是否完成过滤会给出三种说明：
  一是给出没有通过过滤的变异位点；二是`PASS`表示全部通过了过滤；三是`.` 表示这个位点没有任何过滤

- INFO：以`tag=value`的形式给出位点注释信息；分号`;`分割

- FORMAT：针对样本的注释；冒号`:`分割，并且对应后面各个sample中的信息

```shell
# 大体就是这样
#CHROM 		POS	ID	REF	ALF	QUAL	FILTER	INFO	FORMAT	align.bam
AF086833	60	.	A	T	54	.	DP=43	GT:PL	0/1:51,0,48	
```

> 上面是关于VCF的大体了解，下面具体看看重要的部分

#### 具体--REF/ALT

`REF`表示在`POS`位置的参考等位基因；`ALT`表示在`POS`位置的**所有**变异

开始想象场景：我们用软件发现，在60bp处发现了**一个变异**，是A变成了T，那么应该这么表示：

```shell
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	align.bam
AF086833	60	.	A	T	54	.	DP=43	GT:PL	0/1:51,0,48
```

如果**同一个位置检测到有两种**可能发生的变异呢？

```shell
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	align.bam
AF086833	60	.	A	T,C	43.2	.	DP=95	GT:PL	1/2:102,124,42,108,0,48
```

事情还没完，刚刚两个例子都只是一个变异位点，但实际上有可能**多个位点**发生缺失（**而这才是VCF复杂的开始**），例如：58、59、60碱基（GGA）发生了缺失

```shell
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	align.bam
AF086833	55	.	TGAGGA	TGA	182	.	INDEL;IDV=42	GT:PL	0/1:215,0,255
AF086833	60	.	A	C,T	39.8	.	DP=126	GT:PL	1/2:99,87,47,86,0,49
```

这里注意：虽然我们看到缺失是从58碱基开始发生的，但是记录的是从55碱基。**你会不会好奇**：为什么是55而不是准确的58？为什么要表示成`TGAGGA->TGA`，而不是`GAGGA->GA`或者`GGA->空 ` 呢？这个需要引入一个新词汇”variant normalizatoin“，也就是说，所有的结果展示都是有规定的。

> variant normalization:  2015年Unified representation of genetic variants文章中就论述了这个一致性标记的问题，其中写道：
>
> A genetic variant can be represented in the Variant Call Format (VCF) in multiple different ways. Inconsistent representation of variants between variant callers and analyses will magnify discrepancies between them and complicate variant filtering and duplicate removal.
>
> 于是制定了VCF的标准化，来方便交流，规定以下几点：
>
> 用尽可能少的字母来表示变异位点；
>
> 等位基因长度不为0；
>
> 变异向左”贪婪比对“ （也就是说：一直向左比对，直到不匹配为止，然后以最左边的碱基位置表示变异的起始位置）
>
> https://academic.oup.com/bioinformatics/article/31/13/2202/196142

因此，这里写成`TGAGGA->TGA`就是表示缺失位点`GGA`向左最远可以匹配到第55号T碱基处

#### 具体--FORMAT

假设得到的VCF中9-11列内容如下:

```shell
#FORMAT		sample1			sample2
GT:PL		0/1:51,0,48		1/1:34,32,0
```

表示的意思就是：在样本1中发现的变异中含有`GT=0/1`和`PL=51,0,48`这样的信息，样本2中的变异含有`GT=1/1`和`PL=34,32,0`这样的信息

看到两个名词`GT`和`PL`，那么它们是什么意思呢？

我们可以回过头去Header部分看一看，将会找到如下内容：

```shell
##FORMAT=<ID=GT, Number=1, Type=String, Description="Genotype">
##FORMAT=<ID=PL, Number=G, Type=Integer, Description="List of Phred-scaled genotype likelihoods">
```

还感觉看不太懂？不着急，往下接着学

#### 具体--Genotypes

尽管我们可以根据`REF和ALT`知道了碱基发生了怎样的变化，但是我们想知道这个变化是只是发生在一个DNA拷贝中，还是两个拷贝都有？需要用一个指标来**量化**这种变异~**基因型（Genotype）**。`GT`就是用来表示样本中这个位点的基因型，其中0表示参考`REF`，1表示变异`ALT`的第一个entry，2表示`ALT`的第二个entry（以此类推）

对于二倍体生物，GT表示了一个样本中的等位基因：

- `0/0` 表示样本是纯合子，并且和参考的等位基因一样
- `0/1`表示样本是杂合子，有一个参考的等位基因，一个变异的等位基因
- `1/2` 样本是杂合子，两个都是变异的等位基因
- `1/1` 样本是纯合子，且两个都是变异的第一个等位基因
- `2/2`样本是纯合子，且两个都是变异的第二个等位基因（以此类推）

当然，如果不是二倍体，命名原理也是一样：单倍体（Haploid）只有一个GT值；多倍体有多个GT值

#### 具体--Genotype likelihoods

直白地说就是”基因型可能性“，就是用来衡量不同基因型可能发生的概率，这是利用p-value统计，因此**0表示可能性最大**，例如：

```shell
GT:PL	0/1:51,0,48
```

其中`PL`这一项有三个数值，分别对应三种可能的基因型（`0/0`，`0/1`，`1/1`）发生的概率：第一个数值51表示基因型为`0/0`的概率是`Phred值51` ，也就是`1x10^6` ；第二个数值0表示基因型为`0/1`的概率是0（和GT判断的一致）；第三个数值48表示基因型为`1/1`的概率是`1x10^5` 

#### 具体--Allele depth and depth of coverage

软件判断是那种基因型，到底是不是发生了变异，是需要一定的统计方法的，主体就是之前比对的结果BAM文件，其中包含了reads的比对信息，这里就是根据reads比对的数量进行判断

- **AD**（DepthPerAlleleBySample):  **unfiltered** allele depth 就是有多少reads出现了某个等位基因（其中也包含了没有经过variant caller过滤的reads），但是不包括没意义的reads（就是那些统计结果不过关的，没法说服软件相信这个等位基因）
- **DP**（Coverage）： **filtered** depth, at the sample level 只有通过variant caller软件过滤后的reads才能计算入内，但是DP也纳入了那些经过过滤但没有意义的reads（uninformative reads），这一点又和AD不同

#### 具体--Genotype Quality

GQ就是用Phred值来表示GT判断的准确性，它和PL相似，但是取值不同。PL最小值0表示最准确，**GQ一般取PL的第二个小的值**（除非第二小的PL大于99）。**在GATK分析中，GQ最大就限制在99**，因为超过99其实是没有什么意义的，并且还多占字符。因此，如果GATK中发现PL值中第二个小的值比99还要大，软件就将GQ标为99。用GQ值就可以得到，第一位和第二位之间到底差了多少，因此可以快速判断分析的准不准，选择第一个靠不靠谱

#### 做一个小总结--VCF帮助判断基因型

> 现在来看看NA12878基因（1: 899282）的统计情况

```shell
1   899282  rs28548431  C   T   [CLIPPED] GT:AD:DP:GQ:PL    0/1:1,3:4:26:103,0,26
```

这个位点的`GT=0/1`，可以判断基因型是`C/T`；`GQ=26`表示排名第二基因型的可能性是0.0025，结果不是很好，因为虽然判断的第一位基因型的`PL` 为0很可靠，但是毕竟相差不多，很难推翻第二位（如果GQ值再大一些，我们就更有信心说明判断的`C/T`基因型是正确的）。当然，这里GQ的原因很有可能是取样太少，只有4条reads在这个位点作为参考（`DP=4`），这四条中有1条带参考的碱基信息，另外3条与参考不一致，存在变异（`AD=1,3`）

**因此，重点的结论来啦！**尽管我们相信，这个位点确实存在变异，但是假阳性依然存在，也就意味着基因型判断结果不一定是杂合子`C/T`，有一定的可能是变异纯合`T/T`(`PL(1/1)=26`)，但一定不可能是参考纯合`C/C` （`PL(0/0)=103`）

参考：

VCF short summary：http://www.htslib.org/doc/vcf.html

VCF Poster：http://vcftools.sourceforge.net/VCF-poster.pdf

VCF 说明书: http://samtools.github.io/hts-specs/

GATK解释VCF：https://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it

Difference between QUAL and GQ annotations  https://software.broadinstitute.org/gatk/documentation/article.php?id=4860

---

### 如何得到VCF？

VCF一般是利用"variant caller"或者"SNP caller"的工具对BAM比对文件操作得到的

#### bowtie2+samtools+bcftools

> 需要用到bowtie2的测试数据集 
> 之前要生成VCF或者它的二进制BCF，需要用samtools mpileup，后来bcftools将mpileup加入了自己的功能中，避免了使用mpileup+bcftools call pipeline时版本冲突报错的问题，直接一步到位

```shell
# 下载bowtie2
cd ~/test
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.3/bowtie2-2.3.4.3-linux-x86_64.zip 
unzip bowtie2-2.3.4.3-linux-x86_64.zip 
cd bowtie2-2.3.4.3-linux-x86_64/example/reads

# bowtie2 比对 + samtools排序
wkd=~/test/bowtie2-2.3.4.3-linux-x86_64

$wkd/bowtie2 -x $wkd/example/index/lambda_virus -1 reads_1.fq -2 reads_2.fq | samtools sort -@ 5 -o lamda.bam -

bcftools mpileup -f $wkd/example/reference/lambda_virus.fa lamda.bam |bcftools call -vm > lamda.vcf 
```

#### GATK

```shell
# From https://www.biostars.org/p/307422/
1. trim reads
2. bwa mem align to genome
3. mark duplicates (e.g. picardtools =》MarkDuplicates)
4. use HaplotypeCaller to generate gvcf
5. CombineGVCFs
6. GenotypeGVCFs on the combined gvcf
7. filter your vcf however you want
8. You can do base recalibration iteratively now if you want with the filtered vcf.
```



参考：http://www.bio-info-trainee.com/3577.html

bcftools 说明书 https://samtools.github.io/bcftools/bcftools.html

拓展阅读：GATK的深度学习Convolutional Neural Networks（CNN）vs. Google DeepVariant vs. GATK VQSR  https://www.biostars.org/p/301751/ [其中 Andrew认为GATK是寻找变异的金标准，但是程序有点繁琐，好用不好入门。有推荐使用samtools+bcftools的，并且它们在鉴定SNVs有优势，但InDel方面就欠缺一些]

https://gatkforums.broadinstitute.org/gatk/discussion/10996/deep-learning-in-gatk4

【**！放在公众号！**】小福利：2018.11最新的GATK Worksheet 链接：https://share.weiyun.com/5vefLMG 密码：yf5322

---

### VCF基本操作

> 得到了VCF文件只是第一步，下面还需要一定的技巧获得想要的数据

#### 利用bcftools

##### 下载测试数据及索引

```shell
# 测试文件是19号染色体：400-500kb
wget http://data.biostarhandbook.com/variant/subset_hg19.vcf.gz
wget http://data.biostarhandbook.com/variant/subset_hg19.vcf.gz.tbi
```



参考：

#### 利用GATK

> 我们想要找到可靠的结果，一般会采用两种以上的软件进行分析，对于找变异也是如此。
>
> 想象几个场景：我们用不同的软件对同一个样本进行分析，然后想找到它们共同包含的变异信息；对于不同的样本，我们想知道它们有哪些相同变异，又有哪些不同变异；我们自己的样本中有没有比较特殊的变异...

##### GATK SelectVariants

GATK的这个功能就是用来选择VCF文件中符合指定条件的变异，它的基本用法：

```shell
# 输入参数
-V : 原始的VCF文件 [Required参数]
-select ： 指定特定条件(e.g. AF<0.25; DP>1000)
-sn/sf/se : 指定筛选的样本名称/文件/正则条件
# 输出参数
-o : 筛选后的VCF文件
# 筛选设置
-ef : 去掉不符合筛选条件的变异
-env : 去掉没有变异的行
-selectType : 选择特定种类的变异(SNP, InDel ...)
-xl_{sn/sf/se}: 去掉特定的样本
```

输入参数中，`sn（sampleName）`表示单个样本名，`sf(sampleFile)`表示包含多个样本名的文件，`se(sampleExpressions)`表示用正则表示出样本名。

需要注意的是：SelectVariants功能只是标出了符合条件的变异，并不会默认去除不符合条件的变异。
如果只想保存符合条件的变异，使用`ef(excludeFilteredVariants)` ；
如果从多个样本中选取变异信息，但其实它们都存在共同的没有变异的位点，就可以用`env(excludeNonVariants)`去除；
利用`ef与env`，可以去除不是这个人种的变异位点

例如：我们想看看自己的样本中有多少变异存在于1000Genome中

```shell
java -jar GATK.jar \
	-R reference.fasta \
	-T SelectVariants \
	-V 1000G_ALL_Genotypes.vcf.gz \
	-sf samples.list \ # 将文件名都存在这个list中
	-ef \ #--excludeFiltered 保留PASS与.的位点
	-env \ #去除样本中没有在1000Genome中出现的位点
	-o 1000G_.vcf \
# 最后生成的文件中只包含过滤没有fail并且至少一个样本中有变异的位点
```

> 更详细的参数说明参考官方：https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php

##### GATK CombineVariants

上面是说提取符合特定条件的子集，许多时候还需要合并多个VCF文件。当合并多个文件时，如果包含同一个位点，但位点信息不一样，就需要考虑设置合并参数。基本用法：

```shell
# 输入参数
-V : 多个vcf文件(多个文件可以添加tag描述，用于区分)
# 可选参数
-genotypeMergeOptions ：#关于合并同一个样本的不同基因型
-filteredAreUncalled ： #忽略filtered变异位点
-filteredRecordsMergeType ： #如何对待filter status不一样的位点
-priority #对于多个vcf文件，设置它们的基因型优先级(需要配合之前添加-V时的tag描述，逗号分隔，来表示不同的文件)
# 输出参数
-o : 
```

`genotypeMergeOptions`的设置：如果目的是`Merge two separate callsets`，就可以用`UNIQUIFY`；如果目的是`Get the union of calls made on the same samples`，就可以用`PRIORITIZE` 

`filteredRecordsMergeType`的设置：`KEEP_IF_ANY_UNFILTERED`是至少有一个filter status不一致就不合并这个位点；`KEEP_IF_ALL_UNFILTERED`是当所有的status不一致才不合并这个位点；`KEEP_UNCONDITIONAL`不管filter status，只要是位点信息，就合并

例如，我们有两个vcf文件需要合并

```shell
java -jar GATK.jar \
	-R reference.fa \
	-T CombineVariants \
	-V:hua file1.vcf \ # 如果文件比较多，可以写到list中，然后只用-V file.list
	-V:dou file2.vcf \
	-filteredAreUncalled \
	-o output.vcf
# 结果文件中，包含了两个输入文件中的变异，每个变异位点都有一个"set=="栏，标识位点来源，可以来自file1可以来自file2，或者来自二者
```

> 官方说明：https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_CombineVariants.php

##### Combine联合Select

例如，我们想找到某个病人特有变异

```shell
java -jar GATK.jar \
	-R reference.fa \
	-T CombineVariants \
	-V:1000G 1000G_ALL_Sites.vcf \ #ESP6500 SNP文件
	-V:ESP ESP6500.snps.vcf \ #1000G数据库文件
	-V:D2D D2D.vcf \ #假设使我们得到的vcf文件
	-o patient.D2D.combined.vcf # 合并的文件
java -jar GATK.jar \
	-R reference.fa \
	-T SelectVariants \
	-V patient.D2D.combined.vcf \ # 加载合并好的vcf文件
	-select "set==D2D" \ # 选择标记的D2D位点
	-o patient.D2D.private.vcf #仅在该病人中有，在ESP6500和1000G中没有的变异位点
```

### 变异位点的注释

> 我们得到变异位点，但仅仅是知道了它们在基因组上的位置信息和相关的碱基信息。那么还存在许许多多的疑问没有解决：
>
> 这个位点是在基因上吗？是内含子还是外显子区域？这个突变对基因功能产生了什么影响？对于转录翻译有没有影响？除了研究的样本，还有没有其他样本也出现了这个变异？有的话是什么人种，又是什么病例？
>
> 这些问题都要靠**变异注释**来解决

一般来说，**变异注释分为**：突变频率注释、变异的蛋白功能危害注释、剪切位点突变危害注释、突变相关的疾病注释

#### 突变频率注释

> 做这个内容的数据库有许多，其中比较重要的有dbSNP、1000人基因组项目（1000 Genome）、ExAC、gnomeAD

- **dbSNP**（The single-nucleotide polymorphism database）：http://www.ncbi.nlm.nih.gov/SNP/  NCBI与人类基因组研究所合作建立，包含了SNP、短重复序列、微卫星标记等来源、检测方法、基因型信息、上下游序列、人群分布频率等

- **1000G （千人基因组项目）** 研究时限：2008-2015年；汇集30个人种、3904个样本WGS和WES测序结果。目前已被ANNOVAR收纳为变异位点在正常人群中进行突变频率注释的数据库，实际分析中也应该将1000G的不同人群作为control组进行疾病关联分析

  它的构建总共分为**4个阶段：**
  **一、Pilot phase** 【A map of human genome variation from population-scale sequencing *Nature 467, 1061–1073 (28 October 2010)*】

  ![image.png](https://upload-images.jianshu.io/upload_images/9376801-a2a71d99625073e2.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

  **二、Phase one** 【An integrated map of genetic variation from 1,092 human genomes *Nature 491, 56–65 (01 November 2012)*】

  **三、Phase two**

  **四、Phase three** 【A global reference for human genetic variation *Nature 526, 68–74 (01 October 2015)*； An integrated map of structural variation in 2,504 human genomes *Nature 526, 75–81 (01 October 2015)*】

  参考：http://www.internationalgenome.org/about

  含有两个**下载镜像**：ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/

  ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/

  其中**Phase3所有的样本信息**下载：ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/

- **ExAC**（Exome Aggregation Consortium）：整合了60706个人的WES测序数据及相关遗传信息，包含超过1000万种基因变异信息 http://exac.broadinstitute.org/ 。包括了AFR（African）、AMR（Admixed American）、EAS（East Asian）、FIN（Finnish）、NFE（NON-finnish European）、SAS（South Asia）等种群的突变频率（AF）信息http://exac.broadinstitute.org/faq

  ![image.png](https://upload-images.jianshu.io/upload_images/9376801-36ae5e7cad2c4635.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

- **gnomeAD**：（Genome Aggregation Database）博得研究所支持建立，包含了千人基因组、ESP数据库以及绝大部分的ExAC数据库。目前有125,748个外显子数据和15,708个基因组数据 http://gnomad.broadinstitute.org/

#### 变异的蛋白功能危害注释

- **PROVEAN**：（Protein Variation Effect Analyzer）http://provean.jcvi.org/index.php 用来预测SNP或者InDel是否影响蛋白质的生物功能，不仅可以对CDS区域的非同义突变进行预测，还可以对CDS区域的非移码InDel对蛋白功能的影响进行预测，并将结果大致分为：危害、可以容忍、无害
- **SIFT**：（Sorting Intolerant From Tolerant）https://sift.bii.a-star.edu.sg/ 根据氨基酸在蛋白序列中的保守程度来预测氨基酸的变化对蛋白功能造成的影响。其中保守程度是比对进化关系较近的蛋白序列得到，分值（SIFT-score）表示突变对蛋白序列的影响，**分值越小越严重** ，一般认为：SIFT值小于0.05为有害（D：Deleterious），大于0.05表示容忍（T：Tolerance）
- **Polyphen2_HAVR**: (Polymorphism Phenotyping v2)  http://genetics.bwh.harvard.edu/pph2/dokuwiki/downloads 根据HumanVar数据库预测突变对蛋白的影响，来诊断孟德尔遗传病。**分值**表示SNP导致蛋白结构或功能改变的可能性，**越大越严重** 
- **Polyphen2_HDIV**： 根据HumanDiv数据库预测**，分值越大越严重**
- **LRT**：也是基于序列保守性进行预测（像SIFT和Polyphen）http://www.genetics.wustl.edu/jflab/lrt_query.html 。对每一个测试的密码子，LRT将来自31个物种的氨基酸进行比对来预测突变的危害。结果的**有害突变（D：Deleterious）**表示：突变来自高度保守的密码子；突变氨基酸在其他比对的真核哺乳动物中不存在。**中性突变（N: Neutral）**表示：突变发生在非高度保守的密码子；突变的氨基酸至少在一个进行比对的真核哺乳动物中发现

#### 剪切位点突变危害注释

> 如果突变发生在剪切位点附近，我们可以判断它对剪切的危害。可以用的软件有：DbscSNV、Spidex、MaxEntScan

- **DbscSNV**：属于VEP（Variant Effect Predictor）插件，http://asia.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#plugins_existing

  由AdaBoost与Random Forest开发，它根据突变前后分值的变化来预测剪切位点的突变危害性

- **Spidex**：http://www.openbioinformatics.org/annovar/spidex_download_form.php 基于深度学习，因此预测的剪切变异可能距离常规的剪切位点比较远（这一点和DbscSNV不同）

- **MaxEntScan**：对**5'**剪切位点附近的6bp内含子与3bp的编码区(http://genes.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html)以及**3‘** 剪切位点附近20bp的内含子与3bp的编码区内突变进行预测(http://genes.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq_acc.html)，按照突变前后分值变化来得到结论（认为有危害：突变后比突变前分值降低15%以上）

#### 突变相关的疾病注释

- **OMIM**：（Online Mendelian Inheritance in Man）https://www.omim.org/在线人类孟德尔遗传信息数据库，包含了遗传性的基因疾病信息与表型信息，目前收录了16000多个基因词条和5400多表型词条

- **HGMD**：（The Human Gene Mutation Database）1996年创立的人类基因突变数据库，目前包括240,269个变异，覆盖9976个基因。收集的突变包含了SNP、InDel、CNV、SV、基因重组等，可以说是遗传病变异检测金标准数据库。有两个版本，一个是免费的学术public版，但更新慢（http://www.hgmd.cf.ac.uk/ac/index.php）；另一个是收费可试用的Professional版（https://www.qiagenbioinformatics.com/products/human-gene-mutation-database/），包含的变异数量也更多

- **ClinVar**：2013年创立，是一个已报道突变与疾病表型关联数据库，https://www.ncbi.nlm.nih.gov/clinvar/。数据主要来源是OMIM、dbSNP、locus specific database等开源数据库，对变异位点的审核比较缺乏，因此会包含报道中冲突的致病位点

  ![HGMD and ClinVar: Avoiding the Knowledge Blind Spot](https://upload-images.jianshu.io/upload_images/9376801-da9e3993df8d3ff2.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
