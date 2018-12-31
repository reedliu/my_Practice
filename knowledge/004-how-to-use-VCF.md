# 变异信息那些事

> 刘小泽写于18.12.31
>
> 一般我们会从WES的上游得到SNP、InDel等信息，这些重要的信息都保存在VCF中，那么怎么对这些变异进行提取、评估与解释呢？一起来学习一下

### VCF是什么？

https://www.jianshu.com/p/957efb50108f

### 如何得到VCF？



### VCF基本操作

#### 利用vcftools



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

- **PROVEAN**：（Protein Variation Effect Analyzer）用来预测SNP或者InDel是否影响蛋白质的生物功能，不仅可以对CDS区域的非同义突变进行预测，还可以对CDS区域的非移码InDel对蛋白功能的影响进行预测，并将结果大致分为：危害、可以容忍、无害
- **SIFT**：（Sorting Intolerant From Tolerant）根据氨基酸在蛋白序列中的保守程度来预测氨基酸的变化对蛋白功能造成的影响。其中保守程度是比对进化关系较近的蛋白序列得到，分值（SIFT-score）表示突变对蛋白序列的影响，**分值越小越严重** 
- **Polyphen2_HAVR**: (Polymorphism Phenotyping v2) 根据HumanVar数据库预测突变对蛋白的影响，来诊断孟德尔遗传病。**分值**表示SNP导致蛋白结构或功能改变的可能性，**越大越严重**
- **Polyphen2_HDIV**： 根据HumanDiv数据库预测**，分值越大越严重**
- **LRT**：也是基于序列保守性进行预测（像SIFT和Polyphen）。对每一个测试的密码子，LRT将来自31个物种的氨基酸进行比对来预测突变的危害。结果的**有害突变（D：Deleterious）**表示：突变来自高度保守的密码子；突变氨基酸在其他比对的真核哺乳动物中不存在。**中性突变（Neutral）**表示：突变发生在非高度保守的密码子；突变的氨基酸至少在一个进行比对的真核哺乳动物中发现

#### 剪切位点突变危害注释

> 如果突变发生在剪切位点附近，我们可以判断它对剪切的危害。可以用的软件有：DbscSNV、Spidex、MaxEntScan

- **DbscSNV**：由AdaBoost与Random Forest开发，它根据突变前后分值的变化来预测剪切位点的突变危害性
- **Spidex**：基于深度学习，因此预测的剪切变异可能距离常规的剪切位点比较远（这一点和DbscSNV不同）
- **MaxEntScan**：