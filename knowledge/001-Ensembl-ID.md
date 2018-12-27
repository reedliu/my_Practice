# 关于Ensembl ID

我们常见的Ensembl ID比如：

> ENSG00000279928.1
>
> ENSG00000279457.2
>
> 详细说明：https://asia.ensembl.org/Help/Faq?id=488

- 其中，ENS是固定字符，表示它是属于Ensembl数据库的。默认物种是人，如果是小鼠就要用`ENSMUS`开头，关于物种代码：http://www.ensembl.org/info/genome/stable_ids/index.html
- G表示：这个ID指向一个基因；E指向Exon；FM指向 protein family；GT指向gene tree；P指向protein；R指向regulary feature；T指向transcript
- 后面11位数字部分如`00000279928` 表示基因真正的编号
- 小数点后不一定每个都有，表示这个ID在数据库中做了几次变更，比如`.1`做了1次变更，在分析时需要去除

