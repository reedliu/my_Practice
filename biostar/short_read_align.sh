mkdir -p ~/tmp/bwa_test/ref

###################
# 下载基因组(19Kb)[安装entrez-direct]
###################
conda install -c bioconda entrez-direct
efetch -db=nuccore -format=fasta -id=AF086833 > ~/tmp/bwa_test/ref/1976.fa

# 构建索引
ref=~/tmp/bwa_test/ref/1976.fa
bwa index $ref

# ls ~/tmp/bwa_test/ref

###################
# 下载测试数据
###################
mkdir -p ~/tmp/bwa_test/raw && cd ~/tmp/bwa_test/raw

# 获取全部的埃博拉病毒项目的测序数据
esearch -db sra -query PRJNA257197 | efetch -format runinfo > runinfo.csv

# 挑选SRR1972739下载
cat runinfo.csv| grep SRR1972739 | cut -d "," -f 10 | xargs -n1 wget {}

# 解压成fq文件 [为了测试，选取前1万条reads]
fastq-dump -X10000 --split-files SRR1972739

###################
# 比对
###################
mkdir -p ~/tmp/bwa_test/align && cd ~/tmp/bwa_test/align

r1=~/tmp/bwa_test/raw/SRR1972739_1.fastq
r2=~/tmp/bwa_test/raw/SRR1972739_2.fastq
ref=~/tmp/bwa_test/ref/1976.fa

baw mem $ref $r1 $r2 > bwa.sam

# 比对结束看一下tag
cat bwa.sam | cut -f 12-20 | head