#!/bin/bash
 
######################################################################
#This software is designed to simply the process of sequence alignment and help
#researchers get the blast results fast.
#Author: Yunze Liu
#E-mail: jieandze1314@gmail.com
#Website: www.jieandze1314.com
#Copyright 2019
#License: Creative Commons Attribution-ShareAlike 4.0
#http://creativecommons.org/licenses/by-sa/4.0/legalcode
######################################################################
 
#Set Script Name variable
SCRIPT=`basename ${BASH_SOURCE[0]}`
 
#Initialize variables to default values.
OPT_A=A
OPT_E=E
OPT_T=T
OPT_C=C
 
#Set fonts for Help
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`
 
#Help function
function HELP {
  echo -e "***CREATED BY Yunze Liu (jieandze1314@gmail.com)***"
  echo -e \\n"Help documentation for ${BOLD}${SCRIPT}.${NORM}"\\n
  echo "Command line switches are optional. The following switches are recognized."
  echo "${REV}-t${NORM}  --Sets the type of blast database ${BOLD}ncbidb/owndb${NORM}. Default is ${BOLD}ncbidb${NORM}."
  echo "${REV}-T${NORM}  --If sets -t owndb,you can build your own database with ${BOLD}Target fasta sequence${NORM}."
  echo "${REV}-c${NORM}  --Sets the cpu number for option ${BOLD}c${NORM}. Default is ${BOLD}5${NORM}."
  echo "${REV}-e${NORM}  --Sets ${BOLD}yes/no${NORM} for extraction. Default is ${BOLD}no${NORM}."
  echo "${REV}-S${NORM}  --If sets -e TRUE, you can set ${BOLD}start${NORM} position. Default is ${BOLD}0${NORM}."
  echo "${REV}-E${NORM}  --If sets -e TRUE, you can set ${BOLD}end${NORM} position. "
  echo "${REV}-L${NORM}  --If sets -e TRUE, you can set ${BOLD}length${NORM} for each line. Default is ${BOLD}100${NORM}."
  echo -e "${REV}-h${NORM}  --Displays this help message."\\n
  echo -e "Example 1: ${BOLD}$SCRIPT -t ncbidb -c 5 -e yes -S 0 -E 600 -L 100 ${NORM}"\\n
  echo -e "Example 2: ${BOLD}$SCRIPT -t owndb -T /home/ncbi/database/owndb/own.fasta -c 3 -e no  ${NORM}"\\n
  exit 1
}
 
#Check the number of arguments. If none are passed, print help and exit.
NUMARGS=$#
echo -e \\n"Number of arguments: $NUMARGS"
if [ $NUMARGS -eq 0 ]; then
  HELP
fi
 
### Start getopts code ###
 
#Parse command line flags
type='ncbidb'
target=''
cpu=5
extract=yes
start=0
end=''
length=100

 
while getopts :t:T:c:e:S:E:L:h FLAG; do
  case $FLAG in
    t)  #set option "t"
      type=$OPTARG
      ;;
    T)  #set option "T"
      target=$OPTARG
      ;;      
    c)  #set option "c"
      cpu=$OPTARG
      ;;
    e)  #set option "e"
      extract=$OPTARG
      ;;
    S)  #set option "S"
      start=$OPTARG
      ;;
    E)  #set option "E"
      end=$OPTARG
      ;;  
    L)  #set option "L"
      length=$OPTARG
      ;; 
    h)  #show help
      HELP
      ;;
    \?) #unrecognized option - show help
      echo -e \\n"Option -${BOLD}$OPTARG${NORM} not allowed."
      HELP
      ;;
  esac
done
 
shift $((OPTIND-1))  #This tells getopts to move on to the next argument.
 
### End getopts code ###


cat<<EOF >rename_seq.sh
ls * | fgrep '(' &>/dev/null
RETVAL1=\$?
if [ \$RETVAL1 -eq 0 ];then
	rename \( _ *
	rename \) _ *
fi

ls * | egrep '^[0-9]' &>/dev/null
RETVAL2=\$?
if [ \$RETVAL2 -eq 0 ];then
	ls * | egrep '^[0-9]'| xargs -n1 -i mv {} _{} 
fi

file * | fgrep -v "ASCII text" | fgrep -v "shell" &>/dev/null
RETVAL3=\$?
if [ \$RETVAL3 -eq 0 ];then
	mkdir wrong_files 
	file *| fgrep -v  "ASCII text" | fgrep -v "executable" | awk -F ":" '{print \$1}'| xargs -n1 -i mv {} wrong_files
	echo "wrong files found and moved to wrong_files"
fi
file *.seq | fgrep -w  "with no line terminators" &>/dev/null
RETVAL4=\$?
if [ \$RETVAL4 -eq 0 ];then
	file *.seq | fgrep -w  "with no line terminators" | awk -F ":" '{print \$1}' | xargs -n1 rename .seq .fasta
	for file in *.fasta;do name=\${file%.*};sed -i '1i>'\$name \$file;echo >> \$file ;done
fi
rename .seq .fasta *.seq
cat *.fasta > all.fasta
EOF

sh rename_seq.sh



cat<<EOF >extract_seq.pl
#!/usr/bin/perl -w
use strict;
if (scalar @ARGV!=4){
	die "Usage: This is used to extract specific sequence from fasta 
	perl <fasta> <start pos> <end pos> <length per line>";
}
open FA,"<\$ARGV[0]";
$/=">";<FA>;
while (<FA>){
	chomp;
	my (\$id,\$seq)= (split /\n/,\$_,2)[0,1];
	\$seq=~ s/\n//g;
	my \$subseq=substr(\$seq,\$ARGV[1],\$ARGV[2]);
	\$subseq=~ s/(\w{\$ARGV[3]})/\$1\n/g;
	print ">\$id\n\$subseq\n";

}
EOF




if [ "$extract" == "yes" ]; then
	perl extract_seq.pl all.fasta $start $end $length > all.edit.fasta
	if [ "$type" == "owndb" ]; then
    dir=`dirname $target`
    makeblastdb -in $target -dbtype nucl -parse_seqids -out ${dir}/owndb
    blastn -query all.edit.fasta -db ${dir}/owndb \
      -outfmt 7 \
      -out "all_edit_owndb.csv" -num_threads $cpu
	elif [ "$type" == "ncbidb" ]; then
    blastn -query all.edit.fasta -db nt \
      -outfmt '7 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
      -out "all_edit_ncbidb.csv" -evalue 1e-5 \
      -perc_identity 99 -num_alignments 10 -num_threads $cpu
  else
    echo "** Please check the argument! **"
	fi
elif [ "$extract" == "no" ]; then
	if [ "$type" == "owndb" ]; then
    dir=`dirname $target`
    makeblastdb -in $target -dbtype nucl -parse_seqids -out ${dir}/owndb
    blastn -query all.fasta -db ${dir}/owndb \
      -outfmt 7 \
      -out "all_owndb.csv" -num_threads $cpu
  elif [ "$type" == "ncbidb" ]; then
    blastn -query all.fasta -db nt \
      -outfmt '7 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
      -out "all_ncbidb.csv" -evalue 1e-5 \
      -perc_identity 99 -num_alignments 10 -num_threads $cpu
  else
    echo "** Please check the argument! **"
  fi
else
  echo "** Please check the argument! **"
fi


cat<<EOF >extract_owndb_result.r
options(stringsAsFactors = FALSE)
x <- read.csv('all_edit_owndb.csv',header = F)
col <- unlist(c(x[4,]))
col[1] <- substring(col[1],11,nchar(col[1]))
names(col) <- NULL
x2 <- read.csv('all_edit_owndb.csv',header = F,sep = "\t",comment.char = "#")
colnames(x2) <- col
write.csv(x2,"RESULT_owndb_extract.csv",row.names = FALSE)
EOF

cat<<EOF >extract_ncbidb_result.r
options(stringsAsFactors = FALSE)
x <- read.csv('all_edit_ncbidb.csv',header = F)
col <- unlist(c(x[4,]))
col[1] <- substring(col[1],11,nchar(col[1]))
names(col) <- NULL
x2 <- read.csv('all_edit_ncbidb.csv',header = F,sep = "\t",comment.char = "#")
colnames(x2) <- col
write.csv(x2,"RESULT_ncbidb_extract.csv",row.names = FALSE)
EOF

cat<<EOF >owndb_result.r
options(stringsAsFactors = FALSE)
x <- read.csv('all_owndb.csv',header = F)
col <- unlist(c(x[4,]))
col[1] <- substring(col[1],11,nchar(col[1]))
names(col) <- NULL
x2 <- read.csv('all_owndb.csv',header = F,sep = "\t",comment.char = "#")
colnames(x2) <- col
write.csv(x2,"RESULT_owndb.csv",row.names = FALSE)
EOF

cat<<EOF >ncbidb_result.r
options(stringsAsFactors = FALSE)
x <- read.csv('all_ncbidb.csv,header = F)
col <- unlist(c(x[4,]))
col[1] <- substring(col[1],11,nchar(col[1]))
names(col) <- NULL
x2 <- read.csv('all_ncbidb.csv',header = F,sep = "\t",comment.char = "#")
colnames(x2) <- col
write.csv(x2,"RESULT_ncbidb.csv",row.names = FALSE)
EOF



if [ "$extract" == "yes" ]; then
  if [ "$type" == "owndb" ]; then
	  Rscript extract_owndb_result.r
  else 
    Rscript extract_ncbidb_result.r
  fi
else 
  if [ "$type" == "owndb" ]; then
    Rscript owndb_result.r
  else 
    Rscript ncbidb_result.r
  fi
fi

