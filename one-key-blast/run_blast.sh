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
  echo -e \\n"Help documentation for ${BOLD}${SCRIPT}.${NORM}"\\n
  echo -e "${REV}Basic usage:${NORM} ${BOLD}$SCRIPT file.ext${NORM}"\\n
  echo "Command line switches are optional. The following switches are recognized."
  echo "${REV}-e${NORM}  --Sets TRUE/FALSE for extraction of fasta ${BOLD}e${NORM}. Default is ${BOLD}FALSE${NORM}."
  echo "${REV}-S${NORM}  --If sets -e yes, you can set ${BOLD}start${NORM} position. Default is ${BOLD}0${NORM}."
  echo "${REV}-E${NORM}  --If sets -e yes, you can set ${BOLD}end${NORM} position. "
  echo "${REV}-L${NORM}  --If sets -e yes, you can set ${BOLD}length${NORM} for each line. Default is ${BOLD}100${NORM}."
  echo "${REV}-t${NORM}  --Sets the type of blast database ${BOLD}ncbidb/owndb${NORM}. Default is ${BOLD}ncbidb${NORM}."
  echo "${REV}-c${NORM}  --Sets the cpu number for option ${BOLD}c${NORM}. Default is ${BOLD}5${NORM}."
  echo -e "${REV}-h${NORM}  --Displays this help message."\\n
  echo -e "Example: ${BOLD}$SCRIPT -t ncbidb -c 5 -e TRUE -S 0 -E 600 -L 100 ${NORM}"\\n
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
cpu=5
extract=FALSE
start=0
end=''
length=100
 
while getopts :t:c:e:S:E:L:h FLAG; do
  case $FLAG in
    t)  #set option "t"
      type=$OPTARG
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
      #在这里如果你不想打印完整的帮助信息，只想显示简单的错误信息，去掉上面的两行，同时使用下面的两行。
      #echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
      #exit 2
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

if [ "$extract" == "TRUE" ]; then
	perl extract_seq.pl all.fasta $start $end $length > all.edit.fasta
	if [ "$type" == "ncbidb" ]; then
		blastn -query all.edit.fasta -db nt -outfmt 10 \
			-out "all.edit.ncbidb.csv" -evalue 1e-5 \
			-perc_identity 99 -num_alignments 10 -num_threads $cpu
	else
		blastn -query all.edit.fasta -db owndb -outfmt 10 \
			-out "all.edit.owndb.csv" -num_threads $cpu
	fi
else
	if [ "$type" == "ncbidb" ]; then
		blastn -query all.fasta -db nt -outfmt 10 \
			-out "all.ncbidb.csv" -evalue 1e-5 \
			-perc_identity 99 -num_alignments 10 -num_threads $cpu
	else
		blastn -query all.fasta -db owndb -outfmt 10 \
			-out "all.owndb.csv" -num_threads $cpu
	fi

fi



cat<<EOF >tidy_result.r
options(stringsAsFactors = FALSE)
if(file.exists("all.edit.ncbidb.csv")){
	x <- read.csv("all.edit.ncbidb.csv")
}
elif(file.exists("all.edit.owndb.csv")){
	x <- read.csv("all.edit.owndb.csv")
}
elif(file.exists("all.ncbidb.csv")){
	x <- read.csv("all.ncbidb.csv")
}
else(x <- read.csv("all.owndb.csv"))

col <- substring(x[3,1],11,nchar(x[3,1]))
colnames <- unlist(strsplit(col,", "))
x2 <- na.omit(x)
colnames(x2) <- colnames
write.csv(x2,"output.csv",row.names = FALSE)
EOF

if [ "$execute" == "TRUE" ]; then
	Rscript tidy_result.r
fi

