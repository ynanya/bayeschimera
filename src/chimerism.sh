#! /bin/bash		

#$ -S /bin/bash		
#$ -o log/		
#$ -e log/	
#$ -cwd


usage_exit() {
	echo   "" 1>&2
	echo   "" 1>&2
        echo "Usage: chimerism.sh [OPTIONS] PATH_TO_RECIPIENT_BAM_FILE  PATH_TO_CHIMERA_BAM_FILE " 1>&2
	echo "" 1>&2
	echo "OPTIONS: " 1>&2
	echo "[-d dir] : path to SNPs list file."  1>&2
	echo "[-h]: show this message. " 1>&2
	echo "" 1>&2	
        exit 1
}


	
while getopts d:h OPT
do
    case $OPT in
        d)  OPT_FLAG_d=1; OPT_VALUE_D=$OPTARG
            ;;
        h)  usage_exit
            ;;
        \?) usage_exit
            ;;
    esac
done



shift $(($OPTIND - 1))


## path to SNPs file 
SNPs_path=data/snps_list.txt
if [[ -n "${OPT_FLAG_d+UNDEF}" ]];then
  	echo "opt=-d, SNPs_path=${OPT_VALUE_D}"
fi

if [ ! -e $SNPs_path ]; then
  	echo   "$SNPs_path does not exist." 1>&2	
        exit 1
fi





arg1=$1 # path to bam file of purely recipient genome
arg1=`echo ${arg1} | sed -e "s/[\r\n]\+//g"`

arg2=$2 # path to bam file of post transplant DNA
arg2=`echo ${arg2} | sed -e "s/[\r\n]\+//g"`


filename1=`basename $arg1 | sed 's/.markdup.bam$//'`
filename2=`basename $arg2 | sed 's/.markdup.bam$//'`

echo $filename1
echo $filename2


if [ ! -e $arg1 ]; then
  	echo $filename1 does not exist
	exit 1
fi
if [ ! -e $arg2 ]; then
  	echo $filename2 does not exist
	exit 1
fi

if [ -d tmp/${filename1}_${filename2} ]; then
  	rm -r  tmp/${filename1}_${filename2}
fi

mkdir tmp/${filename1}_${filename2}
work_dir=tmp/${filename1}_${filename2}



export R=/usr/local/package/r/current3.6/bin/R


#  execute pileup 
echo  "R --vanilla --slave --args  $filename1 $arg1 $work_dir $SNPs_path< src/pileup.R"
R --vanilla --slave --args  $filename1 $arg1 $work_dir $SNPs_path< src/pileup.R 
echo  "R --vanilla --slave --args  $filename2 $arg2 $work_dir $SNPs_path<  src/pileup.R" 
R --vanilla --slave --args  $filename2 $arg2 $work_dir $SNPs_path<  src/pileup.R 


$R --vanilla --slave --args tmp/${filename1}_${filename2}/${filename1}calls2.txt  tmp/${filename1}_${filename2}/${filename2}calls2.txt  < src/stan.R 

