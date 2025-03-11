#!/bin/bash
# by Sanzhen Liu
# 6/11/2022

version=0.01;

default_f1suffix='.R1.pair.fq';
default_f2suffix='.R2.pair.fq';
default_cpus=4;
default_parser="/homes/liu3zhen/scripts/sam/samparser.bwa.pl";
default_parser_para="-e 60 -m 2 100 --tail 60 100 --gap 5 --insert 100 750";

RED='\033[0;31m'
NC='\033[0m' # No Color
usage() {
	echo -e "${RED}Prerequirement${NC}: bwa, samtools"
	echo -e "${RED}Usage${NC}: $0 -f <fastq> -r <ref> [other options]" >&2
	echo "   -f: fastq file; required" >&2
	echo "   -r: bwa indexed database; required" >&2
	echo "   -1: suffix of first pair of fastq ($default_f1suffix)" >&2
	echo "   -2: suffix of second pair of fastq ($default_f2suffix)" >&2
	echo "   -c: number of cpus ($default_cpus)" >&2
	echo "   -m: modules to load" >&2
	echo "   -p: bwa alignment parser ($default_parser)" >&2;
	echo "   -a: alignment parser criteria ($default_parser_para)" >&2;
	echo "   -h: help information" >&2
}

while getopts ":f:r:1:2:c:p:a:m:vh" opt; do
case $opt in
	f) fq=$OPTARG;;
	r) ref=$OPTARG;;
	1) f1suffix=$OPTARG;;
	2) f2suffix=$OPTARG;;
	c) cpus=$OPTARG;;
	p) parser=$OPTARG;;
	a) para=$OPTARG;;
	m) modules+=($OPTARG);;
	v) echo $version; exit;;
	h) usage; exit;;
\?) echo "Invalid options: -$OPTARG." >&2; exit;;
:) echo "Option -$OPTARG requires an argument" >&2; exit;;
esac
done

###############################################
# modules
###############################################
cmd_check () {
	input_cmd=$1
	which $input_cmd &>/dev/null
	if [ $? -eq 1 ]; then
		echo -e "${RED}ERROR${NC}: $input_cmd not available." >&2
		exit;
	fi
}

file_check () {
	input_file=$1
	if [ ! -f $input_file ]; then
		echo -e "${RED}ERROR${NC}: $input_file does not exit." >&2
		exit;
	fi
}

###############################################
### check required parameters
###############################################
if [ -z $fq ] || [ -z $ref ]; then
	echo -e "${RED}ERROR${NC}: Required parameters: -f; -r." >&2
	usage;
	exit;
fi

file_check $fq  # check input data

if [ ! -f ${ref}.bwt ]; then
	echo -e "${RED}ERROR${NC}: BWA databae $ref does not exit." >&2
	exit;
fi

if [ -z $f1suffix ]; then
	f1suffix=$default_f1suffix;
fi

if [ -z $f2suffix ]; then
	f2suffix=$default_f2suffix;
fi

if [ -z $cpus ]; then
	cpus=$default_cpus;
fi

if [ -z $parser ]; then
	parser=$default_parser;
fi

if [ -z $para ]; then
	para=$default_parser_para;
fi

for module in "${modules[@]}"; do
	module load $module;
	if [ $? -eq 0 ]; then 
		echo "module "$module" loaded";
	fi
done

###############################################
# input fastq
###############################################
fq2=$(echo $fq | sed "s/$f1suffix/$f2suffix/g");
echo -e "${RED}input data:${NC}" >&2
echo "    read 1: "$fq  >&2
echo "    read 2: "$fq2  >&2

# if fq are gzip files
fq_extension="${fq##*.}"
if [ $fq_extension == "gz" ]; then
	new_fq=`echo $fq | sed 's/.*\///g' | sed 's/.gz//g' | sed 's/^/./g'`;
	new_fq2=`echo $fq2 | sed 's/.*\///g' | sed 's/.gz//g' | sed 's/^/./g'`;
	gunzip -c $fq > $new_fq
	gunzip -c $fq2 > $new_fq2
	fq=$new_fq; fq2=$new_fq2
	new_f1suffix=`echo $fq1suffix | sed 's/.gz//g'`;
	#new_f2suffix=`echo $fq2suffix | sed 's/.gz//g'`
	fq1suffix=$new_fq1suffix;
	#fq2suffix=$new_fq2suffix;
fi



###############################################
# check requirements:
###############################################
cmd_check bwa;
cmd_check samtools;
file_check $parser;

###############################################
# run and output:
###############################################
out=$(echo $fq | sed 's/.*\///g' | sed 's/^\.//g' | sed "s/$f1suffix//g");
sam=${out}.sam

echo -e "${RED}Output prefix:${NC}" >&2
echo "    $out"  >&2

echo -e "${red}reference db:${nc}" >&2
echo "    $ref" >&2

### aln
group_info="@RG\tID:${out}\tSM:${out}"
bwa mem -t $cpus -R "$group_info" $ref $fq $fq2 1>${out}.sam 2>${out}.aln.log
if [ $? -eq 1 ]; then
	echo -e "${RED}ERROR${NC}: BWA alignment failed." >&2
	rm $new_fq; rm $new_fq2;
	exit;
fi

# cleanup
if [ $fq_extension == "gz" ]; then
	rm $new_fq; rm $new_fq2
fi


### filter
perl $parser -i ${out}.sam \
	$para \
	1> ${out}.parse.sam  \
	2>${out}.parse.log

if [ $? -eq 1 ]; then
	echo -e "${RED}ERROR${NC}: Alignment parsing failed." >&2
	exit;
fi

### convert SAM to BAM:
samtools view -@ $cpus -bS ${out}.parse.sam | samtools sort -@ $cpus -o ${out}.bam
if [ $? -eq 1 ]; then
	echo -e "${RED}ERROR${NC}: SAMtools sam2bam conversion failed." >&2
	exit;
fi

### Index sorted BAM:
samtools index -@ $cpus ${out}.bam
if [ $? -eq 1 ]; then
	echo -e "${RED}ERROR${NC}: SAMtools index failed." >&2
	exit;
fi


