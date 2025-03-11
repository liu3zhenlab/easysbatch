#!/bin/bash
# by Sanzhen Liu
# 9/12/2021

version=0.01;

default_f1suffix='.R1.pair.fq';
default_f2suffix='.R2.pair.fq';
default_cpus=4;
script_dir=`echo $0 | sed 's/[^\/]*\/[^\/]*$//g'` # remove file and the direct subdirectory
default_parser=$script_dir"/utils/samparser.bwa.pl";
default_aggregator=$script_dir"/utils/alignment.log.aggregate.pl";
parser_para="-e 60 -m 5 100 --tail 5 100 --gap 10 --insert 100 600";
cleanup=0

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
	echo "   -p: parser script for filtering bwa alignments ($default_parser)" >&2;
	echo "   -a: aggregator script for merging parser log ($default_aggregator)" >&2;
	echo "   -d: deletion of SAM output if specified" >&2;
	echo "   -v: version information" >&2;
	echo "   -h: help information" >&2
}

while getopts ":f:r:1:2:c:p:m:vdh" opt; do
case $opt in
	f) fq=$OPTARG;;
	r) ref=$OPTARG;;
	1) f1suffix=$OPTARG;;
	2) f2suffix=$OPTARG;;
	c) cpus=$OPTARG;;
	p) parser=$OPTARG;;
	m) modules+=($OPTARG);;
	a) aggregator=$OPTARG;;
	d) cleanup=1;;
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

if [ -z $aggregator ]; then
	aggregator=$default_aggregator;
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
fq_extension="${fq##*.}"  # suffix
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
bwaout=${out}.sam
group_info="@RG\tID:${out}\tSM:${out}"
echo $group_info
bwa mem -t $cpus -R "$group_info" $ref $fq $fq2 1>${out}.sam 2>${out}.aln.log
if [ $? -eq 1 ]; then
	echo -e "${RED}ERROR${NC}: BWA alignment failed." >&2
	if [ $fq_extension == "gz" ]; then
		rm $new_fq; rm $new_fq2;
	fi
	exit;
fi

# cleanup
if [ $fq_extension == "gz" ]; then
	rm $new_fq; rm $new_fq2
fi

# split alignments
naln=`wc -l ${bwaout} | sed 's/ .*//g'`;
nlines=`expr $naln / $cpus`;
split -d -l $nlines ${bwaout} ${bwaout}

# tmp
tmp_split_parser=${bwaout}.split.parse.tmp.sh
echo "#!/bin/bash" > $tmp_split_parser
for sam in ${bwaout}[0-9]*[0-9]; do
	echo -e perl ${parser} -i ${sam} $parser_para 1\>${sam}.parse 2\>${sam}.parse.log >> $tmp_split_parser
done

### filter
xargs --arg-file=$tmp_split_parser --max-proc=$cpus --replace --verbose /bin/sh -c "{}";

if [ $? -eq 1 ]; then
	echo -e "${RED}ERROR${NC}: Alignment parsing failed." >&2
	exit;
fi

rm $tmp_split_parser

### merge
cat ${bwaout}[0-9]*.parse > ${out}.parse.sam
perl $aggregator ${bwaout} ${bwaout}[0-9]*.parse.log > ${out}.parse.log
rm ${bwaout}[0-9]*

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
else
	if [ $cleanup -eq 1 ]; then
		rm ${out}.sam
		rm ${out}.parse.sam
	fi
fi

