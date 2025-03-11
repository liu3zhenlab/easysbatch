#!/bin/bash
cmd="sh /homes/liu3zhen/scripts2/easysbatch/pipelines/fungi_fq_bwa_bam.parallel.sh"
bwaDB=/homes/liu3zhen/references/fungi/magnaporthe/B71Ref2/genome/bwa/B71Ref2
indir=../1_trim
fq_suffix=".R1.pair.fq.gz"
ncpu=16
mem_gb=2
prefix=1c_aln

$cmd -v # version

perl /homes/liu3zhen/scripts2/easysbatch/esbatch \
	--indir $indir \
	--cmd $cmd \
	--varPara $bwaDB \
	--inpattern $fq_suffix \
	--threads $ncpu \
	--mem $mem_gb \
	--prefix $prefix \
	--opt4in "-f" \
	--fixPara "--cpus "$ncpu \
	--opt4var "-r" \
	--submit

