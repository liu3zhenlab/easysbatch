#!/bin/bash
cmd="sh /homes/liu3zhen/scripts2/easysbatch/pipelines/sam2bam.sh"
indir=.
pattern=".sam$"
ncpu=16
mem_gb=4
module="SAMtools"

$cmd -v # version

perl /homes/liu3zhen/scripts2/easysbatch/esbatch \
	--indir $indir \
	--cmd "$cmd" \
	--inpattern $pattern \
	--threads $ncpu \
	--mem $mem_gb \
	--opt4in "-s" \
	--fixPara "-m "$module" -c "$ncpu \
	--submit

