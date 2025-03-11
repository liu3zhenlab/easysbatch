#!/bin/bash
# script to align 

star_cmd=STAR
star_db=/homes/liu3zhen/references/tomato/genome/star
protocol=star_pe
fq_dir=../31_trim
fq_pattern=R1.pair.fq$
prefix=321o_star

perl /homes/liu3zhen/scripts2/easysbatch/esbatch \
	--prefix $prefix \
	--indir $fq_dir \
	--inpattern $fq_pattern \
	--cmd $star_cmd \
	--preset $protocol \
	--in2 's/R1.pair.fq/R2.pair.fq/g' \
	--varPara $star_db \
	--opt4var "--genomeDir" \
	--time "3-00:00:00" \
	--threads 16

#--presetDB /homes/liu3zhen/scripts2/easysbatch/lib/parameters/mo_star0_se.para \

