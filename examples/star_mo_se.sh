#!/bin/bash
# script to align 

star_cmd=STAR
star_db=/homes/liu3zhen/references/fungi/magnaporthe/B71Ref2/genome/star
protocol=mo_star0_se
fq_dir=../2-trim
fq_pattern=.fq$
prefix=starRun

perl /homes/liu3zhen/scripts2/easysbatch/esbatch \
	--prefix $prefix \
	--indir $fq_dir \
	--inpattern $fq_pattern \
	--cmd $star_cmd \
	--preset $protocol \
	--varPara $star_db \
	--opt4var "--genomeDir" \
	--time "3-00:00:00" \
	--threads 16

#--presetDB /homes/liu3zhen/scripts2/easysbatch/lib/parameters/mo_star0_se.para \

