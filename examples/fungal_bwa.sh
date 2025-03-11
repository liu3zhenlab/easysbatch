perl ~/scripts2/easysbatch/esbatch \
	--indir "../10a_fq" \
	--cmd "sh /homes/liu3zhen/scripts2/easysbatch/pipelines/fungi_fq_bwa_bam.parallel.sh" \
	--varPara "/homes/liu3zhen/references/fungi/magnaporthe/B71Ref2/genome/bwa/B71Ref2" \
	--inpattern ".R1.pair.fq" \
	--threads 16 \
	--mem 1 \
	--prefix 1c_aln \
	--opt4in "-f" \
	--fixPara "-c 16" \
	--opt4var "-r" \
	--submit

