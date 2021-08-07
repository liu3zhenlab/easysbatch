# easysbatch
scripts and parameter documents for slurm array submissions

## Background
To submit slurm jobs in an array manner is a good way to management running job on the slurm system. However, to submit array jobs take time to organize codes and input files. The design of the script "esbatch" is to make the procedure easier. In addition, we think it is necessary to document parameters for regular analyses so those parameters can be repeatedly and consistently used. Therefore, we develop this easy-sbatch (esbatch) script and prepare parameter documents for regular analyses.  

## Exampels
### case 1: gzip files
**task**: is to gzip \*fastq files in the directory of data  
*One way* is to generate a list of input files and then run esbatch:
```
ls -1 data/*fastq > fqlist
esbatch --inlist --cmd gzip --submit
```
*The other way* is to generate the file list and run esbatch at the same time:
```
esbatch --indir data --inpattern fastq$ --cmd gzip --submit
```
`--submit` is specified for the slurm submission. Without `--submit`, only submission script is generated. By default, "aj" is the prefix and the submission script `aj.sbatch` is generated. The code in `aj.sbatch` can be double checked or modified and then submitted manually through the following script.
```
sbatch aj.sbatch
```
### case 2: star alignment
**task**: to align multiple single-end RNA-Seq datasets to a reference genome and count reads per gene  
We use `--preset mo_star0_se` to apply preset parameters stored in the parameter file [mo_star0_se.para](lib/mo_star0_se.para). The file can be supplied by using `--presetDB <path>/mo_star0_se.para`. If `--presetDB` is not specified, the program will search `mo_star0_se.para` in the `lib` subdirectory  under the directory containing the main script `esbatch`.

Below shows the content of [mo_star0_se.para](lib/mo_star0_se.para), in which information starts form preset and ends with "===".
```
preset: mo_star0_se
mem: 16g
time: 1-00:00:00
threads: 4
inposition: first
opt4in: --readFilesIn 
opt4out: --outFileNamePrefix
in2out: s/[fastq$|fq$]//
fixPara: --alignIntronMax 3000 --outSAMattrIHstart 0 --outSAMmultNmax 1 ...
note: STAR alignment for fungal RNA-Seq of fungal strains identical or highly similar to the reference strain
===
```
**Note**:
`in2out` specifies the replacement in an input file to generate an output filename or prefix.

```
#!/bin/bash
star_cmd=STAR
star_db=<path-to-STAR-indexed database>
protocol=mo_star0_se
fq_dir=<path-to-fastq>
fq_pattern=.fq$
esbatch \
    --indir $fq_dir \
    --inpattern $fq_pattern \
    --cmd $star_cmd \
    --preset $protocol \
    --varPara $star_db \
    --opt4var "--genomeDir" \
    --submit
```
In addition, the parameter `--varPara` allows to add additional parameters that are not with fixed values. In this case, the indexed reference is add here. `--opt4var` is the option for this paramter. In the command, `--genomeDir <path-to-STAR-indexed database>` will be shown.

**Notes**: For paired-end (PE) read alignments, `--in2` can be used to add the 2nd read file. `--in2` could be a file or the replacement code (e.g., s/R1.pair/R2.pair/) to generate the 2nd input by making a change from the 1st input.  

## Full usages of esbatch
```
Usage: esbatch --inlist <path_to_input_files> --cmd <command and parater> [options]
      --inlist* <file>    a file lists input files with one file per row, including path
                          required if --inpattern is not input
      --indir <str>       path to input data if --inlist is not provided (.)
      --inpattern* <str>  regular expression pattern to filter files in --indir
                          required if --inlist is not input
      --cmd* <str>        command (e.g., gzip); required
      --presetDB <file>   file to store preset parameters for one or multiple analyses
                          see below for an example for one analysis: mo_star0_se
      --preset <str>      name of a certain analysis in the --presetDB file
                          for the example case, mo_star0_se should be used
      --mem <num>         Gb memory requested (8g)
      --time <time>       running time requested; (1-00:00:00)
      --threads <num>     number of threads (4)
      --prefix <str>      prefix for output (aj)
      --inposition <str>  the position of input parameters relative to other parameters in the command line
                          default is at the second to the last
                          only "first" and "last" are the options.
      --opt4in <str>      the option for input used in --cmd program (e.g., -i), none by default 
      --in2 <str>         the 2nd input; assuming the same directory as --indir or files in --inlist
                          if in s///g format, the replacement will be applied on 1st input to generate 2nd input
      --opt4in2 <str>     the option for 2nd input (--in2), none by default
      --outdir <str>      path of output data; to be created if not existed. (.) 
      --opt4out <str>     the option for output in --cmd program (e.g., -o or >), none by default
      --fixPara <str>     string for parameters with fixed values for --cmd, none by default
      --varPara <str>     variable parameter; multiple invocations allowed
      --opt4var <str>     the option for --varPara; equal number invocations as --varPara
                          if multiple invocations exist, first --varPara is paired with first --opt4var, and so on.
      --in2out <str>      the replacement to change input name to output name via the program of "sed"
                          (e.g., s/fq$/sam/)
      --submit            submit slurm job if specified
      --version           version information
      --help              help information
```
## BUG REPORT
Please report any bugs or suggestion on github or by email to Sanzhen Liu (liu3zhen@gmail.com).

## LICENSE
The script and parameter documents are distributed under MIT licence.
