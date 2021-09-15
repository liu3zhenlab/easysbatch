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
*The other way* is to input the pattern to recognize files and automatically generate the file list for gzip:
```
esbatch --indir data --inpattern fastq$ --cmd gzip --submit
```
`--submit` is specified for the slurm submission. Without `--submit`, only submission script is generated. By default, "aj" is the prefix and the submission script `aj.sbatch` is generated. The code in `aj.sbatch` can be double checked or modified and then submitted manually through the following script.
```
sbatch aj.sbatch
```
### case 2: star alignment
**task**: to align multiple single-end RNA-Seq datasets to a reference genome and count reads per gene  
We use `--preset zm_star0_se.json` to apply preset parameters stored in the parameter file [zm_star0_se.json](libs/zm_star0_se.json).  

Below shows the content of [zm_star0_se.json](libs/zm_star0_se.json).

```
{
  "_comment": "preset: zm_star0_se",
  "_comment": "version: 0.1",
  "_comment": "STAR alignment for maize RNA-Seq to a low-polymorphism reference",
  "_comment": "tested on STAR version 2.7.9a_2021-06-25",
  "mem": "36g",
  "time": "1-00:00:00",
  "threads": "8",
  "inposition": "first",
  "opt4in": "--readFilesIn", 
  "opt4out": "--outFileNamePrefix",
  "in2out": "s/[fastq$|fq$]//"
  "fixPara": {
    "_comment": "STAR parameters",
    "--alignIntronMax": "30000",
    "--outSAMattrIHstart": "0",
    "--outSAMmultNmax": "1",
    "--outSAMstrandField": "intronMotif",
    "--outFilterIntronMotifs": "RemoveNoncanonicalUnannotated",
    "--outSAMtype": "BAM SortedByCoordinate", 
    "--quantMode": "GeneCounts",
    "--outFilterMismatchNmax": "1",
    "--outFilterMismatchNoverLmax": "0.01",
    "--outFilterMatchNmin": "60",
    "--outSJfilterReads": "Unique",
    "--outFilterMultimapNmax": "1",
    "--outFilterMultimapScoreRange": "2",
    "--outSAMmapqUnique": "60",
    "--outFilterMatchNminOverLread": "0.98"
  }
}

```

**Note**:
`in2out` specifies the replacement in an input file to generate an output filename or prefix.

```
#!/bin/bash
star_cmd=STAR
star_db=<path-to-STAR-indexed database>
protocol=zm_star0_se.json
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
      --preset <str>      preset .json file to store parameters; if the path is not included, it will be searched automatically.
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
