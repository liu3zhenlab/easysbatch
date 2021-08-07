# easysbatch
scripts and parameter documents for slurm array submissions

## Background
To submit slurm jobs in an array manner is a good way to management running job on the slurm system. However, to submit array jobs take time to organize codes and input files. The design of the script "esbatch" is to make the procedure easier. In addition, we think it is necessary to document parameters for regular analyses so those parameters can be repeatedly and consistently used. Therefore, we develop this easy-sbatch (esbatch) script and prepare parameter documents for regular analyses.  

## Usages of esbatch
```
Usage: perl ./esbatch --inlist <path_to_input_files> --cmd <command and parater> [options]
    [Options]
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
