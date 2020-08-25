import sys
from snakemake import shell
try:
    from snakemake import shell
except ImportError:
    print("MISSING DEPENDENCY: Snakemake is not installed. You can install with the command \"conda install -c bioconda -c conda-forge snakemake\" or \"pip3 install snakemake\"")
    sys.exit()



def _workflow(workdir, snakefile, configfile, cores, memory, unlock, rerun_incomplete, keep_going):

    cmd = 'snakemake -s {snakefile} --config wd={workdir} memory={memory} ' \
          '--cores {cores} --configfile {configfile} '
    
    if rerun_incomplete:
        cmd += '--rerun-incomplete '
    if keep_going:
        cmd += '--keep-going '
    if unlock:
        cmd += '--unlock '
    
    cmd = cmd.format(snakefile=snakefile, configfile=configfile, workdir=workdir, memory=memory, cores=cores)
    print('COMMAND:', cmd)
    shell(cmd)
