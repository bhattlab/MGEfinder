from snakemake import shell


def _workflow(workdir, snakefile, configfile, cores, memory, unlock, rerun_incomplete, keep_going):

    if unlock == True:
        cmd = 'snakemake -s {snakefile} --configfile {configfile} --config wd={workdir} --unlock'.format(
              snakefile=snakefile, configfile=configfile, workdir=workdir
        )
    cmd = 'snakemake -s {snakefile} --config wd={workdir} memory={memory} ' \
          '--cores {cores} --configfile {configfile}'
    if rerun_incomplete:
        cmd += '--rerun-incomplete '
    if keep_going:
        cmd += '--keep-going'
    
    cmd = cmd.format(snakefile=snakefile, configfile=configfile, workdir=workdir, memory=memory, cores=cores)
    print('COMMAND:', cmd)
    shell(cmd)
