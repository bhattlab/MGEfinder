from os.path import dirname, join
from snakemake import snakemake


def _workflow(workdir, snakefile, config, cores, memory, unlock, rerun_incomplete, keep_going):

    snakemake(snakefile, configfile=config, config={'wd': workdir, 'memory': memory}, cores=cores, unlock=unlock,
              force_incomplete=rerun_incomplete, keepgoing=keep_going)