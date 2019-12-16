import sys
from snakemake.shell import shell

class DependencyChecker:

    def __init__(self, tool, version=None):
        self.tool = tool
        self.version = version

    def check(self, cmd='{tool} --version 2>&1', extract_version= lambda x: x):
        self.check_exists()
        output = shell(cmd.format(tool=self.tool), read=True).decode('utf-8').strip()
        version = extract_version(output)
        print("Current version of {tool}: {version}".format(tool=self.tool, version=version))
        print("Expected version of {tool}: {version}".format(tool=self.tool, version=self.version))
        if version != self.version:
            print("WARNING: Current version of {tool} does not match expectaction. Proceed with caution".format(tool=self.tool))


    def check_exists(self):
        try:
            path = shell('which {tool}'.format(tool=self.tool), read=True)
        except:
            print("MISSING DEPENDENCY: {tool} is not installed. Please install before continuing.".format(tool=self.tool))
            sys.exit()
            


def check_dependencies():

    print("#### CHECKING DEPENDENCIES ####")
    snakemake_checker = DependencyChecker('snakemake', '3.13.3')
    snakemake_checker.check()

    emboss_checker = DependencyChecker('einverted', 'EMBOSS:6.6.0.0')
    emboss_checker.check()

    bowtie2_checker = DependencyChecker('bowtie2', '2.3.5')
    bowtie2_checker.check(extract_version=lambda x: x.split()[2])

    samtools_checker = DependencyChecker('samtools', '1.9')
    samtools_checker.check(extract_version=lambda x: x.split()[1])

    cdhit_checker = DependencyChecker('cd-hit', '4.8.1')
    cdhit_checker.check(cmd='{tool} || exit 0', extract_version=lambda x: x.split()[3])
    print("###############################")
