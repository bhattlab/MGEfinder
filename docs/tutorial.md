[Back to main page](../README.md)  

# Tutorial

This tutorial will analyze ten *E. faecium* isolates using the *MGEfinder* command-line tool.
In our own analysis, we found that *E. faecium* isolates had especially high numbers of novel insertions.

## Downloading the test dataset
The test dataset includes about 2.4 GB of compressed data. It can be downloaded using the following commands:

    wget https://mdurrant.s3-us-west-1.amazonaws.com/mustache/test_workdir.tar.gz
    tar -zxvf test_workdir.tar.gz
    rm test_workdir.tar.gz


This will result in a new "working directory" with the test data inside called `test_workdir`. The complete *MGEfinder* workflow 
several includes different steps that can be run using their individual commands. But all of these steps have been combined into 
a single command called `workflow`, which can be run on a properly formatted working directory.


    test_workdir/
    ├── 00.assembly/
    │   ├── ERR1036032.fna
    │   ├── ERR1036049.fna
    │   ├── ERR1036051.fna
    │   ├── ERR1078777.fna
    │   ├── ERR1078789.fna
    │   ├── ERR1195862.fna
    │   ├── ERR1541798.fna
    │   ├── ERR1541854.fna
    │   └── ERR1541922.fna
    ├── 00.bam/
    │   ├── ERR1036032.efae_GCF_900639545.bam
    │   ├── ERR1036032.efae_GCF_900639545.bam.bai
    │   ├── ERR1036049.efae_GCF_900639545.bam
    │   ├── ERR1036049.efae_GCF_900639545.bam.bai
    │   ├── ERR1036051.efae_GCF_900639545.bam
    │   ├── ERR1036051.efae_GCF_900639545.bam.bai
    │   ├── ERR1078777.efae_GCF_900639545.bam
    │   ├── ERR1078777.efae_GCF_900639545.bam.bai
    │   ├── ERR1078789.efae_GCF_900639545.bam
    │   ├── ERR1078789.efae_GCF_900639545.bam.bai
    │   ├── ERR1195862.efae_GCF_900639545.bam
    │   ├── ERR1195862.efae_GCF_900639545.bam.bai
    │   ├── ERR1541798.efae_GCF_900639545.bam
    │   ├── ERR1541798.efae_GCF_900639545.bam.bai
    │   ├── ERR1541854.efae_GCF_900639545.bam
    │   ├── ERR1541854.efae_GCF_900639545.bam.bai
    │   ├── ERR1541922.efae_GCF_900639545.bam
    │   ├── ERR1541922.efae_GCF_900639545.bam.bai
    │   ├── ERR1541932.efae_GCF_900639545.bam
    │   └── ERR1541932.efae_GCF_900639545.bam.bai
    └── 00.genome/
        └── efae_GCF_900639545.fna
        

## Generating the test dataset
We have provided these files to limit computation time to complete this tutorial. Here, we show you how these files 
were generated and how the working directory is organized so that you can analyze your own raw data in the future.

Skip to the next section if this is not of interest to you.

First, you can download the SRA run files using the [sra-tools](https://github.com/ncbi/sra-tools) command: 

    fastq-dump --gzip --split-files ERR1036032
    
These FASTQ files were then deduplicated using the [SuperDeduper](https://github.com/ibest/HTStream) command from the 
HTStream toolset with the following command:
    
    hts_SuperDeduper -g -1 ERR1036032_1.fastq.gz -2 ERR1036032_2.fastq.gz -p ERR1036032.nodup
 
And then the FASTQ files were trimmed using the 
[Trime Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) package and the command:

    trim_galore --fastqc --paired ERR1036032.nodup_R1.fastq.gz ERR1036032.nodup_R2.fastq.gz
    
The reference genome was indexed for alignment with [BWA MEM](http://bio-bwa.sourceforge.net/) and the alignment was 
performed using the commands:

    bwa index efae_GCF_900639545.fna
    bwa mem efae_GCF_900639545.fna ERR1036032.nodup_R1_val_1.fq.gz ERR1036032.nodup_R2_val_2.fq.gz > ERR1036032.efae_GCF_900639545.sam

There is a built-in command in *mgefinder*, called `mgefinder formatbam` which takes this SAM file and prepares it for 
analysis by *mgefinder*:

    mgefinder formatbam ERR1036032.efae_GCF_900639545.sam ERR1036032.efae_GCF_900639545.bam

This step generates the `ERR1036032.efae_GCF_900639545.bam` and `ERR1036032.efae_GCF_900639545.bam.bai` files that we have already 
provided to you for this tutorial.

We are working to remove this step so that you can use BAM files without this, but for the time being it is best to use 
this step when possible.

Finally, to generate the assembled contig file `ERR1036032.fna`, you can use the 
[SPAdes](http://cab.spbu.ru/software/spades/) command:

    spades.py -1 ERR1036032.nodup_R1_val_1.fq.gz -2 ERR1036032.nodup_R2_val_2.fq.gz -o ERR1036032
    
The assembled contigs are available in the `ERR1036032/contigs.fasta` file.

We repeated these steps for all ten isolates. We then organized these files into a working directory. It is imperative
when you make your own working directories, you use this directory structure and file names:

    workdir/
        ├── 00.assembly/
        │   ├── <sample1>.fna
        │   ├── <sample1>.fna
        │   └── <sample1>.fna
        ├── 00.bam/
        │   ├── <sample1>.<genome>.bam
        │   ├── <sample1>.<genome>.bam.bai
        │   ├── <sample2>.<genome>.bam
        │   ├── <sample2>.<genome>.bam.bai
        │   ├── <sample3>.<genome>.bam
        │   └── <sample3>.<genome>.bam.bai
        └── 00.genome/
            └── <genome>.fna

Be sure to keep the file naming consistent with what is shown above, including directory names and suffixes. 
`<genome>` refers to the name of the genome to which the sample was aligned. Make sure that the only periods in your file names are the ones shown above. Files in the directory can be symlinks to the actual files. 

## Analyzing the test dataset using *mgefinder workflow*

This working directory can be analyzed using a single command:

    mgefinder workflow test_workdir/
    
This will begin running a snakemake workflow that will execute all steps in the *MGEfinder* pipeline on the
contents of the working directory.

We tested this workflow on a Google Cloud Virtual machine with 4 vCPUs and 26 GB memory. By default, `mgefinder
workflow` uses one core. This default setting takes 7 min 37 seconds to analyze these 10 isolates. The number of
cores can be increased:

    mgefinder workflow --cores 4 test_workdir/

Which takes 3 min 3 seconds on the same machine. Certain steps can be done in parallel, but there are bottlenecks.
If you have hundreds of isolates and available memory, we recommend increasing the memory using `--memory` from
the default of `--memory 16000` to `--memory 50000`, where the parameter is the amount of available memory in
megabases.

This workflow generates three more directories in the `test_workdir` working directory:

    01.mgefinder
    02.database
    03.results
    
The `01.mgefinder` directory contains intermediate *MGEfinder* files that may be of interest to individuals with
more advanced knowledge of how the tool works. This includes the output of the `find`, `pair`, and `inferseq`
commands.

The `02.database` directory contains a dynamically-constructed database of identified elements that were 
identified when running *MGEfinder*. These are necessary for the `inferseq-database` commands. This file
may be of value to the user, depending on the use case.

Most of the files of value can be found in the `03.results` directory. In the case of this test dataset,
the files that were generated are named:

    workdir/
        └── efae_GCF_900639545/
            ├── 01.clusterseq.efae_GCF_900639545.tsv
            ├── 02.genotype.efae_GCF_900639545.tsv
            ├── 03.summarize.efae_GCF_900639545.clusters.tsv
            ├── 03.summarize.efae_GCF_900639545.groups.tsv
            ├── 04.makefasta.efae_GCF_900639545.all_seqs.fna
            └── 04.makefasta.efae_GCF_900639545.repr_seqs.fna

The file `01.clusterseq.efae_GCF_900639545.tsv` a complete list of all inferred sequences
for all ten files. This file includes the `sample`, `pair_id` (refers to results of the `pair` output files),
`method` (sequence inference method, such as `inferred_assembly_with_full_context`), `loc` (the location
of where the sequence was found, either in the reference genome, assembly, or dynamically-constructed database),
the `inferred_seq_length`, the `seqid` (an identifier for each unique sequence), `cluster` (an identifier for the 
sequence cluster), `group` (an identifier for the cluster group), and the `inferred_seq` (nucleotide sequence of
the element).

The file `02.genotype.efae_GCF_900639545.tsv` includes insertion genotypes for all of the isolates. This
includes the `sample`, the position of the insertion within the reference genome, the `seqid`, `cluster`, and
`group` that indicates the identity of the inserted sequence. The final column indicates the `conf` (confidence) 
level of the insertion. This refers primarily to how the idnetity of the sequence was inferred. The highest confidence
inferred sequence is `IAwFC`, which indicates that the insertion was found in the expected location in the
assembly. The degree of confidence that you should accept depends largely on the use case. If, for example,
you are performing a re-sequencing experiment, it may be sufficient to use the `IDB` column, which includes
sequences inferred from the reference and dynamically-constructed database.

The files `03.summarize.efae_GCF_900639545.clusters.tsv` and `03.summarize.efae_GCF_900639545.groups.tsv` provide
summary statistics for each of the sequence clusters and groups, respectively. This can be used to perform further
QC on the identified elements, and to stratify elements by their transposability. Some important columns include 
`num_unique_sites_all`, which refers to the number of unique sites in the reference genome where sequence cluster was
found at least once in at least one isolate. and `num_unique_sites_unambig`, which refers to the number of unique sites
where the cluster is the unambiguously assigned element (often two clusters will be assigned to a single site when the
exact element cluster cannot be identified). See the manual for more details on these output files.

The files `04.makefasta.efae_GCF_900639545.all_seqs.fna` and `04.makefasta.efae_GCF_900639545.repr_seqs.fna` include FASTA
files of all identified elements found in the `genotype` and `summarize` files. The `*.all_seqs.fna` includes all unique
sequences identified (even if they differ by a single base pair), and `*.repr_seqs.fna` includes the representative
sequence for each element cluster.

## Next steps
This tutorial produced a list of candidate integrative mobile genetic elements for these ten E. faecium isolates,
and genotyped insertions with respect to an E. faecium reference genome. The next steps taken will depend largely on 
the type of analysis that you want to perform. If you want to search these candidate insertions for antibiotic 
resistance genes, you can upload the `04.makefasta.efae_GCF_900639545.all_seqs.fna` FASTA file to 
to [ResFinder](https://cge.cbs.dtu.dk/services/ResFinder/). If you wish to identify transposases, we recommend
the HMMs provided by [ISEScan](https://www.ncbi.nlm.nih.gov/pubmed/29077810). If you wish to identify prophage
elements, you can upload the `04.makefasta.efae_GCF_900639545.repr_seqs.fna` file to [PHASTER](http://phaster.ca/).

Please refer to the [manuscript](https://doi.org/10.1101/527788) for more ideas on how you can analyze your data.

## Other considerations
We hope we have made the limitations of this approach clear, and we urge you to carefully consider these limitations 
when analyzing your own data.

Some additional details to keep in mind:
* *MGEfinder* is not designed to identify very small indels. By default it identifies insertions as small as 70 base
pairs, and under certain conditions it can reliably identify insertions as short as 30 base pairs. By default, the upper
limit on the insertion size is 200 kbp.
* It is possible that some identified elements are in fact assembly errors, and we urge you to consider this when making
conclusions about specific insertions.
* This workflow is not configured to handle resequencing experiments by default, please consider refer to the manual
if you are analyzing resequencing data.

Good luck, and if you have any questions please submit an issue [here](https://github.com/bhattlab/MGEfinder/issues).

[NEXT: Detailed User Manual](manual.md)
