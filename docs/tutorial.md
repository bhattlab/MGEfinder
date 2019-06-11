[Back to main page](../README.md)  

# Tutorial

This tutorial will analyze a single *E. coli* isolate using the *mustache* command-line tool.

We have selected the SRA sample with the Run ID SRR3180793. This is an E. coli isolate sequenced as part of the CDC's
PulseNet program. In our own analysis, we found that this isolate had an especially high number of insertions.

## Downloading the test dataset
The test dataset includes about 460 MB of compressed data. It can be downloaded using the following commands:

    wget https://s3-us-west-1.amazonaws.com/mdurrant/mustache/mustache_test_dataset.tar.gz
    tar -zxvf mustache_test_dataset.tar.gz
    rm mustache_test_dataset.tar.gz
    cd mustache_test_dataset

In the test dataset are 5 files that will be used in this tutorial:
* `CFT073.fna` - the E. coli reference genome that we will use for this analysis.
* `query_database.fna` - A database of genetic elements that we can use as a query.
* `SRR3180793.CFT073.bam` - The SRR3180793 reads aligned to the `CFT073.fna` genome.
* `SRR3180793.CFT073.bam.bai` - The BAM index produced using the `samtools index` command. 
* `SRR3180793.contigs.fna` - The contigs generated from genomic assembly of the SRR3180793 FASTQ files.

## Generating the test dataset
We have provided these files to limit computation time to complete this tutorial. Here, we show you how these files 
were generated so that you can analyze your own raw data in the future.

Skip to the next section if this is not of interest to you.

First, you can download the SRA run file using the [sra-tools](https://github.com/ncbi/sra-tools) command: 

    fastq-dump --gzip --split-files SRR3180793
    
These FASTQ files were then deduplicated using the [SuperDeduper](https://github.com/ibest/HTStream) command from the 
HTStream toolset with the following command:
    
    hts_SuperDeduper -g -1 SRR3180793_1.fastq.gz -2 SRR3180793_2.fastq.gz -p SRR3180793.nodup
 
And then the FASTQ files were trimmed using the 
[Trime Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) package and the command:

    trim_galore --fastqc --paired SRR3180793.nodup_R1.fastq.gz SRR3180793.nodup_R2.fastq.gz
    
The reference genome was indexed for alignment with [BWA MEM](http://bio-bwa.sourceforge.net/) and the alignment was 
performed using the commands:

    bwa index CFT073.fna
    bwa mem CFT073.fna SRR3180793.nodup_R1_val_1.fq.gz SRR3180793.nodup_R2_val_2.fq.gz > SRR3180793.CFT073.sam

There is a built-in command in *mustache*, called `mustache formatbam` which takes this SAM file and prepares it for 
analysis by *mustache*:

    mustache formatbam SRR3180793.CFT073.sam SRR3180793.CFT073.bam

This step generates the `SRR3180793.CFT073.bam` and `SRR3180793.CFT073.bam.bai` files that we have already 
provided to you for this tutorial.

Finally, to generate the assembled contig file `SRR3180793.contigs.fna`, you can use the 
[SPAdes](http://cab.spbu.ru/software/spades/) command:

    spades.py -1 SRR3180793.nodup_R1_val_1.fq.gz -2 SRR3180793.nodup_R2_val_2.fq.gz -o SRR3180793
    
The assembled contigs are available in the `SRR3180793/contigs.fasta` file.

The query database `query_database.fna` was generated from our analysis of thousands of randomly downloaded *E. coli* 
isolates available in the SRA database. A description of how this file was generated is available in our 
[preprint](https://www.biorxiv.org/content/10.1101/527788v1).


## Analyzing the test dataset using *mustache*

We present a simple tutorial designed to explain key parts of the *mustache* workflow.

### `mustache findflanks`
Once a properly formatted BAM file is available, you can run the `findflanks` command in *mustache*. This command is
described in greater detail [here](manual.md). You can run this command on the test dataset as:
  
    mustache findflanks SRR3180793.CFT073.bam

This will generate a file called `mustache.findflanks.tsv`. Alternatively, you can specify the filename using:

    mustache findflanks SRR3180793.CFT073.bam -o SRR3180793.findflanks.tsv
    
This step should take about 3.5 minutes to complete for this sample.
The output file includes the consensus flanks generated from analyzing the clipped-ends of reads aligned to the 
CFT073 E. coli genome, and information about the reads that support the existence of each flank. For more information 
about this file, please read the [manual](manual.md).

### `mustache pairflanks`
The next step is to pair these candidate flanks together, giving us a candidate pair that represents the 5' and 3' 
flanks of the insertion. This is performed using the `pairflanks` command in `mustache` as follows:

    mustache pairflanks mustache.findflanks.tsv SRR3180793.CFT073.bam CFT073.fna

This step should take about 1 minute to complete. As you can see, this requires as input the `mustache.findflanks.tsv` 
file (the output of the `findflanks` command), the file `SRR3180793.CFT073.bam` (the alignment of the reads to the 
reference genome), and the reference genome used `CFT073.fna`.

The output of this command is a file called `mustache.pairflanks.tsv`. This includes all of the candidate flank pairs 
generated from analyzing the SRR3180793 sample, details about the reads supporting the insertion, and other important
information that may be of interest to you in your own analysis. Please see details about this output file in the 
[manual](manual.md).

### `mustache inferseq-reference`
This next step is used to infer the identity of the candidate insertions that were identified in the `pairflanks` step.
If you suspect that the inserted element already exists within the reference genome, as is often the case in a 
adaptive laboratory evolution experiments, then be sure to use this command.

While *mustache* will automatically index the reference genome if it is not already indexed, we recommend that you index
the reference genome beforehand with the command:

    bowtie2-build -o 0 -q CFT073.fna CFT073.fna

You can then run the `inferseq-reference` command as:

    mustache inferseq-reference mustache.pairflanks.tsv CFT073.fna
    
This step should take about 5 seconds to complete. This command requires as input the file generated by the `pairflanks`
command and the reference genome. The output of this command is in the `mustache.inferseq_reference.tsv` file. This file
includes the elements that were inferred from the reference genome and are a possible source of the insertion. These 
inferred sequences can be mapped back to the candidate insertions in the `mustache.pairflanks.tsv` file using the 
`pair_id` column. Please see details about this output file in the [manual](manual.md).

### `mustache inferseq-assembly`
This next step is used to infer the identity of the candidate insertions that were identified in the previous step.
This step uses the assembled contigs of the isolate itself to infer the identity of the inserted elements. This is 
useful when you suspect that the reference genome and the isolate are more distantly related, or that novel genetic material may
be the source of the insertions.

While *mustache* will automatically index the assembled contigs if it is not already indexed, we recommend that you index
the assembled contigs beforehand with the command:

    bowtie2-build -o 0 -q SRR3180793.contigs.fna SRR3180793.contigs.fna

You can then run the `inferseq-assembly` command as:

    mustache inferseq-assembly mustache.pairflanks.tsv SRR3180793.CFT073.bam SRR3180793.contigs.fna CFT073.fna
    
This step should take about 30 seconds to complete. This command requires as input the file generated by the `pairflanks`
command, the file `SRR3180793.CFT073.bam` (the alignment of the reads to the reference genome), the assembled contigs, 
and the reference genome. The output of this command is in the `mustache.inferseq_assembly.tsv` file. This file
includes the elements that were inferred from the assembled contigs and are a possible source of the insertion. These 
inferred sequences can be mapped back to the candidate insertions in the `mustache.pairflanks.tsv` file using the 
`pair_id` column. Please see details about this output file in the [manual](manual.md).

### `mustache inferseq-overlap`
This next step is used to infer the identity of the candidate insertions that were identified in the previous step.
This step attempts to overlap the two candidate flanks to determine the identity of the full inserted element. This is 
useful for inferring the identity of insertions that are relatively small, depending heavily on the length of the reads
in the library.

You can run the `inferseq-overlap` command as:

    mustache inferseq-overlap mustache.pairflanks.tsv
    
This step should take about 5 seconds to complete. This command only requires as input the file generated by the 
`pairflanks` command. The output of this command is in the `mustache.inferseq_overlap.tsv` file. Again, this approach is
limited by the length of the reads. These inferred sequences can be mapped back to the candidate insertions in the 
`mustache.pairflanks.tsv` file using the `pair_id` column. Please see details about this output file in the 
[manual](manual.md).


### `mustache inferseq-database`
This next step is used to infer the identity of the candidate insertions from a query database of previously identified
insertions. This is useful when the insertions that you are looking for or already known, or you want to look for
insertions that were identified when analyzing a different isolate. The database used in tutorial 
`query_database.fna` was generated by running `mustache` on thousands of E. coli isolates, and compiling the identified
insertions into a single fasta file. This is called a 'dynamically-constructed database' because it was generated from
from using *mustache* itself on thousands of isolates, not from a curated online database.

While *mustache* will automatically index the query database if it is not already indexed, we recommend that you index
the query database beforehand with the command:

    bowtie2-build -o 0 -q query_database.fna query_database.fna

You can run the `inferseq-database` command as:

    mustache inferseq-database mustache.pairflanks.tsv query_database.fna
    
This step should take about 10 seconds to complete. This command requires as input the file generated by the 
`pairflanks` command, and the query database. The output of this command is in the `mustache.inferseq_database.tsv` 
file. Please see details about this output file in the [manual](manual.md). This is arguably the most sensitive 
approach, provided that the database contains all of the inserted elements of interest. 

[Detailed user manual](manual.md)

## Next steps
This tutorial produced a list of candidate insertions for the E. coli isolate SRR3180793 (`mustache.pairflanks.tsv`), as
well as the inferred identity of many of the insertions using a variety of inference techniques. Of the 149 candidate 
insertions identified in the `mustache.pairflanks.tsv` file, a sequence was inferred for 137 using at least one of the
inference approaches described here. Those candidate insertions whose identity could not be inferred may include large 
genomic inversions, or insertions that were too complicated to be properly inferred from short-read sequencing using
*mustache*.

The next steps taken will depend largely on the type of analysis that you want to perform. If you want to search these
candidate insertions for antibiotic resistance genes, you can create a FASTA file of the inferred sites and upload the file
to [ResFinder](https://cge.cbs.dtu.dk/services/ResFinder/). Other more complicated analyses may require additional
processing of the data. 

A given candidate insertion may have multiple inferred sequences that belong to it. This ambiguity makes it difficult to
make a final decision about the identity of a given insertion. A description of the approach that we took in our own
analysis can be described in detail in our [preprint](https://www.biorxiv.org/content/10.1101/527788v1). Briefly,
we chose to cluster all inferred sequences at 90% similarity across 85% of each sequence using CD-HIT-EST. This meant 
that sequences that were highly similar to each other were assumed to be the same insertion. We then prioritized the 
sequences that were inferred by the different inference methods described here. The `inferseq-assembly` step was given
top priority, as these elements come directly from the organism's assembly, and are most likely very close to the true
identity of the element. The `inferseq-overlap` sequences were given second priority. The `inferseq-database` and the
`inferseq-reference` sequences were given last priority, meaning the inferred sequences here were only used if no
`inferseq-assembly` or `inferseq-overlap` sequences could be inferred.

Even still, ambiguities may exist. We tried to resolve some of these remaining ambiguities, and a description of the 
approach we took can be found in our [preprint](https://www.biorxiv.org/content/10.1101/527788v1).

## Other considerations
We hope we have made the limitations of this approach clear, and we urge you to carefully consider these limitations 
when analyzing your own data.

Some additional details to keep in mind:
* *Mustache* is not designed to identify small insertions. If an insertion is smaller than your library's read length,
we urge you to ignore them as they cannot be reliably identified. We recommend a lower size cutoff of 300 bp when
analyzing sequences.
* Similarly, *mustache* has not been evaluated with respect to very large insertions. The larger an insertion,
the less likely that it will be able to assemble, and it will not be identified in the `inferseq-assembly` step. We
used an upper cutoff of 10 kilobase pairs in our own analysis.

Good luck, and if you have any questions please submit an issue [here](https://github.com/durrantmm/mustache/issues).

[NEXT: Detailed User Manual](manual.md)