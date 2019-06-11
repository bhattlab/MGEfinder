[Back to main page](../README.md)  

# User Manual 
This user manual is designed as a reference for users who have questions about *mustache*. It includes details about
how to execute important commands, the parameters available, brief descriptions of how the command works, and 
descriptions of the output files.

## Primary *mustache* commands
We now describe each of the primary *mustache* commands. These commands are at the center of the *mustache* workflow,
and should not be ignored when doing a complete analysis.

### `findflanks`
![alt text](img/findflanks.png)

This command finds insertion sites and reconstructs the flanks of inserted sequence.

#### `findflanks`: Input and parameters
The `findflanks` command takes a BAM file as input. A local alignment algorithm, such as 
[BWA MEM](http://bio-bwa.sourceforge.net/) should be used to generate this BAM file.

You can run the command as

    mustache findflanks BAMFILE
    

Additional parameters include:

    --min_softclip_length, -minlen
    --min_softclip_count, -mincount
    --min_count_consensus, -mcc
    --min_alignment_quality, -minq
    --min_softclip_ratio, -minratio
    --min_alignment_inner_length, -minial
    --min_distance_to_mate, -mindist
    --max_indel_ratio, -maxir

The parameter `min_softclip_length` takes an integer. 
For a softclipped site to be considered, there must be at least one softclipped read of this length. 
The default used is 8.
 
The parameter `min_softclip_count` takes an integer.
It requires that for a clipped site to be considered it must be supported by at least this number of reads.
The default used is 4.

The parameter `min_count_consensus` takes an integer.
When generating a consensus sequence for the insertion flank, it requires that each site must have at least this many 
reads supporting it.
The default used is 2.

The parameter `min_alignment_quality` takes an integer.
This determines the minimum mapping quality of each alignment in order for it to be considered.
The default used is 20.

The parameter `min_softclip_ratio` takes a float between 0 and 1.0.
It sets a minimum for the ratio of clipped reads to all reads at a given site.
This filters out very small indels that could disrupt the flank-pairing step.
The default used is 0.15.

The parameter `min_alignment_inner_length` takes an integer.
It filters out reads that are softclipped on both ends, and where the inner aligned portion is shorter than this length.
Filtering out such reads removes a potential false positives.
The default is 21, this parameter should not be changed.

The parameter `min_distance_to_mate` takes an integer.
It filters out clipped-sites that have no other candidate flank pair nearby. 
In the `pairflanks` step, we only pair flanks with each other if they are within 20 bp of 
each other. Filtering out sites without a candidate pair can dramatically speed up the 
analysis. However, if you want to keep all candidate insertion sites, you can increase this
number to some very high integer. The default is 22.

#### `findflanks`: Description of implementation
The `findflanks` algorithm works by identifying candidate insertion sites by searching for clipped-end sites in locally 
aligned reads. To generate a consensus sequence of the candidate flank, we use a trie-based approach intended to filter 
out spurious reads when building a consensus sequence.  This also allows our algorithm to distinguish between multiple 
insertions at a single site, which may be observed in a metagenomic sample, for example (analysis of metagenomic samples
is in development).

Briefly, we build a sequence Trie of all of the reads found clipped at a given site. We then generate all of the 
unique, complete sequences by traversing the trie. We then cluster these sequences together by sequence similarity, 
and process the different clusters separately.

We then generate a consensus sequence for a given insertion flank by traversing down the trie, choosing the base with 
the most support at each step. Evidence for a given base is calculated as the sum of all base-pair quality scores 
originally reported in the FASTQ file.

The consensus sequence ends when the number of reads supporting a position drops below `--min_count_consensus`.

#### `findflanks`: Output file format
By default, `findflanks` will write to a file named `mustache.findflanks.tsv` (The name can be changed with the 
`--output_file` parameter). The columns of this file are described as follows:

1. `flank_id` - An arbitrary identifier used internally by *mustache*.
2. `contig` - The name of the contig where the flank was identified.
3. `pos` - The exact base pair position of where the clipped-ends begin, the presumed start of the insertion. 
This location is 0-based with respect to the reference genome used.
4. `orient` - The orientation of the flank. `5p` refers to the 5' end of an insertion, which suggests this flank runs 
from the 5' to the 3' direction with respect to the reference genome. `3p` refers to the 3' end of an insertion, which 
suggests the flank runs from the 3' to the 5' direction with respect to the reference genome.
5. `softclip_count_5p` - The total number of reads clipped from the 5' to 3' direction at this site.
6. `softclip_count_3p` - The total number of reads clipped from the 3' to 5' direction at this site.
7. `runthrough_count` - The total number of reads that run through the insertion site without being clipped.
8. `small_insertion_count_5p` - The total number of reads that include a small insertion from the 5' to 3' direction at
this site. If this number is high, it indicates that the detected insertion is most likely a small indel. `findflanks` 
uses this to filter out false positives.
9. `small_insertion_count_3p` - Same as above, but for small insertions running from the 3' to the 5' direction.
10. `deletion_count` - The number of reads with a deletion at the insertion site, also used to filter out false 
positives.
11. `upstream_deletion_count` - The number of reads that contain a deletion one base pair upstream of this insertion 
site. Also
useful for removing false positives
12. `downstream_deletion_count` - Same as above, but for reads with deletions one base pair downstream of this insertion
site.
13. `total_count` - The sum of the `softclip_count_5p`, `softclip_count_3p`, `runthrough_count`, 
`small_insertion_count_5p`, `small_insertion_count_3p`, and `deletion_count` columns.
14. `consensus_softclip_count` - The number of clipped reads at this site that specifically support the 
consensus flank sequence on this line. This will be used as the `softclip_count_5p` and `softclip_count_3p` columns
in the file generated by running the `pairflanks` command.
15. `consensus_seq` - The consensus sequence for the insertion flank identified at this site.
    
This intermediate file is not particularly valuable to users on its own, and it is primarily intended to be used
internally by *mustache* in the next step.

### `pairflanks`
This command pairs insertion flanks with each other to represent the 5' and 3' flanks of a candidate insertion.

#### `pairflanks`: Input and parameters
You can run the `pairflanks` command as

    mustache pairflanks FLANKSFILE BAMFILE REFERENCE_GENOME

Where `FLANKSFILE` is the output of the `findflanks` command, `BAMFILE` is the binary sequence alignment file used as 
input for `findflanks`, and `REFERENCE_GENOME` is the reference genome that `BAMFILE` was aligned to.

Additional parameters include:

    --max_direct_repeat_length, -maxdr
    --min_alignment_quality, -minq
    --min_alignment_inner_length, -minial
    --max_junction_spanning_prop, -maxjsp

The parameter `max_direct_repeat_length` takes an integer. This specifies the maximum distance that oppositely-oriented
insertion flanks can be from each other in order to consider pairing them together. Since insertions often cause direct
repeats at the insertion site, the location of the insertion flanks are often separated by several base pairs. 
The default for this parameter is 20. This means that our algorithm by default will NOT identify true insertions that 
create direct repeats that exceed 20 base pairs.

The parameter `min_alignment_quality` takes an integer. This determines the minimum mapping quality of each alignment
in order for it to be considered. The default used is 20. This should be the same as the `--min_alignment_quality`
parameter in `findflanks`.

The parameter `min_alignment_inner_length` takes an integer.
It filters out reads that are softclipped on both ends, and where the inner aligned portion is shorter than this length.
Filtering out such reads removes a potential false positives.
The default is 21, this parameter should not be changed, and it should be the same as the value of 
`--min_alignment_inner_length` used in `findflanks.`

The parameter `max_junction_spanning_prop` takes an integer. This is used to filter out low-confidence insertion pairs,
or insertions that are occurring in duplicated regions. We expect that very few reads will fully span the insertion
junction, as that would indicate that the insertion site does not contain the insertion somewhere in the isolate. 
If the number of reads spanning the insertion junction without being clipped exceeds this proportion of the total reads
at the site, then it will be ignored. By default, this parameter is 0.15.

#### `pairflanks`: Description of implementation
The `pairflanks` command uses a variety of techniques to pair flanks with each other. It first filters insertion flanks 
by the `--max_direct_repeat_length` parameter. It then does pairwise comparisons between all nearby, oppositely-oriented 
insertion flanks to determine if they share inverted repeats at their termini, a common feature of many prokaryotic 
insertion sequences. It then prioritizes all pairs by the following parameters, in order: 
1. The length of the inverted repeat identified, if any. Pairs with longer shared inverted repeats receive higher 
priority.
2. The difference in the number of softclipped reads that support the insertion flanks. 
Pairs that have a similar number of reads supporting each flank are given higher priority.
3. The difference in the length of the recovered insertion flanks. 
If the two insertion flanks have very similar lengths, they are given higher priority than pairs with disparate flank 
lengths.

It then filters sites with too many junction-spanning reads (see `max_junction_spanning_prop`), and infers the identity
of the direct repeat created by the insertion, if any. We realize that these are somewhat arbitrary filters, and we 
welcome suggestions from the users on how this might be improved.

#### `pairflanks`: Output file format
By default, `pairflanks` will write to a file named `mustache.pairflanks.tsv` (The name can be changed with the 
`--output_file` parameter). The columns of this file are described as follows:

1. `pair_id` - An arbitrary identifier used to identify the candidate insertion. This is an important identifier that
is used to relate the output of the `pairflanks` command to the output of the various `inferseq` commands.
2. `contig` - The name of the contig where the pair was identified.
3. `pos_5p` - The exact base pair position of where the clipped-ends of the 5' insertion flank begin, the presumed 
start of the insertion. This location is 0-based with respect to the reference genome used.
3. `pos_3p` - The exact base pair position of where the clipped-ends of the 3' insertion flank begin, the presumed 
end of the insertion. This location is 0-based with respect to the reference genome used.
4. `softclip_count_5p` - The number of clipped reads at this site that support the 5' insertion flank. This is taken 
directly from the `consensus_softclip_count` of the `findflanks` output file.
5. `softclip_count_4p` - The number of clipped reads at this site that support the 3' insertion flank. This is also 
taken directly from the `consensus_softclip_count` of the `findflanks` output file.
6.  `total_count_5p` - The total number of reads found at the insertion site for the 5' flank.
7. `total_count_3p` - The total number of reads found at the inseriton site for the 3' flank.
8. `spanning_count` - The total reads that are found to span the 5' and 3' insertion sites by 10 bp in both directions.
9. `has_IR` - A boolean True/False indicating whether or not the flank pair contains inverted repeats at their end, a 
common feature of many inserted elements.
10. `IR_length` - The length of the detected inverted repeat.
11. `IR_5p` - The inverted repeat in the 5' flank.
12. `IR_3p` - The inverted repeat in the 3' flank.
13. `seq_5p` - The consensus sequence of the 5' flank.
14. `seq_3p` - The consensus sequence of the 3' flank.
15. `direct_repeat_reference` - The sequence of the direct repeat created, using the reference genome to infer this
direct repeat. 
16. `direct_repeat_reads_consensus` - The sequence of the direct repeat created, but using information from the reads 
themselves to determine the a more accurate direct repeat. This can be useful if the direct repeat of the isolate
actually contains some base pair substitutions that do not exist in the reference genome.

This file will likely have intrinsic value to the user, as they indicate candidate insertions. However, these candidates
must be further investigated before any final conclusions about their identity can be made. In some cases, for example,
these may actually represent genomic inversions, not insertions. To infer the identity of these insertions, several
inference approaches are implemented and described in the following sections.

### `inferseq-reference`
![alt text](img/inferseqreference.png)

This command infers the identity of insertions by aligning the flanks of candidate pairs to a reference genome.

#### `inferseq-reference`: Input and parameters
The `inferseq-reference` command takes the output of the `pairflanks` command and a reference genome as input.

While *mustache* will automatically index the reference genome if it is not already indexed, we recommend that you index 
the reference genome beforehand with the command:

    bowtie2-build -o 0 -q INFERSEQ_REFERENCE INFERSEQ_REFERENCE

Where INFERSEQ_REFERENCE is the reference genome of interest.

[BWA MEM](http://bio-bwa.sourceforge.net/) should be used to generate this BAM file.

You can then run the command as

    mustache inferseq-reference PAIRSFILE INFERSEQ_REFERENCE
    
Additional parameters include:

    --min_perc_identity, -minident
    --max_internal_softclip_prop, -maxclip
    --max_inferseq_size, -maxsize
    --min_inferseq_size, -minsize
    --keep-intermediate/--no-keep-intermediate


The parameter `min_perc_identity` takes an float between 0 and 1. When aligning candidate insertion flanks to the
reference genome, it will only consider alignments that exceed this percentage identity with the reference. The default
for this parameter is 0.95. 

The parameter `max_internal_softclip_prop` takes an float between 0 and 1. This is an additional candidate insertion
flank alignment filter. If the aligned flanks are internally clipped by a proportion of their total length that exceeds
this number, then the alignment is excluded. For example, imagine the alignment:

    -------->          <----------
    flank1                  flank2

Imagine the vertical bar `|` indicating that one of the alignments is clipped at a specific site, such as:
    
    -------->         <---|-------
    flank1                  flank2

If the proportion of the clipped end of the flank2 alignment exceeds 0.05 (by default), then this alignment will be
excluded.

The parameter `max_inferseq_size` excludes inferred sequences that are larger than `max_inferseq_size`. The default for
this parameter is 500 kilobase pairs, but we urge users to be cautious when working with very large insertions, as these
may be false positives.

The parameter `min_inferseq_size` excludes inferred sequences that are smaller than `min_inferseq_size`. The default for
this parameter is 1 base pair, but in reality *mustache* is not well suited to identify insertions that are smaller than
the read length of the sequencing library.

The `--keep-intermediate/--no-keep-intermediate` will determine whether or not the intermediate alignment file will
be deleted or kept, which can be useful for debugging purposes.

#### `inferseq-reference`: Description of implementation
A schematic of how this step is implemented is shown above, but more details are given here. First, if the reference
genome is not already bowtie2 indexed, *mustache* will index the genome. Next, it will align the individual flanks 
to the reference genome. It will then filter the alignments according to the `min_perc_identity` parameter, and it will
also filter alignments that are clipped to on their outer edge. It will then iterate through the sorted alignments and
pair flanks with each other. This approach will actually take the minimum non-overlapping alignment pairs, and ignore
all other combinations of alignments, as described in the schematic below:

![alt text](img/minimum_nonoverlapping_windows.png)

Next, candidate pairs will be filtered according to the `max_internal_softclip_prop` parameter, and the 
inferred sequence size filters (`max_inferseq_size` and `min_inferseq_size`). Finally using the sum of the alignment
scores of both flanks in each pair, we filter out all pairs that do not have the maximum alignment score. This leaves
us with a set of final inferred sequences of equal quality.

#### `inferseq-reference`: Output file format
By default, `inferseq-reference` will write to a file named `mustache.inferseq_reference.tsv` (The name can be changed 
with the `--output_file` parameter). The columns of this file are described as follows:

1. `pair_id` - An identifier used to map the inferred sequence in this file to the candidate insertion pairs described
in the `pairflanks` file. That is, all lines in the `mustache.inferseq_reference.tsv` with the value `20` in the 
`pair_id` column represent an inferred sequence identified by the `inferseq_reference` command that corresponds to the
candidate insertion pair identified in the `pairflanks` file that is labeled with `20` in the `pair_id` column.
2. `method` - The method used to infer the sequence, which will always be `inferred_reference` when running the 
`inferseq-reference` command.
3. `loc` - The location of the inferred sequence as it was found in the reference genome. This is in the format
`CONTIG:START-END` where `CONTIG` is the contig where the inferred sequence was located, `START` is the base pair
start of the inferred sequence, and `END` is the base pair end of the inferred sequence.
4. `inferred_seq_length` - The length of the inferred sequence.
5. `inferred_seq` - The identity of the inferred sequence in full.

This file is of value to the user, as it infers the full identity of the candidate insertions by referring to a
reference genome. In certain situations this may be sufficient, such as when the sequenced isolate differs only
slightly from the reference genome (like in a resequencing experiment). But if the inserted element may not exist in
an available reference genome, alternative inference approaches may be of interest (see below).

It should also be noted that the reference genome used here does not necessarily have to be the reference genome that
the sequencing reads were aligned to.

### `inferseq-assembly`
![alt text](img/inferseqassembly.png)

This command infers the identity of insertions by aligning the flanks of candidate pairs to a sequence assembly of the
isolate of interest. We recommend using [SPAdes](http://cab.spbu.ru/software/spades/)  to assemble the genome of the 
isolate of interest.

#### `inferseq-assembly`: Input and parameters
The `inferseq-assembly` command takes as input the output of the `pairflanks` command, the BAM file of the sequencing 
reads aligned to the reference genome, an assembly file in FASTA format of the sequencing reads derived from the isolate 
of interest, and the reference genome used when initially aligning the isolates reads. In this step, it is imperative
that the sequence assembly comes from the isolate of interest.

While mustache will automatically index the sequence assembly if it is not already indexed, we recommend that you index 
the assembly beforehand with the command: 

    bowtie2-build -o 0 -q INFERSEQ_ASSEMBLY INFERSEQ_ASSEMBLY

You can then run the command as

    mustache inferseq-assembly PAIRSFILE BAMFILE INFERSEQ_ASSEMBLY INFERSEQ_REFERENCE
    
Additional parameters include:

    --min_perc_identity, -minident
    --max_internal_softclip_prop, -maxclip
    --max_inferseq_size, -maxsize
    --min_inferseq_size, -minsize
    --keep-intermediate/--no-keep-intermediate


These additional parameters are identical in function to the parameters described in the 
"`inferseq-reference`: Input and parameters" section above. Please refer to this section for more details. 

#### `inferseq-assembly`: Description of implementation
Details about how this step is implemented are shown in Figures a and b above. The sequence of steps is very similar to
those described in the "`inferseq-reference`: Description of implementation" section above, with important differences. 
Two alignments to the assembly are performed. In the first step, 25 base pairs of reference genome are appended to the 
beginning of the flanks before alignment (See Figure a above). This allows us to determine if the insertion has 
assembled within the expected sequence context, which gives us greater confidence that it the inferred identity of the 
insertion is of high quality. Then the flanks are aligned without the appended context sequence to infer the identity of
the insertion. Both of these types of inferred sequences are then reported in the output file.

#### `inferseq-assembly`: Output file format
By default, `inferseq-assembly` will write to a file named `mustache.inferseq_assembly.tsv` (The name can be changed 
with the `--output_file` parameter). The columns of this file are the same as those found in the
"`inferseq-reference`: Output file format" section above, but with the following differences:

2. `method` - The method used to infer the sequence. With the `inferseq-assembly` command, this will be one of three
different values: 1) `inferred_assembly_with_full_context`, meaning that the inserted was inferred from the sequence
assembly with the 25 base pairs of context sequence appended on both sides, the highest quality inferred sequence. 2)
`inferred_assembly_with_half_context` means that the sequence was inferred in context, but one side of the inserted
element was found at the end of a contig, meaning that the sequence context of one side is still ambiguous. 3) 
`inferred_assembly_without_context` indicates the the sequence was inferred from the assembly, but without the
additional sequence context. This is common when the inserted sequence is a large repeated element and it cannot be
properly assembled.
3. `loc` - The location of the inferred sequence as it was found in the sequence assembly. This is in the format
`CONTIG:START-END` where `CONTIG` is the contig where the inferred sequence was located, `START` is the base pair
start of the inferred sequence, and `END` is the base pair end of the inferred sequence.

This file is of value to the user, as it infers the full identity of the candidate insertions by referring to sequence
assembly. This is of value when the sequenced isolate and the reference genome used are assumed to be quite different
from each other. It is limited, however, by the quality of the genome assembly. Large mobile genetic elements that 
cannot be properly assembled will be missed.

### `inferseq-overlap`
In certain situations, the identity of the inserted sequence can be inferred by attempting to overlap the flanks.
For example:

    flank1 ----------------->
                        ||||
                       <-------------------flank 2

In the situation depicted above, `flank1` and `flank2` overlap, and they can be merged into a single insertion.

This command performs this merging operation, if possible.

#### `inferseq-overlap`: Input and parameters
The `inferseq-overlap` command takes as input only the output of the `pairflanks` command.

You can run the command as

    mustache inferseq-overlap PAIRSFILE
    
Additional parameters include:

    --min_overlap_score, -minscore
    --min_overlap_perc_identity, -minopi

The parameter `min_overlap_score` determines the minimum alignment score necessary for two flanks to be merged. 
Matched positions increase the score by 1, and mismatch positions decrease the score by 1. 

The parameter `min_overlap_perc_identity` is the minimum sequence identity between the overlapping portion of the two
flanks that is required for them to merge.

#### `inferseq-overlap`: Description of implementation
The merging procedure taken here is quite simple. Flanks are scanned across each other one base pair at a time.
The overlapping portions are checked for matching base pairs, with matches increasing the overlap score by one, and 
mismatches decreasing the score by one. The parameters `min_overlap_score` and `min_overlap_perc_identity` are used
to determine if the flanks are merged into a single sequence.


#### `inferseq-overlap`: Output file format
By default, `inferseq-overlap` will write to a file named `mustache.inferseq_overlap.tsv` (The name can be changed 
with the `--output_file` parameter). The columns of this file are the same as those found in the
"`inferseq-reference`: Output file format" section above, but with the following differences:

2. `method` - The method used to infer the sequence which will always be `inferrred_overlap` when using the 
`inferseq-overlap` command.
3. `loc` - The portions of flank1 and flank2 that are kept to form the final merged sequence. This is in the form
`seq_5p:START_5p-END_5p;seq_3p:START_3p-END_3p`.

This is valuable to the user primarily for filtering out insertions that are smaller in size, as this step is 
fundamentally limited by the read length of the library. *Mustache* is not well-suited to identify smaller
insertions, and this step can be used to filter out these small insertions so they do not interfere with downstream
analyses.

### `inferseq-database`
![alt text](img/inferseqdatabase.png)

This command infers the identity of insertions by aligning the flanks of candidate pairs to a database of query
sequences. Ideally, this database is filtered to remove redundant sequences, which may dramatically reduce the amount
of time required to run this step.

#### `inferseq-database`: Input and parameters
The `inferseq-database` command takes as input the output of the `pairflanks` command, and a FASTA file of the query
sequences of interest.

While mustache will automatically index the query sequence database if it is not already indexed, we recommend that you
index the database beforehand with the command: 

    bowtie2-build -o 0 -q INFERSEQ_DATABASE INFERSEQ_DATABASE

You can then run the command as:

    mustache inferseq-database PAIRSFILE INFERSEQ_DATABASE
    
Additional parameters include:

    --min_perc_identity, -minident
    --max_internal_softclip_prop, -maxclip
    --max_edge_distance, -maxedgedist
    --keep-intermediate / --no-keep-intermediate
  

The parameter `min_perc_identity` takes an float between 0 and 1. When aligning candidate insertion flanks to the
query sequence database, it will only consider alignments that exceed this percentage identity with the reference. 
The default for this parameter is 0.90, which is more lenient than the cutoff used in the `inferseq-reference` and
`inferseq-assembly` commands. 

The parameter `max_internal_softclip_prop` takes an float between 0 and 1. This is an additional candidate insertion
flank alignment filter. If the aligned flanks are internally clipped by a proportion of their total length that exceeds
this number, then the alignment is excluded. See "`inferseq-reference`: Input and parameters" for more details on how 
this parameter works.

The parameter `max_edge_distance` takes an integer. This determines how close the aligned flanks must be to the edge
of the sequence in order to keep the alignments (See Figure d above). The default for this parameter is 10.

The `--keep-intermediate/--no-keep-intermediate` will determine whether or not the intermediate alignment file will
be deleted or kept, which can be useful for debugging purposes.

#### `inferseq-database`: Description of implementation
This sequence inference approach assumes that the inserted sequences of interest exist within the query database itself.
It aligns flank pairs to this database, and determines whether or not they align to a given a sequence. Aligned flanks
that do not meet the `min_perc_identity` cutoff are removed, and only inferred sequences where both flank pairs align
to the edge of the sequence of interest are kept, where the aligned flanks are allowed to be within `max_edge_distance`
base pairs of the end of the query sequence.

#### `inferseq-database`: Output file format
By default, `inferseq-database` will write to a file named `mustache.inferseq_database.tsv` (The name can be changed 
with the `--output_file` parameter). The columns of this file are the same as those found in the
"`inferseq-reference`: Output file format" section above, but with the following differences:

2. `method` - The method used to infer the sequence which will always be `inferrred_daatabase` when using the 
`inferseq-database` command.

All sequences in the database that are matched equal matches for a given flank pair will be returned.

This is a valuable approach to taken when information about the inserted elements is already known beforehand. This
database can include several insertion sequences that are relevant to the species under study, for example. It is 
quite sensitive, provided that read coverage is high enough (we recommend above 40x) and the inserted sequence of
interest exists in the query database.

## Additional/ Experimental *mustache* Commands
### `formatbam`
This command is intended to prepare SAM files for analysis by *mustache* as input. While this step is not strictly 
necessary under all circumstances, some features of *mustache* require it (such as the experimental `extendpairs` 
command). We recommend that all SAM/BAM files are processed using this command prior to further analysis if it is 
computationally feasible.

#### `formatbam`: Input and parameters
This command requires only the alignment file SAM/BAM of interest. This is a file thatyou wish to analyze using the
`findflanks` command in *mustache*. You can run this command as follows:

    mustache formatbam IN_SAM OUT_BAM
    
Where `IN_SAM` is the input SAM/BAM file of interest, and `OUT_BAM` is the name of the output file.

Additional parameters include `--single-end`, which is used when the file is a single-end and not paired-end, and
`--keep-tmp-files` which will keep intermediate files without rather than deleting them.
    
#### `formatbam`: Description of implementation
This command first removes secondary alignments from the alignment file. It then formats the file for use by *mustache*,
in particular it stores information about each reads mate so that it is immediately accessible, a step that is 
necessary when using the experimental `extendpairs` command. It will then sort and index the BAM files.

#### `formatbam`: Output file format

The output of this command is just a BAM file that is ready to be analyzed by *mustache*.

### `recall`
The `recall` command is used to collect information about potential insertion sites from a BAM file. For example,
if you have many insertion sites of interest, you may want to determine if these sites were actually deleted in the
reference genome of some of the isolates. This can be useful when comparing isolates to each other. This command was
used in our own GWAS of antibiotic resistance in *E. coli* to determine if some of the insertion sites were deleted in 
some of the isolates. This can allow you to then filter out these unreliable insertion sites.

#### `recall`: Input and parameters
The `recall` command takes as input the output of the `pairflanks` command, and the original `BAM` file of the reads
aligned to the reference genome. The `pairflanks` output file actually only requires three columns: 
`contig`, `pos_5p`, and `pos_3p`. The `recall` command will analyze the sites specified in these columns.

This command can be executed as:
   
    mustache recall PAIRSFILE BAMFILE
    
Additional parameters include:

    --min_alignment_quality, -minq
    --min_alignment_inner_length, -minial
    
Which should be the same as the values used when originally running the `findflanks` command.
By default `min_alignment_quality` is 20 and `min_alignment_inner_length` is 21.

#### `recall`: Description of implementation

Recall works very similarly to the findflanks algorithm, but it does not attempt to reconstruct consensus flanks.
It's only goal is to gather read information about the sites present in the `PAIRSFILE`. This can then be used in other
analyses.

#### `recall`: Output file format
The output format is quite similar to the output of the `findflanks` command. Each row is a specific site (note that
it splits each pair into its individual insertion sites). It includes columns about read counts, such as
`softclip_count_5p`, `softclip_count_3p`, `runthrough_count`, `small_insertion_count_5p`, `small_insertion_count_3p`,
`deletion_count`, `upstream_deletion_count`, `downstream_deletion_count`, and `total_count`. See the
section "`findflanks`: Output file format" for more details on what each of these columns means.


### `extendpairs`

This is an experimental command in *mustache*, so use at your own risk. This command is used to perform a local
assembly on the candidate insertions in order to extend the lengths of the consensus flanks. This effectively
increases the upper limit of the consensus flanks such that it is closer to the length of the library's fragment length,
rather than the read length. This allows one to determine more information about an insertion by a more direct 
measurement, rather than by inference.

This step requires that the user install a working version of [AMOS](http://amos.sourceforge.net/wiki/index.php/AMOS) 
sequence assembly software on their machine. This is not included in the default installation of *mustache*.

#### `extendpairs`: Input and parameters
This command requires as input the output of the `pairflanks` command, and the BAM of the read alignment to the 
reference genome. This BAM file MUST have beem processed by the `formatbam` command. 
The command can then be executed as:

    mustache extendpairs PAIRSFILE BAMFILE
    
Additional parameters include:

    --threads, -t
    
Which specifies the number of threads used to perform the local assemblies. This will increase speed at the cost of
using additional CPUs. This is 1 by default.

#### `extendpairs`: Description of implementation
In brief, `extendpairs` iterates through all of the candidate insertions, and performs local assemblies of all of the
reads found at the insertion site. It is able to draw on information about reads that did not align to the reference
genome near the insertion site, but whose read mate did align. The local assembly of these reads is then used as a 
reference for a candidate flank alignment, and the extended sequence is determined from this alignment.


#### `extendpairs`: Output file format
The output of this file is the same as the output of the `pairflanks` command, with two additional columns added
at the end:

`extended_5p` is a boolean True/False column that indicates whether or not the local assembly was able to extend the
5' flank.

`extended_3p` is a boolean True/False column that indicates whether or not the local assembly was able to extend the
3' flank.


Again, this is an experimental command, but it may be of use to the user.

[Back to main page](../README.md)
