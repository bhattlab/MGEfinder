![mustache image](https://github.com/bhattlab/MGEfinder/blob/master/docs/img/mustache.png)

# *MGEfinder* - A toolbox for identifying mobile genetic element (MGE) insertions from short-read sequencing data of bacterial isolates.

The command-line tool *MGEfinder* identifies large insertions and genotypes them with respect to a reference genome.

It is designed to work with haploid genomes, and has been tested extensively on bacteria.

It can identify mobile genetic elements and their sites of insertion using an *ab initio* approach.

Follow the links below to learn more.

## Table of Contents
1. [How it works](https://github.com/bhattlab/MGEfinder/wiki/How-it-works)
2. [Install *MGEfinder*](https://github.com/bhattlab/MGEfinder/wiki/Installation)
3. [Step-by-step tutorial](https://github.com/bhattlab/MGEfinder/wiki/Tutorial)
4. [Detailed user manual](https://github.com/bhattlab/MGEfinder/wiki/User-manual)

## *IMPORTANT UPDATES*

### Database-only mode (when you don't have assemblies)
We have added a new feature in MGEfinder v1.0.5 that allows you to run an MGEfinder workflow without isolate assemblies. This means that you no longer need to assemble each isolate if you already know what elements you are looking for. You just need a working directory that includes this file:

    myWorkdir/00.database/database.fna

Instead of the `myWorkdir/00.assembly` directory. The `00.database/database.fna` file is a FASTA file of the elements that you are searching for. For example, if you just want to look at IS6110 insertions in *M. tuberculosis*, you can place a copy of that IS element in the `00.database/database.fna` file. You can then run the workflow with the command

    mgefinder workflow database myWorkdir
    
While the original *de novo* workflow that requires assemblies can be run with:
  
    mgefinder workflow denovo myWorkdir

Let me know if you have questions, and I hope you find it useful!

### Sensitive mode
This is a flag you can add to the `mgefinder workflow` commands to potentially increase sensitivity. These settings have not been validated, but they should make it easier to identify certain integrative mobile element insertions, in particular insertions that create a direct repeat between 20 bp and 50 bp in length. This should help to increase sensitivity to detect insertions of elements that insert via a tRNA-targeting tyrosine integrase, for example.

You can use it with

    mgefinder workflow database --sensitive myWorkdir
    
or

    mgefinder workflow denovo --sensitive myWorkdir

Let me know if you have any questions by submitting an [issues](https://github.com/bhattlab/MGEfinder/issues).

## Special note on site-specific integrative mobile elements
While *MGEfinder* can detect a wide variety of site-specific integrative mobile elements, it is best suited for transposable elements. For example, certain types of tRNA-targeting mobile elements will be missed by *MGEfinder* by default because of their unique integration mechanism, whereby they replace the target sequence with a new sequence to repair the tRNA. But *MGEfinder* can still find these integrative elements if you adjust the parameters properly. If you use the `--sensitive` flag when running `mgefinder workflow denovo` or `mgefinder workflow database`, you should be able to better detect these elements. If you would like more advice on how to do this, please open an [issue](https://github.com/bhattlab/MGEfinder/issues) with your request.

## Publication & Presentation
Durrant, M. G., Li, M. M., Siranosian, B. A., Montgomery, S. B. & Bhatt, A. S. [A Bioinformatic Analysis of Integrative Mobile Genetic Elements Highlights Their Role in Bacterial Adaptation](https://www.cell.com/cell-host-microbe/fulltext/S1931-3128(19)30546-3). *Cell Host & Microbe* 0, (2019)

Please also consider viewing my [video presentation](https://jrnlclub.org/research-films/mgefinder) of this paper, provided through [JRNLclub](https://jrnlclub.org/research-films/mgefinder).

## Questions / Comments
Please submit any questions or comments to our [issues handler](https://github.com/bhattlab/MGEfinder/issues). 
