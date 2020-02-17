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

## Special note on site-specific integrative mobile elements
While *MGEfinder* can detect a wide variety of site-specific integrative mobile elements, it is best suited for transposable elements. For example, certain types of tRNA-targeting mobile elements will be missed by *MGEfinder* by default because of their unique integration mechanism, whereby they replace the target sequence with a new sequence to repair the tRNA. But *MGEfinder* can still find these integrative elements if you adjust the parameters properly. If you would like advice on how to do this, please open an [issue](https://github.com/bhattlab/MGEfinder/issues) with your request.

## Publication
Durrant, M. G., Li, M. M., Siranosian, B. A., Montgomery, S. B. & Bhatt, A. S. [A Bioinformatic Analysis of Integrative Mobile Genetic Elements Highlights Their Role in Bacterial Adaptation](https://www.cell.com/cell-host-microbe/fulltext/S1931-3128(19)30546-3). *Cell Host & Microbe* 0, (2019)

## Questions / Comments
Please submit any questions or comments to our [issues handler](https://github.com/durrantmm/mgefinder/issues). 
