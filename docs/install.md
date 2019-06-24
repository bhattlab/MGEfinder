[Back to main page](../README.md)  

# Installing *mgefinder*
*MGEfinder* is currently only available on linux-64 systems.

First, install [miniconda3](https://conda.io/en/master/miniconda.html). This is an environment management system that 
should keep everything organized.

Once installed, we recommend installing *MGEfinder* as follows:

1) Activate your conda environment or create a new conda environment:
    
        conda create -n mgefinder
        source activate mgefinder
    
    This environment activation step must be repeated whenever using *MGEfinder*.

2) Now install *mgefinder* with the command

        conda install -c bioconda -c conda-forge -c anaconda -c mdurrant mgefinder 
    
Once complete, you can check to see if *mgefinder* installed properly by typing

    mgefinder
   
This can then be called from anywhere on the file system while in the `mgefinder` conda environment.

[NEXT: Step-by-step tutorial](tutorial.md)



