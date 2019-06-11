[Back to main page](../README.md)  

# Installing *mustache*
First, install [miniconda3](https://conda.io/en/master/miniconda.html). This is an environment management system that 
should keep everything organized.

Once installed, you have two options for downloading mustache:

1) To download the current development version:

        git clone https://github.com/durrantmm/mustache.git
        cd mustache

2) To download the most recent *mustache* release:

        wget https://github.com/durrantmm/mustache/archive/v1.0.0.tar.gz
        tar -zxvf v1.0.0.tar.gz
        rm v1.0.0.tar.gz
        cd mustache-1.0.0
    
Install the mustache conda environment and other dependencies with

    bash install.sh

Now activate the environment with
    
    source activate mustache
    
This activation step must be repeated whenever using *mustache*.

Now install *mustache* with the command

    pip install .
    
Once complete, you can check to see if *mustache* installed properly by simply typing

    mustache
   
This can then be called from anywhere on the file system while in the `mustache` conda environment.


[NEXT: Step-by-step tutorial](tutorial.md)



