[Back to main page](../README.md)  

# Installing *mgefinder*
First, install [miniconda3](https://conda.io/en/master/miniconda.html). This is an environment management system that 
should keep everything organized.

Once installed, you have two options for downloading mgefinder:

1) To download the current development version:

        git clone https://github.com/durrantmm/mgefinder.git
        cd mgefinder

2) To download the most recent *mgefinder* release:

        wget https://github.com/durrantmm/mgefinder/archive/v1.0.0.tar.gz
        tar -zxvf v1.0.0.tar.gz
        rm v1.0.0.tar.gz
        cd mgefinder-1.0.0
    
Install the mgefinder conda environment and other dependencies with

    bash install.sh

Now activate the environment with
    
    source activate mgefinder
    
This activation step must be repeated whenever using *mgefinder*.

Now install *mgefinder* with the command

    pip install .
    
Once complete, you can check to see if *mgefinder* installed properly by simply typing

    mgefinder
   
This can then be called from anywhere on the file system while in the `mgefinder` conda environment.


[NEXT: Step-by-step tutorial](tutorial.md)



