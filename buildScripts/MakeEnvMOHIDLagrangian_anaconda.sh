#!/bin/bash

#### title MOHID-Lagrangian python installer

echo '----------------------------------------------------------------------'
echo '      Python environment installer for MOHID Lagrangian               '
echo '---------------------------------------------------------------------- '

echo Installing python packages in a new environment MOHID-Lagrangian...

### Make the changes to fit your setup ###

# Just for HPC computer. If you have miniconda or conda in your path skip this.  
spack load miniconda3@24.3.0

#--------------------------------------------------------
conda create --name MOHID-Lagrangian python=3.11
source activate MOHID-Lagrangian
pip3 install -r requirements.txt 

conda list







