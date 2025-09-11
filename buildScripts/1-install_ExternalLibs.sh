#!/bin/bash

### Make the changes to fit your setup ###


inteldir=/home/software/spack/opt/spack/linux-ubuntu24.04-sapphirerapids/gcc-13.3.0/intel-oneapi-compilers-2023.2.1-cz7grxfezyuoifpdwhgn6dcdbezbydez
source $inteldir/setvars.sh

LagrangianMaster=~/lagrangian/MOHID-Lagrangian

root_dir=$PWD

cd $LagrangianMaster/ExternalLibs/

./MakeLibraries.sh -intel 


cd $root_dir
cp -r $LagrangianMaster/ExternalLibs $root_dir
cp -r $LagrangianMaster/src $root_dir


