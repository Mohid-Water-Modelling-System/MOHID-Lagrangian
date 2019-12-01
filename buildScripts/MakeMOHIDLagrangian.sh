#!/bin/bash

dirout=build

cd ..

# "dirout" is created to store results or it is cleaned if it already exists
if [ -e $dirout ]; then
  rm -f -r $dirout
fi
mkdir $dirout

echo Build directory is $dirout

# building initial libs
cd $dirout
echo Running cmake for MOHID Lagrangian $1
if [[ $1 == "Intel" ]]; then
cmake -Wno-dev -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_Fortran_COMPILER_ID="Intel" ..
else
echo Running cmake for MOHID Lagrangian
cmake -Wno-dev ..
fi
make