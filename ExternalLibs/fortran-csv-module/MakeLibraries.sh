#!/bin/bash

dirout=build

# "dirout" is created to store results or it is cleaned if it already exists
if [ -e $dirout ]; then
  rm -f -r $dirout
fi
mkdir $dirout

# building initial libs
cd $dirout
echo ------------------------------------
echo Running cmake for fortran-cvs-module library
cmake -Wno-dev .. -DCMAKE_BUILD_TYPE=Release
echo ------------------------------------
echo Running make for fortran-cvs-module library
make

cd ..
cd ..
