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
echo Running cmake for MOHID Lagrangian
cmake -Wno-dev ..
echo Running make for MOHID Lagrangian
make