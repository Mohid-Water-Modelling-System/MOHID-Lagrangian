#!/bin/bash

# "name" and "dirout" are named according to the testcase

name=build
dirout=${name}_linux

# "dirout" is created to store results or it is cleaned if it already exists
if [ -e $dirout ]; then
  rm -f -r $dirout
fi
mkdir $dirout

# building initial libs
cd $dirout
echo ------------------------------------
echo Running cmake for small base libraries
cmake -Wno-dev .. -DCMAKE_BUILD_TYPE=Release
echo ------------------------------------
echo Running make for nice and small libraries
make

cd ..
cd ..
