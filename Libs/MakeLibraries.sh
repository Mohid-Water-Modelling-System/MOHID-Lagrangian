#!/bin/bash

# "name" and "dirout" are named according to the testcase

name=Libs_build
dirout=${name}_linux

# "dirout" is created to store results or it is cleaned if it already exists
if [ -e $dirout ]; then
  rm -f -r $dirout
fi
mkdir $dirout

echo Build directory is $dirout

# building initial libs
cd $dirout
diroutpwd=$PWD
echo Running cmake for nice and small libraries
cmake -Wno-dev ..
echo Running make for nice and small libraries
make

# building complex libs
cd ../fox
diroutfox=build
# "diroutfox" is created to store results or it is cleaned if it already exists
if [ -e $diroutfox ]; then
  rm -f -r $diroutfox
fi
mkdir $diroutfox
cd $diroutfox
echo Running cmake for fox library
cmake -Wno-dev ..
echo Running make for fox library
make
echo coppying files from built fox library to $dirout
cp -R modules/. ${diroutpwd}/modules/
cp -R lib/. ${diroutpwd}/lib/