#!/bin/bash

# "name" and "dirout" are named according to the testcase

name=Libs_build
dirout=${name}_linux

# "dirout" is created to store results or it is cleaned if it already exists
if [ -e $dirout ]; then
  rm -f -r $dirout
fi
mkdir $dirout

echo ------------------------------------
echo Build directory is $dirout

# building initial libs
cd $dirout
diroutpwd=$PWD
echo ------------------------------------
echo Running cmake for nice and small libraries
cmake -Wno-dev ..
echo Running make for nice and small libraries
make

# building complex libs
cd ../netcdf-fortran
diroutnetcdff=build
# "diroutnetcdff" is created to store results or it is cleaned if it already exists
if [ -e $diroutnetcdff ]; then
  rm -f -r $diroutnetcdff
fi
mkdir $diroutnetcdff
cd $diroutnetcdff
echo ------------------------------------
echo Running cmake for netCDF-Fortran library, hope you installed netCDF in your system
cmake -Wno-dev .. -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=OFF -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTING=OFF -DENABLE_TESTS=OFF
echo Running make for netCDF-Fortran library
make
echo copying files from built netCDF-Fortran library to $dirout
find . -iname "*.mod" -exec cp -R '{}' ${diroutpwd}/modules/ \;
find . -iname "*.a" -exec cp -R '{}' ${diroutpwd}/lib/ \;

cd ../../fox
diroutfox=build
# "diroutfox" is created to store results or it is cleaned if it already exists
if [ -e $diroutfox ]; then
  rm -f -r $diroutfox
fi
mkdir $diroutfox
cd $diroutfox
echo ------------------------------------
echo Running cmake for fox library
cmake -Wno-dev ..
echo Running make for fox library
make
echo copying files from built fox library to $dirout
cp -R modules/. ${diroutpwd}/modules/
cp -R lib/. ${diroutpwd}/lib/

echo ------------------------------------
echo All done, have a nice day