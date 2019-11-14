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


#pyv="$(python -V 2>&1)"
#ver=${pyv:7:1}
#if [ $ver != "3" ];
#then 
#  echo 'Python 3 is required to use MOHIDLagrangian Software'
#  echo 'Please, install it' 
#  exit 1
#fi

#echo 'Installing required python packages'
#if ! type conda &> /dev/null;
#then
#     conda install netcdf4 xarray numba vtk
#elif ! type pip &> /dev/null;
#then
#     pip install netcdf4 xarray numba vtk
#else;
#then
#    echo 'You do not have <conda> or <pip> installed'
#    echo 'Please install it before continue'
#    exit 1
#fi
