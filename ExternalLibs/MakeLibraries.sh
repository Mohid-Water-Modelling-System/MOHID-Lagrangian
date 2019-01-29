#!/bin/bash

curr_dir=$PWD

base_dir=Base
datetime_dir=datetime-fortran
csv_dir=fortran-csv-module
fox_dir=fox
zlib_dir=zlib-1.2.11
hdf5_dir=HDF5_1.8.17_x64/Linux
netcdfc_dir=Netcdf_4.4.1.1/Linux
netcdff_dir=netcdf-fortran-4.4.4

cd $base_dir
./MakeLibraries.sh
cd ..

cd $datetime_dir
./MakeLibraries.sh
cd ..

cd $csv_dir
./MakeLibraries.sh
cd ..

cd $fox_dir
./MakeLibraries.sh
cd ..

cd $zlib_dir
mkdir build
./configure --prefix=build
make install
cd ..

cd $hdf5_dir
./configure --with-zlib=$curr_dir/$zlib_dir/build --enable-fortran --enable-fortran2003
make install
cd ..

cd $netcdfc_dir
mkdir build
./configure --prefix=$curr_dir/$netcdfc_dir/build
make install
cd ..

cd $netcdff_dir
mkdir build
./configure --prefix=$curr_dir/$netcdff_dir/build
make install
cd ..

echo ------------------------------------
echo All done, have a nice day
wait 