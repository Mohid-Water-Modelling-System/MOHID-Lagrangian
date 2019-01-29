#!/bin/bash

curr_dir=$PWD

base_dir=Base
datetime_dir=datetime-fortran
csv_dir=fortran-csv-module
fox_dir=fox

zlib_dir=zlib-1.2.11
hdf5_dir=HDF5_1.8.17_x64/Linux

echo we are at $curr_dir
wait

cd $hdf5_dir
./configure --with-zlib=$curr_dir/$zlib_dir/build --enable-fortran --enable-fortran2003
make install

echo ------------------------------------
echo All done, have a nice day
wait 