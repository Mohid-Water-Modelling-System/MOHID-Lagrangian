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
proj4Base_dir=Proj4/Linux
proj4C_dir=$proj4Base_dir/proj-4.9.3
proj4F_dir=$proj4Base_dir/proj4-fortran
mohid_dir=MOHID

cd $base_dir
./MakeLibraries.sh
cd $curr_dir

cd $datetime_dir
./MakeLibraries.sh
cd $curr_dir

cd $csv_dir
./MakeLibraries.sh
cd $curr_dir

cd $fox_dir
./MakeLibraries.sh
cd $curr_dir

cd $zlib_dir
mkdir build
./configure --prefix=build
make install
cd $curr_dir

cd $hdf5_dir
./configure --with-zlib=$curr_dir/$zlib_dir/build --enable-fortran --enable-fortran2003
make install
cd $curr_dir

cd $netcdfc_dir
mkdir build
./configure --prefix=$curr_dir/$netcdfc_dir/build
make install
cd $curr_dir

cd $netcdff_dir
mkdir build
./configure --prefix=$curr_dir/$netcdff_dir/build
make install
cd $curr_dir

cd $proj4C_dir
mkdir build
./configure --prefix=$curr_dir/$proj4C_dir/build
make install
cd $curr_dir

cd $proj4F_dir
mkdir build
./bootstrap
./configure --with-proj4=$curr_dir/$proj4C_dir/build --prefix=$curr_dir/$proj4F_dir/build CFLAGS=f2cFortran
make install
cd $curr_dir

cd $mohid_dir
./MakeLibraries.sh
cd $curr_dir


echo ------------------------------------
echo All done, have a nice day
wait 