#!/bin/bash

root_libs_dir=$PWD

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
cd $root_libs_dir

cd $datetime_dir
./MakeLibraries.sh
cd $root_libs_dir

cd $csv_dir
./MakeLibraries.sh
cd $root_libs_dir

cd $fox_dir
./MakeLibraries.sh
cd $root_libs_dir

cd $zlib_dir
mkdir build
./configure --prefix=$root_libs_dir/$zlib_dir/build || exit
make install || exit
cd $root_libs_dir

cd $hdf5_dir
./configure --with-zlib=$root_libs_dir/$zlib_dir/build --enable-fortran --enable-fortran2003 --disable-shared || exit
make install || exit
cd $root_libs_dir

cd $netcdfc_dir
mkdir build
export CPPFLAGS="-I$hdf5_dir/hdf5/include -I$zlib_dir/build/include" 
export LDFLAGS="-L$hdf5_dir/hdf5/lib -L$zlib_dir/build/lib"  
./configure --enable-netcdf-4 --prefix=$pwd/build || exit
make install || exit
cd $root_libs_dir

cd $netcdff_dir
mkdir build
./configure --prefix=$root_libs_dir/$netcdff_dir/build || exit
make install || exit
cd $root_libs_dir

cd $proj4C_dir
mkdir build
./configure --prefix=$root_libs_dir/$proj4C_dir/build || exit
make install || exit
cd $root_libs_dir

cd $proj4F_dir
mkdir build
./bootstrap
./configure --with-proj4=$root_libs_dir/$proj4C_dir/build --prefix=$root_libs_dir/$proj4F_dir/build CFLAGS=f2cFortran || exit
make install || exit
cd $root_libs_dir

cd $mohid_dir
./MakeLibraries.sh
cd $root_libs_dir


echo ------------------------------------
echo All done, have a nice day
wait 