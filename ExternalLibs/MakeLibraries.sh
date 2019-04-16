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


if [ $# -eq 1 ]
then
        if [ $1='-intel' ]
        then
                export CFLAGS='-O3 -xHost -ip -fPIC'
                export CXXFLAGS='-O3 -xHost -ip'
                export FCFLAGS='-O3 -xHost -ip'
                export FFLAGS='-O3 -xHost -ip'
                export CC=/opt/intel/bin/icc
                export FC=/opt/intel/bin/ifort
                if [ $# -eq 2 ]
                then
                        export CC=$2
                fi
                if [ $# -eq 3 ]
                then
                        export FC=$3
                fi
        fi
	
	if [ $1 = '-gfortran' ]
	then
		export CFLAGS='-03 -fPIC'
		export CXXFLAGS='-03'
		export FCFLAGS='-03'
		export FFLAGS='-03'
		export CC=/usr/bin/gcc
		export FC=/usr/bin/gfortran
		if [ $# -eq 2 ]
		then
			export CC=$2
		fi
		if [ $# -eq 3 ] 
		then
			export FC=$3
		fi
	fi
	
fi

#---------------------------DO NOT CHANGE------------------------------------------
cd $base_dir
rm -rf /build
./MakeLibraries.sh
cd $root_libs_dir

cd $datetime_dir
rm -rf /build
./MakeLibraries.sh
cd $root_libs_dir

cd $csv_dir
rm -rf /build
./MakeLibraries.sh
cd $root_libs_dir

cd $fox_dir
rm -rf /build
./MakeLibraries.sh
cd $root_libs_dir

cd $zlib_dir
make clean
mkdir build
./configure --prefix=$root_libs_dir/$zlib_dir/build || exit
make install || exit
cd $root_libs_dir

cd $hdf5_dir
make clean
./configure --with-zlib=$root_libs_dir/$zlib_dir/build --enable-fortran --enable-fortran2003 --disable-shared || exit
make install || exit
cd $root_libs_dir

cd $netcdfc_dir
make clean
mkdir build
export CPPFLAGS="-I$root_libs_dir/$hdf5_dir/hdf5/include -I$root_libs_dir/$zlib_dir/build/include" 
export LDFLAGS="-L$root_libs_dir/$hdf5_dir/hdf5/lib -L$root_libs_dir/$zlib_dir/build/lib"  
./configure --enable-netcdf-4 --prefix=$root_libs_dir/$netcdfc_dir/build || exit
make install || exit
cd $root_libs_dir

cd $netcdff_dir
make clean
export CPPFLAGS="-I$root_libs_dir/$hdf5_dir/hdf5/include -I$root_libs_dir/$zlib_dir/build/include"
export LDFLAGS="-L$root_libs_dir/$hdf5_dir/hdf5/lib -L$root_libs_dir/$zlib_dir/build/lib"
export NETCDF=$root_libs_dir/$netcdfc_dir/build
export PATH=$PATH:$NETCDF/bin:$NETCDF/lib:$NETCDF/include
export NETCDF_ROOT=$NETCDF
export NETCDF4_ROOT=$NETCDF
export NETCDF_LIB=$NETCDF/lib
export NETCDF_INC=$NETCDF/include
export NETCDF_GF_ROOT=$NETCDF
export NETCDF4_GF_ROOT=$NETCDF
export NETCDF_GF_LIB=$NETCDF/lib
export NETCDF_GF_INC=$NETCDF/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDF_LIB
export CPPFLAGS="$CPPFLAGS -I$NETCDF_INC"
export LDFLAGS="$LDFLAGS -L$NETCDF_LIB"
mkdir build
./configure --prefix=$root_libs_dir/$netcdff_dir/build || exit
make install || exit
cd $root_libs_dir

cd $proj4C_dir
make clean
mkdir build
./configure --prefix=$root_libs_dir/$proj4C_dir/build  || exit
make install || exit
ln -sf $root_libs_dir/$proj4C_dir/build/lib/libproj.so $root_libs_dir/$proj4C_dir/build/lib/libproj4.so
cd $root_libs_dir

cd $proj4F_dir
make clean
mkdir build
cp cgfortran.h cfortran.h
if [ $# -eq 1 ]
then
        if [ $1='-intel' ]
        then
        	cp cintelfortran.h cfortran.h
        
	fi
fi
./bootstrap
./configure --with-proj4=$root_libs_dir/$proj4C_dir/build --prefix=$root_libs_dir/$proj4F_dir/build || exit
make install || exit
cd $root_libs_dir

cd $mohid_dir
rm -rf /build
./MakeLibraries.sh
cd $root_libs_dir


echo ------------------------------------
echo All done, have a nice day
wait 
