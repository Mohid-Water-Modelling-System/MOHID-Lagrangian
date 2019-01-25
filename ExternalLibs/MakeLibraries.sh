#!/bin/bash

base_dir=Base
datetime_dir=datetime-fortran
csv_dir=fortran-csv-module
fox_dir=fox

cd $base_dir
MakeLibraries.sh
cd ..

cd $datetime_dir
MakeLibraries.sh
cd ..

cd $csv_dir
MakeLibraries.sh
cd ..

cd $fox_dir
MakeLibraries.sh
cd ..

echo ------------------------------------
echo All done, have a nice day
wait 