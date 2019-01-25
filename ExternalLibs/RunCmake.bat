@echo off
cls

set base_dir=Base
set datetime_dir=datetime-fortran
set csv_dir=fortran-csv-module
set fox_dir=fox

cd %base_dir%
call RunCmake.bat
cd ..

cd %datetime_dir%
call RunCmake.bat
cd ..

cd %csv_dir%
call RunCmake.bat
cd ..

cd %fox_dir%
call RunCmake.bat
cd ..

echo Congrats, the solutions are generated. You should now build them
pause
