@echo off
cls

set build_dir_global=Libs_build_win

set LOCAL_NCF4=C:/Program Files/netCDF 4.6.1
set dir_netcdff=netcdf-fortran
set build_dir_netcdff=build_win

set dir_fox=fox
set build_dir_fox=build_win

echo ----------------------------------------------------------------------
echo       Generating build solutions for nice and small libraries
echo ----------------------------------------------------------------------

rem "build_dir" is created to store results or it is cleaned if it already exists
if exist %build_dir_global% del /Q %build_dir_global%\*.*
if not exist %build_dir_global% mkdir %build_dir_global%
cd %build_dir_global%
rem run cmake to build the libs
cmake .. -DCMAKE_GENERATOR_PLATFORM=x64 -DCMAKE_BUILD_TYPE="RELEASE"
cd ..
if not "%ERRORLEVEL%" == "0" goto fail

echo ----------------------------------------------------------------------
echo            Generating build solution for Fox xml parser
echo ----------------------------------------------------------------------

cd %dir_fox%
if exist %build_dir_fox% del /Q %build_dir_fox%\*.*
if not exist %build_dir_fox% mkdir %build_dir_fox%
cd %build_dir_fox%
rem run cmake for Fox
cmake .. -DCMAKE_GENERATOR_PLATFORM=x64 -DCMAKE_BUILD_TYPE="RELEASE" 
cd ..
cd ..

echo Congrats, the solutions are generated. You should now build them
pause