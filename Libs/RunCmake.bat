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
echo   Generating build solution for NetCDF-Fortran. Cross your fingers.
echo   If you run into problems check that the paths to the NetCDF install 
echo   are correct.
echo ----------------------------------------------------------------------

cd %dir_netcdff%
if exist %build_dir_netcdff% del /Q %build_dir_netcdff%\*.*
if not exist %build_dir_netcdff% mkdir %build_dir_netcdff%
cd %build_dir_netcdff%
rem run cmake for NetCDF Fortran
cmake .. -DCMAKE_GENERATOR_PLATFORM=x64 -DCMAKE_BUILD_TYPE="RELEASE" -DBUILD_SHARED_LIBS="OFF" -DCMAKE_PREFIX_PATH:FILEPATH="%LOCAL_NCF4%" 
cd ..
cd ..

echo ----------------------------------------------------------------------
echo   Generating build solution for Fox xml parser. Cross your fingers.
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
