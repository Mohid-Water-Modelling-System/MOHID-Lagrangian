@echo off
cls

set build_dir_global=Libs_build_win

set LOCAL_NCF4=C:/Program Files/netCDF 4.6.1
set dir_netcdff=netcdf-fortran
set build_dir_netcdff=build_win

set dir_fox=fox
set build_dir_fox=build_win

echo ----------------------------------------------------------------------
echo   Generating build solution for NetCDF-Fortran. Cross your fingers.
echo   If you run into problems check that the paths to the NetCDF install 
echo   are correct.
echo ----------------------------------------------------------------------

if exist %build_dir_netcdff% del /Q %build_dir_netcdff%\*.*
if not exist %build_dir_netcdff% mkdir %build_dir_netcdff%
cd %build_dir_netcdff%
rem run cmake for NetCDF Fortran
cmake .. -DCMAKE_GENERATOR_PLATFORM=x64 -DCMAKE_BUILD_TYPE="RELEASE" -DBUILD_DAP:BOOL="0" -DENABLE_NETCDF4:BOOL="0" -DBUILD_PARALLEL:STRING="ON" -DUSE_NETCDF4:BOOL="0" -DBUILD_SHARED_LIBS="OFF" -DBUILD_TESTING:BOOL="1" -DBUILD_EXAMPLES:BOOL="1" -DENABLE_TESTS:BOOL="1" -DCMAKE_PREFIX_PATH:FILEPATH="%LOCAL_NCF4%" 
cd ..

echo Congrats, the solutions are generated. You should now build them
pause
