@echo off
cls

set build_dir=build

cd ..

echo ----------------------------------------------------------------------
echo       Generating build solution for MOHID Lagrangian
echo ----------------------------------------------------------------------

rem "build_dir" is created to store results or it is cleaned if it already exists
if exist %build_dir% del /Q %build_dir%\*.*
if not exist %build_dir% mkdir %build_dir%
cd %build_dir%
copy CMakeLists.txt %build_dir%\
rem run cmake to build the libs
rem cmake .. -DCMAKE_GENERATOR_PLATFORM=x64 -DCMAKE_BUILD_TYPE="RELEASE" -DCMAKE_Fortran_COMPILER="C:\\PROGRA~2\\Intel\\oneAPI\\compiler\\latest\\bin\\ifx.exe"
cmake .. -DCMAKE_GENERATOR_PLATFORM=x64 -DCMAKE_BUILD_TYPE="RELEASE" -DCMAKE_GENERATOR_TOOLSET="v143,fortran=ifx"

cd ..

echo Congrats, the requested solution is generated
pause
