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
rem run cmake to build the libs
cmake .. -DCMAKE_GENERATOR_PLATFORM=x64 -DCMAKE_BUILD_TYPE="RELEASE"
cd ..

echo Congrats, the requested solution is generated
pause
