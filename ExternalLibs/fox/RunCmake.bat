@echo off

set build_dir_fox=build

echo ----------------------------------------------------------------------
echo           Generating build solution for Fox xml parser
echo ----------------------------------------------------------------------

rem "build_dir" is created to store results or it is cleaned if it already exists
if exist %build_dir_fox% del /Q %build_dir_fox%\*.*
if not exist %build_dir_fox% mkdir %build_dir_fox%
cd %build_dir_fox%
rem run cmake for Fox
cmake .. -DCMAKE_GENERATOR_PLATFORM=x64 -DCMAKE_BUILD_TYPE="RELEASE" -DCMAKE_GENERATOR_TOOLSET="v143,fortran=ifx"
cd ..
if not "%ERRORLEVEL%" == "0" goto fail

