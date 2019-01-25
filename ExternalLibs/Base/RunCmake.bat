
set build_dir_global=build_win

echo ----------------------------------------------------------------------
echo       Generating build solutions for small base libraries
echo ----------------------------------------------------------------------

rem "build_dir" is created to store results or it is cleaned if it already exists
if exist %build_dir_global% del /Q %build_dir_global%\*.*
if not exist %build_dir_global% mkdir %build_dir_global%
cd %build_dir_global%
rem run cmake to build the libs
cmake .. -DCMAKE_GENERATOR_PLATFORM=x64 -DCMAKE_BUILD_TYPE="RELEASE" -Wno-dev
cd ..
if not "%ERRORLEVEL%" == "0" goto fail

