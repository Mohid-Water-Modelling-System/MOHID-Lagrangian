
set build_dir_fox=build_win

echo ----------------------------------------------------------------------
echo           Generating build solution for Fox xml parser
echo ----------------------------------------------------------------------

if exist %build_dir_fox% del /Q %build_dir_fox%\*.*
if not exist %build_dir_fox% mkdir %build_dir_fox%
cd %build_dir_fox%
rem run cmake for Fox
cmake .. -DCMAKE_GENERATOR_PLATFORM=x64 -DCMAKE_BUILD_TYPE="RELEASE" 
cd ..
cd ..

