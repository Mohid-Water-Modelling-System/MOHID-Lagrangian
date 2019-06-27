@echo off
cls

rem "name" and "dirout" are named according to the case

set name=Arousa2D
set dirout=%name%_out

rem "executables" are renamed and called from their directory

set tools=../../build/bin/RELEASE
set mohidlagrangian="%tools%/MOHIDLagrangian.exe"

set preprocessorDir=../../src/MOHIDLagrangianPreProcessor
set PreProcessor="%preprocessorDir%/MOHIDLagrangianPreProcessor.py"

rem "dirout" is created to store results or it is cleaned if it already exists
if exist %dirout% del /Q %dirout%\*.*
if not exist %dirout% mkdir %dirout%

copy %name%_Def.xml %dirout% 
ren %dirout%\%name%_Def.xml %name%.xml

rem CODES are executed according the selected parameters of execution in this case

python %PreProcessor% -i %dirout%/%name%.xml -o %dirout%

%mohidlagrangian% -i %dirout%/%name%.xml -o %dirout%
if not "%ERRORLEVEL%" == "0" goto fail

:success
echo All done
goto end

:fail
echo Execution aborted.

:end
pause
