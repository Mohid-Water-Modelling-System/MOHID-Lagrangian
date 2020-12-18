@echo off
cls

call activate MOHID-Lagrangian

rem "name" and "dirout" are named according to the case

set name=Arousa2D_case
set dirout=%name%_out

rem "executables" are renamed and called from their directory

set postProcessorDir=../../src\MOHIDLagrangianPostProcessor
set postProcessor="%postProcessorDir%/MOHIDLagrangianPostprocessor.py"

python -W ignore %postProcessor% -i %name%.xml -o %dirout% 

:success
echo All done
goto end

:fail
echo Execution aborted.

:end
pause
