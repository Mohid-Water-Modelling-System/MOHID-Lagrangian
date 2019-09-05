@echo off
cls

rem "name" and "dirout" are named according to the case

set name=Arousa2D_case
set dirout=%name%_out

rem "executables" are renamed and called from their directory

set tools=../../build/bin/RELEASE
set mohidlagrangian="%tools%/MOHIDLagrangian.exe"

set preprocessorDir=../../src/MOHIDLagrangianPreProcessor
set PreProcessor="%preprocessorDir%/MOHIDLagrangianPreProcessor.py"

set postProcessorDir=../../../src\MOHIDLagrangianPostProcessor
set postProcessor="%postProcessorDir%/MOHIDLagrangianPostprocessor.py"

cd %dirout%
python %postProcessor% %name%.xml concentrations residence_time

:success
echo All done
goto end

:fail
echo Execution aborted.

:end
pause