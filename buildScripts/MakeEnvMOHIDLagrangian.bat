@echo off
title MOHID-Lagrangian python installer

echo ----------------------------------------------------------------------
echo       Python packages installer for MOHID Lagrangian
echo ----------------------------------------------------------------------

where conda >nul 2>nul
if %ERRORLEVEL% NEQ 0 ( 
	echo Anaconda is not installed or is not in the PATH!
	echo Install it or add it to your path.
) else (
	echo Installing python packages in a new enironment MOHID-Lagrangian...
	conda create --name MOHID-Lagrangian python=3.6 --file requirements.txt
	conda activate MOHID-Lagrangian
)

pause