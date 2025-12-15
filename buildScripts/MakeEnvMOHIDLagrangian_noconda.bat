@echo off
title MOHID-Lagrangian python installer

echo ----------------------------------------------------------------------
echo       Python environment installer for MOHID Lagrangian
echo ----------------------------------------------------------------------


echo Installing python packages in a new environment MOHID-Lagrangian...
call python -m venv MOHID-Lagrangian
call MOHID-Lagrangian\Scripts\activate
call pip install -r requirements.txt
call pip list


pause