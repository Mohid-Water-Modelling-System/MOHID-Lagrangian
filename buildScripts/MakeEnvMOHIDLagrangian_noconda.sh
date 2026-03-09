#!/bin/bash
# MOHID-Lagrangian Python installer

echo "----------------------------------------------------------------------"
echo "      Python environment installer for MOHID Lagrangian"
echo "----------------------------------------------------------------------"

echo "Installing python packages in a new environment MOHID-Lagrangian..."
python3 -m venv MOHID-Lagrangian
source MOHID-Lagrangian/bin/activate
pip install -r requirements.txt
pip list

echo "Installation complete. Press any key to exit."
read -n 1 -s -r
