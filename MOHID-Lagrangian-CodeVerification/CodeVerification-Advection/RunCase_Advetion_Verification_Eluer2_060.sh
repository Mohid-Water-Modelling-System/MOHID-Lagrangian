#!/bin/bash

# Resolve paths from this script so it can be launched from any directory.
caseDir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)" || exit 1
repoRoot="$(cd -- "$caseDir/../.." && pwd)" || exit 1
cd -- "$caseDir" || exit 1

clear

# Activate conda environment from shell. Minimum conda version required 4.6
eval "$(conda shell.bash hook)"
conda activate MOHID-Lagrangian

# Increase stacksizes in Linux.
ulimit -s hard
export OMP_STACKSIZE=512M

# "name" and "dirout" are named according to the testcase
name=Advetion_Verification_Eluer2_060
dirout=${name}_out
data=data

# "executables" are renamed and called from their directory
tools="${repoRoot}/build/bin"
mohidlagrangian=${tools}/MOHIDLagrangian

preprocessorDir="${repoRoot}/src/MOHIDLagrangianPreProcessor"
mohidPreprocessor=${preprocessorDir}/MOHIDLagrangianPreProcessor.py

postProcessorDir="${repoRoot}/src/MOHIDLagrangianPostProcessor"
mohidPostprocessor=${postProcessorDir}/MOHIDLagrangianPostProcessor.py

# Check configuration before clearing any previous results.
python analyze_results.py --check-case "${name}.xml"
errcode=$?
if [ $errcode -ne 0 ]; then
  exit "$errcode"
fi

# "dirout" is created to store results or it is cleaned if it already exists
if [ -e "$dirout" ]; then
  rm -f -r "$dirout"
fi
mkdir "$dirout"

cp "${name}.xml" "$dirout/"
cp -r "$data" "$dirout/"

# CODES are executed according the selected parameters of execution in this testcase

python "$mohidPreprocessor" -i "$dirout/${name}.xml" -o "$dirout"
errcode=$?

if [ $errcode -eq 0 ]; then
  "$mohidlagrangian" -i "$dirout/${name}.xml" -o "$dirout"
  errcode=$?
fi

# python -W ignore "$mohidPostprocessor" -i "${name}.xml" -o "$dirout"

if [ $errcode -eq 0 ]; then
  echo All done
else
  echo Execution aborted
fi
read -n1 -r -p "Press any key to continue..." key
echo
