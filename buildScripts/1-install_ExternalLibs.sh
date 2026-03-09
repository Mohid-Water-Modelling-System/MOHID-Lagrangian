#!/bin/bash

### Make the changes to fit your setup ###

#!/bin/bash

#!/bin/bash

# List of folders to remove
folders=("MOHID-Lagrangian" "build" "src" "ExternalLibs")

for folder in "${folders[@]}"; do
    if [ -d "$folder" ]; then
        echo "Removing folder: $folder"
        
        # Try with nohup first
        nohup rm -rf "$folder" 2>/dev/null &
        wait $!  # Wait for the nohup process to finish

        # Check if folder still exists
        if [ -d "$folder" ]; then
            echo "nohup failed, trying normal rm..."
            rm -rf "$folder"
        fi
    else
        echo "Folder $folder does not exist, skipping."
    fi
done

echo "Done."


inteldir=/home/software/spack/opt/spack/linux-ubuntu24.04-sapphirerapids/gcc-13.3.0/intel-oneapi-compilers-2023.2.1-cz7grxfezyuoifpdwhgn6dcdbezbydez
source $inteldir/setvars.sh

LagrangianMaster=~/lagrangian/MOHID-Lagrangian

root_dir=$PWD

cd $LagrangianMaster/ExternalLibs/

./MakeLibraries.sh -intel 


cd $root_dir
cp -r $LagrangianMaster/ExternalLibs $root_dir
cp -r $LagrangianMaster/src $root_dir


