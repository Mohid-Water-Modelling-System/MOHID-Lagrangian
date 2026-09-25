#!/bin/bash

### Make the changes to fit your setup ###


# Directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# MOHID-Lagrangian root directory
LagrangianMaster="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "MOHID-Lagrangian root:"
echo "    $LagrangianMaster"

# -------------------------------------------------------------------------
# Keep the tracked ExternalLib files available after cloning, but prevent
# local modifications from appearing in normal Git status output.
# New files inside ExternalLib remain ignored through .gitignore.
git ls-files -z ExternalLibs/ |   xargs -0 -r git update-index --no-skip-worktree --
# -------------------------------------------------------------------------

if git -C "$LagrangianMaster" rev-parse --is-inside-work-tree >/dev/null 2>&1; then

    git -C "$LagrangianMaster" ls-files -z ExternalLibs/ | \
        xargs -0 -r git -C "$LagrangianMaster" update-index --skip-worktree --

fi

# -------------------------------------------------------------------------
#If a user later wants to modify ExternalLibs and commit those changes, they just need
# to remove the skip-worktree flag first.
#To tack  whole folder:
#git ls-files -z ExternalLibs/ |   xargs -0 -r git update-index --no-skip-worktree --
#To tack  one file:
#git update-index --no-skip-worktree ExternalLibs/path/to/file
# -------------------------------------------------------------------------
#
#Normal installation:
#    skip-worktree ON
#    |
#ExternalLibs local changes are hidden
#
#Developer wants to modify ExternalLibs:
#    --no-skip-worktree
#    |
#git add / commit / push
#    |
#optional: --skip-worktree again
#
# -------------------------------------------------------------------------

folders=("MOHID-Lagrangian" "build" "src" "ExternalLibs")

for folder in "${folders[@]}"; do

    if [ -d "$folder" ]; then

        echo "Removing $folder ..."
        nohup rm -rf "$folder" > /dev/null 2>&1 &

    fi

done

wait

echo "Done."


# -----------------------------------------------------------------------------
# Find and load Spack
# -----------------------------------------------------------------------------

if command -v spack >/dev/null 2>&1; then

    echo "Using existing Spack:"
    echo "    $(command -v spack)"

elif [ -f "/home/software/spack/share/spack/setup-env.sh" ]; then

    echo "Using Spack installed by 0-basics-ExternalLibs.sh:"
    echo "    /home/software/spack"

    source /home/software/spack/share/spack/setup-env.sh

else

    echo "ERROR: Spack was not found."
    echo
    echo "Load the Spack installation available on this server"
    echo "or run 0-basics-ExternalLibs.sh first."

    exit 1

fi


# -----------------------------------------------------------------------------
# Find Intel oneAPI compiler
# -----------------------------------------------------------------------------

inteldir=$(spack location -i --first intel-oneapi-compilers@2023.2.1 2>/dev/null)


if [ -z "$inteldir" ] || [ ! -f "$inteldir/setvars.sh" ]; then

    echo "ERROR: intel-oneapi-compilers@2023.2.1 was not found."
    echo
    echo "Run 0-basics-ExternalLibs.sh first."

    exit 1

fi


echo "Intel oneAPI compiler:"
echo "    $inteldir"


source "$inteldir/setvars.sh"




# Regenerate Proj4 Autotools files
PROJ4_DIR="$LagrangianMaster/ExternalLibs/Proj4/Linux/proj-4.9.3"

if [ -f "$PROJ4_DIR/configure.ac" ]; then

    echo "Regenerating Proj4 build files..."

    cd "$PROJ4_DIR" || exit 1

    autoreconf -fi || exit 1

    rm -f config.status config.log
    find . -name Makefile -type f -delete

fi


cd "$LagrangianMaster/ExternalLibs/" || exit 1

find . -type f -name "*.sh" -exec chmod +x {} +
find . -type f -name "configure" -exec chmod +x {} +
find . -type f -name "mkdirs" -exec chmod +x {} +
find . -type f -name "install-sh" -exec chmod +x {} +

./MakeLibraries.sh -intel




# Return to the buildScripts directory
#cd "$SCRIPT_DIR"

# Optional: copy ExternalLibs and src into buildScripts
#cp -r "$LagrangianMaster/ExternalLibs" "$SCRIPT_DIR"
#cp -r "$LagrangianMaster/src" "$SCRIPT_DIR"