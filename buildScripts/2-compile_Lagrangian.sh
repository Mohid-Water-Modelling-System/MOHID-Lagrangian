#!/bin/bash


#########   choose the OPTION of compilation by the number at left #################

                    # 1 - RELEASE_X64
                    # 2 - DEBUG_X64

OPTION=1



### Make the changes to fit your setup ###

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LagrangianMaster="$(cd "$SCRIPT_DIR/.." && pwd)"


# -----------------------------------------------------------------------------
# Find and load Spack
# -----------------------------------------------------------------------------

# First, check whether Spack is already available on the server.
if command -v spack >/dev/null 2>&1; then

    echo "Using existing Spack:"
    echo "    $(command -v spack)"

# Otherwise, check the location used by 0-basics-ExternalLibs.sh.
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
# Find Intel oneAPI compiler through Spack
# -----------------------------------------------------------------------------

inteldir=$(spack location -i --first intel-oneapi-compilers@2023.2.1 2>/dev/null)

if [ -z "$inteldir" ] || [ ! -f "$inteldir/setvars.sh" ]; then

    echo "ERROR: intel-oneapi-compilers@2023.2.1 was not found."
    exit 1

fi

echo "Intel oneAPI compiler:"
echo "    $inteldir"

source "$inteldir/setvars.sh"





build_dir=$LagrangianMaster/build
output_path=$build_dir/bin


#---------------  Setting text colors --------------------

set -e

# echo colors ##
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color
OK=${GREEN}OK${NC}
NOK=${RED}NOK${NC}
ERROR=${RED}ERROR${NC}
WARNING=${RED}warning${NC}
Orange='\033[93m'
Cyan='\033[96m'

#------------------------------------------------------------

START_OF_COMPILE=`date`


cd $LagrangianMaster/

case $OPTION in
            
    1)  # RELEASE_X64

       
        NAMEEXE=MOHIDLagrangian
    
        find . -name \$NAMEEXE -type f -delete
        
        echo -e " ${Orange} "
        echo '        ---------------------------------------------------------'
        echo '                      Compiling ' $NAMEEXE  
        echo '        ---------------------------------------------------------'
        echo -e " ${NC} "
        
        ./buildScripts/MakeMOHIDLagrangian.sh -intel > $SCRIPT_DIR/compile.log
 
        ;;

    2)  # DEBUG_X64
    
        NAMEEXE=MOHIDLagrangian_debug.exe
    
        #find -name $output_path/$NAMEEXE -type f -delete
      
        find . -name \$NAMEEXE -type f -delete
        
        echo -e " ${Orange} "
        echo '        ---------------------------------------------------------'
        echo '                      Compiling ' $NAMEEXE  
        echo '        ---------------------------------------------------------'
        echo -e " ${NC} "
        
        ./buildScripts/MakeMOHIDLagrangian.sh -intel -debug > $SCRIPT_DIR/compile.log
                
        ;;
esac


if [ "$NAMEEXE" != "MOHIDLagrangian" ]; then
  mv "$output_path/MOHIDLagrangian" "$output_path/$NAMEEXE"
fi


#cp -r $build_dir $LagrangianMaster/buildScripts


END_OF_COMPILE=`date`

  if [ ! -f "$output_path/$NAMEEXE" ]; then
    echo -e "                      ${ERROR} $NAMEEXE File not created!"
    echo 
    exit 1
  else
    echo -e "                      compile $NAMEEXE ${OK}  "
    echo
  fi



echo "=========================================================================="
echo "build started:    $START_OF_COMPILE"
echo "build completed:  $END_OF_COMPILE"
echo
echo "--->                  Executables ready                               <---"
echo
#ls -l $output_path/*.exe --color=auto
echo -e "                      $output_path/$NAMEEXE "
echo
echo "=========================================================================="

exit 0