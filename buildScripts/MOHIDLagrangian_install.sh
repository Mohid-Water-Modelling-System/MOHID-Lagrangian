#!/bin/bash

#==============================================================================
#     MULTI-OS INSTALLER FOR MOHID-LAGRANGIAN (UBUNTU & ROCKY LINUX)
#==============================================================================
# This script installs MOHID-Lagrangian and all its dependencies.
# Compatible with:
#   - Debian/Ubuntu based systems (apt)
#   - Rocky Linux/RHEL/CentOS based systems (dnf/yum)
#==============================================================================

# --- Error handling: the script will stop if a command fails ---
set -e

#==============================================================================
#     CONFIGURATION (MODIFY IF NECESSARY)
#==============================================================================

# 1. INSTALLATION DIRECTORY NAME
#    Defines the folder name inside your home directory.
#    - Default: "lagrangian" -> installs to $HOME/lagrangian
INSTALL_FOLDER_NAME="lagrangian"

# 2. GITHUB BRANCH
#    Define the branch you want to clone.
#    - Default: "master"
TARGET_BRANCH="master"

# 3. COMPILE EXTERNAL LIBRARIES?
#    Set to true to compile ExternalLibs (required for first installation).
#    Set to false to skip this step (if they are already compiled).
COMPILE_LIBS=true

# 4. BUILD TYPE
#    Options: "RELEASE" or "DEBUG"
BUILD_TYPE="RELEASE"

# 5. TOOL VERSIONS (Spack)
INTEL_COMPILER_VERSION="intel-oneapi-compilers@2023.2.1"
MINICONDA_VERSION="miniconda3@24.3.0"

# 6. PROJECT FOLDER NAME (Internal repo folder name)
PROJECT_DIR_NAME="MOHID-Lagrangian"

#==============================================================================
#     PATH SETUP (AUTOMATIC)
#==============================================================================

# Construct the absolute path based on the folder name provided above
LAGRANGIAN_BASE_DIR="$HOME/$INSTALL_FOLDER_NAME"

echo "Configuration Loaded:"
echo "-> Path:   $LAGRANGIAN_BASE_DIR"
echo "-> Branch: $TARGET_BRANCH"
echo "-> Libs:   $COMPILE_LIBS"
echo "--------------------------------------------------------"

#==============================================================================
#     HELPER FUNCTIONS
#==============================================================================

# --- Function to log messages to the terminal ---
log_message() {
    echo ""
    echo "******************************************************************************"
    echo "✓ $(date '+%Y-%m-%d %H:%M:%S') - $1"
    echo "******************************************************************************"
}

# --- Function to check if a Spack package is installed ---
is_spack_installed() {
    spack find "$1" | grep -q "$1"
    return $?
}

# --- Function to install a Spack package if it doesn't exist ---
install_spack_package() {
    if is_spack_installed "$1"; then
        log_message "Spack package '$1' already installed. Skipping..."
    else
        log_message "Installing Spack package '$1'..."
        spack install "$1"
    fi
    spack load "$1"
    log_message "Spack package '$1' loaded."
}

#==============================================================================
#     MAIN INSTALLATION FUNCTIONS
#==============================================================================

# --- 1. Install system dependencies and Spack (OS Adapted) ---
install_base_dependencies() {
    
    # OS Detection
    if [ -f /etc/os-release ]; then
        . /etc/os-release
        OS_ID=$ID
        OS_LIKE=$ID_LIKE
    else
        log_message "WARNING: Cannot detect OS. Assuming dependencies are met."
        OS_ID="unknown"
    fi

    log_message "Detected OS: $OS_ID ($OS_LIKE)"

    if [[ "$OS_ID" == "ubuntu" || "$OS_ID" == "debian" || "$OS_LIKE" == *"debian"* ]]; then
        # --- UBUNTU / DEBIAN ---
        log_message "Installing system dependencies using apt (Ubuntu/Debian)..."
        sudo apt-get update
        sudo apt-get install -y git vim python3 python3-pip python3-venv build-essential

    elif [[ "$OS_ID" == "rocky" || "$OS_ID" == "centos" || "$OS_ID" == "rhel" || "$OS_LIKE" == *"rhel"* || "$OS_LIKE" == *"fedora"* ]]; then
        # --- ROCKY LINUX / RHEL ---
        log_message "Installing system dependencies using dnf (Rocky/RHEL)..."
        
        # Install EPEL release if not present (sometimes needed for specific packages)
        if ! rpm -q epel-release > /dev/null 2>&1; then
             log_message "Installing EPEL repository..."
             sudo dnf install -y epel-release
        fi
        
        sudo dnf update -y
        # 'Development Tools' is the RHEL equivalent of 'build-essential'
        sudo dnf groupinstall -y "Development Tools"
        sudo dnf install -y git vim python3 python3-pip which

    else
        log_message "WARNING: Unsupported OS family. Skipping system package installation."
    fi

    # --- Common Python Dependencies (Clingo) ---
    # Install Clingo if not present
    if ! pip3 list 2>/dev/null | grep -q "clingo"; then
        log_message "Installing Clingo..."
        pip3 install --user clingo
    else
        log_message "Clingo is already installed. Skipping..."
    fi

    # --- Spack Setup ---
    log_message "Configuring Spack..."
    if [ ! -d "$HOME/spack" ]; then
        cd "$HOME"
        git clone -c feature.manyFiles=true https://github.com/spack/spack.git
    else
        log_message "Spack is already cloned. Skipping..."
    fi

    # Add Spack to the shell environment
    if ! grep -q "spack/setup-env.sh" ~/.bashrc; then
        echo ". $HOME/spack/share/spack/setup-env.sh" >> ~/.bashrc
        log_message "Spack added to .bashrc."
    fi
    
    # Source spack for current session
    source "$HOME/spack/share/spack/setup-env.sh"

    log_message "Installing tools via Spack..."
    install_spack_package "$INTEL_COMPILER_VERSION"
    install_spack_package "$MINICONDA_VERSION"
    install_spack_package "cmake"
    install_spack_package "automake"
    install_spack_package "autoconf"
    install_spack_package "libtool"
    install_spack_package "m4"
    install_spack_package "perl"
}

# --- 2. Clone the MOHID-Lagrangian repository ---
clone_repository() {
    log_message "Cloning the MOHID-Lagrangian repository..."
    
    # Ensure base directory exists
    if [ ! -d "$LAGRANGIAN_BASE_DIR" ]; then
        mkdir -p "$LAGRANGIAN_BASE_DIR"
    fi
    
    cd "$LAGRANGIAN_BASE_DIR"

    if [ ! -d "$PROJECT_DIR_NAME" ]; then
        # Uses the variable TARGET_BRANCH set in configuration
        log_message "Cloning branch: $TARGET_BRANCH"
        git clone -b "$TARGET_BRANCH" https://github.com/Mohid-Water-Modelling-System/MOHID-Lagrangian.git
    else
        log_message "Project directory already exists. Skipping clone."
        # Optional: Log current branch for verification
        cd "$PROJECT_DIR_NAME"
        current_branch=$(git rev-parse --abbrev-ref HEAD)
        log_message "Existing repo is on branch: $current_branch"
        cd ..
    fi
}

# --- 3. Apply necessary code modifications ---
apply_modifications() {
    local project_path="$LAGRANGIAN_BASE_DIR/$PROJECT_DIR_NAME"
    log_message "Applying code modifications in $project_path..."

    mkdir -p "$project_path/src/MOHIDLagrangianPostProcessor"
    mkdir -p "$project_path/src/MOHIDLagrangianPreProcessor"
    log_message "Directories for PostProcessor and PreProcessor created."
}


# --- 4. Compile external libraries and the main project ---
compile_project() {
    local project_path="$LAGRANGIAN_BASE_DIR/$PROJECT_DIR_NAME"
    log_message "Starting project compilation..."

    # Load Intel compiler environment variables
    local intel_dir
    intel_dir=$(spack find -p "$INTEL_COMPILER_VERSION" | grep -v -- '---' | tail -n1 | awk '{print $2}')
    
    if [ -z "$intel_dir" ]; then
        log_message "ERROR: Intel compiler path not found via Spack."
        exit 1
    fi

    log_message "Using Intel Compilers from: $intel_dir"
    source "$intel_dir/setvars.sh"

    # --- Integrated logic from "1-install_ExternalLibs.sh" ---
    if [ "$COMPILE_LIBS" = true ]; then
        log_message "Compiling external libraries..."
        cd "$project_path/ExternalLibs/"
        # The MakeLibraries.sh script is already in the repository, we just run it
        # Ensure script is executable
        chmod +x MakeLibraries.sh
        ./MakeLibraries.sh -intel
        log_message "External libraries compiled successfully."
    else
        log_message "SKIP: User configuration set COMPILE_LIBS to false."
    fi

    # --- Integrated logic from "2-compile_Lagrangian.sh" and "MakeMOHIDLagrangian.sh" ---
    log_message "Compiling MOHIDLagrangian (mode: $BUILD_TYPE)..."
    cd "$project_path"
    
    local build_dir="build"
    # Clean previous build
    if [ -d "$build_dir" ]; then
        rm -rf "$build_dir"
    fi
    mkdir "$build_dir"
    cd "$build_dir"

    # Configure with CMake according to BUILD_TYPE
    if [ "$BUILD_TYPE" == "DEBUG" ]; then
        cmake -Wno-dev -DCMAKE_BUILD_TYPE=DEBUG -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_Fortran_COMPILER_ID="Intel" ..
    else # RELEASE by default
        cmake -Wno-dev -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_Fortran_COMPILER_ID="Intel" ..
    fi
    
    # Compile with make
    make

    # Rename executable
    local exe_name_suffix=""
    if [ "$BUILD_TYPE" == "DEBUG" ]; then
        exe_name_suffix="_debug"
        local final_exe_name="MOHIDLagrangian${exe_name_suffix}"
        if [ -f "$project_path/$build_dir/bin/MOHIDLagrangian" ]; then
             mv "$project_path/$build_dir/bin/MOHIDLagrangian" "$project_path/$build_dir/bin/$final_exe_name"
        fi
    else  # RELEASE by default
        local final_exe_name="MOHIDLagrangian"
    fi
    
    if [ -f "$project_path/$build_dir/bin/$final_exe_name" ]; then
        log_message "Compilation finished successfully. Executable: $final_exe_name"
    else
        log_message "ERROR: Compilation failed. Executable not found at $project_path/$build_dir/bin/$final_exe_name"
        exit 1
    fi
}


# --- 5. Configure the Python environment ---
configure_python_environment() {
    local project_path="$LAGRANGIAN_BASE_DIR/$PROJECT_DIR_NAME"
    local env_name="MOHID-Lagrangian"
    log_message "Configuring Python environment '$env_name'..."

    # Load Miniconda from Spack
    source "$HOME/spack/share/spack/setup-env.sh"
    spack load "$MINICONDA_VERSION"
    
    # NEW: Initialize conda for the current bash session
    # This allows the script to recognize 'conda activate'
    eval "$(conda shell.bash hook)"
    
    # NEW: Permanently add conda initialization to .bashrc 
    # so 'conda' command works in new terminals
    conda init bash

    # Create Conda environment
    if conda env list | grep -q "$env_name"; then
        log_message "Conda environment '$env_name' already exists. Skipping creation."
    else
        log_message "Creating Conda environment '$env_name' with Python 3.11..."
        conda create --name "$env_name" python=3.11 -y
    fi

    log_message "Installing Python packages from requirements.txt..."
    # Use conda run to ensure we are using the correct environment's pip
    conda run -n "$env_name" pip install -r "$project_path/buildScripts/requirements.txt"
    
    log_message "Python environment configured."
}

#==============================================================================
#     SCRIPT EXECUTION
#==============================================================================

main() {
    log_message "Starting MOHID-Lagrangian installation..."
    
    install_base_dependencies
    clone_repository
    apply_modifications
    compile_project
    configure_python_environment

    # NEW: Add Spack and Miniconda auto-load to .bashrc for future sessions
    if ! grep -q "spack load $MINICONDA_VERSION" ~/.bashrc; then
        echo "spack load $MINICONDA_VERSION" >> ~/.bashrc
    fi

    log_message "Installation completed successfully!"
    echo ""
    echo "IMPORTANT: Run 'source ~/.bashrc' now to enable conda in this window."
    echo ""
    log_message "To run a test case:"
    echo "1. conda activate MOHID-Lagrangian"
    echo "2. cd $LAGRANGIAN_BASE_DIR/$PROJECT_DIR_NAME/RUN_Cases/Tagus3D_case"
    echo "3. ./RunCase.sh"
    echo ""
}

main
