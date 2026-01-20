#!/bin/bash

# Color configuration for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

# -------------------------------
#  NEW ROOT CHECK
# -------------------------------

if [ "$EUID" -ne 0 ]; then
    echo -e "${RED}============================================================${NC}"
    echo -e "${RED}  NOTICE: You do not have root privileges.                 ${NC}"
    echo -e "${RED}============================================================${NC}"
    echo -e ""
    echo -e "${YELLOW}The script will check if the required packages are installed, NO INSTALLATION WILL BE DONE WITHOUT ROOT PRIVILEGES.${NC}"

    APT_PKGS=( git vim m4 autotools-dev autoconf cmake gfortran wget "automake (1.15)" )
    PY_PKGS=( python3 python3-pip python3-venv clingo )
    SPACK_PKGS=( intel-oneapi-compilers@2023.2.1 miniconda3@24.3.0 cmake automake autoconf libtool m4 perl )

    echo -e "\n${YELLOW}--- Checking current installations ---${NC}"

    # --------------------------
    # CHECK APT PACKAGES
    # --------------------------
    echo -e "\n${YELLOW}>> APT Packages${NC}"
    for pkg in git vim m4 autotools-dev autoconf cmake gfortran wget; do
        if dpkg -l | grep -q "^ii  $pkg "; then
            echo -e "${GREEN}[INSTALLED]${NC} $pkg"
        else
            echo -e "${RED}[NOT INSTALLED]${NC} $pkg"
        fi
    done

    # Special case Automake 1.15
    if dpkg -l | grep -q "automake"; then
        echo -e "${GREEN}[INSTALLED]${NC} automake (check version manually if needed)"
    else
        echo -e "${RED}[NOT INSTALLED]${NC} automake 1.15"
    fi

    # --------------------------
    # CHECK PYTHON PACKAGES
    # --------------------------
    echo -e "\n${YELLOW}>> Python Packages${NC}"
	for pkg in python3 python3-pip python3-venv; do
		if dpkg -l | grep -q "^ii  $pkg "; then
			echo -e "${GREEN}[INSTALLED]${NC} $pkg"
		elif [ "$pkg" = "python3-venv" ] && python3 -m venv --help >/dev/null 2>&1; then
			echo -e "${GREEN}[INSTALLED]${NC} $pkg"
		else
			echo -e "${RED}[NOT INSTALLED]${NC} $pkg"
		fi
	done


    if python3 -m pip show clingo >/dev/null 2>&1; then
        echo -e "${GREEN}[INSTALLED]${NC} clingo (pip)"
    else
        echo -e "${RED}[NOT INSTALLED]${NC} clingo (pip)"
    fi

    # --------------------------
    # CHECK SPACK PACKAGES IF SPACK EXISTS
    # --------------------------
    echo -e "\n${YELLOW}>> Spack Packages${NC}"
    if command -v spack >/dev/null 2>&1; then
        echo -e "${GREEN}Spack detected. Checking packages...${NC}"
        for pkg in "${SPACK_PKGS[@]}"; do
            if spack find "$pkg" >/dev/null 2>&1; then
                echo -e "${GREEN}[INSTALLED]${NC} $pkg"
            else
                echo -e "${RED}[NOT INSTALLED]${NC} $pkg"
            fi
        done
    else
        echo -e "${RED}Spack is NOT installed.${NC}"
    fi

    echo -e "\n${RED}In order to install all required packages, run this script with root privileges.\n${NC}"
    exit 1
fi

# ====================================================================
#  From this point, the original script continues
# ====================================================================

# Define paths and required versions
SOFTWARE_DIR="/home/software"
SPACK_DIR="$SOFTWARE_DIR/spack"
LAGRANGIAN_RES="$SOFTWARE_DIR/lagrangianResources"
AUTOMAKE_URL="http://mirrors.edge.kernel.org/ubuntu/pool/main/a/automake-1.15/automake_1.15.1-3ubuntu2_all.deb"

# Helper function to check root privileges
am_i_root() {
    if [ "$EUID" -ne 0 ]; then
        return 1
    fi
    return 0
}

echo -e "${YELLOW}--- Starting Environment and Dependencies Check for MOHID Lagrangian ---${NC}"

# -----------------------------------------------------------------------------
# 1. SYSTEM PACKAGE INSTALLATION
# -----------------------------------------------------------------------------
echo -e "\n${YELLOW}>> Checking base system packages${NC}"

# Required package list
PACKAGES=("git" "vim" "m4" "autotools-dev" "autoconf" "cmake" "gfortran" "wget")

# Update indices if root
if am_i_root; then
    echo "Updating repositories (apt update)..."
    apt update -qq
fi

# Check installation status
MISSING_PKGS=()
for pkg in "${PACKAGES[@]}"; do
    if ! dpkg -l | grep -q "ii  $pkg "; then
        MISSING_PKGS+=("$pkg")
    else
        echo -e "${GREEN}[OK] $pkg is installed.${NC}"
    fi
done

# Install missing packages or report error
if [ ${#MISSING_PKGS[@]} -ne 0 ]; then
    if am_i_root; then
        echo -e "${YELLOW}Installing missing packages: ${MISSING_PKGS[*]}...${NC}"
        apt install -y "${MISSING_PKGS[@]}"
    else
        echo -e "${RED}[ERROR] Missing required packages: ${MISSING_PKGS[*]}.${NC}"
        echo -e "Administrator privileges are required to continue."
        exit 1
    fi
fi

# Specific check for Automake 1.15
echo -e "\n${YELLOW}>> Checking specific Automake version (1.15)${NC}"
if ! dpkg -l | grep -q "automake"; then
    if am_i_root; then
        echo "Downloading and installing automake 1.15..."
        wget -q "$AUTOMAKE_URL" -O /tmp/automake.deb
        apt install -y /tmp/automake.deb
        rm /tmp/automake.deb
    else
        echo -e "${RED}[ERROR] Automake 1.15 not detected.${NC}"
        echo -e "Administrator privileges are required for manual installation."
        exit 1
    fi
else
    echo -e "${GREEN}[OK] Automake detected.${NC}"
fi

# Spack setup in /home/software
echo -e "\n${YELLOW}>> Checking Spack installation${NC}"
if [ ! -d "$SPACK_DIR" ]; then
    if am_i_root; then
        echo "Creating software directory and installing Spack..."
        mkdir -p "$SOFTWARE_DIR"
        cd "$SOFTWARE_DIR" || exit
        git clone -c feature.manyFiles=true https://github.com/spack/spack.git
        chmod -R 755 "$SOFTWARE_DIR"
    else
        echo -e "${RED}[ERROR] Directory $SPACK_DIR does not exist.${NC}"
        echo -e "Contact the administrator for initial installation."
        exit 1
    fi
else
    echo -e "${GREEN}[OK] Spack found at $SPACK_DIR.${NC}"
fi

# Global environment variable setup
echo -e "\n${YELLOW}>> Checking global environment variables (/etc/environment)${NC}"
ENV_FILE="/etc/environment"
SPACK_ENV_PATH="$SPACK_DIR/share/spack/setup-env.sh"

if grep -q "$SPACK_ENV_PATH" "$ENV_FILE"; then
    echo -e "${GREEN}[OK] Environment configuration detected.${NC}"
else
    if am_i_root; then
        echo "Configuring environment variables for Spack..."
        cp "$ENV_FILE" "${ENV_FILE}.bak"
        echo ". $SPACK_ENV_PATH" >> "$ENV_FILE"

        if ! grep -q "$SPACK_DIR/bin" "$ENV_FILE"; then
             sed -i "s|PATH=\"|PATH=\"$SPACK_DIR/bin:|g" "$ENV_FILE"
        fi
        echo -e "${GREEN}Variables configured. You must restart session or reload environment.${NC}"
    else
        echo -e "${RED}[ERROR] Global environment not configured.${NC}"
        echo -e "Administrator must add Spack setup script to /etc/environment."
        exit 1
    fi
fi

# Load Spack environment temporarily
export PATH="$SPACK_DIR/bin:$PATH"
. "$SPACK_ENV_PATH"

# -----------------------------------------------------------------------------
# 2. SPACK AND PYTHON DEPENDENCIES
# -----------------------------------------------------------------------------
echo -e "\n${YELLOW}>> Checking Python environment and Spack libraries${NC}"

PY_PKGS=("python3" "python3-pip" "python3-venv")
MISSING_PY=()

for pkg in "${PY_PKGS[@]}"; do
    if ! dpkg -l | grep -q "ii  $pkg "; then
        MISSING_PY+=("$pkg")
    fi
done

if [ ${#MISSING_PY[@]} -ne 0 ]; then
    if am_i_root; then
        apt install -y "${MISSING_PY[@]}"
        apt install -y python3.12-venv 2>/dev/null || echo "Info: Using default python3-venv."
    else
        echo -e "${RED}[ERROR] Missing Python packages: ${MISSING_PY[*]}.${NC}"
        exit 1
    fi
fi

# Install Clingo
if ! python3 -m pip show clingo >/dev/null 2>&1; then
    echo "Installing Clingo library..."
    pip3 install --user --break-system-packages clingo || {
        echo -e "${RED}Error installing clingo.${NC}"; exit 1;
    }
else
    echo -e "${GREEN}[OK] Clingo library installed.${NC}"
fi

# Check Spack binary
if ! command -v spack &> /dev/null; then
    echo -e "${RED}[ERROR] 'spack' command not available in current PATH.${NC}"
    exit 1
fi

echo -e "${YELLOW}Checking Spack compilers and libraries...${NC}"

install_spack_pkg() {
    PKG_NAME=$1
    if spack find "$PKG_NAME" >/dev/null 2>&1; then
        echo -e "${GREEN}[OK] Spack: $PKG_NAME already installed.${NC}"
    else
        echo -e "${YELLOW}Installing $PKG_NAME (this may take a while)...${NC}"
        if [ ! -w "$SOFTWARE_DIR" ] && ! am_i_root; then
            echo -e "${RED}[ERROR] No write permissions in $SOFTWARE_DIR.${NC}"
            exit 1
        fi
        spack install "$PKG_NAME"
    fi
}

# Install specific packages
install_spack_pkg "intel-oneapi-compilers@2023.2.1"
install_spack_pkg "miniconda3@24.3.0"

# Build tools
for tool in cmake automake autoconf libtool m4 perl; do
    install_spack_pkg "$tool"
done

echo -e "\n${GREEN}--- Setup Complete ---${NC}"
echo -e "Remember to run 'eval spack load --sh <package>' to load necessary libraries before compiling."
