#!/bin/bash

# =============================================================================
# MOHID-Lagrangian
# Basic system dependencies and Spack installation
# =============================================================================


# -----------------------------------------------------------------------------
# Color configuration for output
# -----------------------------------------------------------------------------

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'


# -----------------------------------------------------------------------------
# Define paths and required versions
# -----------------------------------------------------------------------------

SOFTWARE_DIR="/home/software"
SPACK_DIR="$SOFTWARE_DIR/spack"

SPACK_VERSION="v1.2.2"
EXPECTED_SPACK_VERSION="${SPACK_VERSION#v}"

SPACK_ENV_PATH="$SPACK_DIR/share/spack/setup-env.sh"

LAGRANGIAN_RES="$SOFTWARE_DIR/lagrangianResources"

# NEW: Reference environment based on the working Zeus installation
REQUIRED_UBUNTU_VERSION="24.04"
GNU_COMPILER_MAJOR="13"


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

am_i_root() {
    if [ "$EUID" -ne 0 ]; then
        return 1
    fi

    return 0
}


apt_pkg_installed() {
    dpkg-query -W -f='${Status}' "$1" 2>/dev/null | \
        grep -q "install ok installed"
}


# =============================================================================
# ROOT CHECK
# =============================================================================

if ! am_i_root; then

    echo -e "${RED}============================================================${NC}"
    echo -e "${RED}  NOTICE: You do not have root privileges.                 ${NC}"
    echo -e "${RED}============================================================${NC}"

    echo
    echo -e "${YELLOW}The script will check the required packages.${NC}"
    echo -e "${YELLOW}NO INSTALLATION WILL BE DONE WITHOUT ROOT PRIVILEGES.${NC}"


    # -------------------------------------------------------------------------
    # Check APT packages
    # -------------------------------------------------------------------------

    echo -e "\n${YELLOW}>> APT Packages${NC}"

    APT_PKGS=(
        git
        vim
        m4
        autotools-dev
        autoconf
        automake
        libtool
        cmake
        build-essential
        gcc-13
        g++-13
        gfortran-13
        wget
        python3
        python3-pip
        python3-venv
    )

    for pkg in "${APT_PKGS[@]}"; do

        if apt_pkg_installed "$pkg"; then
            echo -e "${GREEN}[INSTALLED]${NC} $pkg"
        else
            echo -e "${RED}[NOT INSTALLED]${NC} $pkg"
        fi

    done


    # -------------------------------------------------------------------------
    # Check compilers
    # -------------------------------------------------------------------------

    echo -e "\n${YELLOW}>> GNU Compilers${NC}"

    for compiler in gcc-13 g++-13 gfortran-13; do

        if command -v "$compiler" >/dev/null 2>&1; then

            echo -e "${GREEN}[INSTALLED]${NC} $compiler -> $(command -v "$compiler")"

        else

            echo -e "${RED}[NOT INSTALLED]${NC} $compiler"

        fi

    done


    # -------------------------------------------------------------------------
    # Check Spack
    # -------------------------------------------------------------------------

    echo -e "\n${YELLOW}>> Spack${NC}"

    if command -v spack >/dev/null 2>&1; then

        echo -e "${GREEN}[INSTALLED]${NC} Spack"
        echo "Location: $(command -v spack)"
        echo "Version:  $(spack --version)"

    elif [ -x "$SPACK_DIR/bin/spack" ]; then

        echo -e "${GREEN}[INSTALLED]${NC} Spack"
        echo "Location: $SPACK_DIR/bin/spack"
        echo "Version:  $("$SPACK_DIR/bin/spack" --version)"

    else

        echo -e "${RED}[NOT INSTALLED]${NC} Spack"

    fi


    echo
    echo -e "${RED}To install the missing dependencies, run:${NC}"
    echo
    echo "    sudo bash $0"
    echo

    exit 1
fi


# =============================================================================
# START INSTALLATION
# =============================================================================

echo -e "${YELLOW}--- Starting Environment and Dependencies Check for MOHID Lagrangian ---${NC}"


# -----------------------------------------------------------------------------
# NEW: CHECK REFERENCE OPERATING SYSTEM
# -----------------------------------------------------------------------------

echo -e "\n${YELLOW}>> Checking operating system${NC}"

if [ ! -f /etc/os-release ]; then

    echo -e "${RED}[ERROR] /etc/os-release was not found.${NC}"
    exit 1

fi


. /etc/os-release


echo "Detected operating system:"
echo "    $PRETTY_NAME"


if [ "$ID" != "ubuntu" ]; then

    echo
    echo -e "${RED}[ERROR] This installation procedure is prepared for Ubuntu.${NC}"
    exit 1

fi


if [ "$VERSION_ID" != "$REQUIRED_UBUNTU_VERSION" ]; then

    echo
    echo -e "${RED}[ERROR] Ubuntu $VERSION_ID detected.${NC}"
    echo
    echo "Reference MOHID-Lagrangian environment:"
    echo "    Ubuntu $REQUIRED_UBUNTU_VERSION"
    echo
    echo "The working Zeus installation uses Ubuntu 24.04."
    echo
    echo "Installation stopped to avoid using an untested system toolchain."
    echo

    exit 1

fi


echo -e "${GREEN}[OK] Ubuntu $VERSION_ID detected.${NC}"


# -----------------------------------------------------------------------------
# 1. SYSTEM PACKAGE INSTALLATION
# -----------------------------------------------------------------------------

echo -e "\n${YELLOW}>> Checking base system packages${NC}"


PACKAGES=(
    git
    vim
    m4
    autotools-dev
    autoconf
    automake
    libtool
    cmake
    build-essential
    gcc-13
    g++-13
    gfortran-13
    wget
    python3
    python3-pip
    python3-venv
)


echo "Updating repositories (apt update)..."

if ! apt update -qq; then

    echo -e "${RED}[ERROR] apt update failed.${NC}"
    exit 1

fi


MISSING_PKGS=()


for pkg in "${PACKAGES[@]}"; do

    if apt_pkg_installed "$pkg"; then

        echo -e "${GREEN}[OK] $pkg is installed.${NC}"

    else

        MISSING_PKGS+=("$pkg")

    fi

done


if [ ${#MISSING_PKGS[@]} -ne 0 ]; then

    echo
    echo -e "${YELLOW}Installing missing packages:${NC}"
    echo "${MISSING_PKGS[*]}"

    if ! apt install -y "${MISSING_PKGS[@]}"; then

        echo -e "${RED}[ERROR] Failed to install required packages.${NC}"
        exit 1

    fi

fi


# -----------------------------------------------------------------------------
# 2. VERIFY GNU COMPILERS
# -----------------------------------------------------------------------------

echo -e "\n${YELLOW}>> Checking GNU $GNU_COMPILER_MAJOR compiler toolchain${NC}"


for compiler in gcc-13 g++-13 gfortran-13; do

    if ! command -v "$compiler" >/dev/null 2>&1; then

        echo -e "${RED}[ERROR] Required compiler '$compiler' was not found.${NC}"
        exit 1

    fi

    echo -e "${GREEN}[OK]${NC} $compiler -> $(command -v "$compiler")"

done


echo

gcc-13 --version | head -n 1
g++-13 --version | head -n 1
gfortran-13 --version | head -n 1


# NEW: Explicitly select the GNU 13 compiler toolchain
export CC=/usr/bin/gcc-13
export CXX=/usr/bin/g++-13
export FC=/usr/bin/gfortran-13
export F77=/usr/bin/gfortran-13
export F90=/usr/bin/gfortran-13


echo
echo "Selected compiler environment:"
echo "    CC=$CC"
echo "    CXX=$CXX"
echo "    FC=$FC"


# NEW: Report glibc version supplied by Ubuntu
echo
echo "glibc:"
ldd --version | head -n 1


# -----------------------------------------------------------------------------
# 3. CHECK AUTOMAKE
# -----------------------------------------------------------------------------

echo -e "\n${YELLOW}>> Checking Automake${NC}"


if command -v automake >/dev/null 2>&1; then

    echo -e "${GREEN}[OK] $(automake --version | head -n 1)${NC}"

else

    echo -e "${RED}[ERROR] Automake was not installed correctly.${NC}"
    exit 1

fi


# -----------------------------------------------------------------------------
# 4. SPACK INSTALLATION
# -----------------------------------------------------------------------------

echo -e "\n${YELLOW}>> Checking Spack installation${NC}"


if [ -d "$SPACK_DIR" ]; then

    if [ -x "$SPACK_DIR/bin/spack" ]; then

        CURRENT_SPACK_VERSION=$(
            "$SPACK_DIR/bin/spack" --version 2>/dev/null | awk '{print $1}'
        )

        echo "Existing Spack version: $CURRENT_SPACK_VERSION"
        echo "Required Spack version: $EXPECTED_SPACK_VERSION"


        if [ "$CURRENT_SPACK_VERSION" != "$EXPECTED_SPACK_VERSION" ]; then

            BACKUP_DIR="${SPACK_DIR}-backup-$(date +%Y%m%d-%H%M%S)"

            echo
            echo -e "${YELLOW}Different Spack version detected.${NC}"
            echo "Moving existing installation to:"
            echo
            echo "    $BACKUP_DIR"
            echo

            if ! mv "$SPACK_DIR" "$BACKUP_DIR"; then

                echo -e "${RED}[ERROR] Could not back up the existing Spack installation.${NC}"
                exit 1

            fi

        fi

    else

        BACKUP_DIR="${SPACK_DIR}-backup-$(date +%Y%m%d-%H%M%S)"

        echo -e "${YELLOW}Incomplete Spack directory detected.${NC}"
        echo "Moving it to:"
        echo
        echo "    $BACKUP_DIR"

        mv "$SPACK_DIR" "$BACKUP_DIR" || exit 1

    fi

fi


if [ ! -d "$SPACK_DIR" ]; then

    echo
    echo "Installing Spack $SPACK_VERSION..."

    mkdir -p "$SOFTWARE_DIR"

    cd "$SOFTWARE_DIR" || exit 1


    if ! git clone \
        -c feature.manyFiles=true \
        --branch "$SPACK_VERSION" \
        --depth 1 \
        https://github.com/spack/spack.git; then

        echo -e "${RED}[ERROR] Failed to clone Spack $SPACK_VERSION.${NC}"
        exit 1

    fi


    chmod -R a+rX "$SPACK_DIR"

else

    echo -e "${GREEN}[OK] Spack found at $SPACK_DIR.${NC}"

fi


# -----------------------------------------------------------------------------
# 5. VERIFY SPACK VERSION
# -----------------------------------------------------------------------------

echo -e "\n${YELLOW}>> Verifying Spack version${NC}"


CURRENT_SPACK_VERSION=$(
    "$SPACK_DIR/bin/spack" --version | awk '{print $1}'
)


if [ "$CURRENT_SPACK_VERSION" != "$EXPECTED_SPACK_VERSION" ]; then

    echo -e "${RED}[ERROR] Incorrect Spack version.${NC}"
    echo "Expected: $EXPECTED_SPACK_VERSION"
    echo "Found:    $CURRENT_SPACK_VERSION"

    exit 1

fi


echo -e "${GREEN}[OK] Spack version: $CURRENT_SPACK_VERSION${NC}"


# -----------------------------------------------------------------------------
# 6. CONFIGURE SPACK ENVIRONMENT
# -----------------------------------------------------------------------------

echo -e "\n${YELLOW}>> Configuring Spack environment${NC}"


SPACK_PROFILE="/etc/profile.d/spack.sh"


cat > "$SPACK_PROFILE" <<EOF
# Spack environment for MOHID-Lagrangian

export SPACK_ROOT="$SPACK_DIR"

if [ -f "\$SPACK_ROOT/share/spack/setup-env.sh" ]; then
    . "\$SPACK_ROOT/share/spack/setup-env.sh"
fi
EOF


chmod 644 "$SPACK_PROFILE"


echo -e "${GREEN}[OK] Spack environment file:${NC}"
echo "    $SPACK_PROFILE"


# -----------------------------------------------------------------------------
# 7. LOAD SPACK IN CURRENT SCRIPT
# -----------------------------------------------------------------------------

echo -e "\n${YELLOW}>> Loading Spack${NC}"


export SPACK_ROOT="$SPACK_DIR"
export PATH="$SPACK_DIR/bin:$PATH"


if [ ! -f "$SPACK_ENV_PATH" ]; then

    echo -e "${RED}[ERROR] Spack setup script was not found:${NC}"
    echo "    $SPACK_ENV_PATH"

    exit 1

fi


source "$SPACK_ENV_PATH"


if ! command -v spack >/dev/null 2>&1; then

    echo -e "${RED}[ERROR] Spack command is not available.${NC}"
    exit 1

fi


echo -e "${GREEN}[OK] Spack loaded:${NC}"
echo "    $(command -v spack)"

echo -e "${GREEN}[OK] Spack version:${NC}"
echo "    $(spack --version)"


# -----------------------------------------------------------------------------
# 8. BOOTSTRAP SPACK
# -----------------------------------------------------------------------------

echo -e "\n${YELLOW}>> Bootstrapping Spack dependencies${NC}"


if ! spack bootstrap now; then

    echo -e "${RED}[ERROR] Spack bootstrap failed.${NC}"
    exit 1

fi


echo -e "${GREEN}[OK] Spack bootstrap completed.${NC}"


# -----------------------------------------------------------------------------
# 9. DETECT SYSTEM COMPILERS
# -----------------------------------------------------------------------------

echo -e "\n${YELLOW}>> Detecting system compilers${NC}"


# GNU 13 has already been selected above through CC/CXX/FC.
if ! spack compiler find --scope=site /usr/bin; then

    echo -e "${RED}[ERROR] Spack could not configure the system compiler.${NC}"
    exit 1

fi


echo
spack compiler list


# -----------------------------------------------------------------------------
# 10. TEST SPACK CONCRETIZATION
# -----------------------------------------------------------------------------

echo -e "\n${YELLOW}>> Testing Spack concretization${NC}"


if ! spack spec cmake >/dev/null; then

    echo
    echo -e "${RED}[ERROR] Spack cannot concretize a C/C++ package.${NC}"
    echo -e "${RED}Check the compiler configuration before continuing.${NC}"
    echo

    spack compiler list

    exit 1

fi


echo -e "${GREEN}[OK] Spack concretization test passed.${NC}"


# -----------------------------------------------------------------------------
# 11. INSTALL SPACK PACKAGES
# -----------------------------------------------------------------------------

echo -e "\n${YELLOW}>> Checking Spack packages${NC}"


install_spack_pkg() {

    PKG_NAME="$1"

    echo
    echo -e "${YELLOW}Checking $PKG_NAME${NC}"


    if spack location -i --first "$PKG_NAME" >/dev/null 2>&1; then

        echo -e "${GREEN}[OK] $PKG_NAME already installed.${NC}"

    else

        echo -e "${YELLOW}Installing $PKG_NAME...${NC}"

        if ! spack install "$PKG_NAME"; then

            echo
            echo -e "${RED}[ERROR] Failed to install:${NC}"
            echo -e "${RED}        $PKG_NAME${NC}"
            echo

            exit 1

        fi

        echo -e "${GREEN}[OK] $PKG_NAME installed successfully.${NC}"

    fi
}


# Intel oneAPI compiler
install_spack_pkg "intel-oneapi-compilers@2023.2.1"


# Miniconda
install_spack_pkg "miniconda3@24.3.0"


# Build tools
for tool in cmake automake autoconf libtool m4 perl; do

    install_spack_pkg "$tool"

done


# -----------------------------------------------------------------------------
# 12. FINAL VERIFICATION
# -----------------------------------------------------------------------------

echo -e "\n${YELLOW}>> Final verification${NC}"


echo
echo "Operating system:"
echo "    $PRETTY_NAME"

echo
echo "GNU compiler:"
echo "    $(gcc-13 --version | head -n 1)"

echo
echo "GNU C++ compiler:"
echo "    $(g++-13 --version | head -n 1)"

echo
echo "GNU Fortran compiler:"
echo "    $(gfortran-13 --version | head -n 1)"

echo
echo "glibc:"
echo "    $(ldd --version | head -n 1)"

echo
echo "Spack:"
echo "    $(command -v spack)"

echo
echo "Spack version:"
echo "    $(spack --version)"

echo
echo "Intel oneAPI compiler:"
spack location -i --first intel-oneapi-compilers@2023.2.1

echo
echo "Miniconda:"
spack location -i --first miniconda3@24.3.0


echo
echo -e "${GREEN}============================================================${NC}"
echo -e "${GREEN} MOHID-Lagrangian dependency setup completed successfully  ${NC}"
echo -e "${GREEN}============================================================${NC}"

echo
echo "Spack installation:"
echo "    $SPACK_DIR"

echo
echo "To load Spack manually:"
echo
echo "    source $SPACK_ENV_PATH"
echo