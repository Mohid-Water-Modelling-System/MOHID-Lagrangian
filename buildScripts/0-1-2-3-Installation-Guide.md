# MOHID-Lagrangian Installation Guide — Sequential Build Scripts

This document describes the installation procedure based on the four sequential build scripts located under `buildScripts/`:

```text
0-basics-ExternalLibs.sh
1-install_ExternalLibs.sh
2-compile_Lagrangian.sh
3-MakeEnvMOHIDLagrangian_noconda.sh
```

The currently tested reference environment is **Ubuntu 24.04** with **GNU 13**, **Spack 1.2.2**, and **Intel oneAPI Compilers 2023.2.1**.

The four scripts divide the installation into four stages:

```text
0 - Prepare the operating system, Spack, and compiler environment
        ↓
1 - Compile the external libraries
        ↓
2 - Compile MOHID-Lagrangian
        ↓
3 - Create the Python virtual environment
```

---

## 1. Reference Environment

The validated reference environment is:

| Component | Version / Setting |
|---|---|
| Operating system | Ubuntu 24.04 |
| GNU C compiler | GCC 13 |
| GNU C++ compiler | G++ 13 |
| GNU Fortran compiler | GFortran 13 |
| glibc | Ubuntu 24.04 system glibc |
| Spack | 1.2.2 |
| Intel oneAPI Compilers | 2023.2.1 |
| Miniconda | 24.3.0 |
| MOHID-Lagrangian compiler | Intel classic `icc`, `icpc`, `ifort` |
| Python environment | Python `venv` in script `3-` |

The installation procedure is currently validated for **Ubuntu 24.04**. Older or newer Ubuntu releases should be treated as separate configurations and tested with their compatible system packages before being considered supported.

---

## 2. Expected Project Structure

The scripts assume the following project structure:

```text
MOHID-Lagrangian/
├── ExternalLibs/
├── src/
├── build/
├── buildScripts/
│   ├── 0-basics-ExternalLibs.sh
│   ├── 1-install_ExternalLibs.sh
│   ├── 2-compile_Lagrangian.sh
│   ├── 3-MakeEnvMOHIDLagrangian_noconda.sh
│   └── requirements.txt
└── RUN_Cases/
```

Scripts `1-` and `2-` determine the MOHID-Lagrangian project root relative to their own location:

```bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LagrangianMaster="$(cd "$SCRIPT_DIR/.." && pwd)"
```

Therefore, the repository can be located in different user directories without hard-coding the complete MOHID-Lagrangian path.

---

# Part I — `0-basics-ExternalLibs.sh`

## 3. Purpose

`0-basics-ExternalLibs.sh` prepares the base operating-system and software environment required to compile MOHID-Lagrangian and its external libraries.

The script:

1. checks the operating system;
2. verifies or installs required Ubuntu packages;
3. selects GNU 13;
4. installs and configures Spack 1.2.2;
5. detects the system compilers in Spack;
6. bootstraps Spack;
7. installs Intel oneAPI Compilers 2023.2.1;
8. installs Miniconda 24.3.0;
9. installs required build tools through Spack;
10. performs a final environment verification.

---

## 4. System-Level Installation Paths

The script uses:

```bash
SOFTWARE_DIR="/home/software"
SPACK_DIR="$SOFTWARE_DIR/spack"
SPACK_ENV_PATH="$SPACK_DIR/share/spack/setup-env.sh"
LAGRANGIAN_RES="$SOFTWARE_DIR/lagrangianResources"
```

Spack is therefore installed under:

```text
/home/software/spack
```

Spack-installed package paths must not be hard-coded because Spack generates architecture- and hash-dependent installation directories.

For example, an Intel installation may look like:

```text
/home/software/spack/opt/spack/linux-skylake/intel-oneapi-compilers-2023.2.1-<spack-hash>
```

The exact directory name can differ between systems.

---

## 5. Operating-System Check

The reference operating system is defined as:

```bash
REQUIRED_UBUNTU_VERSION="24.04"
```

The script reads `/etc/os-release` and verifies that the operating system is Ubuntu 24.04.

If another Ubuntu version is detected, the installation stops to avoid silently using an unvalidated system toolchain.

---

## 6. Required Ubuntu Packages

The script verifies or installs:

```text
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
```

`libtool` is required because the Proj4 preparation performed later by `1-install_ExternalLibs.sh` uses `autoreconf -fi`, which requires `libtoolize`.

---

## 7. Root and Non-Root Behavior

If the script is executed without root privileges, it checks the required packages and reports what is already installed.

It then instructs the user to run:

```bash
sudo bash 0-basics-ExternalLibs.sh
```

System package installation and the `/home/software` Spack installation require root privileges.

A typical execution from `buildScripts/` is:

```bash
cd MOHID-Lagrangian/buildScripts
sudo bash 0-basics-ExternalLibs.sh
```

---

## 8. GNU 13 Toolchain

The required GNU compiler version is:

```bash
GNU_COMPILER_MAJOR="13"
```

The script verifies `gcc-13`, `g++-13`, and `gfortran-13`, then selects them:

```bash
export CC=/usr/bin/gcc-13
export CXX=/usr/bin/g++-13
export FC=/usr/bin/gfortran-13
export F77=/usr/bin/gfortran-13
export F90=/usr/bin/gfortran-13
```

These compilers establish the tested system compiler environment for Spack. The final MOHID-Lagrangian build itself uses the Intel compiler environment loaded later.

---

## 9. glibc

The script reports the system glibc version with:

```bash
ldd --version | head -n 1
```

glibc is not installed or replaced manually.

---

## 10. Spack 1.2.2 Installation

The required version is:

```bash
SPACK_VERSION="v1.2.2"
```

Spack is cloned with:

```bash
git clone \
    -c feature.manyFiles=true \
    --branch "$SPACK_VERSION" \
    --depth 1 \
    https://github.com/spack/spack.git
```

The script also creates `/etc/profile.d/spack.sh` to make the Spack environment available system-wide.

---

## 11. Spack Bootstrap and Compiler Detection

After loading Spack:

```bash
spack bootstrap now
spack compiler find --scope=site /usr/bin
spack compiler list
spack spec cmake
```

These commands bootstrap Spack, register the selected system compilers, and verify that package concretization works.

---

## 12. Packages Installed Through Spack

The main fixed package versions are:

```text
intel-oneapi-compilers@2023.2.1
miniconda3@24.3.0
```

The script also installs:

```text
cmake
automake
autoconf
libtool
m4
perl
```

The Intel path is obtained dynamically with:

```bash
spack location -i --first intel-oneapi-compilers@2023.2.1
```

This avoids hard-coded architecture names and Spack hashes.

---

# Part II — `1-install_ExternalLibs.sh`

## 13. Purpose

`1-install_ExternalLibs.sh` prepares and compiles the external libraries. It detects the project root, loads Spack and Intel oneAPI, regenerates Proj4 4.9.3 Autotools files, sets execution permissions, and runs the Intel external-library build.

---

## 14. Relative Project Paths

The script uses:

```bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LagrangianMaster="$(cd "$SCRIPT_DIR/.." && pwd)"
```

Therefore:

```text
SCRIPT_DIR        → MOHID-Lagrangian/buildScripts
LagrangianMaster → MOHID-Lagrangian
```

No user-specific project path is required.

---

## 15. Spack and Intel Detection

The script first uses an already loaded Spack installation when available. Otherwise it checks:

```text
/home/software/spack/share/spack/setup-env.sh
```

Intel oneAPI is then found dynamically with:

```bash
inteldir=$(spack location -i --first intel-oneapi-compilers@2023.2.1 2>/dev/null)
```

and activated with:

```bash
source "$inteldir/setvars.sh"
```

---

## 16. Proj4 4.9.3 Preparation

The script works only inside:

```text
ExternalLibs/Proj4/Linux/proj-4.9.3
```

and executes:

```bash
autoreconf -fi
rm -f config.status config.log
find . -name Makefile -type f -delete
```

This removes stale generated build files and prevents old Autotools references such as `automake-1.15` from blocking the build.

---

## 17. Execution Permissions

Before compilation:

```bash
find . -type f -name "*.sh" -exec chmod +x {} +
find . -type f -name "configure" -exec chmod +x {} +
find . -type f -name "mkdirs" -exec chmod +x {} +
find . -type f -name "install-sh" -exec chmod +x {} +
```

---

## 18. External Library Compilation

The script enters:

```text
MOHID-Lagrangian/ExternalLibs
```

and runs:

```bash
./MakeLibraries.sh -intel
```

Run the script as a normal user:

```bash
cd MOHID-Lagrangian/buildScripts
./1-install_ExternalLibs.sh
```

### Important note

The current script contains an initial cleanup block that checks for directories named `MOHID-Lagrangian`, `build`, `src`, and `ExternalLibs` relative to the shell working directory. Therefore, the documented execution location is `buildScripts/`.

---

# Part III — `2-compile_Lagrangian.sh`

## 19. Purpose

`2-compile_Lagrangian.sh` compiles the main MOHID-Lagrangian executable.

It determines the project root relative to `buildScripts/`, loads Spack and Intel oneAPI, selects release or debug compilation, calls `MakeMOHIDLagrangian.sh`, writes `compile.log`, and verifies the final executable.

---

## 20. Build Mode

The build mode is selected with:

```bash
OPTION=1
```

Available values:

```text
1 - RELEASE_X64
2 - DEBUG_X64
```

Release produces:

```text
build/bin/MOHIDLagrangian
```

Debug produces:

```text
build/bin/MOHIDLagrangian_debug.exe
```

---

## 21. Build Paths

The script defines:

```bash
build_dir=$LagrangianMaster/build
output_path=$build_dir/bin
```

The project root is determined using the same `SCRIPT_DIR`/`LagrangianMaster` logic as script `1-`.

---

## 22. Release Compilation

For release:

```bash
./buildScripts/MakeMOHIDLagrangian.sh -intel > $SCRIPT_DIR/compile.log
```

The expected executable is:

```text
MOHID-Lagrangian/build/bin/MOHIDLagrangian
```

---

## 23. Debug Compilation

For debug:

```bash
./buildScripts/MakeMOHIDLagrangian.sh -intel -debug > $SCRIPT_DIR/compile.log
```

The resulting executable is renamed to:

```text
MOHID-Lagrangian/build/bin/MOHIDLagrangian_debug.exe
```

---

## 24. Running Script `2-`

```bash
cd MOHID-Lagrangian/buildScripts
./2-compile_Lagrangian.sh
```

The compilation log is:

```text
MOHID-Lagrangian/buildScripts/compile.log
```

---

# Part IV — `3-MakeEnvMOHIDLagrangian_noconda.sh`

## 25. Purpose

`3-MakeEnvMOHIDLagrangian_noconda.sh` creates a Python virtual environment using the system `python3` command.

Unlike the general installer, this script does **not** use Conda.

It executes:

```bash
python3 -m venv MOHID-Lagrangian
source MOHID-Lagrangian/bin/activate
pip install -r requirements.txt
pip list
```

---

## 26. Running Script `3-`

The current script uses paths relative to the directory from which it is executed. Therefore, run it from the directory containing `requirements.txt`, normally:

```bash
cd MOHID-Lagrangian/buildScripts
./3-MakeEnvMOHIDLagrangian_noconda.sh
```

The virtual environment is then created under:

```text
MOHID-Lagrangian/buildScripts/MOHID-Lagrangian/
```

After installation, the script waits for a key press before exiting.

---

# Part V — Complete Installation Procedure

## 27. Recommended Command Sequence

From the repository:

```bash
cd MOHID-Lagrangian/buildScripts
```

### Step 0 — Prepare the system

First inspect the current environment:

```bash
./0-basics-ExternalLibs.sh
```

If installation is required:

```bash
sudo bash 0-basics-ExternalLibs.sh
```

### Step 1 — Compile ExternalLibs

```bash
./1-install_ExternalLibs.sh
```

### Step 2 — Compile MOHID-Lagrangian

```bash
./2-compile_Lagrangian.sh
```

### Step 3 — Create the Python virtual environment

```bash
./3-MakeEnvMOHIDLagrangian_noconda.sh
```

---

## 28. Complete Workflow Summary

```text
0-basics-ExternalLibs.sh
        │
        ├── Verify Ubuntu 24.04
        ├── Install system dependencies
        ├── Select GNU 13
        ├── Install Spack 1.2.2
        ├── Bootstrap Spack
        ├── Register compilers
        ├── Install Intel oneAPI 2023.2.1
        ├── Install Miniconda 24.3.0
        └── Install build tools
        │
        ▼
1-install_ExternalLibs.sh
        │
        ├── Detect project root
        ├── Load Spack
        ├── Detect Intel dynamically
        ├── Regenerate Proj4 4.9.3
        ├── Set execution permissions
        └── Compile ExternalLibs
        │
        ▼
2-compile_Lagrangian.sh
        │
        ├── Detect project root
        ├── Load Spack and Intel
        ├── Select RELEASE or DEBUG
        ├── Compile MOHID-Lagrangian
        ├── Write compile.log
        └── Verify executable
        │
        ▼
3-MakeEnvMOHIDLagrangian_noconda.sh
        │
        ├── Create Python venv
        ├── Activate environment
        ├── Install requirements.txt
        └── Display installed packages
        │
        ▼
Installation ready
```

---

## 29. Important Outputs

### Spack

```text
/home/software/spack
```

### Intel oneAPI

Determine dynamically:

```bash
spack location -i --first intel-oneapi-compilers@2023.2.1
```

### Release executable

```text
MOHID-Lagrangian/build/bin/MOHIDLagrangian
```

### Debug executable

```text
MOHID-Lagrangian/build/bin/MOHIDLagrangian_debug.exe
```

### Compilation log

```text
MOHID-Lagrangian/buildScripts/compile.log
```

### Python virtual environment

With the current script `3-` executed from `buildScripts/`:

```text
MOHID-Lagrangian/buildScripts/MOHID-Lagrangian/
```

---

## 30. Important Notes

- The validated operating system is **Ubuntu 24.04**.
- The validated GNU compiler stack is **GNU 13**.
- Spack is pinned to **1.2.2**.
- Intel oneAPI Compilers are pinned to **2023.2.1**.
- Miniconda is pinned to **24.3.0** by script `0-`, although script `3-` uses Python `venv`, not Conda.
- Do not hard-code Spack-generated package installation paths or hashes.
- Proj4 regeneration is restricted to `ExternalLibs/Proj4/Linux/proj-4.9.3`.
- Scripts `1-` and `2-` derive the project root relative to `buildScripts/`.
- Script `3-` currently depends on the shell working directory for `requirements.txt` and the virtual-environment location.
- Script `0-` is the only stage intended to use root privileges when system packages or `/home/software` must be installed.
- Scripts `1-`, `2-`, and `3-` should normally run as the regular user.
