# MOHID-Lagrangian: Quick Start Guide / Guia de Início Rápido / Guía de Inicio Rápido

## 🇺🇸 English: Quick Start Guide

### **Script Overview**
This script automates the full deployment of the MOHID-Lagrangian model on Ubuntu or Rocky Linux systems. It handles system dependencies, compiler setups via Spack, source code cloning, and Python environment configuration.

### **Configuration Options (Variables)**
Modify these variables at the top of the `MOHIDLagrangian_install.sh` script to customize the process:

| Variable | Default | Description |
| :--- | :--- | :--- |
| `USE_SUDO` | `false` | Set to `true` to allow automatic installation of system packages. If `false`, it will display the required commands for manual installation by an admin. |
| `INSTALL_FOLDER_NAME` | `"lagrangian"` | The directory name inside your `$HOME` folder where everything will be installed. |
| `TARGET_BRANCH` | `"master"` | The specific GitHub branch of the repository you want to clone. |
| `COMPILE_LIBS` | `true` | Whether to compile the mandatory external libraries. Must be `true` for the first installation. |
| `BUILD_TYPE` | `"RELEASE"` | Compilation mode. Use `"RELEASE"` for maximum performance or `"DEBUG"` for development and troubleshooting. |
| `INTEL_VERSION` | `intel-oneapi-compilers` | The version of the Intel compiler to be managed and installed via Spack. |
| `MINICONDA_VERSION` | `miniconda3` | The Miniconda version to be installed via Spack for Python environment management. |

### **How to Execute**
1.  **Grant Execution Permissions**: 
    ```bash
    chmod +x MOHIDLagrangian_install.sh
    ```
2.  **Run the Script**: 
    ```bash
    ./MOHIDLagrangian_install.sh
    ```
3.  **Finalize**: Run `source ~/.bashrc` after completion to update your current terminal session paths.

---

## 🇵🇹 Português: Guia de Início Rápido

### **Visão Geral do Script**
Este script automatiza a instalação completa do modelo MOHID-Lagrangian em sistemas Ubuntu ou Rocky Linux. Ele gere as dependências do sistema, a configuração de compiladores via Spack, o download do código-fonte e a configuração do ambiente Python.

### **Opções de Configuração (Variáveis)**
Você pode ajustar estas variáveis no cabeçalho do script `MOHIDLagrangian_install.sh`:

| Variável | Padrão | Descrição |
| :--- | :--- | :--- |
| `USE_SUDO` | `false` | Defina como `true` para permitir a instalação automática de pacotes do sistema. Se `false`, o script listará os comandos para um administrador. |
| `INSTALL_FOLDER_NAME` | `"lagrangian"` | Nome da pasta no seu diretório pessoal onde tudo será instalado. |
| `TARGET_BRANCH` | `"master"` | A ramificação (branch) do GitHub a ser clonada. |
| `COMPILE_LIBS` | `true` | Define se deve compilar as bibliotecas externas. Deve ser `true` na primeira instalação. |
| `BUILD_TYPE` | `"RELEASE"` | Modo de compilação. Use `"RELEASE"` para desempenho ou `"DEBUG"` para testes. |

### **Como Executar**
1.  **Dar Permissões**: 
    ```bash
    chmod +x MOHIDLagrangian_install.sh
    ```
2.  **Executar**: 
    ```bash
    ./MOHIDLagrangian_install.sh
    ```
3.  **Finalizar**: Execute `source ~/.bashrc` após a conclusão para atualizar a sessão atual do terminal.

---

## 🇪🇸 Castellano: Guía de Inicio Rápido

### **Descripción del Script**
Este script automatiza la implementación completa del modelo MOHID-Lagrangian en sistemas Ubuntu o Rocky Linux. Gestiona las dependencias del sistema, la configuración de compiladores mediante Spack, la descarga del código fuente y la configuración del entorno de Python.

### **Opciones de Configuración (Variables)**
Puedes modificar estas variables en la parte superior del script `MOHIDLagrangian_install.sh`:

| Variable | Por defecto | Descripción |
| :--- | :--- | :--- |
| `USE_SUDO` | `false` | `true` para instalar paquetes del sistema automáticamente. Si es `false`, mostrará los comandos para ejecución manual. |
| `INSTALL_FOLDER_NAME` | `"lagrangian"` | Nombre de la carpeta en tu "Home" donde se instalará todo. |
| `TARGET_BRANCH` | `"master"` | La rama de GitHub que se desea clonar. |
| `COMPILE_LIBS` | `true` | Indica si se deben compilar las librerías externas (necesario en la primera instalación). |
| `BUILD_TYPE` | `"RELEASE"` | `"RELEASE"` para rendimiento óptimo o `"DEBUG"` para desarrollo/depuración. |

### **Cómo Ejecutar**
1.  **Dar Permisos de Ejecución**: 
    ```bash
    chmod +x MOHIDLagrangian_install.sh
    ```
2.  **Ejecutar el Script**: 
    ```bash
    ./MOHIDLagrangian_install.sh
    ```
3.  **Finalizar**: Ejecuta `source ~/.bashrc` al terminar para cargar las nuevas rutas en tu terminal.
