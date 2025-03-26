# DDG_Tutorial_Fast
A fast C/C++ code for DDG_Tutorial

## Overview

This is a fast C/C++ version for the numerical simulation of flexible structures using the discrete differential geometry method. The simple MATLAB version for educational purposes can be found [here](https://github.com/weicheng-huang-mechanics/DDG_Tutorial).

## Prerequisites

- [Ubuntu 18.04 or above](https://ubuntu.com/tutorials/install-ubuntu-desktop#1-overview)
- C++ dependencies
     
   ```bash
   sudo apt-get install libblas-dev liblapack-dev
   sudo apt-get install gfortran
   sudo apt-get install freeglut3-dev
   https://gitlab.com/libeigen/eigen/-/releases/3.4.0
   MKL
   ```

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/weicheng-huang-mechanics/DDG_Tutorial_Fast.git
   cd DDG_Tutorial_Fast
   ```
   
2. Install C++ dependencies

- **Note**: Some of these packages are installed to the system library for convenience. You may want to install locally to e.g., `~/.local` to avoid conflicts with system libraries. Add the `cmake` flag: `-D CMAKE_INSTALL_PREFIX=~/.local`. Then `sudo` is not required to install. You'll need to ensure subsequent builds know where to find the build libraries.
- Update the package list:
  ```bash
  sudo apt update
  ```
    
- [Eigen 3.4.0](http://eigen.tuxfamily.org/index.php?title=Main_Page)
  - Eigen is a C++ template library for linear algebra.
  - Install via APT:
    ```bash
    sudo apt update
    sudo apt install libeigen3-dev
    ```
  - (Optional) Verify installation
    ```bash
    dpkg -s libeigen3-dev | grep Version
    ```

- [LLVM](https://releases.llvm.org/download.html)
  - LLVM is a collection of tools for building compilers and optimizing code.
  - Install via APT:
    ```bash
    sudo apt-get install llvm
    ```
  - (Optional) Verify installation
    ```bash
    llvm-config --version
    ```
    
- [GMP](https://gmplib.org/)
  - GMP is a free library for arbitrary precision arithmetic, operating on signed integers, rational numbers, and floating-point numbers.
  - Install via APT:
    ```bash
    sudo apt install libgmp-dev
    ```
  - (Optional) Verify installation
    ```bash
    dpkg -l | grep libgmp
    ```

- [Intel oneAPI Math Kernel Library (oneMKL)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html?operatingsystem=linux&distributions=webdownload&options=online)
  - Necessary for access to Pardiso, which is used as a sparse matrix solver.
  - Intel MKL is also used as the BLAS / LAPACK backend for Eigen.
  - Install via APT following the [official instruction](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html?operatingsystem=linux&linux-install=apt) Be sure to walk through both "Prerequisites for First-Time Users" and "Install with APT".
  - Check the installation version:
    By default, the installation directory path should be at `/opt/intel/oneapi/mkl`.
    Look for the folder named `mkl`, for example `/opt/intel/oneapi/mkl/2025.0`.
  - Set the MKL environment variable:
    ```bash
    export MKL_DIR=/opt/intel/oneapi/mkl/2025.0     # for newer versions
    ```
  - Add the above corresponding environment variable to your `.bashrc` file.
    ```bash
    nano ~/.bashrc
    ```
  - Reload the `.bashrc` file to apply the changes:
    ```bash
    source ~/.bashrc
    ```
  - (Optional) Verify the MKL installation:
    ```bash
    echo $MKL_DIR
    ```

- [OpenGL / GLUT](https://www.opengl.org/)
  - OpenGL / GLUT is used for rendering the knot through a simple graphic.
  - Install via APT
    ```bash
    sudo apt-get install libglu1-mesa-dev freeglut3-dev mesa-common-dev`
    ```

- [Lapack](https://www.netlib.org/lapack/) (*included in MKL*)

3. Configure the simulation engine


4. To simulate the bilayer robot with customized setting parameters, run

