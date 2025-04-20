# DDG_Tutorial_Fast
A fast C/C++ code for DDG_Tutorial

## Overview

This is a fast C/C++ implementation for the numerical simulation of flexible structures using the discrete differential geometry method. A simplified MATLAB version for educational purposes is available [here](https://github.com/weicheng-huang-mechanics/DDG_Tutorial).

## Dependencies

There are some dependencies required prior to compilation. All the codes are tested in the Ubuntun (linux) systems. For other operating systems, the users should be able to modify the commands below appropriately .

- [Eigen 3.4.0](http://eigen.tuxfamily.org/index.php?title=Main_Page)
  - Eigen is used for various linear algebra operations.
  - DDG_Tutorial_Fast is built with Eigen version 3.4.0 which can be downloaded [here](https://gitlab.com/libeigen/eigen/-/releases/3.4.0). After downloading the source code, install through cmake as follows.
    ```bash
    cd eigen-3.4.0 && mkdir build && cd build
    cmake ..
    sudo make install

- [Intel oneAPI Math Kernel Library (oneMKL)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html?operatingsystem=linux&distributions=webdownload&options=online)
  - Necessary for access to Pardiso, which is used as a sparse matrix solver.
  - Intel MKL is also used as the BLAS / LAPACK backend for Eigen.
  - **Ubuntu**: Follow the below steps.
    ```bash
    cd /tmp
    wget https://registrationcenter-download.intel.com/akdlm/irc_nas/18483/l_onemkl_p_2022.0.2.136.sh

    # This runs an installer, simply follow the instructions.
    sudo sh ./l_onemkl_p_2022.0.2.136.sh
    ```
  - Add one of the following to your .bashrc so that cmake can find the MKL library. Change the directory accordingly if your MKL version is different.
   Note that older versions require setting `MKLROOT` while newer versions require `MKL_DIR`.
   You can find out which one from the cmake error message.
    ```bash
    export MKLROOT=/opt/intel/oneapi/mkl/2022.0.2   # for older versions
    export MKL_DIR=/opt/intel/oneapi/mkl/2024.2     # for newer versions
    ```

- [OpenGL / GLUT](https://www.opengl.org/)
  - OpenGL / GLUT is used for barebones rendering through simple line graphics.
  - Simply install through apt package manager:
    - **Ubuntu**: `sudo apt-get install libglu1-mesa-dev freeglut3-dev mesa-common-dev`

- Lapack (*included in MKL*)


- Some potential required C++ dependencies (usually pre-installed in the system)
     
   ```bash
   sudo apt-get install liblapack-dev
   sudo apt-get install gfortran
   ```

## Compiling and Running Examples in C++

``DDG_Tutorial_Fast`` contains five examples for users to play around. Here, we utilize [CMAKE](https://cmake.org/) to compile the source files. In addition, a cmake flag ``TARGET_BUILD`` is used to specified which example is complied. The flag's value can be ``2d_curve``, ``2d_surface``, ``3d_curve``, ``3d_surface``, and ``hollow_net``. An example for ``2d_curve`` is presented below.


   ```bash
   cp inputdata/2d_curve/option.txt option.txt
   mkdir build && cd build
   cmake .. -DTARGET_BUILD=2d_curve
   make -j4
   cd ..
   ```
Afterwards, simply run the simulation using the ``simDER`` script.

```bash
./simDER option.txt
```
Note that the simulation parameters are contained in the ``option.txt`` and the user can modifed the values of them to see different physical effects in the simulation. 
The geometrical characterization is included in the ``inputdata`` folder. For example, the geometry of ``2d_curve`` is contained in the folder ``inputdata/2d_curve/inputdata``.

### Citation
If our work has helped your research, please cite the following paper.
```
@article{huang2025tutorial,
  title={A tutorial on simulating nonlinear behaviors of flexible structures with the discrete differential geometry (DDG) method},
  author={Huang, Weicheng and Hao, Zhuonan and Li, Jiahao and Tong, Dezhong and Guo, Kexin and Zhang, Yingchao and Gao, Huajian and Hsia, K Jimmy and Liu, Mingchao},
  journal={arXiv preprint arXiv:2504.11417},
  year={2025}
}
```


