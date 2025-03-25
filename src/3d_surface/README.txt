TO run to code, you may need a Linux system and install the following things:

1, eigen3, https://eigen.tuxfamily.org/dox-devel/index.html

2, openGL, sudo apt-get install freeglut3-dev

3, lapack, sudo apt-get install libatlas-base-dev

4, gfortran, sudo apt-get install gfortran

5, pardiso, https://www.pardiso-project.org/. The process to get a licence may be a little bit complex. Ask me if you need help.

6, make:
g++ -I /usr/local/include/eigen-3.3.7/ main.cpp world.cpp setInput.cpp timeStepper.cpp inertialForce.cpp externalGravityForce.cpp dampingForce.cpp elasticStretchingForce.cpp elasticBendingForce.cpp elasticPlate.cpp -lGL -lglut -lGLU -lpardiso600-GNU800-X86-64 -llapack -lgfortran -fopenmp -lpthread -lm -Ofast -o simDER 

7, run the code: 
export OMP_NUM_THREADS=1; ./simDER option.txt
