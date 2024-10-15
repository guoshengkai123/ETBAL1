# ETBAL1
Robust Projection Parameter Calibration in Cryo-ET with L1-norm Optimization

Prerequisites: 
-  Linux or Unix-like operating systems
-  GCC
	sudo apt install gcc
-  CMake
  	sudo apt install cmake
-  Ceres-Solver2.1 (C++ version)
	sudo apt-get install libeigen3-dev libatlas-base-dev libsuitesparse-dev
 	git clone https://ceres-solver.googlesource.com/ceres-solver
	tar zxf ceres-solver-2.1.0.tar.gz
	mkdir ceres-bin
  	cd ceres-bin
	cmake ../ceres-solver-2.1.0
	make -j3
	sudo make install

To run:
	bash run.sh

Note: Before running the program, you need to place the input data in the `build/bin/changedata_output` directory and adjust the corresponding parameters in the `main` function.
