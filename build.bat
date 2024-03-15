@REM @echo off
@echo Cleaning up directory:
set "batch_folder=%~dp0"
cd %batch_folder%
del *.exp *.lib *.o *.dll *.exe *.a

@REM Compile CUDA code into a static library:
nvcc -c cuda_main.cu -o cuda_main.o
ar rcs lib_cuda_main.a cuda_main.o
@REM g++ -c cuda_main.cpp -o cuda_main.o
@REM ar rcs lib_cuda_main.a cuda_main.o

@REM Compile your C++ code
g++ -c main.cpp -o main.o
@REM Link the static library and the object file
g++ main.o -L. -l_cuda_main -o main