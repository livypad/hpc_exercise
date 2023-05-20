# Read Me

UCAS 课程高性能计算系统课程的随附在slides的练习题和题解。

采用cmake进行构建。在mac和wsl上均有限地测试过。请提前安装OpenMP、MPI和Eigen3的依赖。目前依赖 [GoogleTest](https://github.com/google/googletest/tree/main)框架进行测试，GoogleTest目前由cmake自动下载管理。注意，对于使用MPI多进程的程序，由于MPI只能在运行前通过`mpirun -n <arg0> [PROGRAM]...`的形式，确定进程数目。因此，如果需要测试多进程程序的正确性，请单独调用`mpirun`来测试。
