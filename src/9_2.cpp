#include "common.h"
/*
 * 构造MPI数据结构
 * 传输 {int m[3]; float a[2]; char c[5];} x[10];
 */
auto Struct92Datatype(Struct92 *x) -> MPI_Datatype {
  MPI_Datatype struct_datatype;
  int block_lengths[3] = {3, 2, 5};
  MPI_Aint displacements[3];
  MPI_Datatype types[3] = {MPI_INT, MPI_FLOAT, MPI_CHAR};

  MPI_Aint base_addr;
  MPI_Get_address(&x[0], &base_addr);
  MPI_Get_address(&x[0].m_, &displacements[0]);
  MPI_Get_address(&x[0].a_, &displacements[1]);
  MPI_Get_address(&x[0].c_, &displacements[2]);

  for (MPI_Aint &displacement : displacements) {
    displacement -= base_addr;
  }

  MPI_Type_create_struct(3, block_lengths, displacements, types, &struct_datatype);
  MPI_Type_commit(&struct_datatype);
  return struct_datatype;
}

/* 直接把 Struct92 视为连续的byte数组
 * 用MPI_Type_contiguous来传输byte数组
 *
 * WARNING（我不确定）
 * 但是，sizeof()是编译期确定，假设数组和传输的程序是分开编译的
 * 那么，sizeof()的结果可能是不一样的（eg 结构体对齐策略），会有问题
 */
auto Struct92Datatype2() -> MPI_Datatype {
  MPI_Datatype struct_datatype;
  MPI_Type_contiguous(sizeof(Struct92), MPI_BYTE, &struct_datatype);
  MPI_Type_commit(&struct_datatype);
  return struct_datatype;
}
