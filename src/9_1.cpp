#include "common.h"

/*
 * m 行 n 列矩阵，分为 sub_m 行 sub_n 的块
 * 要求构造3种数据类型，要求
 * 1. A00 打包传输出去
 * 2.（A00和A20）或者是 A00 传输出去
 * 3.（A00和A20）直接打包传输出去
 */

auto A00Datatype(int m, int sub_m, int sub_n, MPI_Datatype old_type) -> MPI_Datatype {
  MPI_Datatype new_datatype;
  MPI_Type_vector(sub_n, sub_m, m, old_type, &new_datatype);
  MPI_Type_commit(&new_datatype);
  return new_datatype;
}
auto A2x0Datatype(int m, int sub_m, int sub_n, MPI_Datatype old_type) -> MPI_Datatype {
  auto old_sub_mat_datatype = A00Datatype(m, sub_m, sub_n, old_type);

  int old_type_size;
  MPI_Type_size(old_type, &old_type_size);

  MPI_Datatype new_sub_mat_datatype;
  MPI_Type_create_resized(old_sub_mat_datatype, 0, 2 * m * sub_n * old_type_size, &new_sub_mat_datatype);

  MPI_Type_commit(&new_sub_mat_datatype);
  return new_sub_mat_datatype;
}

auto A0020Datatype(int m, int sub_m, int sub_n, MPI_Datatype old_type) -> MPI_Datatype {
  MPI_Datatype new_datatype;

  auto block_lengths = new int[sub_n * 2];
  auto displacements = new int[sub_n * 2];

  int old_type_size;
  MPI_Type_size(old_type, &old_type_size);

  for (int i = 0; i < sub_n; i++) {
    block_lengths[i] = sub_m;
    displacements[i] = m * i;
  }
  for (int i = sub_n; i < 2 * sub_n; i++) {
    block_lengths[i] = sub_m;
    displacements[i] = sub_n * m + m * i;  // i的起点是sub_n，因此隐含的包含了从2倍sub_n开始的偏移
  }

  MPI_Type_indexed(2 * sub_n, block_lengths, displacements, old_type, &new_datatype);
  MPI_Type_commit(&new_datatype);

  return new_datatype;
}
