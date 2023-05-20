#include "common.h"
/*
 * 构造一个数据类型，传输分块矩阵 nxn 的下三角矩阵块(行优先，mxm)
 */

auto TriMatrixDatatype(int n, int m, MPI_Datatype old_datatype) -> MPI_Datatype {
  int old_type_size;
  MPI_Type_size(old_datatype, &old_type_size);

  MPI_Datatype new_datatype;
  MPI_Datatype tmp_datatype;
  auto tmp_block_len = new int[m];
  auto tmp_displacements = new int[m];
  for (int i = 0; i < m; i++) {
    tmp_block_len[i] = i + 1;
    tmp_displacements[i] = n * i;
  }
  MPI_Type_indexed(m, tmp_block_len, tmp_displacements, old_datatype, &tmp_datatype);
  MPI_Type_commit(&tmp_datatype);

  MPI_Type_create_resized(tmp_datatype, 0, (m * n + m) * old_type_size, &new_datatype);
  MPI_Type_free(&tmp_datatype);
  MPI_Type_commit(&new_datatype);
  return new_datatype;
}
