#include "common.h"

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
