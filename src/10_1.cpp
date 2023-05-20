#include "common.h"

auto MyAll2All(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
               MPI_Datatype recvtype, MPI_Comm comm) -> int {
  int rank;
  int size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int send_tp_size;
  int recv_tp_size;
  MPI_Type_size(sendtype, &send_tp_size);
  MPI_Type_size(recvtype, &recv_tp_size);

  for (int j = 0; j < sendcount; j++) {
    std::memcpy(reinterpret_cast<char *>(recvbuf) + (sendcount * rank + j) * recv_tp_size,
                reinterpret_cast<const char *>(sendbuf) + (sendcount * rank + j) * send_tp_size, send_tp_size);
  }

  for (int i = 0; i < size; i++) {
    // 按从0->size-1的顺序发送
    if (rank == i) {  // send
      for (int j = 0; j < size; j++) {
        if (j == rank) {
          continue;
        }
        MPI_Send(reinterpret_cast<const char *>(sendbuf) + (sendcount * j) * send_tp_size, sendcount, sendtype, j, rank,
                 comm);
      }
    } else {  // recv
      MPI_Recv(reinterpret_cast<char *>(recvbuf) + (sendcount * i) * recv_tp_size, recvcount, recvtype, i, i, comm,
               MPI_STATUS_IGNORE);
    }
  }
  return 1;
}
