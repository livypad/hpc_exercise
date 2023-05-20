#include <memory>
#include "gtest/gtest.h"

#include "common.h"

TEST(TEST10_1, RightTest) {
  int rank;
  int np;

  MPI_Init(nullptr, nullptr);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  int elem_len = 1;
  int len = elem_len * np;

  auto send_buffer = new int[len];
  auto recv_buffer = new int[len];
  auto check_buffer = new int[len];

  for (int i = 0; i < len; i++) {
    send_buffer[i] = i + rank * len;
    recv_buffer[i] = -1;
    check_buffer[i] = -1;
  }

  MPI_Alltoall(send_buffer, elem_len, MPI_INT, check_buffer, elem_len, MPI_INT, MPI_COMM_WORLD);
  MyAll2All(send_buffer, elem_len, MPI_INT, recv_buffer, elem_len, MPI_INT, MPI_COMM_WORLD);

  for (int i = 0; i < len; i++) {
    EXPECT_EQ(recv_buffer[i], check_buffer[i]);
    LOG_DEBUG("Thread %d, recv[%d] = %d,check[%d] = %d", rank, i, recv_buffer[i],i, check_buffer[i]);
  }

  MPI_Finalize();
}
