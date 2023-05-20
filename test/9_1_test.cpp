#include <memory>
#include "gtest/gtest.h"

#include "common.h"

TEST(TEST9_1, RightTest) {
  const int m = 8;
  const int n = 16;

  const int sub_m = 2;
  const int sub_n = 4;

  int rank;
  int size;
  MPI_Init(nullptr, nullptr);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size < 2) {
    printf("This test needs at least 2 processes\n");
    GTEST_SKIP();
  }

  auto a = new int[m * n];

  if (rank == 1) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        a[i * m + j] = i * m + j + 1;
      }
    }
  } else if (rank == 0) {
    std::memset(a, 0, m * n * sizeof(int));
  }

  auto a00_datatype = A00Datatype(m, sub_m, sub_n, MPI_INT);

  if (rank == 1) {
    MPI_Send(a, 1, a00_datatype, 0, 0, MPI_COMM_WORLD);
  } else if (rank == 0) {
    MPI_Recv(a, 1, a00_datatype, 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        if (i < sub_n && j < sub_m) {
          EXPECT_EQ(a[i * m + j], i * m + j + 1);
        } else {
          EXPECT_EQ(a[i * m + j], 0);
        }
      }
    }
  }

  auto a2x_0_datatype = A2x0Datatype(m, sub_m, sub_n, MPI_INT);
  // 只传输 A00，sent count=1
  if (rank == 1) {
    MPI_Send(a, 1, a2x_0_datatype, 0, 0, MPI_COMM_WORLD);
  } else if (rank == 0) {
    std::memset(a, 0, m * n * sizeof(int));

    MPI_Recv(a, 1, a2x_0_datatype, 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        //        LOG_DEBUG("Thread %d, a[%d][%d] = %d", rank, i, j, a[i * m + j]);
        if (j < sub_m && i < sub_n) {
          EXPECT_EQ(a[i * m + j], i * m + j + 1);
        } else {
          EXPECT_EQ(a[i * m + j], 0);
        }
      }
    }
  }

  // 传输 A00 A20，sent count=2
  // 如果需要传输 A00 A20 A40，sent count=3即可
  if (rank == 1) {
    MPI_Send(a, 2, a2x_0_datatype, 0, 0, MPI_COMM_WORLD);
  } else if (rank == 0) {
    std::memset(a, 0, m * n * sizeof(int));
    MPI_Recv(a, 2, a2x_0_datatype, 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        if (j < sub_m && (i < sub_n || (i >= 2 * sub_n && i < 3 * sub_n))) {
          EXPECT_EQ(a[i * m + j], i * m + j + 1);
        } else {
          EXPECT_EQ(a[i * m + j], 0);
        }
      }
    }
  }

  auto a0020_datatype = A0020Datatype(m, sub_m, sub_n, MPI_INT);
  if (rank == 1) {
    MPI_Send(a, 1, a0020_datatype, 0, 0, MPI_COMM_WORLD);
  } else if (rank == 0) {
    std::memset(a, 0, m * n * sizeof(int));
    MPI_Recv(a, 1, a0020_datatype, 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        if (j < sub_m && (i < sub_n || (i >= 2 * sub_n && i < 3 * sub_n))) {
          EXPECT_EQ(a[i * m + j], i * m + j + 1);
        } else {
          EXPECT_EQ(a[i * m + j], 0);
        }
      }
    }
  }

  MPI_Type_free(&a00_datatype);
  MPI_Type_free(&a2x_0_datatype);
  delete[] a;
  MPI_Finalize();
}
