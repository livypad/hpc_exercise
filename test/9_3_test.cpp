#include <memory>
#include "gtest/gtest.h"

#include "common.h"

inline auto CheckRange(const int n, const int m, const int count, const int i, const int j) -> bool {
  int mat_rank = i / m;
  if (mat_rank < count) {
    return j >= m * mat_rank && j <= m * mat_rank + i % m;
  }
  return false;
}
TEST(TEST9_3, RightTest) {
  int rank;
  int size;
  MPI_Init(nullptr, nullptr);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size < 2) {
    printf("This test needs at least 2 processes\n");
    GTEST_SKIP();
  }

  const int n = 10;
  const int m = 2;
  const int send_count = 3;
  EXPECT_LE(send_count, n / m);

  auto a = new int[n * n];
  if (rank == 1) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        a[i * n + j] = i * n + j + 1;
      }
    }
  } else if (rank == 0) {
    std::memset(a, 0, n * n * sizeof(int));
  }

  auto tri_mat_datatype = TriMatrixDatatype(n, m, MPI_INT);

  if (rank == 1) {
    MPI_Send(a, send_count, tri_mat_datatype, 0, 0, MPI_COMM_WORLD);
  } else if (rank == 0) {
    MPI_Recv(a, send_count, tri_mat_datatype, 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (CheckRange(n, m, send_count, i, j)) {
          EXPECT_EQ(a[i * n + j], i * n + j + 1);
        } else {
          EXPECT_EQ(a[i * n + j], 0);
        }
      }
    }
  }

  delete[] a;

  MPI_Type_free(&tri_mat_datatype);
  MPI_Finalize();
}