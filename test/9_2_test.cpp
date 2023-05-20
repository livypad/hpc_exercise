#include <memory>
#include "gtest/gtest.h"

#include "common.h"
/*
 * 构造MPI数据结构
 * 传输 {int m[3]; float a[2]; char c[5];} x[10];
 */
TEST(TEST9_2, RightTest) {
  int rank;
  int size;
  MPI_Init(nullptr, nullptr);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size < 2) {
    printf("This test needs at least 2 processes\n");
    GTEST_SKIP();
  }

  auto x = new Struct92[10];

  // Initialize x
  if (rank == 1) {
    for (int i = 0; i < 10; i++) {
      for (int j = 0; j < 3; j++) {
        x[i].m_[j] = i * 10 + j;
      }
      for (int j = 0; j < 2; j++) {
        x[i].a_[j] = static_cast<float>(0.1 * (i * 10 + j));
      }
      for (int j = 0; j < 5; j++) {
        x[i].c_[j] = static_cast<char>('A' + i * 5 + j);
      }
    }
  } else if (rank == 0) {
    for (int i = 0; i < 10; i++) {
      std::memset(x->m_, 0, 3 * sizeof(int));
      std::memset(x->a_, 0, 2 * sizeof(float));
      std::memset(x->c_, 0, 5 * sizeof(char));
    }
  }

  auto struct9_2_datatype = Struct92Datatype(x);

  if (rank == 1) {
    MPI_Send(x, 10, struct9_2_datatype, 0, 0, MPI_COMM_WORLD);
  } else if (rank == 0) {
    MPI_Recv(x, 10, struct9_2_datatype, 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < 10; i++) {
      for (int j = 0; j < 3; j++) {
        EXPECT_EQ(x[i].m_[j], i * 10 + j);
      }
      for (int j = 0; j < 2; j++) {
        EXPECT_DOUBLE_EQ(x[i].a_[j], static_cast<float>(0.1 * (i * 10 + j)));
      }
      for (int j = 0; j < 5; j++) {
        EXPECT_EQ(x[i].c_[j], static_cast<char>('A' + i * 5 + j));
      }
    }
  }

  MPI_Type_free(&struct9_2_datatype);
  MPI_Finalize();

  delete[] x;
}