#include "CepGen/Utils/Algebra.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

using namespace cepgen;

int main(int argc, char* argv[]) {
  ArgumentsParser(argc, argv).parse();

  {  // test matrix/vector coordinates retrieval
    const auto A = Matrix{{1, 2}, {3, 4}};
    CG_DEBUG("main") << "A =\n" << A << ".";
    CG_TEST_EQUAL(A(0, 1), 2, "coordinates");
    const auto diag = Vector{1, 4};
    CG_TEST_EQUAL(A.diagonal(), diag, "diagonals");
  }

  const auto B = Matrix{{1, 2, 3}, {4, 5, 6}};
  CG_DEBUG("main") << "B =\n" << B << ".";

  {  // test transposition
    auto Bt = B.transposed();
    CG_TEST_EQUAL(B.numRows(), Bt.numColumns(), "transposed dim.1");
    CG_TEST_EQUAL(B.numColumns(), Bt.numRows(), "transposed dim.2");
    const auto &b01 = B(0, 1), &bt10 = Bt(1, 0), &b11 = B(1, 1), &bt11 = Bt(1, 1);
    CG_TEST_EQUAL(b01, bt10, "transposed coord.1");
    CG_TEST_EQUAL(b11, bt11, "transposed coord.2");
  }

  {  // test matrix/vector multiplication
    const auto v = Vector{7, 8, 9}, res = Vector{50, 122};
    CG_TEST_EQUAL(B * v, res, "matrix-vector mult.");
  }

  {  // test matrix/matrix multiplication
    const auto C = Matrix{{7, 8, 9}, {10, 11, 12}, {13, 14, 15}}, res = Matrix{{66, 72, 78}, {156, 171, 186}};
    CG_TEST_EQUAL(B * C, res, "matrix-matrix mult.");
  }

  const auto D = Matrix{// 4x4
                        {0.18, 0.60, 0.57, 0.96},
                        {0.41, 0.24, 0.99, 0.58},
                        {0.14, 0.30, 0.97, 0.66},
                        {0.51, 0.13, 0.19, 0.85}};

  {  // test linear equations solving
    const auto w = Vector{1., 2., 3., 4.};
    const auto res = D * (D % w);
    // D%w should be close enough to Vector{ -4.05205, -12.6056, 1.66091, 8.69377 }
    for (size_t i = 0; i < res.numRows(); ++i)
      CG_TEST_EQUIV(res(i), w(i), "lin.alg.coord." + std::to_string(i));
  }

  {  // test matrix inversion
    const auto Dinv = D.inverted();
    const auto IdD = Matrix::identity(D.numRows());
    const auto ZeD = Matrix::zero(D.numColumns());
    CG_TEST_EQUAL((Dinv * D - IdD).truncate(), ZeD, "D*D^{-1}");
    CG_TEST_EQUAL((D * Dinv - IdD).truncate(), ZeD, "D^{-1}*D");
  }

  return 0;
}
