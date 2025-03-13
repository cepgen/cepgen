/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2025  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <string>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Algebra.h"

using namespace cepgen;

typedef std::unique_ptr<gsl_permutation, void (*)(gsl_permutation*)> Permutation;

Matrix::Matrix(size_t num_rows, size_t num_cols) : gsl_mat_(gsl_matrix_alloc(num_rows, num_cols), gsl_matrix_free) {}

Matrix::Matrix(const std::initializer_list<Vector>& mat)
    : gsl_mat_(gsl_matrix_alloc(mat.size(), mat.begin()->size()), gsl_matrix_free) {
  if (mat.size() == 0)
    return;

  size_t ix = 0;
  for (const auto& line : mat) {
    if (line.size() != numColumns())
      throw CG_FATAL("Matrix") << "Invalid initialiser; all lines must have the same number of elements!";
    for (size_t iy = 0; iy < line.size(); ++iy)
      operator()(ix, iy) = line(iy);
    ++ix;
  }
}

Matrix::Matrix(const Matrix& oth) : gsl_mat_(gsl_matrix_alloc(oth.numRows(), oth.numColumns()), gsl_matrix_free) {
  if (gsl_matrix_memcpy(gsl_mat_.get(), oth.gsl_mat_.get()) != 0)
    throw CG_FATAL("Matrix") << "Failed to clone a matrix!";
}

Matrix& Matrix::operator=(const Matrix& oth) {
  gsl_mat_.reset(gsl_matrix_alloc(oth.numRows(), oth.numColumns()));
  if (gsl_matrix_memcpy(gsl_mat_.get(), oth.gsl_mat_.get()) != 0)
    throw CG_FATAL("Matrix") << "Failed to clone a matrix!";
  return *this;
}

Matrix::operator Vector() const {
  if (numRows() == 1)
    return row(0);
  if (numColumns() == 1)
    return column(0);
  throw CG_FATAL("Matrix:Vector") << "Only 1xN matrices can be converted to vectors.";
}

Matrix Matrix::zero(size_t num_rows, size_t num_cols) {
  if (num_cols == 0ull)
    num_cols = num_rows;
  Matrix out(num_rows, num_cols);
  gsl_matrix_set_zero(out.gsl_mat_.get());
  return out;
}

Matrix Matrix::uniform(size_t num_rows, size_t num_cols, double value) {
  Matrix out(num_rows, num_cols);
  gsl_matrix_set_all(out.gsl_mat_.get(), value);
  return out;
}

Matrix Matrix::identity(size_t n) {
  Matrix out(n, n);
  gsl_matrix_set_identity(out.gsl_mat_.get());
  return out;
}

Matrix Matrix::diagonal(const Vector& vec) {
  auto out = Matrix::uniform(vec.size(), vec.size(), 0.);
  for (size_t i = 0; i < vec.size(); ++i)
    out(i, i) = vec(i);
  return out;
}

size_t Matrix::numColumns() const { return gsl_mat_->size2; }

size_t Matrix::numRows() const { return gsl_mat_->size1; }

Matrix Matrix::subset(size_t min_y, size_t min_x, size_t max_y, size_t max_x) const {
  const size_t num_x = max_x - min_y, num_y = max_y - min_y;
  Matrix out(num_y, num_x);
  gsl_matrix_const_submatrix(out.gsl_mat_.get(), min_y, min_x, num_y, num_x);
  return out;
}

bool Matrix::operator==(const Matrix& oth) const { return gsl_matrix_equal(gsl_mat_.get(), oth.gsl_mat_.get()) == 1; }

Matrix& Matrix::operator*=(double val) {
  gsl_matrix_scale(gsl_mat_.get(), val);
  return *this;
}

Matrix& Matrix::operator*=(const Vector& vec) {
  *this = *this * vec;
  return *this;
}

Matrix& Matrix::operator*=(const Matrix& mat) {
  *this = *this * mat;
  return *this;
}

Matrix& Matrix::operator/=(double val) { return operator*=(1. / val); }

Matrix Matrix::operator-() const { return zero(numRows(), numColumns()) - *this; }

Matrix& Matrix::operator+=(const Matrix& oth) {
  gsl_matrix_add(gsl_mat_.get(), oth.gsl_mat_.get());
  return *this;
}

Matrix& Matrix::operator-=(const Matrix& oth) {
  gsl_matrix_sub(gsl_mat_.get(), oth.gsl_mat_.get());
  return *this;
}

Vector Matrix::operator%(const Vector& vec) const {
  Matrix cpy(*this);
  Permutation perm(gsl_permutation_alloc(vec.size()), gsl_permutation_free);
  gsl_vector_view vec_v = gsl_matrix_column(vec.gsl_mat_.get(), 0);
  int s;
  gsl_linalg_LU_decomp(cpy.gsl_mat_.get(), perm.get(), &s);
  std::vector out(vec.size(), 0.);
  gsl_vector_view res_v = gsl_vector_view_array(out.data(), out.size());
  gsl_linalg_LU_solve(cpy.gsl_mat_.get(), perm.get(), &vec_v.vector, &res_v.vector);
  return Vector(static_cast<VectorRef&>(res_v));
}

double& Matrix::operator()(size_t ix, size_t iy) {
  if (ix >= numRows() || iy >= numColumns())
    throw CG_ERROR("Matrix:operator()") << "Invalid coordinates for " << numColumns() << "*" << numRows()
                                        << " matrix: (" << ix << ", " << iy << ").";
  return *gsl_matrix_ptr(gsl_mat_.get(), ix, iy);
}

double Matrix::operator()(size_t ix, size_t iy) const {
  if (ix >= numRows() || iy >= numColumns())
    throw CG_ERROR("Matrix:operator()") << "Invalid coordinates for " << numColumns() << "*" << numRows()
                                        << " matrix: (" << ix << ", " << iy << ").";
  return gsl_matrix_get(gsl_mat_.get(), ix, iy);
}

Matrix::Indices Matrix::imin() const {
  size_t imax, jmax;
  gsl_matrix_min_index(gsl_mat_.get(), &imax, &jmax);
  return std::make_pair(imax, jmax);
}

Matrix::Indices Matrix::imax() const {
  size_t imax, jmax;
  gsl_matrix_max_index(gsl_mat_.get(), &imax, &jmax);
  return std::make_pair(imax, jmax);
}

double Matrix::min() const { return gsl_matrix_min(gsl_mat_.get()); }

double Matrix::max() const { return gsl_matrix_max(gsl_mat_.get()); }

bool Matrix::null() const { return gsl_matrix_isnull(gsl_mat_.get()) == 1; }

bool Matrix::positive() const { return gsl_matrix_ispos(gsl_mat_.get()) == 1; }

bool Matrix::negative() const { return gsl_matrix_isneg(gsl_mat_.get()) == 1; }

bool Matrix::nonNegative() const { return gsl_matrix_isnonneg(gsl_mat_.get()) == 1; }

Matrix& Matrix::truncate(double min) {
  for (size_t iy = 0; iy < numRows(); ++iy)
    for (size_t ix = 0; ix < numColumns(); ++ix) {
      auto& val = operator()(iy, ix);
      if (val < min)
        val = 0.;
    }
  return *this;
}

Matrix& Matrix::transpose() {
  //gsl_matrix_transpose(gsl_mat_.get()); // only works for square matrices in GSL
  auto transposed = Matrix(numColumns(), numRows());
  for (size_t ix = 0; ix < numRows(); ++ix)
    for (size_t iy = 0; iy < numColumns(); ++iy)
      transposed(iy, ix) = operator()(ix, iy);
  *this = transposed;
  return *this;
}

Matrix Matrix::transposed() const { return Matrix(*this).transpose(); }

Matrix& Matrix::invert() {
  if (numRows() != numColumns()) {
    CG_WARNING("Matrix:inverted") << "Matrix inversion only works on square matrices.";
    return *this;
  }
  Matrix cpy(*this);
  Permutation perm(gsl_permutation_alloc(numRows()), gsl_permutation_free);
  int s;
  gsl_linalg_LU_decomp(cpy.gsl_mat_.get(), perm.get(), &s);
  gsl_linalg_LU_invert(cpy.gsl_mat_.get(), perm.get(), gsl_mat_.get());
  return *this;
}

Matrix Matrix::inverted() const { return Matrix(*this).invert(); }

Vector Matrix::column(size_t i) const { return gsl_matrix_const_column(gsl_mat_.get(), i); }

VectorRef Matrix::column(size_t i) { return static_cast<VectorRef>(gsl_matrix_column(gsl_mat_.get(), i)); }

Vector Matrix::row(size_t i) const { return gsl_matrix_const_row(gsl_mat_.get(), i); }

VectorRef Matrix::row(size_t i) { return static_cast<VectorRef>(gsl_matrix_row(gsl_mat_.get(), i)); }

Vector Matrix::diagonal() const { return gsl_matrix_const_diagonal(gsl_mat_.get()); }

VectorRef Matrix::diagonal() { return static_cast<VectorRef>(gsl_matrix_diagonal(gsl_mat_.get())); }

//--- friend operators

Matrix operator/(const Matrix& lhs, double val) { return lhs * (1. / val); }

namespace cepgen {
  Matrix operator+(const Matrix& lhs, const Matrix& rhs) {
    Matrix out(lhs);
    out += rhs;
    return out;
  }

  Matrix operator-(const Matrix& lhs, const Matrix& rhs) {
    Matrix out(lhs);
    out -= rhs;
    return out;
  }

  Matrix operator*(const Matrix& lhs, double val) {
    Matrix out(lhs);
    out *= val;
    return out;
  }

  Matrix operator*(double val, const Matrix& lhs) {
    Matrix out(lhs);
    out *= val;
    return out;
  }

  Vector operator*(const Matrix& mat, const Vector& vec) {
    gsl_vector_view vec_v = gsl_matrix_column(vec.gsl_mat_.get(), 0);
    std::vector out(mat.numRows(), 0.);
    gsl_vector_view res_v = gsl_vector_view_array(out.data(), out.size());
    gsl_blas_dgemv(CblasNoTrans, 1., mat.gsl_mat_.get(), &vec_v.vector, 0., &res_v.vector);
    return Vector(static_cast<VectorRef&>(res_v));
  }

  Matrix operator*(const Matrix& mat1, const Matrix& mat2) {
    Matrix out(mat1.numRows(), mat2.numColumns());
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., mat1.gsl_mat_.get(), mat2.gsl_mat_.get(), 0., out.gsl_mat_.get());
    return out;
  }

  std::ostream& operator<<(std::ostream& os, const Matrix& mat) {
    os << "(";
    std::string sep;
    for (size_t i = 0; i < mat.numRows(); ++i)
      os << sep << mat.row(i), sep = "\n ";
    return os << ")";
  }
}  // namespace cepgen

//---------------------------------------------------------------------------

VectorRef::VectorRef(const gsl_vector_view& view) : gsl_vector_view(view) {}

VectorRef& VectorRef::operator=(const Vector& vec) {
  for (size_t i = 0; i < vec.size(); ++i)
    operator()(i) = vec(i);
  return *this;
}

VectorRef::operator Vector() const { return Vector(*this); }

bool VectorRef::operator==(const Vector& vec) const { return Vector(*this) == vec; }

double& VectorRef::operator()(size_t i) { return *gsl_vector_ptr(&vector, i); }

double VectorRef::operator()(size_t i) const { return gsl_vector_get(&vector, i); }

std::ostream& operator<<(std::ostream& os, const VectorRef& ref) { return os << "Ref" << Vector(ref); }

//---------------------------------------------------------------------------

Vector::Vector(size_t num_coord, double def) : Matrix(uniform(num_coord, 1, def)) {}

Vector::Vector(const std::initializer_list<double>& vec) : Matrix(vec.size(), 1) {
  size_t i = 0;
  for (const auto& elem : vec)
    operator()(i++) = elem;
}

Vector::Vector(const VectorRef& ref) : Matrix(ref.vector.size, 1) {
  for (size_t i = 0; i < ref.vector.size; ++i)
    Matrix::operator()(i, 0) = ref(i);
}

Vector::Vector(const gsl_vector_const_view& vec) : Matrix(vec.vector.size, 1) {
  for (size_t i = 0; i < vec.vector.size; ++i)
    Matrix::operator()(i, 0) = gsl_vector_get(&vec.vector, i);
}

size_t Vector::size() const { return numRows(); }

Vector Vector::subset(size_t min, size_t max) const { return Vector(Matrix::subset(min, 0, max, 0).column(0)); }

double& Vector::operator()(size_t i) { return Matrix::operator()(i, 0); }

double Vector::operator()(size_t i) const { return Matrix::operator()(i, 0); }

double Vector::dot(const Vector& oth) const {
  if (size() != oth.size()) {
    CG_WARNING("Vector:dot") << "Scalar product of two vectors only defined for same-length vectors.";
    return 0.;
  }
  double out = 0.;
  for (size_t i = 0; i < size(); ++i)
    out += operator()(i) * oth(i);
  return out;
}

Vector Vector::cross(const Vector& oth) const {
  if (size() != oth.size())
    throw CG_FATAL("Vector:cross") << "Vector/cross product of two vectors only defined for same-length vectors.";
  if (oth.size() != 3)
    throw CG_FATAL("Vector:cross") << "Vector product only implemented for 3-vectors.";
  //FIXME extend to N-dimensional vectors
  return Vector{(operator()(1) * oth(2) - operator()(2) * oth(1)),
                (operator()(2) * oth(0) - operator()(0) * oth(2)),
                (operator()(0) * oth(1) - operator()(1) * oth(0))};
}

namespace cepgen {
  std::ostream& operator<<(std::ostream& os, const Vector& vec) {
    os << "(";
    std::string sep;
    for (size_t i = 0; i < vec.size(); ++i)
      os << sep << vec(i), sep = ", ";
    return os << ")";
  }
}  // namespace cepgen
