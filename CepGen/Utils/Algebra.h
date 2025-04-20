/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2024  Laurent Forthomme
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

#ifndef CepGen_Utils_Algebra_h
#define CepGen_Utils_Algebra_h

#include <gsl/gsl_matrix_double.h>

#include <initializer_list>
#include <memory>

namespace cepgen {
  class Vector;
  class VectorRef;
  /// A \f$n\times m\f$ matrix object
  class Matrix {
  public:
    /// Object constructor
    /// \param[in] num_rows number of (horizontal) rows for the matrix
    /// \param[in] num_cols number of (vertical) columns for the matrix
    explicit Matrix(size_t num_rows, size_t num_cols);
    /// Object constructor (from vectors)
    /// \param[in] vectors set of (vector) rows
    Matrix(const std::initializer_list<Vector>& vectors);

    Matrix(const Matrix&);             ///< Copy constructor
    Matrix& operator=(const Matrix&);  ///< Assignment operator
    operator Vector() const;           ///< Implicit conversion to vector

    /// Build a zeroed matrix
    /// \param[in] num_rows number of (horizontal) rows for the matrix
    /// \param[in] num_cols number of (vertical) columns for the matrix
    static Matrix zero(size_t num_rows, size_t num_cols = 0ull);
    /// Build a uniform matrix
    /// \param[in] num_rows number of (horizontal) rows for the matrix
    /// \param[in] num_cols number of (vertical) columns for the matrix
    /// \param[in] value uniform value for all the matrix components
    static Matrix uniform(size_t num_rows, size_t num_cols, double value = 1.);
    static Matrix identity(size_t);         ///< Build a (square) identity matrix
    static Matrix diagonal(const Vector&);  ///< Build a (square) diagonal matrix from its diagonal vector

    size_t numColumns() const;  ///< Number of (vertical) columns
    size_t numRows() const;     ///< Number of (horizontal) rows

    bool operator==(const Matrix&) const;                                           ///< Equality operator
    inline bool operator!=(const Matrix& oth) const { return !(operator==(oth)); }  ///< Inequality operator

    /// Extract a subset of the matrix as a new object
    /// \param[in] min_y first vertical index for the subset
    /// \param[in] min_x first horizontal index for the subset
    /// \param[in] max_y last vertical index for the subset
    /// \param[in] max_x last horizontal index for the subset
    Matrix subset(size_t min_y, size_t min_x, size_t max_y = 0ull, size_t max_x = 0ull) const;

    Matrix& operator*=(double);               ///< Multiplication by a scalar operator
    Matrix& operator*=(const Vector&);        ///< Multiplication by a vector operator
    Matrix& operator*=(const Matrix&);        ///< Multiplication by a matrix operator
    Matrix& operator/=(double);               ///< Division by a scalar operator
    Matrix operator-() const;                 ///< Unary inverse operator
    Matrix& operator+=(const Matrix&);        ///< Addition of another matrix
    Matrix& operator-=(const Matrix&);        ///< Subtraction of another matrix
    double& operator()(size_t, size_t);       ///< Component access operator
    double operator()(size_t, size_t) const;  ///< Component retrieval operator

    Vector operator%(const Vector&) const;  ///< Solving operator (from LU decomposition)

    friend Matrix operator*(double, const Matrix&);         ///< Multiplication of a matrix by a scalar
    friend Matrix operator*(const Matrix&, double);         ///< Multiplication of a matrix by a scalar
    friend Vector operator*(const Matrix&, const Vector&);  ///< Multiplication of a matrix by a vector
    friend Matrix operator*(const Matrix&, const Matrix&);  ///< Multiplication of a matrix by another matrix
    friend Matrix operator/(const Matrix&, double);         ///< Division of a matrix by a scalar
    friend Matrix operator+(const Matrix&, const Matrix&);  ///< Addition of two matrices
    friend Matrix operator-(const Matrix&, const Matrix&);  ///< Subtraction of two matrices

    using Indices = std::pair<size_t, size_t>;
    Indices imin() const;  ///< Index (row, column) of the minimum matrix element
    Indices imax() const;  ///< Index (row, column) of the maximum matrix element

    double min() const;  ///< Minimum matrix element
    double max() const;  ///< Maximum matrix element

    bool null() const;         ///< Is the matrix uniformly null?
    bool positive() const;     ///< Is the matrix positive-defined?
    bool negative() const;     ///< Is the matrix negative-defined?
    bool nonNegative() const;  ///< Is the matrix non-negative-defined?

    friend std::ostream& operator<<(std::ostream&, const Matrix&);  ///< Printout of matrix components

    Matrix& truncate(double min = 1.e-14);  ///< Truncate (specify minimum non-zero value) for all matrix components

    Matrix& transpose();        ///< Transpose the matrix
    Matrix transposed() const;  ///< Return a transposition of this matrix
    Matrix& invert();           ///< Invert the matrix
    Matrix inverted() const;    ///< Return the inverse of this matrix (LU decomposition)

    VectorRef column(size_t);     ///< Return whole column of the matrix
    Vector column(size_t) const;  ///< Return whole column of the matrix
    VectorRef row(size_t);        ///< Return whole row of the matrix
    Vector row(size_t) const;     ///< Return whole row of the matrix
    VectorRef diagonal();         ///< Return the diagonal components of the matrix
    Vector diagonal() const;      ///< Return the diagonal components of the matrix

  private:
    std::unique_ptr<gsl_matrix, void (*)(gsl_matrix*)> gsl_mat_;
  };

  class VectorRef : public gsl_vector_view {
  public:
    explicit VectorRef(const gsl_vector_view&);

    VectorRef& operator=(const Vector&);                                          ///< Assignment operator
    operator Vector() const;                                                      ///< Conversion operator to a vector
    bool operator==(const Vector&) const;                                         ///< Equality with a vector operator
    inline bool operator!=(const Vector& vec) const { return !operator==(vec); }  ///< Inequality with a vector operator
    double& operator()(size_t);                                                   ///< Component access operator
    double operator()(size_t) const;                                              ///< Component retrieval operator
    friend std::ostream& operator<<(std::ostream&, const VectorRef&);             ///< Printout operator

  private:
    friend class Vector;
  };

  /// Specialisation of an \f$m\times 1\f$ matrix
  class Vector final : public Matrix {
  public:
    /// Object constructor
    /// \param[in] num_coord vector multiplicity
    /// \param[in] def default value for all components
    explicit Vector(size_t num_coord, double def = 0.);
    Vector(const std::initializer_list<double>&);  ///< Build a vector from a {...} list of double precision floats
    Vector(const VectorRef&);                      ///< Build a vector from a GSL vector view
    Vector(const gsl_vector_const_view&);          ///< Build a vector from a constant GSL vector view

    size_t size() const;  ///< Vector multiplicity (number of lines)

    /// Extract a subset of the vector
    /// \param[in] min first index for the subset
    /// \param[in] max last index for the subset (last element if not specified)
    Vector subset(size_t min, size_t max = 0ull) const;

    double& operator()(size_t);         ///< Component access operator
    double operator()(size_t) const;    ///< Component retrieval operator
    double dot(const Vector&) const;    ///< Scalar product of two vectors
    Vector cross(const Vector&) const;  ///< Vector product of two vectors

    friend std::ostream& operator<<(std::ostream&, const Vector&);  ///< Printout of vector components
  };
}  // namespace cepgen

#endif
