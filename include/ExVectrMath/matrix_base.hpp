#ifndef EXVECTRMATH_MATRIXBASE_H
#define EXVECTRMATH_MATRIXBASE_H

#include "stdint.h"
#include "math.h"
#include <initializer_list>

namespace VCTR
{
    namespace Math
    {

        /**
         * General Matrix class used for matrix calculations. Can be any type (float, double, int etc.) and any size.
         * Recommended to used the from this derived forms for 3D, 4D etc matrices and specifically Vector and Quaternion for higher performance and more features.
         */
        template <typename TYPE, size_t ROWS, size_t COLS>
        class Matrix
        {
        public:
            /**
             * These are the raw internal values.
             * Access via matrix(row, column) is safer but adds slight overhead.
             * Can be accessed alternatively with matrix[row][column] with (presumably) 0 overhead, but no out of bounds checking.
             */
            TYPE r[ROWS][COLS];

            /**
             * Below are constructors.
             */

            /**
             * Standard constructor that fills all values with 0.
             */
            Matrix();

            /**
             * Contructor that creates a matrix with the diagonal filled with given values and the rest are 0.
             * Similar to identity matrix multplied by the given value.
             * If the entire array should be filled then use the static function Matrix::fill(value);
             * @param value What value to fill diagonal with.
             */
            Matrix(TYPE value);

            /**
             * @brief Construct a new Matrix object using initializer list.
             * @example Matrix<float, 3, 3> mat = {1, 2, 3, 4, 5, 6, 7, 8, 9};
             * @param array
             */
            Matrix(std::initializer_list<TYPE> array);

            /**
             * Conversion constructor. Converts given matrix into type that this matrix is (e.g float -> int)
             * @param matrix Which matrix to copy.
             */
            template <typename TYPE2>
            Matrix(const Matrix<TYPE2, ROWS, COLS> &matrix);

            /**
             * Conversion constructor. Copies given matrix but only a part of it or all of it depending on size.
             */
            template <typename TYPE2, size_t ROWS2, size_t COLS2>
            Matrix(const Matrix<TYPE2, ROWS2, COLS2> &mat, size_t rowStart = 0, size_t colStart = 0);

            /**
             * Below are generators for often needed matrices.
             */

            /**
             * Generates a matrix where all diagonal values are set to given value and the rest are set to 0.
             * @param eyeVal Value the diagonal values should have.
             * @returns generated matrix.
             */
            static Matrix<TYPE, ROWS, COLS> eye(TYPE eyeVal = 1);

            /**
             * Generates a matrix that is filled with the given value.
             * @param fillVal Value to fill array with.
             * @returns generated matrix.
             */
            static Matrix<TYPE, ROWS, COLS> fill(TYPE fillVal);

            /**
             * Below are generial functions to get the size of array or copy it.
             */

            /**
             * @returns number of values a matrix has or ROWS * COLS
             */
            size_t getNumValues() const;

            /**
             * @returns a copy of the matrix.
             */
            Matrix<TYPE, ROWS, COLS> copy() const;

            /**
             * Special functions to calculate transpose, determinant, inverse etc.
             */

            /**
             * @returns transposed matrix.
             */
            Matrix<TYPE, COLS, ROWS> transpose() const;

            /**
             * @note NEEDS FIXING
             * @returns matrix determinant.
             */
            TYPE determinant();

            /**
             * @brief calculates the magnitude of a vector matrix
             *
             * @return TYPE
             */
            TYPE magnitude() const;

            /**
             * @brief normalizes a Nx1 vector matrix.
             *
             * @return Matrix<TYPE, ROWS, COLS>
             */
            Matrix<TYPE, ROWS, COLS> normalize() const;

            /**
             * Used Gauß-Jordan method from https://www.geeksforgeeks.org/finding-inverse-of-a-matrix-using-gauss-jordan-method/
             * @returns inverse of matrix.
             */
            Matrix<TYPE, ROWS, COLS> inverse() const;

            /**
             * @brief Get the Diagonal vector of a square matrix
             *
             * @return Matrix<TYPE, ROWS, COLS>
             */
            Matrix<TYPE, ROWS, 1> diagonal() const;

            /**
             * @brief Computes the cross product of two 3x1 matricies (3D vectors)
             * @param vecB 3D vector
             *
             * @return Matrix<TYPE, 3, 1> (3D vector)
             */
            template<typename TYPE2>
            Matrix<TYPE, 3, 1> cross(const Matrix<TYPE2, 3, 1> &vecB) const;

            /**
             * @brief Calculates this matrix to power given.
             * @note Only powers > 0!
             *
             * @return Matrix<TYPE, ROWS, COLS>
             */
            Matrix<TYPE, ROWS, COLS> pow(size_t power) const;

            /**
             * @brief Calculates e to power of this matrix.
             * @note Approximate calculation. Use increase param iterations for slower speed but better result.
             * @param iterations Precision of calculation. Higher is better but slower. Defaults to 3.
             *
             * @return Matrix<TYPE, ROWS, COLS>
             */
            Matrix<TYPE, ROWS, COLS> ePow(size_t iterations = 3) const;

            /**
             * Below are access operators to access matrix values.
             */

            /**
             * Another way to access the internal data.
             * Low overhead and direct access to internal data, but no protection against going out of bounds.
             */
            TYPE *operator[](size_t row);

            /**
             * Another way to access the internal data. Return value is const.
             * Low overhead and direct access to internal data, but no protection against going out of bounds.
             */
            const TYPE *operator[](size_t row) const;

            /**
             * This operator accesses the first row until the end then goes to the beginning of the next row.
             */
            TYPE &operator()(uint32_t index);

            /**
             * Return value is constant.
             * This operator accesses the first row until the end then goes to the beginning of the next row.
             */
            const TYPE &operator()(uint32_t index) const;

            /**
             * This operator is to access the matrix values while keeping you inside the bounds of data.
             */
            TYPE &operator()(size_t row, size_t column);

            /**
             * Return value is constant.
             * This operator is to access the matrix values while keeping you inside the bounds of data.
             */
            const TYPE &operator()(size_t row, size_t column) const;

            /**
             * Returns part of the matrix.
             * E.g. Matrix<3, 3>::getBlock<2,2>(1, 1); gets the positions marked x:
             * o o o
             * o x x
             * 0 x x
             */
            template <size_t ROWS2, size_t COLS2>
            Matrix<TYPE, ROWS2, COLS2> block(size_t rowStart, size_t colStart) const;

            /**
             * Sets part of the matrix using the given ones values.
             * E.g. Matrix<3, 3>::block<2,2>(Matrix<2, 2>, 1, 1); sets the positions marked x:
             * o o o
             * o x x
             * 0 x x
             * using Matrix<2, 2>'s values.
             */
            template <typename TYPE2, size_t ROWS2, size_t COLS2>
            Matrix<TYPE, ROWS, COLS> &block(const Matrix<TYPE2, ROWS2, COLS2> &mat, size_t rowStart, size_t colStart);

            /**
             * Uses given function to print itsself.
             * @param printf Function to receive a const char*
             */
            void printTo(void (*printf)(const char *, ...)) const;

            /**
             * Below are math operations.
             */

            /**
             * Negates all values in the matrix.
             */
            Matrix<TYPE, ROWS, COLS> operator-() const;

            /**
             * Adds 2 matrices together.
             */
            Matrix<TYPE, ROWS, COLS> operator+(const Matrix<TYPE, ROWS, COLS> &right) const;

            /**
             * Subtracts right matrix from left.
             */
            Matrix<TYPE, ROWS, COLS> operator-(const Matrix<TYPE, ROWS, COLS> &right) const;

            /**
             * Generates an eye matrix from given value on left und then subtracts right matrix from generated value.
             */
            // friend Matrix<TYPE, ROWS, COLS> operator - (const TYPE& val, const Matrix<TYPE, ROWS, COLS>& right);

            /**
             * Multiplies a scaler with the matrix.
             */
            Matrix<TYPE, ROWS, COLS> operator*(const TYPE &scaler) const;

            /**
             * Divides a scaler with the matrix.
             */
            Matrix<TYPE, ROWS, COLS> operator/(const TYPE &scaler) const;

            /**
             * Multiplies a scaler with the matrix just like matrix*float but allows to be swaped to float*matrix.
             */
            // friend Matrix<TYPE, ROWS, COLS> operator * (const TYPE& scaler, const Matrix<TYPE, ROWS, COLS>& right);

            /**
             * Does matrix multiplication with given matrices. Recommended to use auto so save from having to determine matrix size und writing out.
             * E.g: auto C = A*B;
             * This is as normal matrix multiplication not commutative meaning: A*B will result in something different than B*A.
             */
            template <size_t RIGHTCOLS>
            Matrix<TYPE, ROWS, RIGHTCOLS> operator*(const Matrix<TYPE, COLS, RIGHTCOLS> &right) const;

            /**
             * Copies value from one matrix into the other even with differing internal value types.
             */
            template <typename TYPE2>
            Matrix<TYPE, ROWS, COLS> operator=(const Matrix<TYPE2, ROWS, COLS> &right);

            /**
             * Places array into matrix
             */
            // template<typename TYPE2>
            // Matrix<TYPE, ROWS, COLS> operator = (const TYPE2* array);

            /**
             * @brief If the two matricies have the exact same values
             *
             * @tparam TYPE2
             * @param right
             * @return true
             * @return false
             */
            template <typename TYPE2>
            bool operator==(const Matrix<TYPE2, ROWS, COLS> &right);

            /**
             * @brief If any of the values in the two matricies are different
             *
             * @tparam TYPE2
             * @param right
             * @return true
             * @return false
             */
            template <typename TYPE2>
            bool operator!=(const Matrix<TYPE2, ROWS, COLS> &right);
        };

        /**
         * Below are the implementations of all above functions.
         */

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, ROWS, COLS>::Matrix()
        {

            static_assert((ROWS > 0 && COLS > 0), "Matrix must be at least 1X1. Number of Rows or Columns or both are 0.");

            for (size_t row = 0; row < ROWS; row++)
            {
                for (size_t col = 0; col < COLS; col++)
                {

                    r[row][col] = 0;
                }
            }
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, ROWS, COLS>::Matrix(TYPE val)
        {

            static_assert((ROWS > 0 && COLS > 0), "Matrix must be at least 1X1. Number of Rows or Columns or both are 0.");

            for (size_t row = 0; row < ROWS; row++)
            {
                for (size_t col = 0; col < COLS; col++)
                {

                    if (row == col)
                        r[row][col] = val;
                    else
                        r[row][col] = 0;
                }
            }
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, ROWS, COLS>::Matrix(std::initializer_list<TYPE> array)
        {

            if (array.size() != ROWS * COLS)
                return;

            auto it = array.begin();
            for (size_t row = 0; row < ROWS; row++)
            {
                for (size_t col = 0; col < COLS; col++)
                {

                    r[row][col] = *it;
                    it++;
                }
            }
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        template <typename TYPE2>
        Matrix<TYPE, ROWS, COLS>::Matrix(const Matrix<TYPE2, ROWS, COLS> &matrix)
        {

            for (size_t row = 0; row < ROWS; row++)
            {
                for (size_t col = 0; col < COLS; col++)
                {

                    r[row][col] = matrix.r[row][col];
                }
            }
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        template <typename TYPE2, size_t ROWS2, size_t COLS2>
        Matrix<TYPE, ROWS, COLS>::Matrix(const Matrix<TYPE2, ROWS2, COLS2> &mat, size_t rowStart, size_t colStart) : Matrix()
        {
            block(mat, rowStart, colStart);
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, ROWS, COLS> Matrix<TYPE, ROWS, COLS>::eye(TYPE eyeVal)
        {

            Matrix<TYPE, ROWS, COLS> m;

            for (size_t i = 0; i < ROWS && i < COLS; i++)
                m.r[i][i] = eyeVal;

            return m;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, ROWS, COLS> Matrix<TYPE, ROWS, COLS>::fill(TYPE fillVal)
        {

            Matrix<TYPE, ROWS, COLS> m;

            for (size_t row = 0; row < ROWS; row++)
            {
                for (size_t col = 0; col < COLS; col++)
                {

                    m.r[row][col] = fillVal;
                }
            }

            return m;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        size_t Matrix<TYPE, ROWS, COLS>::getNumValues() const
        {
            return (uint32_t)ROWS * COLS;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, ROWS, COLS> Matrix<TYPE, ROWS, COLS>::copy() const
        {
            return *this;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, COLS, ROWS> Matrix<TYPE, ROWS, COLS>::transpose() const
        {

            Matrix<TYPE, COLS, ROWS> m;

            for (size_t row = 0; row < ROWS; row++)
            {
                for (size_t col = 0; col < COLS; col++)
                {

                    m.r[col][row] = r[row][col];
                }
            }

            return m;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, ROWS, COLS> Matrix<TYPE, ROWS, COLS>::inverse() const
        {

            static_assert((ROWS == COLS), "Matrix must be square (NxN) in order to get inverse. Pseudoinverse might solve the issue.");

            TYPE temp;

            Matrix<TYPE, ROWS, COLS * 2> matrix;
            matrix.block(*this, 0, 0);

            // Create the augmented matrix
            // Add the identity matrix
            // of order at the end of original matrix.
            for (int i = 0; i < ROWS; i++)
            {
                matrix[i][i + COLS] = 1;
            }

            // Interchange the row of matrix,
            // interchanging of row will start from the last row
            for (int i = ROWS - 1; i > 0; i--)
            {

                // Swapping each and every element of the two rows
                if (matrix[i - 1][0] < matrix[i][0])
                {
                    for (int j = 0; j < 2 * ROWS; j++)
                    {
                        temp = matrix[i][j];
                        matrix[i][j] = matrix[i - 1][j];
                        matrix[i - 1][j] = temp;
                    }
                }

                // Directly swapping the rows using pointers saves
                // time

                /*if (matrix[i - 1][0] < matrix[i][0]) {
                    TYPE* temp = matrix[i];
                    matrix[i] = matrix[i - 1];
                    matrix[i - 1] = temp;
                }*/
            }

            // Replace a row by sum of itself and a
            // constant multiple of another row of the matrix
            for (int i = 0; i < ROWS; i++)
            {

                for (int j = 0; j < ROWS; j++)
                {

                    if (j != i)
                    {

                        temp = matrix[j][i] / matrix[i][i];
                        for (int k = 0; k < 2 * ROWS; k++)
                        {

                            matrix[j][k] -= matrix[i][k] * temp;
                        }
                    }
                }
            }

            // Multiply each row by a nonzero integer.
            // Divide row element by the diagonal element
            for (int i = 0; i < ROWS; i++)
            {

                temp = matrix[i][i];
                for (int j = 0; j < 2 * COLS; j++)
                {

                    matrix[i][j] = matrix[i][j] / temp;
                }
            }

            Matrix<TYPE, ROWS, COLS> result;
            for (int i = 0; i < ROWS; i++)
            {

                for (int j = 0; j < COLS; j++)
                {

                    result[i][j] = matrix[i][j + COLS];
                }
            }

            return result;
            
        }

        /*template<>
        float Matrix<float, 1, 1>::determinant() {
            return r[0][0];
        }


        template<>
        double Matrix<double, 1, 1>::determinant() {
            return r[0][0];
        }


        template<>
        int32_t Matrix<int32_t, 1, 1>::determinant() {
            return r[0][0];
        }*/

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, ROWS, 1> Matrix<TYPE, ROWS, COLS>::diagonal() const
        {

            static_assert((COLS == ROWS), "Matrix must be square (NxN) to have a diagonal vector");

            Matrix<TYPE, ROWS, 1> diag;

            for (size_t i = 0; i < ROWS; i++)
            {

                diag.r[i][0] = this->r[i][i];
            }

            return diag;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        template <typename TYPE2>
        Matrix<TYPE, 3, 1> Matrix<TYPE, ROWS, COLS>::cross(const Matrix<TYPE2, 3, 1> &vecB) const
        {

            static_assert((ROWS == 3), "Cross product is only defined for 3x1 vectors.");

            Matrix<TYPE, 3, 1> result;
            result.r[0][0] = this->r[1][0] * vecB.r[2][0] - this->r[2][0] * vecB.r[1][0];
            result.r[1][0] = this->r[2][0] * vecB.r[0][0] - this->r[0][0] * vecB.r[2][0];
            result.r[2][0] = this->r[0][0] * vecB.r[1][0] - this->r[1][0] * vecB.r[0][0];

            return result;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, ROWS, COLS> Matrix<TYPE, ROWS, COLS>::pow(size_t power) const
        {

            static_assert((COLS == ROWS), "Matrix must be square (NxN) to calculate power!");

            Matrix<TYPE, ROWS, COLS> m = eye();

            for (size_t i = 0; i < power; i++)
            {
                m = m * (*this);
            }

            return m;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, ROWS, COLS> Matrix<TYPE, ROWS, COLS>::ePow(size_t iterations) const
        {

            static_assert((COLS == ROWS), "Matrix must be square (NxN) to calculate e to power!");

            Matrix<TYPE, ROWS, COLS> m;

            for (size_t i = 0; i < iterations; i++)
            {

                size_t factorial = 1;
                for (size_t j = 0; j < i; j++)
                    factorial *= j;

                m = m + this->pow(i) / factorial;
            }

            return m;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        TYPE Matrix<TYPE, ROWS, COLS>::determinant()
        {

            static_assert((COLS == ROWS), "Matrix must be square (NxN) to have a determinant");

            if (COLS == 2 && ROWS == 2) // Do the 2x2 case explicitly
                return r[0][0] * r[1][1] - r[1][0] * r[0][1];

            if (COLS == 1 && ROWS == 1) // Do the 1x1 case explicitly
                return r[0][0];

            TYPE determinant = 0;
            Matrix<TYPE, ROWS - 1, COLS - 1> temp;

            for (unsigned int j1 = 0; j1 < ROWS; ++j1)
            {
                for (unsigned int i = 1; i < ROWS; ++i)
                {
                    unsigned int j2 = 0;
                    for (unsigned int j = 0; j < ROWS; ++j)
                    {
                        if (j == j1)
                            continue;
                        temp.r[i - 1][j2] = r[i][j];
                        j2++;
                    }
                }
                determinant += ::pow(-1.0, j1 + 2.0) * r[0][j1] * temp.determinant();
            }

            return determinant;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        TYPE Matrix<TYPE, ROWS, COLS>::magnitude() const
        {

            static_assert((COLS == 1), "Matrix must be vector (Nx1) to have a magnitude");

            TYPE buf = 0;

            for (size_t i = 0; i < ROWS; i++)
            {

                buf += r[i][0] * r[i][0];
            }

            return TYPE(sqrt(buf));
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, ROWS, COLS> Matrix<TYPE, ROWS, COLS>::normalize() const
        {

            static_assert((COLS == 1), "Matrix must be vector (Nx1) to have be normalized");

            Matrix<TYPE, ROWS, COLS> copy = *this;

            TYPE mag = copy.magnitude();

            copy = copy / mag;

            return copy;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        TYPE *Matrix<TYPE, ROWS, COLS>::operator[](size_t row)
        {
            return r[row];
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        const TYPE *Matrix<TYPE, ROWS, COLS>::operator[](size_t row) const
        {
            return r[row];
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        TYPE &Matrix<TYPE, ROWS, COLS>::operator()(uint32_t index)
        {
            return r[(index / COLS) % ROWS][index % COLS];
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        const TYPE &Matrix<TYPE, ROWS, COLS>::operator()(uint32_t index) const
        {
            return r[(index / COLS) % ROWS][index % COLS];
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        TYPE &Matrix<TYPE, ROWS, COLS>::operator()(size_t row, size_t column)
        {
            return r[row % ROWS][column % COLS];
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        const TYPE &Matrix<TYPE, ROWS, COLS>::operator()(size_t row, size_t column) const
        {
            return r[row % ROWS][column % COLS];
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        template <size_t ROWS2, size_t COLS2>
        Matrix<TYPE, ROWS2, COLS2> Matrix<TYPE, ROWS, COLS>::block(size_t rowStart, size_t colStart) const
        {
            Matrix<TYPE, ROWS2, COLS2> blockMat;

            for (size_t rw = rowStart; rw < ROWS && rw < ROWS2 + rowStart; rw++)
            {
                for (size_t cw = colStart; cw < COLS && cw < COLS2 + colStart; cw++)
                {
                    blockMat.r[rw - rowStart][cw - colStart] = r[rw][cw];
                }
            }

            return blockMat;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        template <typename TYPE2, size_t ROWS2, size_t COLS2>
        Matrix<TYPE, ROWS, COLS> &Matrix<TYPE, ROWS, COLS>::block(const Matrix<TYPE2, ROWS2, COLS2> &mat, size_t rowStart, size_t colStart)
        {

            for (size_t rw = rowStart; rw < ROWS && rw - rowStart < ROWS2; rw++)
            {

                for (size_t cw = colStart; cw < COLS && cw - colStart < COLS2; cw++)
                {

                    this->r[rw][cw] = mat.r[rw - rowStart][cw - colStart];
                }
            }

            return *this;
        }

        /**
         * Uses given function to print itself. Given print function is expected to work like a normal printf function.
         * @param printf Function pointer to receive a const char* for format and parameters.
         */
        template <typename TYPE, size_t ROWS, size_t COLS>
        void Matrix<TYPE, ROWS, COLS>::printTo(void (*printf)(const char *, ...)) const
        {

            printf("Matrix: \n[");

            for (size_t row = 0; row < ROWS; row++)
            {
                printf(" ");
                for (size_t col = 0; col < COLS; col++)
                {

                    printf(" %f, ", float((*this)[row][col]));
                }
                printf("\n");
            }

            printf("]\n");
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, ROWS, COLS> Matrix<TYPE, ROWS, COLS>::operator-() const
        {
            Matrix<TYPE, ROWS, COLS> m;

            for (size_t row = 0; row < ROWS; row++)
            {
                for (size_t col = 0; col < COLS; col++)
                {

                    m.r[row][col] = -r[row][col];
                }
            }

            return m;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, ROWS, COLS> Matrix<TYPE, ROWS, COLS>::operator+(const Matrix<TYPE, ROWS, COLS> &right) const
        {

            Matrix<TYPE, ROWS, COLS> m;

            for (size_t row = 0; row < ROWS; row++)
            {
                for (size_t col = 0; col < COLS; col++)
                {

                    m.r[row][col] = r[row][col] + right.r[row][col];
                }
            }

            return m;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, ROWS, COLS> Matrix<TYPE, ROWS, COLS>::operator-(const Matrix<TYPE, ROWS, COLS> &right) const
        {

            Matrix<TYPE, ROWS, COLS> m;

            for (size_t row = 0; row < ROWS; row++)
            {
                for (size_t col = 0; col < COLS; col++)
                {

                    m.r[row][col] = r[row][col] - right.r[row][col];
                }
            }

            return m;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, ROWS, COLS> operator-(const TYPE &val, const Matrix<TYPE, ROWS, COLS> &right)
        {

            Matrix<TYPE, ROWS, COLS> m;

            for (size_t row = 0; row < ROWS; row++)
            {
                for (size_t col = 0; col < COLS; col++)
                {

                    if (row == col)
                        m.r[row][col] = val - right.r[row][col];
                    else
                        m.r[row][col] = -right.r[row][col];
                }
            }

            return m;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, ROWS, COLS> Matrix<TYPE, ROWS, COLS>::operator*(const TYPE &scaler) const
        {

            Matrix<TYPE, ROWS, COLS> m;

            for (size_t row = 0; row < ROWS; row++)
            {
                for (size_t col = 0; col < COLS; col++)
                {

                    m.r[row][col] = this->r[row][col] * scaler;
                }
            }

            return m;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, ROWS, COLS> Matrix<TYPE, ROWS, COLS>::operator/(const TYPE &scaler) const
        {

            Matrix<TYPE, ROWS, COLS> m;

            for (size_t row = 0; row < ROWS; row++)
            {
                for (size_t col = 0; col < COLS; col++)
                {

                    m.r[row][col] = this->r[row][col] / scaler;
                }
            }

            return m;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, ROWS, COLS> operator/(const TYPE &scaler, const Matrix<TYPE, ROWS, COLS> &right)
        {

            Matrix<TYPE, ROWS, COLS> m;

            for (size_t row = 0; row < ROWS; row++)
            {
                for (size_t col = 0; col < COLS; col++)
                {

                    m.r[row][col] = scaler / right.r[row][col];
                }
            }

            return m;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        Matrix<TYPE, ROWS, COLS> operator*(const TYPE &scaler, const Matrix<TYPE, ROWS, COLS> &right)
        {

            Matrix<TYPE, ROWS, COLS> m;

            for (size_t row = 0; row < ROWS; row++)
            {
                for (size_t col = 0; col < COLS; col++)
                {

                    m.r[row][col] = scaler * right.r[row][col];
                }
            }

            return m;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        template <size_t RIGHTCOLS>
        Matrix<TYPE, ROWS, RIGHTCOLS> Matrix<TYPE, ROWS, COLS>::operator*(const Matrix<TYPE, COLS, RIGHTCOLS> &right) const
        {

            Matrix<TYPE, ROWS, RIGHTCOLS> m;

            for (size_t lr = 0; lr < ROWS; lr++)
            {
                for (size_t rc = 0; rc < RIGHTCOLS; rc++)
                {

                    for (size_t i = 0; i < COLS; i++)
                    {

                        m.r[lr][rc] = m.r[lr][rc] + this->r[lr][i] * right.r[i][rc];
                    }
                }
            }

            return m;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        template <typename TYPE2>
        Matrix<TYPE, ROWS, COLS> Matrix<TYPE, ROWS, COLS>::operator=(const Matrix<TYPE2, ROWS, COLS> &right)
        {

            for (size_t r = 0; r < ROWS; r++)
            {
                for (size_t c = 0; c < COLS; c++)
                {

                    this->r[r][c] = right.r[r][c];
                }
            }

            return *this;
        }

        /*
        template<typename TYPE, size_t ROWS, size_t COLS>
        template<typename TYPE2>
        Matrix<TYPE, ROWS, COLS> Matrix<TYPE, ROWS, COLS>::operator = (const TYPE2* array) {

            for (size_t i = 0; i < getNumValues(); i++) (*this)(i) = static_cast<TYPE>(array[i]);

            return *this;

        }*/

        template <typename TYPE, size_t ROWS, size_t COLS>
        template <typename TYPE2>
        bool Matrix<TYPE, ROWS, COLS>::operator==(const Matrix<TYPE2, ROWS, COLS> &right)
        {

            for (size_t r = 0; r < ROWS; r++)
            {
                for (size_t c = 0; c < COLS; c++)
                {

                    if (this->r[r][c] != right.r[r][c])
                    {
                        return false;
                    }
                }
            }

            return true;
        }

        template <typename TYPE, size_t ROWS, size_t COLS>
        template <typename TYPE2>
        bool Matrix<TYPE, ROWS, COLS>::operator!=(const Matrix<TYPE2, ROWS, COLS> &right)
        {

            return !((*this) == right);
        }

    }
} // Namespace end

#endif