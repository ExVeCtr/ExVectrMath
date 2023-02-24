#ifndef EXVECTRMATH_MATRIXTYPES_H
#define EXVECTRMATH_MATRIXTYPES_H

#include "stdint.h"

#include "matrix_base.hpp"

namespace VCTR
{
    namespace Math
    {

        /**
         * Creates a matrix with given size and floats as internal type.
         */
        template <size_t ROWS, size_t COLS>
        using Matrix_F = Matrix<float, ROWS, COLS>;

        /**
         * Creates a matrix with given size and doubles as internal type.
         */
        template <size_t ROWS, size_t COLS>
        using Matrix_D = Matrix<double, ROWS, COLS>;

    }
} // Namespace end

#endif