#ifndef EXVECTRMATH_SPECIALISED_TYPES_H
#define EXVECTRMATH_SPECIALISED_TYPES_H

#include "stdint.h"

#include "matrix_base.hpp"
#include "matrix_vector.hpp"

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

        /*template<>
        class Vector<float, 3> : public Vector<float, 3> {
        private:


            Vector();

            Vector(const float &value);

            /**
             * @brief Construct a new Vector object
             *
             * @tparam TYPE2
             * @param vectorMatrix
             */
            /*template <typename TYPE2>
            Vector(const Matrix<TYPE2, 3, 1> &vectorMatrix);

            Vector<float, 4> to4DVec(const float& forthVal = 0);


        };

        template<>
        Vector<float, 3>::Vector() : Vector<float, 3>::Vector() {}

        template<>
        Vector<float, 3>::Vector(const float &value) : Vector<float, 3>::Vector(value) {}

        template<>
        template <typename TYPE2>
        Vector<float, 3>::Vector(const Matrix<TYPE2, 3, 1> &vectorMatrix) : Vector<float, 3>::Vector(vectorMatrix) {}

        template<>
        Vector<float, 4> Vector<float, 3>::to4DVec(const float& forthVal = 0) {
            Vector<TYPE, 4> ret = *this;
            ret(3) = forthVal;
            return ret;
        }*/

    }

} // Namespace end

#endif