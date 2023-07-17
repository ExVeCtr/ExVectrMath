#ifndef EXVECTRMATH_MATRIXVECTOR_H
#define EXVECTRMATH_MATRIXVECTOR_H

#include "matrix_base.hpp"

namespace VCTR
{
    namespace Math
    {

        /**
         * Nx1 Matrix class for vectors.
         */
        template <typename TYPE, size_t ROWS>
        class Vector : public Matrix<TYPE, ROWS, 1>
        {
        public:
            Vector();

            /**
             * @brief Construct a new Vector object and sets all elements to the given value.
             * @param value 
             */
            Vector(const TYPE &value);

            /**
             * @brief Construct a new vector object using initializer list.
             * @example Vector<float, 3> vec = {1.0f, 2.0f, 3.0f};
             * @param array 
             */
            Vector(std::initializer_list<TYPE> array);

            /**
             * @brief Construct a new Vector object
             *
             * @tparam TYPE2
             * @param vectorMatrix
             */
            template <typename TYPE2>
            Vector(const Matrix<TYPE2, ROWS, 1> &vectorMatrix);

            /**
             * @brief Normalizes the vector
             *
             * @return TYPE
             */
            Vector<TYPE, ROWS> normalize() const;

            /**
             * @brief Gets the angle between the vectors
             *
             * @tparam TYPE2
             * @param vecB
             * @return TYPE
             */
            template <typename TYPE2>
            TYPE getAngleTo(const Vector<TYPE2, ROWS> &vecB) const;

            /**
             * @brief Get the Projection of the vector onto the given vector parameter (vecB).
             *
             * @tparam TYPE2
             * @param vecB Vector to project onto.
             * @return TYPE
             */
            template <typename TYPE2>
            Vector<TYPE, ROWS> getProjectionOn(const Vector<TYPE2, ROWS> &vecB) const;

            /**
             * @brief Calculates the dot product between the two vectors.
             *
             * @tparam TYPE2
             * @param vecB
             * @return TYPE
             */
            template <typename TYPE2>
            TYPE operator*(const Vector<TYPE2, ROWS> &vecB) const;

        };

        template <typename TYPE, size_t ROWS>
        Vector<TYPE, ROWS>::Vector()
        {
            for (size_t i = 0; i < ROWS; i++)
                this->r[i][0] = 0.0f;
        }

        template <typename TYPE, size_t ROWS>
        Vector<TYPE, ROWS>::Vector(const TYPE &value)
        {
            for (size_t i = 0; i < ROWS; i++)
                this->r[i][0] = value;
        }

        template <typename TYPE, size_t ROWS>
        Vector<TYPE, ROWS>::Vector(std::initializer_list<TYPE> array)
        {
            size_t i = 0;
            for (auto &element : array)
            {
                this->r[i][0] = element;
                i++;
            }
        }

        template <typename TYPE, size_t ROWS>
        template <typename TYPE2>
        Vector<TYPE, ROWS>::Vector(const Matrix<TYPE2, ROWS, 1> &vectorMatrix)
        {
            for (size_t i = 0; i < ROWS; i++)
                this->r[i][0] = vectorMatrix[i][0];
        }

        template <typename TYPE, size_t ROWS>
        Vector<TYPE, ROWS> Vector<TYPE, ROWS>::normalize() const
        {

            auto copy = *this;

            TYPE mag = copy.magnitude();

            copy = copy / mag;

            return copy;

        }

        template <typename TYPE, size_t ROWS>
        template <typename TYPE2>
        TYPE Vector<TYPE, ROWS>::getAngleTo(const Vector<TYPE2, ROWS> &vecB) const
        {
            return acos(((*this) * vecB) / ((*this).magnitude() * vecB.magnitude()));
        }

        template <typename TYPE, size_t ROWS>
        template <typename TYPE2>
        Vector<TYPE, ROWS> Vector<TYPE, ROWS>::getProjectionOn(const Vector<TYPE2, ROWS> &vecB) const
        {
            return ((*this) * vecB) / (vecB * vecB) * vecB;
        }

        template <typename TYPE, size_t ROWS>
        template <typename TYPE2>
        TYPE Vector<TYPE, ROWS>::operator*(const Vector<TYPE2, ROWS> &vecB) const
        {
            TYPE val = 0;
            for (size_t i = 0; i < ROWS; i++)
                val += this->r[i][0] * vecB[i][0];
            return val;
        }

        /*
        template<typename TYPE, size_t ROWS>
        TYPE Vector<TYPE, ROWS>::magnitude(size_t startIndex, size_t endIndex) {

            TYPE temp = 0;

            for (size_t i = startIndex; i <= endIndex; i++) {

                temp = this->r[i][0]*this->r[i][0];

            }

            return sqrt(temp);

        }


        template<typename TYPE, size_t ROWS>
        Vector<TYPE, ROWS> Vector<TYPE, ROWS>::normalize(size_t startIndex, size_t endIndex) {

            TYPE mag = magnitude(startIndex, endIndex);

            Vector<TYPE, ROWS> normalisedVector;

            for (size_t i = startIndex; i <= endIndex; i++) {

                this->r[i][0] = this->r[i][0]/mag;

            }

        }
        */

        using Vector_F = Vector<float, 3>; // For float vector types

        using Vector_D = Vector<double, 3>; // For double vector types

    }
} // Namespace end

#endif