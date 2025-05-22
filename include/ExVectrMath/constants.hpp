#ifndef EXVECTRMATH_CONSTANTS_H
#define EXVECTRMATH_CONSTANTS_H

#include "stdint.h"

#include "matrix_base.hpp"
#include "matrix_vector.hpp"

namespace VCTR
{
    namespace Math
    {

        /// @brief Gravity on earth in m/s/s
        const float GRAVITY = 9.807f;
        /// @brief Gravity vector on earth in m/s/s. Points up as the force is excerted on an object from what it lays on.
        const Math::Vector_F GRAVITY_3F({0.0f, 0.0f, GRAVITY});

        

    }

    constexpr float DEGREES = M_PI / 180.0f;

} // Namespace end

#endif