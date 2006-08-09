/*
 * mathfunctions.h
 *
 * Copyright (C) 2006 by Universitaet Stuttgart (VIS). Alle Rechte vorbehalten.
 */

#ifndef VISLIB_MATHFUNCTIONS_H_INCLUDED
#define VISLIB_MATHFUNCTIONS_H_INCLUDED
#if (_MSC_VER > 1000)
#pragma once
#endif /* (_MSC_VER > 1000) */


#include <cmath>
#include <cstdlib>
#include <limits>

// TODO

/** Epsilon value for floating point comparison. */
#define FLOAT_EPSILON (1e-5f)

/** Epsilon value for double precision floating point comparison. */
#define DOUBLE_EPSILON (1e-7)


namespace vislib {
namespace math {

    /**
     * Answer the next power of two which is greater or equal to 'n'.
     *
     * This function has logarithmic runtime complexity.
     *
     * Notes: DO NOT INSTANTIATE THE TEMPLATE FOR OTHER DATA TYPES THAN INTEGRAL
     *        NUMBERS (Should not compile for floating point types because of
     *        missing bit shift operators)!
     *
     * @param n A number.
     *
     * @return The smallest power of two with ('n' <= result).
     */
    template<class T> T CalcNextPowerOfTwo(const T n) {
        ASSERT(std::numeric_limits<T>::is_integer);
        
        T retval = static_cast<T>(1);

        while (retval < n) {
            retval <<= 1;
        }

        return retval;
    }


    /**
     * Clamp 'n' to ['minVal', 'maxVal'].
     *
     * @param n      A number.
     * @param minVal The minimum valid number.
     * @param maxVal The maximum valid number, must be greater than 'minVal'.
     *
     * @return The clamped value of 'n'.
     */
    template<class T> inline T Clamp(const T n, const T minVal, 
            const T maxVal) {
        return ((n >= minVal) ? ((n <= maxVal) ? n : maxVal) : minVal);
    }


    /**
     * Answer whether 'm' and 'n' are nearly equal.
     *
     * @param m       A floating point number.
     * @param n       A floating point number.
     * @param epsilon The epsilon value used for comparion. Defaults to
     *                FLOAT_EPSILON.
     *
     * @return true, if 'm' and 'n' are nearly equal, false otherwise.
     */
    inline bool FltEq(const float m, const float n, 
            const float epsilon = FLOAT_EPSILON) {
        return (::fabs(m - n) < epsilon);
    }


    /**
     * Answer whether 'm' and 'n' are nearly equal.
     *
     * @param m       A floating point number.
     * @param n       A floating point number.
     * @param epsilon The epsilon value used for comparion. Defaults to
     *                DOUBLE_EPSILON.
     *
     * @return true, if 'm' and 'n' are nearly equal, false otherwise.
     */
    inline bool FltEq(const double m, const double n, 
            const double epsilon = DOUBLE_EPSILON) {
        return (::fabs(m - n) < epsilon);
    }


    /**
     * Answer the maximum of 'n' and 'm'.
     *
     * @param n An number.
     * @param m Another number.
     *
     * @return The maximum of 'n' and 'm'.
     */
    template<class T> inline T Max(const T n, const T m) {
        return (n > m) ? n : m;
    }


    /**
     * Answer the minimum of 'n' and 'm'.
     *
     * @param n An number.
     * @param m Another number.
     *
     * @return The minimum of 'n' and 'm'.
     */
    template<class T> inline T Min(const T n, const T m) {
        return (n < m) ? n : m;
    }


    /**
     * Compute the square of 'n'.
     * 
     * @param n A number.
     *
     * @return The square of 'n'.
     */
    template<class T> inline T Sqr(const T n) {
        return (n * n);
    }


    /**
     * Answer the square root of 'n'.
     *
     * @param n A number.
     *
     * @return The square root of 'n'.
     */
    inline float Sqrt(const float n) {
        return ::sqrtf(n);
    }


    /**
     * Answer the square root of 'n'.
     *
     * @param n A number.
     *
     * @return The square root of 'n'.
     */
    inline double Sqrt(const double n) {
        return ::sqrt(n);
    }

} /* end namespace math */
} /* end namespace vislib */

#endif /* VISLIB_MATHFUNCTIONS_H_INCLUDED */
