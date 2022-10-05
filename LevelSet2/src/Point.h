//
//  Point.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/11/21.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_Point_h
#define LevelSetMethod2_Point_h

#include "DimensionTypes.h"
#include <array>
#include <cmath>

namespace lsm
{
        /// definition of Point<T, D>
        template<typename T, typename D>
        using Point = std::array<T, D::Dimension_>;

        template<typename D>
        using IntPoint = Point<int, D>;

        template<typename D>
        using DoublePoint = Point<double, D>;
        
        template<typename D>
        using FloatPoint = Point<float, D>;

        template<typename T>
        using Point2d = Point<T, TwoDimension>;

        template<typename T>
        using Point3d = Point<T, ThreeDimension>;

        using IntPoint2d    = IntPoint<TwoDimension>;
        using DoublePoint2d = DoublePoint<TwoDimension>;
        using FloatPoint2d  = FloatPoint<TwoDimension>;

        using IntPoint3d    = IntPoint<ThreeDimension>;
        using DoublePoint3d = DoublePoint<ThreeDimension>;
        using FloatPoint3d  = FloatPoint<ThreeDimension>;

        /// definition of operator==
        template<typename T>
        inline bool operator==(
                const Point<T, TwoDimension>&   lhs,
                const Point<T, TwoDimension>&   rhs)
        {
                return lhs[0] == rhs[0] && lhs[1] == rhs[1];
        }

        template<typename T>
        inline bool operator==(
                const Point<T, ThreeDimension>& lhs,
                const Point<T, ThreeDimension>& rhs)
        {
                return lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2];
        }

        // test ok
        /// definition of operator+= with TwoDimension
        template<typename T>
        inline Point<T, TwoDimension>& operator+=(
                Point<T, TwoDimension>&         lhs,
                const Point<T, TwoDimension>&   rhs)
        {
                lhs[0] += rhs[0];
                lhs[1] += rhs[1];
                return lhs;
        }

        // test ok
        /// definition of operator+= with ThreeDimension
        template<typename T>
        inline Point<T, ThreeDimension>& operator+=(
                Point<T, ThreeDimension>&       lhs,
                const Point<T, ThreeDimension>& rhs)
        {
                lhs[0] += rhs[0];
                lhs[1] += rhs[1];
                lhs[2] += rhs[2];
                return lhs;
        }

        // test ok
        /// definition of operator-= with TwoDimension
        template<typename T>
        inline Point<T, TwoDimension>& operator-=(
                Point<T, TwoDimension>&         lhs,
                const Point<T, TwoDimension>&   rhs)
        {
                lhs[0] -= rhs[0];
                lhs[1] -= rhs[1];
                return lhs;
        }

        // test ok
        /// definition of operator-= with ThreeDimension
        template<typename T>
        inline Point<T, ThreeDimension>& operator-=(
                Point<T, ThreeDimension>&       lhs,
                const Point<T, ThreeDimension>& rhs)
        {
                lhs[0] -= rhs[0];
                lhs[1] -= rhs[1];
                lhs[2] -= rhs[2];
                return lhs;
        }

        // test ok
        /// definition of operator+ with TwoDimension
        template<typename T>
        inline const Point<T, TwoDimension> operator+(
                Point<T, TwoDimension>          lhs,
                const Point<T, TwoDimension>&   rhs)
        {
                return lhs += rhs;
        }

        /// definition of operator+ with ThreeDimension
        template<typename T>
        inline const Point<T, ThreeDimension> operator+(
                Point<T, ThreeDimension>        lhs,
                const Point<T, ThreeDimension>& rhs)
        {
                return lhs += rhs;
        }
        
        /// definition of operator* with ThreeDimension
        template<typename T>
        inline const Point<T, ThreeDimension> operator*(
                const T& t,
                const Point<T, ThreeDimension>& rhs)
        {
                return {{t * rhs[0], t * rhs[1], t * rhs[2]}};
        }

        // test ok
        /// definition of operator- with TwoDimension
        template<typename T>
        inline const Point<T, TwoDimension> operator-(
                Point<T, TwoDimension>          lhs,
                const Point<T, TwoDimension>&   rhs)
        {
                return lhs -= rhs;
        }

        /// definition of operator- with ThreeDimension
        template<typename T>
        inline const Point<T, ThreeDimension> operator-(
                Point<T, ThreeDimension>        lhs,
                const Point<T, ThreeDimension>& rhs)
        {
                return lhs -= rhs;
        }

        /// cast for int <-> double
        template<typename T, typename U>
        inline const Point2d<T> cast(const Point2d<U>& p)
        {
                return {{static_cast<T>(p[0]), static_cast<T>(p[1])}};
        }

        /// cast for int <-> double
        template<typename T, typename U>
        inline const Point3d<T> cast(const Point3d<U>& p)
        {
                return {{static_cast<T>(p[0]), static_cast<T>(p[1]), static_cast<T>(p[2])}};
        }

        // test ok
        template<typename T>
        inline T inner_product(
                const Point<T, TwoDimension>&   lhs,
                const Point<T, TwoDimension>&   rhs)
        {
                return lhs[0] * rhs[0] + lhs[1] * rhs[1];
        }

        // test ok
        template<typename T>
        inline T inner_product(
                const Point<T, ThreeDimension>& lhs,
                const Point<T, ThreeDimension>& rhs)
        {
                return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
        }

        // test ok
        template<typename T, typename D>
        inline double norm(const Point<T, D>& p)
        {
                return std::sqrt(inner_product(p, p));
        }

        template<typename T>
        inline T outer_product(
                const Point<T, TwoDimension>&   a,
                const Point<T, TwoDimension>&   b)
        {
                return a[0] * b[1] - a[1] * b[0];
        }

        template<typename T>
        inline const Point<T, ThreeDimension> outer_product(
                const Point<T, ThreeDimension>& a,
                const Point<T, ThreeDimension>& b)
        {
                return {{a[1] * b[2] - a[2] * b[1],
                        -a[0] * b[2] + a[2] * b[0],
                         a[0] * b[1] - a[1] * b[0]}};
        }
} // namespace lsm

#endif // LevelSetMethod2_Point_h
