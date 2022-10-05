#ifndef LevelSetMethod2_NeighboringPoints_h
#define LevelSetMethod2_NeighboringPoints_h

#include "Point.h"
#include "Utilities.h"
#include "Concepts.h"

#include <boost/concept_check.hpp>
#include <vector>

namespace lsm
{
        /// defines neighboring points
        /**
         *      @tparam D the dimension type
         */
        template<typename D>
        struct NeighboringPoints
        {
                BOOST_CLASS_REQUIRE(D, lsm, DimensionConcept);
                NeighboringPoints();
                IntPoint<D> points_[power(3, D::Dimension_)];

                const IntPoint<D>& operator()(int x, int y) const
                {
                        return points_[(1 + x) + 3 * (1 + y)];
                }

                const IntPoint<D>& operator()(int x, int y, int z) const
                {
                        return points_[(1 + x) + 3 * (1 + y) + 9 * (1 + z)];
                }
        }; // struct NeighboringPoints

        /// explicit specialization with TwoDimension
        template<>
        inline NeighboringPoints<TwoDimension>::NeighboringPoints()
               : points_{
                        {{-1, -1}}, // 0
                        {{ 0, -1}}, // 1
                        {{ 1, -1}}, // 2
                        {{-1,  0}}, // 3
                        {{ 0,  0}}, // 4
                        {{ 1,  0}}, // 5
                        {{-1,  1}}, // 6
                        {{ 0,  1}}, // 7
                        {{ 1,  1}}} // 8
        {}

        /// explicit specialization with ThreeDimension
        template<>
        inline NeighboringPoints<ThreeDimension>::NeighboringPoints()
                : points_{
                        {{-1, -1, -1}}, // 0
                        {{ 0, -1, -1}}, // 1
                        {{ 1, -1, -1}}, // 2
                        {{-1,  0, -1}}, // 3
                        {{ 0,  0, -1}}, // 4
                        {{ 1,  0, -1}}, // 5
                        {{-1,  1, -1}}, // 6
                        {{ 0,  1, -1}}, // 7
                        {{ 1,  1, -1}}, // 8

                        {{-1, -1,  0}}, // 9
                        {{ 0, -1,  0}}, // 10
                        {{ 1, -1,  0}}, // 11
                        {{-1,  0,  0}}, // 12
                        {{ 0,  0,  0}}, // 13
                        {{ 1,  0,  0}}, // 14
                        {{-1,  1,  0}}, // 15
                        {{ 0,  1,  0}}, // 16
                        {{ 1,  1,  0}}, // 17     
                
                        {{-1, -1,  1}}, // 18
                        {{ 0, -1,  1}}, // 19
                        {{ 1, -1,  1}}, // 20
                        {{-1,  0,  1}}, // 21
                        {{ 0,  0,  1}}, // 22
                        {{ 1,  0,  1}}, // 23
                        {{-1,  1,  1}}, // 24
                        {{ 0,  1,  1}}, // 25
                        {{ 1,  1,  1}}} // 26     
        {}

        const NeighboringPoints<TwoDimension> NeighboringPoints2d;
        const NeighboringPoints<ThreeDimension> NeighboringPoints3d;
} // namespace lsm
#endif // LevelSetMethod2_NeighborPoints_h
