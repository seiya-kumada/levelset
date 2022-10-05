//
//  InsideEstimator.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/11/21.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_InsideEstimator_h
#define LevelSetMethod2_InsideEstimator_h

#include "Concepts.h"
#include "Grid.h"
#include "DimensionTypes.h"
#include <boost/concept_check.hpp>

namespace lsm
{
        /// judges whether a given point is inside a given grid or not
        /**
         *      @tparam D the dimension type
         */
        template<typename D>
        class InsideEstimator
        {
                BOOST_CLASS_REQUIRE(D, lsm, DimensionConcept);

        public:
                explicit InsideEstimator(
                        const Grid<D>& grid
                )
                        : grid_(grid) {}
                
                InsideEstimator() = default;
                ~InsideEstimator() = default;

                void set_grid(const Grid<D>& grid)
                {
                        grid_ = grid;
                }
                
                const Grid<D>& get_grid() const
                {
                        return grid_;
                }
                
                bool operator()(const IntPoint<D>& p) const;

        private:
                InsideEstimator(const InsideEstimator&) = delete;
                InsideEstimator& operator=(const InsideEstimator&) = delete;
                
                InsideEstimator(InsideEstimator&&) = delete;
                InsideEstimator& operator=(InsideEstimator&&) = delete;

                Grid<D> grid_;
        }; // class InsideEstimator

        /// explicit specialization with TwoDimension
        template<>
        inline bool InsideEstimator<TwoDimension>::operator()(const IntPoint2d& p) const
        {
                return (grid_.left_ < p[0]) && (p[0] < grid_.right_)
                     && (grid_.top_ < p[1]) && (p[1] < grid_.bottom_);
        }

        /// explicit specialization with ThreeDimension
        template<>
        inline bool InsideEstimator<ThreeDimension>::operator()(const IntPoint3d& p) const
        {
                return (grid_.left_ < p[0]) && (p[0] < grid_.right_)
                     && (grid_.top_ < p[1]) && (p[1] < grid_.bottom_)
                   && (grid_.front_ < p[2]) && (p[2] < grid_.back_);
        }
} // namespace lsm

#endif // LevelSetMethod2_InsideEstimator_h
