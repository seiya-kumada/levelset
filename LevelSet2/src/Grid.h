//
//  Grid.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/11/23.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_Grid_h
#define LevelSetMethod2_Grid_h

#include "DimensionTypes.h"

namespace lsm
{
        template<typename D>
        struct Grid;

        /// explicit specialization with TwoDimension
        template<>
        struct Grid<TwoDimension>
        {
                int left_;
                int top_;
                int right_;
                int bottom_;
        }; // struct Grid<TwoDimension>

        /// explicit specialization with ThreeDimension
        template<>
        struct Grid<ThreeDimension>
        {
                int left_;
                int top_;
                int front_;

                int right_;
                int bottom_;
                int back_;
        }; // struct Grid<ThreeDimension>
        
        typedef Grid<TwoDimension>   Grid2d;
        typedef Grid<ThreeDimension> Grid3d;
} // namespace lsm

#endif // LevelSetMethod2_Grid_h
