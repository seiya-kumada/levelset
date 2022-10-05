//
//  Indexer.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/11/22.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_Indexer_h
#define LevelSetMethod2_Indexer_h

#include "SpaceSize.h"
#include "Point.h"

namespace lsm
{
        /// converts a give point to a correspoinding index
        template<typename D>
        struct Indexer;

        /// explicit specialization with TwoDimension
        template<>
        struct Indexer<TwoDimension>
        {
                int width_;
                explicit Indexer(const SpaceSize<TwoDimension>& size)
                        : width_{size.width_} {}
                
                int operator()(const IntPoint<TwoDimension>& p) const
                {
                        return p[0] + width_ * p[1];
                }
        }; // struct Indexer<TwoDimension>

        /// explicit specialization with ThreeDimension
        template<>
        struct Indexer<ThreeDimension>
        {
                int width_;
                int area_;
                
                explicit Indexer(const SpaceSize<ThreeDimension>& size)
                        : width_{size.width_}
                        , area_{width_ * size.height_} {}
                
                int operator()(const IntPoint<ThreeDimension>& p) const
                {
                        return p[0] + width_ * p[1] + area_ * p[2];
                }
        }; // struct Indexer<ThreeDimension>
        
        typedef Indexer<TwoDimension>   Indexer2d;
        typedef Indexer<ThreeDimension> Indexer3d;
        
} // namespace lsm
#endif
