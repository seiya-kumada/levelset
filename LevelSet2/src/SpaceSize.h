//
//  SpaceSize.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/11/22.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_Size_h
#define LevelSetMethod2_Size_h

#include "DimensionTypes.h"

namespace lsm
{
        template<typename D>
        struct SpaceSize;

        template<>
        struct SpaceSize<TwoDimension>
        {
                int width_;
                int height_;
                int total_;
                
                SpaceSize() = default;
                
                SpaceSize(int w, int h)
                        : width_{w}
                        , height_{h}
                        , total_{w * h} {}
        }; // struct SpaceSize<TwoDimension>

        template<>
        struct SpaceSize<ThreeDimension>
        {
                int width_;
                int height_;
                int depth_;
                int total_;
                
                SpaceSize() = default;
                
                SpaceSize(int w, int h, int d)
                        : width_{w}
                        , height_{h}
                        , depth_{d}
                        , total_{w * h * d} {}
        }; // struct SpaceSize<ThreeDimension>
        
        typedef SpaceSize<TwoDimension>   SpaceSize2d;
        typedef SpaceSize<ThreeDimension> SpaceSize3d;
        
} // namespace lsm

#endif // LevelSetMethod2_Size_h
