//
//  Front.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/11/22.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_Front_h
#define LevelSetMethod2_Front_h

#include "Point.h"
#include "DimensionTypes.h"
#include <vector>

namespace lsm
{
        template<typename D>
        using Front = std::vector<IntPoint<D>>;

        using Front2d = Front<TwoDimension>;
        using Front3d = Front<ThreeDimension>;
} // namespace lsm
#endif // LevelSetMethod2_Front_h
