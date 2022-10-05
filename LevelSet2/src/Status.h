//
//  Status.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/12/16.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_Status_h
#define LevelSetMethod2_Status_h

namespace lsm
{
        enum class Status {
                Farway,
                Band,
                ResetBand,
                Front,
        }; // enum class Status
} // namespace lsm

#endif // LevelSetMethod2_Status_h
