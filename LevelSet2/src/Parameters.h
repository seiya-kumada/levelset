//
//  Parameters.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/11/23.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_Parameters_h
#define LevelSetMethod2_Parameters_h

namespace lsm
{
        struct Parameters
        {
                /// width of a narrow band
                int wband_;
                
                /// width to be reset
                int wreset_;
                
                /// time step
                double time_step_;
                
                /// gain of curvature
                double gain_;
                
                /// constant speed
                double constant_speed_;
                
                double speed_threshold_;
                
        }; // Parameters
} // namespace lsm
#endif // LevelSetMethod2_Parameters_h
