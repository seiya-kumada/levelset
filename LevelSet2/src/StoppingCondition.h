//
//  StoppingCondition.h
//  LevelSetMethod2
//
//  Created by kumada on 2013/01/18.
//  Copyright (c) 2013年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_StoppingCondition_h
#define LevelSetMethod2_StoppingCondition_h

namespace lsm
{
        constexpr double Epsilon = 1.0e-08;
        struct StoppingCondition
        {
                static constexpr int Num_ {3};
                std::array<double, Num_> total_speeds_ {{0, 0, 0}};
                std::size_t counter_ {0};
                
                
                void add_total_speed(double speed)
                {
                        total_speeds_[counter_ % Num_] = speed;
                        ++counter_;
                }
                
                bool is_satisfied() const
                {
                        return std::abs(total_speeds_[0] - total_speeds_[1]) < Epsilon
                            && std::abs(total_speeds_[1] - total_speeds_[2]) < Epsilon;
                }
                
                std::size_t get_counter() const
                {
                        return counter_;
                }
        };

        
}
#endif // LevelSetMethod2_StoppingCondition_h
