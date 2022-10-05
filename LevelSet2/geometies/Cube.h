//
//  Cube.h
//  LevelSetMethod2
//
//  Created by kumada on 2013/01/25.
//  Copyright (c) 2013年 kumada. All rights reserved.
//

#ifndef __LevelSetMethod2__Cube__
#define __LevelSetMethod2__Cube__

#include "Geometry.h"

namespace lsm
{
        class Cube : public Geometry
        {
        public:
                Cube(bool is_solid, const IntPoint3d& center, int size);
                virtual ~Cube() {};

        private:
                IntPoint3d center_;
                int        size_;
                IntPoint3d left_top_front_;
                IntPoint3d right_bottom_back_;
        
                bool is_inside(const IntPoint3d& p) const;
                virtual void assign_value(const IntPoint3d& p, std::uint8_t& v) const override;
                virtual void draw_geometry(const double matrix[16]) const override;
        }; // class Cube
} // namespace lsm


#endif /* defined(__LevelSetMethod2__Cube__) */
