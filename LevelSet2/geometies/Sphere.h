//
//  Sphere.h
//  LevelSetMethod2
//
//  Created by kumada on 2013/01/24.
//  Copyright (c) 2013年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_Sphere_h
#define LevelSetMethod2_Sphere_h

#include "Geometry.h"

namespace lsm
{
        class Sphere : public Geometry
        {
        public:
                Sphere(bool is_solid, const IntPoint3d& center, int raidus);
                virtual ~Sphere() {};
                
        private:
                IntPoint3d center_;
                int        radius_;
                
                virtual void assign_value(const IntPoint3d& p, std::uint8_t& v) const override;
                virtual void draw_geometry(const double matrix[16]) const override;
        }; // class Sphere
} // namespace lsm

#endif // LevelSetMethod2_Sphere_h
