//
//  Geometry.h
//  LevelSetMethod2
//
//  Created by kumada on 2013/01/24.
//  Copyright (c) 2013年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_Geometry_h
#define LevelSetMethod2_Geometry_h

#include "Point.h"

namespace lsm
{
        class Geometry
        {
        public:
                Geometry(bool is_solid)
                        : is_solid_(is_solid) {}
                
                virtual ~Geometry() = 0;
                
                void operator()(const IntPoint3d& p, std::uint8_t& v) const
                {
                        assign_value(p, v);
                }
                
                void draw(const double matrix[16]) const
                {
                        draw_geometry(matrix);
                }
                
        protected:
                bool is_solid_;
                
        private:
                Geometry(const Geometry&) = delete;
                Geometry& operator=(const Geometry&) = delete;
                Geometry(Geometry&&) = delete;
                Geometry& operator=(Geometry&&) = delete;
                
                virtual void assign_value(const IntPoint3d& p, std::uint8_t& v) const = 0;
                virtual void draw_geometry(const double matrix[16]) const = 0;
        }; // class Geometry
} // namespace lsm

#endif // LevelSetMethod2_Geometry_h
