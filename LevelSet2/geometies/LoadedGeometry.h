//
//  LoadedGeometry.h
//  LevelSetMethod2
//
//  Created by kumada on 2013/02/10.
//  Copyright (c) 2013年 kumada. All rights reserved.
//

#ifndef __LevelSetMethod2__LoadedGeometry__
#define __LevelSetMethod2__LoadedGeometry__

#include "Geometry.h"
#include "Point.h"
#include <vector>

namespace lsm
{
        class LoadedGeometry : public Geometry
        {
        public:
                LoadedGeometry(bool is_solid, const std::vector<FloatPoint3d>& object);
                virtual ~LoadedGeometry() {};

        private:
                const std::vector<FloatPoint3d>& mixed_vertices_;
                std::size_t vertex_size_;
                
                virtual void assign_value(const IntPoint3d& p, std::uint8_t& v) const override {}
                virtual void draw_geometry(const double matrix[16]) const override;
        }; // class LoadedGeometry
} // namespace lsm

#endif /* defined(__LevelSetMethod2__LoadedGeometry__) */
