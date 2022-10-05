//
//  GeometryGenerator.h
//  LevelSetMethod2
//
//  Created by kumada on 2013/01/19.
//  Copyright (c) 2013年 kumada. All rights reserved.
//

#ifndef __LevelSetMethod2__GeometryGenerator__
#define __LevelSetMethod2__GeometryGenerator__

#include <vector>
#include "Point.h"
#include "SpaceSize.h"
#include <memory>

namespace lsm
{
        class Geometry;
        class GeometryGenerator
        {
        public:
                GeometryGenerator(const SpaceSize3d& space_size);
                ~GeometryGenerator() = default;
                
                void add_geometry(const std::shared_ptr<Geometry>& geometry);
                void generate(std::vector<std::uint8_t>& space) const;
                
                const std::vector<std::shared_ptr<Geometry>>& get_geometries() const
                {
                        return geometries_;
                }
                
                const SpaceSize3d& get_space_size() const
                {
                        return space_size_;
                }
                
        private:
                GeometryGenerator(const GeometryGenerator&) = delete;
                GeometryGenerator& operator=(const GeometryGenerator&) = delete;
                
                GeometryGenerator(GeometryGenerator&&) = delete;
                GeometryGenerator& operator=(GeometryGenerator&&) = delete;
                
                const SpaceSize3d&                     space_size_;
                std::vector<std::shared_ptr<Geometry>> geometries_;
        }; // GeometryGenerator
} // namespace lsm

#endif /* defined(__LevelSetMethod2__GeometryGenerator__) */
