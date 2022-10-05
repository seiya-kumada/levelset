//
//  GeometryGenerator.cpp
//  LevelSetMethod2
//
//  Created by kumada on 2013/01/19.
//  Copyright (c) 2013年 kumada. All rights reserved.
//

#include "GeometryGenerator.h"
#include <boost/range/algorithm/for_each.hpp>
#include "Geometry.h"
#include <memory>

using namespace lsm;

GeometryGenerator::GeometryGenerator(const SpaceSize3d& space_size)
        : space_size_(space_size) {}

void GeometryGenerator::add_geometry(const std::shared_ptr<Geometry>& geometry)
{
        geometries_.push_back(geometry);
}

void GeometryGenerator::generate(std::vector<std::uint8_t>& space) const
{
        const int width = space_size_.width_;
        const int height = space_size_.height_;
        const int depth = space_size_.depth_;
        const int area = width * height;
        
        for ( int k = 0, ak = area * k; k < depth; ++k, ak += area ) {
                for ( int j = 0, wj = ak + width * j; j < height; ++j, wj += width ) {
                        std::uint8_t* p = &space[wj];
                        for ( int i = 0; i < width; ++i ) {
                                for ( const auto& g : geometries_ ) {
                                        (*g)({{i, j, k}}, p[i]);
                                }
                        }
                }
        }
}