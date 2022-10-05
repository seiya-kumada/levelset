//
//  GeometryLoader.h
//  AssimpTest
//
//  Created by kumada on 2013/01/31.
//  Copyright (c) 2013年 kumada. All rights reserved.
//

#ifndef __AssimpTest__GeometryLoader__
#define __AssimpTest__GeometryLoader__

#include <string>
#include <vector>
#include "Point.h"
#include "Grid.h"
#include "Geometry.h"
#include "SpaceSize.h"
#include <memory>
#include <assimp/Importer.hpp> // C++ importer interface

namespace lsm
{
        /// loads 3D model by means of assimp library.
        class GeometryLoader
        {
        public:
                struct Triangle
                {
                        const FloatPoint3d* vertices_[3];
                        const FloatPoint3d* normal_;
                };
        
                GeometryLoader() = default;
                ~GeometryLoader() = default;
            
                void read_file(const std::string& filename);
                void extract_scene(const aiScene* scene);
                void scale_scene(float max_length);
                void move_to_center(const FloatPoint3d& center);
                
                const std::vector<FloatPoint3d>& get_vertices() const
                {
                        return vertices_;
                }
                
                const std::vector<FloatPoint3d>& get_normals() const
                {
                        return normals_;
                }
                
                const FloatPoint3d& get_left_top_front_corner() const
                {
                        return left_top_front_corner_;
                }
                
                const FloatPoint3d& get_right_bottom_back_corner() const
                {
                        return right_bottom_back_corner_;
                }
                
                /// converts to an input data passed to an instance of the class, LevelSetMethod<D>.
                void convert(const SpaceSize3d& space_size);
                void convert(const SpaceSize3d& space_size, std::vector<std::uint8_t>& data);
                
                /// get an input data passed to an instance of the class, LevelSetMethod<D>
                const std::vector<std::uint8_t>& get_data() const
                {
                        return data_;
                }
                
                /// used for drawing the object 
                const std::vector<std::shared_ptr<Geometry>>& get_geometries() const
                {
                        return geometries_;
                }

        private:
                enum Label0
                {
                        Left = 0,
                        Top,
                        Front,
                };
                
                enum Label1
                {
                        Right = 0,
                        Bottom,
                        Back,
                };
        
                /// assimp library module
                Assimp::Importer          importer_;
                
                /// extracted vertices and normals
                std::vector<FloatPoint3d> vertices_;
                std::vector<FloatPoint3d> normals_;
                
                /// bounding box of an object
                FloatPoint3d              left_top_front_corner_;
                FloatPoint3d              right_bottom_back_corner_;
                
                /// an input data passed to an instance of the class, LevelSetMethod<D>
                std::vector<std::uint8_t> data_;
                
                /// data used for drawing an object
                std::vector<FloatPoint3d> mixed_vertices_;
                
                std::vector<std::shared_ptr<Geometry>> geometries_;
                
                GeometryLoader(const GeometryLoader&) = delete;
                GeometryLoader& operator=(const GeometryLoader&) = delete;
                GeometryLoader(GeometryLoader&&) = delete;
                GeometryLoader& operator=(GeometryLoader&&) = delete;
                
                void decide_range(int label0, int label1, float value)
                {
                        if ( left_top_front_corner_[label0] > value ) {
                                left_top_front_corner_[label0] = value;
                        }
                        if ( right_bottom_back_corner_[label1] < value ) {
                                right_bottom_back_corner_[label1] = value;
                        }
                }
                
                void decide_range(const aiVector3D& p)
                {
                        decide_range(Label0::Left,  Label1::Right,  p.x);
                        decide_range(Label0::Top,   Label1::Bottom, p.y);
                        decide_range(Label0::Front, Label1::Back,   p.z);
                }
                
                float get_length(int label0, int label1) const
                {
                        return right_bottom_back_corner_[label0] - left_top_front_corner_[label1];
                }
                
                const FloatPoint3d calculate_centroid() const;
                
        }; // class GeometryLoader
} // namespace lsm
#endif /* defined(__AssimpTest__GeometryLoader__) */
