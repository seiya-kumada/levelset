//
//  GeometryLoader.cpp
//  AssimpTest
//
//  Created by kumada on 2013/01/31.
//  Copyright (c) 2013年 kumada. All rights reserved.
//

#include "GeometryLoader.h"
#include "Debug.h"
#include <boost/range/algorithm/for_each.hpp>
#include <boost/range/numeric.hpp>
#include "LoadedGeometry.h"
#include <assimp/scene.h>

using namespace lsm;

namespace
{
        inline void initialize_corner_point(FloatPoint3d& p, float initial_value)
        {
                p[0] = initial_value;
                p[1] = initial_value;
                p[2] = initial_value;
        }
        
        inline void scale_corner_point(FloatPoint3d& p, float scale)
        {
                p[0] *= scale;
                p[1] *= scale;
                p[2] *= scale;
        }
        
        inline void move_corner_point(FloatPoint3d& p, float x_distance, float y_distance, float z_distance)
        {
                p[0] += x_distance;
                p[1] += y_distance;
                p[2] += z_distance;
        }
        
        constexpr std::size_t VertexSizeOfTriangle {3};

        // test ok
        void construct_triangles(
                const std::vector<FloatPoint3d>&   vertices,
                const std::vector<FloatPoint3d>&   normals,
                std::vector<GeometryLoader::Triangle>& triangles
        ) {
                const std::size_t vertex_size = vertices.size();
                triangles.reserve(vertex_size / 3);
                for ( std::size_t i = 0; i < vertex_size; i += VertexSizeOfTriangle ) {
                        const FloatPoint3d& v0 = vertices[i + 0];
                        const FloatPoint3d& v1 = vertices[i + 1];
                        const FloatPoint3d& v2 = vertices[i + 2];
                        const GeometryLoader::Triangle triangle { {&v0, &v1, &v2}, &normals[i] };
                        triangles.push_back(std::move(triangle));
                }
        }
        
        constexpr int Offset {2};
        
        // test ok
        void calculate_bounding_box(
                const GeometryLoader::Triangle& triangle,
                Grid3d&                         bounding_box
        ) {
                const FloatPoint3d& v0 = *triangle.vertices_[0];
                const FloatPoint3d& v1 = *triangle.vertices_[1];
                const FloatPoint3d& v2 = *triangle.vertices_[2];
                
                bounding_box.left_   = static_cast<int>(std::min(std::min(v0[0], v1[0]), v2[0])) - Offset;
                bounding_box.right_  = static_cast<int>(std::max(std::max(v0[0], v1[0]), v2[0])) + Offset;
                bounding_box.top_    = static_cast<int>(std::min(std::min(v0[1], v1[1]), v2[1])) - Offset;
                bounding_box.bottom_ = static_cast<int>(std::max(std::max(v0[1], v1[1]), v2[1])) + Offset;
                bounding_box.front_  = static_cast<int>(std::min(std::min(v0[2], v1[2]), v2[2])) - Offset;
                bounding_box.back_   = static_cast<int>(std::max(std::max(v0[2], v1[2]), v2[2])) + Offset;
        }
        
        constexpr float Epsilon = 1.0f;
        
        // test ok
        inline bool is_on_plane(
                const FloatPoint3d&             p,
                const GeometryLoader::Triangle& triangle
        ) {
                const FloatPoint3d& v0 = *triangle.vertices_[0];
                return std::abs(inner_product(cast<float>(p) - v0, *triangle.normal_)) <= Epsilon;
        }
        
        // test ok
        inline bool is_on_triangle(
                const FloatPoint3d&             p,
                const GeometryLoader::Triangle& triangle
        ) {
                const FloatPoint3d& a = *triangle.vertices_[0];
                const FloatPoint3d& b = *triangle.vertices_[1];
                const FloatPoint3d& c = *triangle.vertices_[2];
                const FloatPoint3d& n = *triangle.normal_;
                
                return inner_product(outer_product(b - a, p - a), n) >= 0 &&
                       inner_product(outer_product(c - b, p - b), n) >= 0 &&
                       inner_product(outer_product(a - c, p - c), n) >= 0;
        }
        
        void search_points_on_triangle(
                const GeometryLoader::Triangle& triangle,
                const SpaceSize3d&              space_size,
                std::vector<std::uint8_t>&      data,
                std::vector<FloatPoint3d>&      mixed_vertices
        ) {
                Grid3d bounding_box;
                calculate_bounding_box(triangle, bounding_box);
                
                const int left = bounding_box.left_;
                const int right = bounding_box.right_;
                const int top = bounding_box.top_;
                const int bottom = bounding_box.bottom_;
                const int front = bounding_box.front_;
                const int back = bounding_box.back_;
                
                const int width  = space_size.width_;
                const int height = space_size.height_;
                const int area   = width * height;
                
                for ( int k = front, ak = area * k; k < back; ++k, ak += area ) {
                        for ( int j = top, wj = ak + width * j; j < bottom; ++j, wj += width ) {
                                std::uint8_t* p = &data[wj];
                                for ( int i = left; i < right; ++i ) {
                                        const FloatPoint3d v {{static_cast<float>(i), static_cast<float>(j), static_cast<float>(k)}};
                                        if ( is_on_plane(v, triangle) ) {
                                                if ( is_on_triangle(v, triangle) ) {
                                                        p[i] = 128;
                                                        mixed_vertices.push_back(*triangle.normal_);
                                                        mixed_vertices.push_back(std::move(v));
                                                }
                                        }
                                }
                        }
                }
        }
}

void GeometryLoader::read_file(const std::string& filename)
{
        const aiScene* scene = importer_.ReadFile(filename.c_str(), 0);
        if ( !scene ) {
                throw std::runtime_error(">> unable to read file");
        }
        extract_scene(scene);
}

void GeometryLoader::convert(const SpaceSize3d& space_size)
{
        data_.clear();
        data_.resize(space_size.total_, 255);
        convert(space_size, data_);
}

void GeometryLoader::convert(const SpaceSize3d& space_size, std::vector<std::uint8_t>& data)
{
        std::vector<Triangle> triangles;
        construct_triangles(vertices_, normals_, triangles);
        
        mixed_vertices_.clear();
        for ( const auto& triangle : triangles ) {
                search_points_on_triangle(triangle, space_size, data, mixed_vertices_);
        }
        
        geometries_.clear();
        geometries_.push_back(std::make_shared<LoadedGeometry>(false, mixed_vertices_));
}

// test ok
void GeometryLoader::extract_scene(const aiScene* scene)
{
        if ( scene->mNumMeshes != 1 ) {
                throw std::runtime_error("scene->mNumMeshes must be 1");
        }
                
        const aiMesh* mesh = scene->mMeshes[0];
        const std::size_t vertex_size = mesh->mNumVertices;
        
        vertices_.clear();
        vertices_.reserve(vertex_size / VertexSizeOfTriangle);
        
        normals_.clear();
        normals_.reserve(vertex_size / VertexSizeOfTriangle);
        
        initialize_corner_point(left_top_front_corner_, std::numeric_limits<float>::max());
        initialize_corner_point(right_bottom_back_corner_, std::numeric_limits<float>::min());
        
        for ( std::size_t i = 0 ; i < vertex_size; i += VertexSizeOfTriangle ) {
                const aiVector3D& p0 = mesh->mVertices[i + 0];
                const aiVector3D& p1 = mesh->mVertices[i + 1];
                const aiVector3D& p2 = mesh->mVertices[i + 2];
                
                decide_range(p0);
                decide_range(p1);
                decide_range(p2);
                
                vertices_.push_back({{p0.x, p0.y, p0.z}});
                vertices_.push_back({{p1.x, p1.y, p1.z}});
                vertices_.push_back({{p2.x, p2.y, p2.z}});
                
                const aiVector3D& n0 = mesh->mNormals[i + 0];
                const aiVector3D& n1 = mesh->mNormals[i + 1];
                const aiVector3D& n2 = mesh->mNormals[i + 2];
                
                normals_.push_back({{n0.x, n0.y, n0.z}});
                normals_.push_back({{n1.x, n1.y, n1.z}});
                normals_.push_back({{n2.x, n2.y, n2.z}});
        }
}

// test ok
void GeometryLoader::scale_scene(float standard_length)
{
        const float width  {get_length(Label1::Right,  Label0::Left)};
        const float height {get_length(Label1::Bottom, Label0::Top)};
        const float depth  {get_length(Label1::Back,   Label0::Front)};
        
        const float max_length {std::max(width, std::max(height, depth))};
        const float scale      {standard_length / max_length};
        
        boost::for_each(vertices_,
                [&](FloatPoint3d& p)
                {
                        p[0] *= scale;
                        p[1] *= scale;
                        p[2] *= scale;
                }
        );
        scale_corner_point(left_top_front_corner_, scale);
        scale_corner_point(right_bottom_back_corner_, scale);
}

// test ok
void GeometryLoader::move_to_center(const FloatPoint3d& center)
{
        const FloatPoint3d centroid = calculate_centroid();
        const float x_distance {center[0] - centroid[0]};
        const float y_distance {center[1] - centroid[1]};
        const float z_distance {center[2] - centroid[2]};
        
        for ( std::size_t i = 0, n = vertices_.size(); i < n; ++i ) {
                FloatPoint3d& v = vertices_[i];
                v[0] += x_distance;
                v[1] += y_distance;
                v[2] += z_distance;
        }
        move_corner_point(left_top_front_corner_, x_distance, y_distance, z_distance);
        move_corner_point(right_bottom_back_corner_, x_distance, y_distance, z_distance);
}

// test ok
const FloatPoint3d GeometryLoader::calculate_centroid() const
{
        FloatPoint3d c {{0.0f, 0.0f, 0.0f}};
        for ( const auto& p : vertices_ ) {
                c += p;
        }
        const std::size_t n {vertices_.size()};
        return {{c[0] / n, c[1] / n, c[2] / n}};
}

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
#if(UNIT_TEST_GeometryLoader)
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/random.hpp>
#include <boost/range/algorithm/generate.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <assimp/Importer.hpp> // C++ importer interface
#include <assimp/scene.h> // Output data structure
#include <assimp/postprocess.h> // Post processing flags

namespace
{
        template<typename T>
        std::ostream& operator<<(std::ostream& os, const Point3d<T>& p)
        {
                return os << p[0] << ", " << p[1] << ", " << p[2];
        }
                
        const float a[] = {50, 0, 0};
        const float b[] = {0, 50, 0};
        const float c[] = {0, 0, 50};
        const float d[] = {0, 0, 0};
        
        const FloatPoint3d av {{50,  0,  0}};
        const FloatPoint3d bv {{ 0, 50,  0}};
        const FloatPoint3d cv {{ 0,  0, 50}};
        const FloatPoint3d dv {{ 0,  0,  0}};
        
        const FloatPoint3d an = {{0, 0, 0}};
        const FloatPoint3d bn = {{1, 0, 0}};
        const FloatPoint3d cn = {{0, 1, 0}};
        const FloatPoint3d dn = {{0, 0, 1}};
        
        
        void set_triangle(std::vector<FloatPoint3d>& vertices, const FloatPoint3d& a, const FloatPoint3d& b, const FloatPoint3d& c)
        {
                vertices.push_back(a);
                vertices.push_back(b);
                vertices.push_back(c);
        }
        
        void set_vertices(std::vector<FloatPoint3d>& vertices)
        {
                set_triangle(vertices, av, bv, cv);
                set_triangle(vertices, dv, cv, bv);
                set_triangle(vertices, cv, dv, av);
                set_triangle(vertices, av, dv, bv);
        }
        
        void set_normal(const FloatPoint3d& normal, std::vector<FloatPoint3d>& normals)
        {
                for ( int i = 0; i < 3; ++i ) {
                        normals.push_back(normal);
                }
        }
        
        void set_normals(std::vector<FloatPoint3d>& normals)
        {
                set_normal(an, normals);
                set_normal(bn, normals);
                set_normal(cn, normals);
                set_normal(dn, normals);
        }
        
        void test_construct_triangles()
        {
                std::vector<FloatPoint3d> vertices;
                set_vertices(vertices);
                
                std::vector<FloatPoint3d> normals;
                set_normals(normals);
                
                std::vector<GeometryLoader::Triangle> triangles;
                construct_triangles(vertices, normals, triangles);
                
                BOOST_CHECK_EQUAL(triangles.size(), 4);

                BOOST_CHECK(triangles[0].vertices_[0] == &vertices[0]);
                BOOST_CHECK(triangles[0].vertices_[1] == &vertices[1]);
                BOOST_CHECK(triangles[0].vertices_[2] == &vertices[2]);
                BOOST_CHECK(triangles[0].normal_ == &normals[0]);

                BOOST_CHECK(triangles[1].vertices_[0] == &vertices[3]);
                BOOST_CHECK(triangles[1].vertices_[1] == &vertices[4]);
                BOOST_CHECK(triangles[1].vertices_[2] == &vertices[5]);
                BOOST_CHECK(triangles[1].normal_ == &normals[3]);

                BOOST_CHECK(triangles[2].vertices_[0] == &vertices[6]);
                BOOST_CHECK(triangles[2].vertices_[1] == &vertices[7]);
                BOOST_CHECK(triangles[2].vertices_[2] == &vertices[8]);
                BOOST_CHECK(triangles[2].normal_ == &normals[6]);

                BOOST_CHECK(triangles[3].vertices_[0] == &vertices[9]);
                BOOST_CHECK(triangles[3].vertices_[1] == &vertices[10]);
                BOOST_CHECK(triangles[3].vertices_[2] == &vertices[11]);
                BOOST_CHECK(triangles[3].normal_ == &normals[9]);
        }
        
        void test_extract_scene()
        {
                Assimp::Importer importer;
                const char* filename = "/Users/kumada/Projects/level-set-method/AssimpTest/AssimpTest/src/test.stl";
                const aiScene* scene = importer.ReadFile(filename, 0);
                if ( scene == nullptr ) {
                        std::cout << "unable to load an input file\n";
                        return;
                }
                GeometryLoader loader;
                loader.extract_scene(scene);
                const std::vector<FloatPoint3d>& vertices = loader.get_vertices();
                const std::vector<FloatPoint3d>& normals = loader.get_normals();
                BOOST_CHECK_EQUAL(normals.size(), 6);
                
                BOOST_CHECK(vertices[0] == FloatPoint3d({{1, 0, 0}}));
                BOOST_CHECK(vertices[1] == FloatPoint3d({{2, 0, 0}}));
                BOOST_CHECK(vertices[2] == FloatPoint3d({{3, 0, 0}}));
                
                BOOST_CHECK(vertices[3] == FloatPoint3d({{1, 2, 0}}));
                BOOST_CHECK(vertices[4] == FloatPoint3d({{2, 1, 0}}));
                BOOST_CHECK(vertices[5] == FloatPoint3d({{3, 0, 2}}));
                
                BOOST_CHECK(normals[0] == FloatPoint3d({{1, 1, 0}}));
                BOOST_CHECK(normals[1] == FloatPoint3d({{1, 1, 0}}));
                BOOST_CHECK(normals[2] == FloatPoint3d({{1, 1, 0}}));
                
                BOOST_CHECK(normals[3] == FloatPoint3d({{1, 2, 0}}));
                BOOST_CHECK(normals[4] == FloatPoint3d({{1, 2, 0}}));
                BOOST_CHECK(normals[5] == FloatPoint3d({{1, 2, 0}}));
                
                const FloatPoint3d& ltf = loader.get_left_top_front_corner();
                const FloatPoint3d& rbb = loader.get_right_bottom_back_corner();
                
                BOOST_CHECK(ltf == FloatPoint3d({{1, 0, 0}}));
                BOOST_CHECK(rbb == FloatPoint3d({{3, 2, 2}}));
                
        }

        void test_scale_scene()
        {
                Assimp::Importer importer;
                const char* filename = "/Users/kumada/Projects/level-set-method/AssimpTest/AssimpTest/src/test.stl";
                const aiScene* scene = importer.ReadFile(filename, 0);
                if ( scene == nullptr ) {
                        std::cout << "unable to load an input file\n";
                        return;
                }
                GeometryLoader loader;
                loader.extract_scene(scene);
                loader.scale_scene(100);
                
                const std::vector<FloatPoint3d>& vertices = loader.get_vertices();
                
                BOOST_CHECK(vertices[0] == FloatPoint3d({{50, 0, 0}}));;
                BOOST_CHECK(vertices[1] == FloatPoint3d({{100, 0, 0}}));
                BOOST_CHECK(vertices[2] == FloatPoint3d({{150, 0, 0}}));
                BOOST_CHECK(vertices[3] == FloatPoint3d({{50, 100, 0}}));
                BOOST_CHECK(vertices[4] == FloatPoint3d({{100, 50, 0}}));
                BOOST_CHECK(vertices[5] == FloatPoint3d({{150, 0, 100}}));

                const FloatPoint3d& ltf = loader.get_left_top_front_corner();
                const FloatPoint3d& rbb = loader.get_right_bottom_back_corner();
                
                BOOST_CHECK(ltf == FloatPoint3d({{50, 0, 0}}));
                BOOST_CHECK(rbb == FloatPoint3d({{150, 100, 100}}));
        }

        void test_move_scene()
        {
                Assimp::Importer importer;
                const char* filename = "/Users/kumada/Projects/level-set-method/AssimpTest/AssimpTest/src/test.stl";
                const aiScene* scene = importer.ReadFile(filename, 0);
                if ( scene == nullptr ) {
                        std::cout << "unable to load an input file\n";
                        return;
                }
                GeometryLoader loader;
                loader.extract_scene(scene);
                loader.scale_scene(100);
                FloatPoint3d center {{0, 0, 0}};
                loader.move_to_center(center);
                const std::vector<FloatPoint3d>& vertices = loader.get_vertices();

                FloatPoint3d c {{0.0f, 0.0f, 0.0f}};
                for ( const auto& p : vertices ) {
                        c += p;
                }
                
                constexpr float epsilon = 1.0e-05;
                BOOST_CHECK_CLOSE(c[0], 0.0f, epsilon);
                BOOST_CHECK_CLOSE(c[1], 0.0f, epsilon);
                BOOST_CHECK(std::abs(c[2] - 0.0f) < epsilon);
                
                const FloatPoint3d& ltf = loader.get_left_top_front_corner();
                const FloatPoint3d& rbb = loader.get_right_bottom_back_corner();
                
                BOOST_CHECK(ltf == FloatPoint3d({{-50, -25, -50.0f/3}}));
                BOOST_CHECK(rbb == FloatPoint3d({{50, 75, 250.0f/3}}));
        }

        void test_is_on_plane()
        {
                const FloatPoint3d normal {{1 / std::sqrt(2.0f), 1 / std::sqrt(2.0f), 1 / std::sqrt(2.0f)}};
                GeometryLoader::Triangle triangle;
                triangle.vertices_[0] = &av;
                triangle.vertices_[1] = &bv;
                triangle.vertices_[2] = &cv;
                triangle.normal_ = &normal;
                
                const FloatPoint3d p0 {{0, 0, 0}};
                BOOST_CHECK(!is_on_plane(p0, triangle));

                const FloatPoint3d p1 {{1, 1, 48}};
                BOOST_CHECK(is_on_plane(p1, triangle));
                
                const FloatPoint3d p2 {{2, 3, 48}};
                BOOST_CHECK(!is_on_plane(p2, triangle));
                
                const FloatPoint3d p3 {{10, 1, 39}};
                BOOST_CHECK(is_on_plane(p3, triangle));
        }
        
        void test_is_on_triangle()
        {
                const FloatPoint3d normal {{1 / std::sqrt(2.0f), 1 / std::sqrt(2.0f), 1 / std::sqrt(2.0f)}};
                GeometryLoader::Triangle triangle;
                triangle.vertices_[0] = &av;
                triangle.vertices_[1] = &bv;
                triangle.vertices_[2] = &cv;
                triangle.normal_ = &normal;
                
                const FloatPoint3d p1 {{1, 1, 48}};
                BOOST_CHECK(is_on_triangle(p1, triangle));
                
                const FloatPoint3d p3 {{100, 1, 39}};
                BOOST_CHECK(!is_on_triangle(p3, triangle));
                
                
        }
        
        void test_calculate_bounding_box()
        {
                std::vector<FloatPoint3d> vertices;
                set_vertices(vertices);
                
                std::vector<FloatPoint3d> normals;
                set_normals(normals);
                
                std::vector<GeometryLoader::Triangle> triangles;
                construct_triangles(vertices, normals, triangles);

                Grid3d bounding_box;
                calculate_bounding_box(triangles[0], bounding_box);
                
                BOOST_CHECK(bounding_box.left_ == -Offset);
                BOOST_CHECK(bounding_box.right_ == 50 + Offset);
                BOOST_CHECK(bounding_box.top_ == -Offset);
                BOOST_CHECK(bounding_box.bottom_ == 50 + Offset);
                BOOST_CHECK(bounding_box.front_ == -Offset);
                BOOST_CHECK(bounding_box.back_ == 50 + Offset);
                
        }
        
        void test_point()
        {
                const FloatPoint3d a {{1, 2, 3}};
                const FloatPoint3d b = 3.0f * a;
                BOOST_CHECK(FloatPoint3d({{3, 6, 9}}) == b);
        }        
}

BOOST_AUTO_TEST_CASE(TEST_GeometryLoader)
{
        test_construct_triangles();
        test_extract_scene();
        test_scale_scene();
        test_move_scene();
        test_is_on_plane();
        test_is_on_triangle();
        test_calculate_bounding_box();
        test_point();
        std::cout << "GeometryLoader\n";
}

#endif // UNIT_TEST_GeometryLoader