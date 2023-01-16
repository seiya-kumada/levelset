//
//  LevelSetMethod.cpp
//  LevelSetMethod2
//
//  Created by kumada on 2012/11/23.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#if(UNIT_TEST_LevelSetMethod)
#define BOOST_TEST_DYN_LINK
#include "LevelSetMethod.h"
#include <iomanip>
#include <boost/test/unit_test.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <chrono>

using namespace lsm;

namespace
{
        //_/_/_/ test for initialize_distance_map _/_/_/
         
        template<typename D>
        struct CheckerSelector;
        
        // specilization with TwoDimension
        template<>
        struct CheckerSelector<TwoDimension>
        {
                static const int Count_ = 37;
                static const std::vector<int> Sizes_;
                static const std::vector<std::vector<IntPoint<TwoDimension>>> Points_;
        };
        const std::vector<int> CheckerSelector<TwoDimension>::Sizes_ = {1, 4, 4, 4, 8, 4, 4, 8};
        const std::vector<std::vector<IntPoint<TwoDimension>>> CheckerSelector<TwoDimension>::Points_ = {
                {{{0, 0}}},
                {{{-1, 0}}, {{0, -1}}, {{0, 1}}, {{1, 0}}},
                {{{-1, -1}}, {{-1, 1}}, {{1, -1}}, {{1, 1}}},
                {{{-2, 0}}, {{0, -2}}, {{0, 2}}, {{2, 0}}},
                {{{-2, -1}}, {{-2, 1}}, {{-1, -2}}, {{-1, 2}}, {{1, -2}}, {{1, 2}}, {{2, -1}}, {{2, 1}}},
                {{{-2, -2}}, {{-2, 2}}, {{2, -2}}, {{2, 2}}},
                {{{-3, 0}}, {{0, -3}}, {{0, 3}}, {{3, 0}}},
                {{{-3, -1}}, {{-3, 1}}, {{-1, -3}}, {{-1, 3}}, {{1, -3}}, {{1, 3}}, {{3, -1}}, {{3, 1}}},
        };

      // specilization with ThreeDimension
        template<>
        struct CheckerSelector<ThreeDimension>
        {
                static const int Count_ = 179;
                static const std::vector<int> Sizes_;
                static const std::vector<std::vector<IntPoint<ThreeDimension>>> Points_;
        };
        const std::vector<int> CheckerSelector<ThreeDimension>::Sizes_ = {1, 6, 12, 8, 6, 24, 24, 12, 30, 24, 24, 8};
        const std::vector<std::vector<IntPoint<ThreeDimension>>> CheckerSelector<ThreeDimension>::Points_ = {
                {{{0, 0}}},
                {{{-1, 0, 0}}, {{0, -1, 0}}, {{0, 0, -1}}, {{0, 0, 1}}, {{0, 1, 0}}, {{1, 0, 0}}},
                {{{-1, -1, 0}}, {{-1, 0, -1}}, {{-1, 0, 1}}, {{-1, 1, 0}}, {{0, -1, -1}}, {{0, -1, 1}}, {{0, 1, -1}}, {{0, 1, 1}}, {{1, -1, 0}}, {{1, 0, -1}}, {{1, 0, 1}}, {{1, 1, 0}}},
                {{{-1, -1, -1}}, {{-1, -1, 1}}, {{-1, 1, -1}}, {{-1, 1, 1}}, {{1, -1, -1}}, {{1, -1, 1}}, {{1, 1, -1}}, {{1, 1, 1}}},
                {{{-2, 0, 0}}, {{0, -2, 0}}, {{0, 0, -2}}, {{0, 0, 2}}, {{0, 2, 0}}, {{2, 0, 0}}},
                {{{-2, -1, 0}}, {{-2, 0, -1}}, {{-2, 0, 1}}, {{-2, 1, 0}}, {{-1, -2, 0}}, {{-1, 0, -2}}, {{-1, 0, 2}}, {{-1, 2, 0}}, {{0, -2, -1}}, {{0, -2, 1}}, {{0, -1, -2}}, {{0, -1, 2}}, {{0, 1, -2}}, {{0, 1, 2}}, {{0, 2, -1}}, {{0, 2, 1}}, {{1, -2, 0}}, {{1, 0, -2}}, {{1, 0, 2}}, {{1, 2, 0}}, {{2, -1, 0}}, {{2, 0, -1}}, {{2, 0, 1}}, {{2, 1, 0}}},
                {{{-2, -1, -1}}, {{-2, -1, 1}}, {{-2, 1, -1}}, {{-2, 1, 1}}, {{-1, -2, -1}}, {{-1, -2, 1}}, {{-1, -1, -2}}, {{-1, -1, 2}}, {{-1, 1, -2}}, {{-1, 1, 2}}, {{-1, 2, -1}}, {{-1, 2, 1}}, {{1, -2, -1}}, {{1, -2, 1}}, {{1, -1, -2}}, {{1, -1, 2}}, {{1, 1, -2}}, {{1, 1, 2}}, {{1, 2, -1}}, {{1, 2, 1}}, {{2, -1, -1}}, {{2, -1, 1}}, {{2, 1, -1}}, {{2, 1, 1}}},
                {{{-2, -2, 0}}, {{-2, 0, -2}}, {{-2, 0, 2}}, {{-2, 2, 0}}, {{0, -2, -2}}, {{0, -2, 2}}, {{0, 2, -2}}, {{0, 2, 2}}, {{2, -2, 0}}, {{2, 0, -2}}, {{2, 0, 2}}, {{2, 2, 0}}},
                {{{-3, 0, 0}}, {{-2, -2, -1}}, {{-2, -2, 1}}, {{-2, -1, -2}}, {{-2, -1, 2}}, {{-2, 1, -2}}, {{-2, 1, 2}}, {{-2, 2, -1}}, {{-2, 2, 1}}, {{-1, -2, -2}}, {{-1, -2, 2}}, {{-1, 2, -2}}, {{-1, 2, 2}}, {{0, -3, 0}}, {{0, 0, -3}}, {{0, 0, 3}}, {{0, 3, 0}}, {{1, -2, -2}}, {{1, -2, 2}}, {{1, 2, -2}}, {{1, 2, 2}}, {{2, -2, -1}}, {{2, -2, 1}}, {{2, -1, -2}}, {{2, -1, 2}}, {{2, 1, -2}}, {{2, 1, 2}}, {{2, 2, -1}}, {{2, 2, 1}}, {{3, 0, 0}}},
                {{{-3, -1, 0}}, {{-3, 0, -1}}, {{-3, 0, 1}}, {{-3, 1, 0}}, {{-1, -3, 0}}, {{-1, 0, -3}}, {{-1, 0, 3}}, {{-1, 3, 0}}, {{0, -3, -1}}, {{0, -3, 1}}, {{0, -1, -3}}, {{0, -1, 3}}, {{0, 1, -3}}, {{0, 1, 3}}, {{0, 3, -1}}, {{0, 3, 1}}, {{1, -3, 0}}, {{1, 0, -3}}, {{1, 0, 3}}, {{1, 3, 0}}, {{3, -1, 0}}, {{3, 0, -1}}, {{3, 0, 1}}, {{3, 1, 0}}},
                {{{-3, -1, -1}}, {{-3, -1, 1}}, {{-3, 1, -1}}, {{-3, 1, 1}}, {{-1, -3, -1}}, {{-1, -3, 1}}, {{-1, -1, -3}}, {{-1, -1, 3}}, {{-1, 1, -3}}, {{-1, 1, 3}}, {{-1, 3, -1}}, {{-1, 3, 1}}, {{1, -3, -1}}, {{1, -3, 1}}, {{1, -1, -3}}, {{1, -1, 3}}, {{1, 1, -3}}, {{1, 1, 3}}, {{1, 3, -1}}, {{1, 3, 1}}, {{3, -1, -1}}, {{3, -1, 1}}, {{3, 1, -1}}, {{3, 1, 1}}},
                {{{-2, -2, -2}}, {{-2, -2, 2}}, {{-2, 2, -2}}, {{-2, 2, 2}}, {{2, -2, -2}}, {{2, -2, 2}}, {{2, 2, -2}}, {{2, 2, 2}}},
        };
        
      template<typename D>
        void display_point(const IntPoint<D>& p)
        {
                std::cout << "{";
                for ( int i = 0; i < D::Dimension_; ++i ) {
                        std::cout << p[i];
                        if ( i == D::Dimension_ - 1 ) {
                                std::cout << "}, ";
                        } else {
                                std::cout << ", ";
                        }
                }
        }
        
        
        const SpaceSize<TwoDimension> create_size(TwoDimension)
        {
                return SpaceSize<TwoDimension>{3, 3};
        }

        const SpaceSize<ThreeDimension> create_size(ThreeDimension)
        {
                return SpaceSize<ThreeDimension>{3, 3, 3};
        }
#if(LSM_FAST)
#else
        template<typename D>
        void test_initialize_distance_map()
        {
                Parameters params;
                params.wband_ = 3;
                SpaceSize<D> size = create_size(D());
                LevelSetMethod<D> lms{params, size};
                
                lms.initialize_distance_map();
                
                const auto& map = lms.get_distance_map();
                typedef typename std::multimap<double, IntPoint<D>>::value_type value_type;
                
                BOOST_CHECK(CheckerSelector<D>::Count_ == static_cast<int>(map.size()));
                auto beg = map.begin();
                auto end = map.end();
                int index = 0;
                while ( beg != end ) {
                        double key = beg->first;
//                        std::cout << key << " : ";
                        auto range = map.equal_range(key);
                        const std::ptrdiff_t size = std::distance(range.first, range.second);
                        BOOST_CHECK(CheckerSelector<D>::Sizes_[index] == size);
//                        std::cout << "[" << size << "]: ";
                        int k = 0;
                        boost::for_each(range,
                                [&](const value_type& v)
                                {
                                        const auto& p = v.second;
                                        BOOST_CHECK(CheckerSelector<D>::Points_[index][k] == p);
//                                        display_point<D>(p);
                                        ++k;
                                }
                        );
//                        std::cout << std::endl;
                        beg = range.second;
                        ++index;
                }
        }
#endif
        void test_initialize_over_all_2d()
        {
                Parameters params;
                params.wband_ = 3;
                
                InitialFront<TwoDimension> initial_front;
                initial_front.vertices_[0] = {{10, 15}};
                initial_front.vertices_[1] = {{82, 74}};
                SpaceSize<TwoDimension> size {101, 143};
                LevelSetMethod<TwoDimension> lms{params, size};
                lms.initialize_along_front(initial_front);
                lms.initialize_over_all(initial_front);
                const auto& phi = lms.get_phi();
                const auto& grid = lms.get_grid();
                const auto& statuses = lms.get_statuses();
                const auto& indexer = lms.get_indexer();
                
                const int width = size.width_;
                const int height = size.height_;
                InsideEstimator<TwoDimension> is_inside{grid};
                for ( int j = 0; j < height; ++j ) {
                        for ( int i = 0; i < width; ++i ) {
                                const IntPoint2d p {{i, j}};
                                int index = indexer(p);
                                if ( statuses[index] != Status::Front ) {
                                        if ( is_inside(p) ) {
                                                BOOST_CHECK(phi[index] == -params.wband_);
                                        } else {
                                                BOOST_CHECK(phi[index] == params.wband_);
                                        }
                                }
                        }
                }
        }
        
        void test_initialize_along_front_2d()
        {
                Parameters params;
                params.wband_ = 3;
                
                InitialFront<TwoDimension> initial_front;
                initial_front.vertices_[0] = {{10, 15}};
                initial_front.vertices_[1] = {{82, 74}};
                SpaceSize<TwoDimension> size {101, 143};
                LevelSetMethod<TwoDimension> lms{params, size};
                lms.initialize_along_front(initial_front);
                
                const auto& phi = lms.get_phi();
                const auto& statuses = lms.get_statuses();
                const auto& front = lms.get_front();
                Indexer<TwoDimension> indexer{size};
                const int left = initial_front.vertices_[0][0];
                const int top = initial_front.vertices_[0][1];
                const int right = initial_front.vertices_[1][0];
                const int bottom = initial_front.vertices_[1][1];
                int index;
                std::size_t k = 0;
                for ( int j = top; j < bottom; ++j ) {
                        index = indexer({{left, j}});
                        BOOST_CHECK(phi[index] == 0 && statuses[index] == Status::Front);
                        index = indexer({{right, j}});
                        BOOST_CHECK(phi[index] == 0 && statuses[index] == Status::Front);
                        k += 2;
                }

                for ( int i = left; i < right; ++i ) {
                        index = indexer({{i, top}});
                        BOOST_CHECK(phi[index] == 0 && statuses[index] == Status::Front);
                        index = indexer({{i, bottom}});
                        BOOST_CHECK(phi[index] == 0 && statuses[index] == Status::Front);
                        k += 2;
                }
                BOOST_CHECK(front.size() == k);
        }
        
        void test_initialize_over_all_3d()
        {
                Parameters params;
                params.wband_ = 3;
                
                InitialFront<ThreeDimension> initial_front;
                initial_front.vertices_[0] = {{10, 15, 32}};
                initial_front.vertices_[1] = {{82, 74, 61}};
                SpaceSize<ThreeDimension> size {101, 143, 131};
                LevelSetMethod<ThreeDimension> lms{params, size};
                lms.initialize_along_front(initial_front);
                lms.initialize_over_all(initial_front);
                const auto& phi = lms.get_phi();
                const auto& grid = lms.get_grid();
                const auto& statuses = lms.get_statuses();
                const auto& indexer = lms.get_indexer();
                
                const int width = size.width_;
                const int height = size.height_;
                const int depth = size.depth_;
                
                InsideEstimator<ThreeDimension> is_inside{grid};
                for ( int k = 0; k < depth; ++k ) {
                        for ( int j = 0; j < height; ++j ) {
                                for ( int i = 0; i < width; ++i ) {
                                        const IntPoint3d p {{i, j, k}};
                                        int index = indexer(p);
                                        if ( statuses[index] != Status::Front ) {
                                                if ( is_inside(p) ) {
                                                        BOOST_CHECK(phi[index] == -params.wband_);
                                                } else {
                                                        BOOST_CHECK(phi[index] == params.wband_);
                                                }
                                        }
                                }
                        }
                }
        }
        


        void test_initialize_along_front_3d()
        {
                Parameters params;
                params.wband_ = 3;
                
                InitialFront<ThreeDimension> initial_front;
                initial_front.vertices_[0] = {{10, 15, 32}};
                initial_front.vertices_[1] = {{82, 74, 61}};
                SpaceSize<ThreeDimension> size {101, 143, 131};
                LevelSetMethod<ThreeDimension> lms{params, size};
                lms.initialize_along_front(initial_front);
                
                const auto& phi = lms.get_phi();
                const auto& statuses = lms.get_statuses();
                //const auto& front = lms.get_front();
                Indexer<ThreeDimension> indexer{size};
                const int left = initial_front.vertices_[0][0];
                const int top = initial_front.vertices_[0][1];
                const int front_ = initial_front.vertices_[0][2];
                const int right = initial_front.vertices_[1][0];
                const int bottom = initial_front.vertices_[1][1];
                const int back = initial_front.vertices_[1][2];
                
                int index;
                std::size_t s = 0;
                for ( int j = top; j < bottom; ++j ) {
                        for ( int i = left; i < right; ++i ) {
                                index = indexer({{i, j, front_}});
                                BOOST_CHECK(phi[index] == 0 && statuses[index] == Status::Front);
                                index = indexer({{i, j, back}});
                                BOOST_CHECK(phi[index] == 0 && statuses[index] == Status::Front);
                                s += 2;
                        }
                }

                for ( int k = front_; k < back; ++k ) {
                        for ( int i = left; i < right; ++i ) {
                                index = indexer({{i, top, k}});
                                BOOST_CHECK(phi[index] == 0 && statuses[index] == Status::Front);
                                index = indexer({{i, bottom, k}});
                                BOOST_CHECK(phi[index] == 0 && statuses[index] == Status::Front);
                                s += 2;
                        }
                }

                for ( int j = top; j < bottom; ++j ) {
                        for ( int k = front_; k < back; ++k ) {
                                index = indexer({{left, j, k}});
                                BOOST_CHECK(phi[index] == 0 && statuses[index] == Status::Front);
                                index = indexer({{right, j, k}});
                                BOOST_CHECK(phi[index] == 0 && statuses[index] == Status::Front);
                                s += 2;
                        }
                }
        }
}

namespace lsm
{

struct LevelSetMethodTester
{
        template<typename D>
        static const Grid<D> create_speed_grid(LevelSetMethod<D>& lsm)
        {
                return lsm.create_space_without_edge(D());
        }

        template<typename D>
        static double set_speed_on_front(LevelSetMethod<D>& lsm)
        {
                return lsm.set_speed_on_front();
        }

        template<typename D>
        static std::vector<double>& get_speed(LevelSetMethod<D>& lsm)
        {
                return lsm.speed_;
        }

        template<typename D>
        static void clear_speed_within_narrow_band(LevelSetMethod<D>& lsm, bool resets)
        {
               lsm.clear_speed_within_narrow_band(resets);
        }

        template<typename D>
        static std::vector<IntPoint<D>>& get_narrow_band(LevelSetMethod<D>& lsm)
        {
                return lsm.narrow_band_;
        }

        template<typename D>
        static std::vector<double>& get_phi(LevelSetMethod<D>& lsm)
        {
                return lsm.phi_;
        }


        template<typename D>
        static std::vector<double>& get_dphi(LevelSetMethod<D>& lsm)
        {
                return lsm.dphi_;
        }

        template<typename D>
        static Front<D>& get_front(LevelSetMethod<D>& lsm)
        {
                return lsm.front_;
        }
        
        template<typename D>
        static void copy_nearest_speed_to_narrow_band(LevelSetMethod<D>& lsm, bool resets)
        {
                //const auto start = std::chrono::steady_clock::now();
                lsm.copy_nearest_speed_to_narrow_band(resets);
                //const auto end = std::chrono::steady_clock::now();
                //std::cout << std::chrono::duration<double, std::milli>(end - start).count() << std::endl;
                
        }

        template<typename D>
        static std::vector<Status>& get_statuses(LevelSetMethod<D>& lsm)
        {
                return lsm.statuses_;
        }
        
        template<typename D>
        static void register_to_narrow_band(LevelSetMethod<D>& lsm)
        {
                lsm.register_to_narrow_band(D());
        }
        
        template<typename D>
        static const std::vector<IntPoint<D>>& get_narrow_band(const LevelSetMethod<D>& lsm)
        {
                return lsm.narrow_band_;
        }
        
};
} // namespace lsm

namespace 
{
        void test_create_speed_grid_2d()
        {
                Parameters params;
                params.wband_ = 3;
                SpaceSize<TwoDimension> size {101, 143};
                LevelSetMethod<TwoDimension> lsm{params, size};
                const auto grid = LevelSetMethodTester::create_speed_grid(lsm);
                BOOST_CHECK(grid.left_ == 0);
                BOOST_CHECK(grid.right_ == size.width_ - 1);
                BOOST_CHECK(grid.top_ == 0);
                BOOST_CHECK(grid.bottom_ == size.height_ - 1);
        }

        void test_create_speed_grid_3d()
        {
                Parameters params;
                params.wband_ = 3;
                SpaceSize<ThreeDimension> size {101, 143, 3};
                LevelSetMethod<ThreeDimension> lsm{params, size};
                const auto grid = LevelSetMethodTester::create_speed_grid(lsm);
                BOOST_CHECK(grid.left_ == 0);
                BOOST_CHECK(grid.right_ == size.width_ - 1);
                BOOST_CHECK(grid.top_ == 0);
                BOOST_CHECK(grid.bottom_ == size.height_ - 1);
                BOOST_CHECK(grid.front_ == 0);
                BOOST_CHECK(grid.back_ == size.depth_ - 1);
        }

        const std::vector<std::uint8_t> get_input_gray(
                const SpaceSize<TwoDimension>& size,
                const InitialFront<TwoDimension>& front)
        {
                std::vector<std::uint8_t> gray(size.total_, 1);
                const int left = front.vertices_[0][0];
                const int top = front.vertices_[0][1];
                const int right = front.vertices_[1][0];
                const int bottom = front.vertices_[1][1];

                Indexer<TwoDimension> indexer{size}; 
                for ( int i = left; i <= right; ++i ) {
                        gray[indexer({{i, top}})] = 0;
                        gray[indexer({{i, bottom}})] = 0;
                }
                for ( int i = top; i <= bottom; ++i ) {
                        gray[indexer({{left, i}})] = 0;
                        gray[indexer({{right, i}})] = 0;
                }

                return std::move(gray);
        }

        // !!!
        const std::vector<std::uint8_t> get_input_gray(
                const SpaceSize<ThreeDimension>& size,
                const InitialFront<ThreeDimension>& front)
        {
                std::vector<std::uint8_t> gray(size.total_, 1);
                const int left = front.vertices_[0][0];
                const int top = front.vertices_[0][1];
                const int front_ = front.vertices_[0][2];

                const int right = front.vertices_[1][0];
                const int bottom = front.vertices_[1][1];
                const int back = front.vertices_[1][2];

                Indexer<ThreeDimension> indexer{size}; 
                // xy
                for ( int j = top; j < bottom; ++j ) {
                        for ( int i = left; i < right; ++i ) {
                                gray[indexer({{i, j, front_}})] = 0;
                                gray[indexer({{i, j, back}})] = 0;
                        }
                }
                // yz
                for ( int j = top; j < bottom; ++j ) {
                        for ( int k = front_; k < back; ++k ) {
                                gray[indexer({{left, j, k}})] = 0;
                                gray[indexer({{right, j, k}})] = 0;
                        }
                }
                // xz
                for ( int i = left; i < right; ++i ) {
                        for ( int k = front_; k < back; ++k ) {
                                gray[indexer({{i, top, k}})] = 0;
                                gray[indexer({{i, bottom, k}})] = 0;
                        }
                }
                return std::move(gray);
        }

        template<typename T, typename U>
        void display_buffer(const SpaceSize<TwoDimension>& size, const std::vector<T>& buffer)
        {
                std::cout << "{\n";
                for ( int j = 0, wj = j * size.width_; j < size.height_; ++j, wj += size.width_ ) {
                        const T* p = &buffer[wj];
                        for ( int i = 0; i < size.width_; ++i ) {
                                std::cout << std::setw(10) << static_cast<U>(p[i]) << ", ";
                        }
                        std::cout << std::endl;
                }
                std::cout << "}\n";
        }

        void test_set_speed_on_front_2d()
        {
                Parameters params;
                params.wband_ = 1;
                params.constant_speed_ = 1;
                params.gain_ = 2;
                
                InitialFront<TwoDimension> initial_front;
                const int left = 2;
                const int top = 3;
                const int right = 8;
                const int bottom = 7;
                initial_front.vertices_[0] = {{left, top}};
                initial_front.vertices_[1] = {{right, bottom}};
                
                SpaceSize<TwoDimension> size {11, 11};
                LevelSetMethod<TwoDimension> lsm{params, size};
                lsm.initialize_along_front(initial_front);
                lsm.initialize_over_all(initial_front);

                // input image
                std::vector<std::uint8_t> input = get_input_gray(size, initial_front);
                //display_buffer<std::uint8_t, int>(size, input);

                // set an input image
                auto& gray = lsm.input_object();
                gray.resize(size.total_);
                std::memcpy(&gray[0], &input[0], size.total_);
 
                auto& speed = LevelSetMethodTester::get_speed(lsm);
                speed.resize(size.total_, 0);

                lsm.calculate_speed_factors();
                
                //display_buffer<double, double>(size, speed);
                double fs = LevelSetMethodTester::set_speed_on_front(lsm);
                //std::cout << "fs = " << fs << std::endl;
//                display_buffer<double, double>(size, speed);

                BOOST_CHECK(fs != 0.0); 
                for (int j = 0; j < size.height_; ++j) {
                        int wj = size.width_ * j;
                        const double* p = &speed[wj];
                        for ( int i = 0; i < size.width_; ++i ) {
                                if ( left <= i && i <= right && j == top ) {
                                        BOOST_CHECK(p[i] != 0.0);
                                } else if ( left <= i && i <= right && j == bottom ) {
                                        BOOST_CHECK(p[i] != 0.0);
                                } else if ( i == right && top <= j && j <= bottom ) {
                                        BOOST_CHECK(p[i] != 0.0); 
                                } else if ( i == left && top <= j && j <= bottom ) {
                                        BOOST_CHECK(p[i] != 0.0); 
                                } else {
                                        BOOST_CHECK(p[i] == 0);
                                }
                        }
                }
                      //speed_[i] *= (parameters_.constant_speed_ - parameters_.gain_ * kappa);
        }
        
        void test_set_speed_on_front_3d()
        {
                Parameters params;
                params.wband_ = 1;
                params.constant_speed_ = 1;
                params.gain_ = 2;
                
                InitialFront<ThreeDimension> initial_front;
                const int left = 2;
                const int top = 3;
                const int right = 8;
                const int bottom = 7;
                const int front = 3;
                const int back = 7;

                initial_front.vertices_[0] = {{left, top, front}};
                initial_front.vertices_[1] = {{right, bottom, back}};
                
                SpaceSize<ThreeDimension> size {11, 11, 11};
                LevelSetMethod<ThreeDimension> lsm{params, size};
                lsm.initialize_along_front(initial_front);
                lsm.initialize_over_all(initial_front);
                
              // input image
                std::vector<std::uint8_t> input = get_input_gray(size, initial_front);

                // set an input image
                auto& gray = lsm.input_object();
                gray.resize(size.total_);
                std::memcpy(&gray[0], &input[0], size.total_);

                auto& speed = LevelSetMethodTester::get_speed(lsm);
                speed.resize(size.total_, 0);

                lsm.calculate_speed_factors();
                lsm.initialize_narrow_band();

                double fs = LevelSetMethodTester::set_speed_on_front(lsm);
                BOOST_CHECK(fs != 0.0);

                const int width = size.width_;
                const int height = size.height_;
                const int area = width * height;
                const int depth = size.depth_;

                for ( int k = 0, ak = area * k; k < depth; ++k, ak += area ) { 
                        for ( int j = 0, wj = ak + width * j; j < height; ++j, j += width ) {
                                const double* p = &speed[wj];
                                //std::cout << "p = " << p << std::endl;
                                for ( int i = 0; i < size.width_; ++i ) {
                                        if ( (left <= i && i <= right && top <= j && j <= bottom) 
                                                        && (k == front || k == back) ) {
                                                BOOST_CHECK(p[i] != 0.0);
                                        } else if ( (left <= i && i <= right && front <= k && k <= back) 
                                                        && (j == top || j == bottom) ) {
                                                BOOST_CHECK(p[i] != 0.0);
                                        } else if ( (top <= j && j <= bottom && front <= k && k <= back) 
                                                        && (i == left || i == right) ) { 
                                                BOOST_CHECK(p[i] != 0.0);
                                        } else {
                                                //std::cout << "i = " << i << std::endl;
                                                BOOST_CHECK(p[i] == 0.0);
                                        }
                                }
                        }
                }
        }

        void set_narrow_band(std::vector<IntPoint<TwoDimension>>& narrow_band)
        {
                for ( int i = 1; i < 10; ++i ) {
                       narrow_band.push_back({{i, 2}});
                       narrow_band.push_back({{i, 3}});
                       narrow_band.push_back({{i, 4}});
                       narrow_band.push_back({{i, 6}});
                       narrow_band.push_back({{i, 7}});
                       narrow_band.push_back({{i, 8}});
                }
               narrow_band.push_back({{1, 5}});
               narrow_band.push_back({{2, 5}});
               narrow_band.push_back({{3, 5}});
               narrow_band.push_back({{7, 5}});
               narrow_band.push_back({{8, 5}});
               narrow_band.push_back({{9, 5}});
        }

        void set_narrow_band(std::vector<IntPoint<ThreeDimension>>& narrow_band)
        {
                for ( int k = 2; k <= 8; ++k ) {
                        for ( int j = 2; j <= 8; ++j ) {
                                for ( int i = 1; i <= 9; ++i ) {
                                        if ( j == 5 && k == 5 && (i == 4 || i == 5 || i == 6) ) {
                                                //        
                                        } else {
                                                narrow_band.push_back({{i, j, k}}); 
                                        }
                                }
                        }
                }
        }

        bool is_within_narrow_band(const IntPoint3d& p)
        {
                int i = p[0];
                int j = p[1];
                int k = p[2];
                if ( 2 <= k && k <= 8 ) {
                        if ( 2 <= j && j <= 8 ) {
                                if ( 1 <= i && i <= 9 ) {
                                        if ( j == 5 && k == 5 && (i == 4 || i == 5 || i == 6) ) {
                                                        
                                        } else {
                                                return true; 
                                        }
                                }
                        }
                }
                return false;
        }

        void check_buffer(const std::vector<double>& buffer, const SpaceSize<TwoDimension>& size)
        {
                const int w = size.width_;
                const int h = size.height_;
                
                for ( int j = 0, wj = w * j; j < h; ++j, wj += w ) {
                        const double* p = &buffer[wj];
                        for ( int i = 0; i < w; ++i ) {
                                if ( 1 <= i && i <= 9 ) {
                                        if ( j == 2 || j == 3 || j == 4 || j == 6 || j == 7 || j == 8 ) {
                                                BOOST_CHECK(p[i] == 0.0); 
                                        }
                                } else if ( j == 5 ) {
                                        if ( i == 1 || i == 2 || i == 3 || i == 7 || i == 8 || i == 9 ) {
                                                BOOST_CHECK(p[i] == 0.0); 
                                        }
                                } else {
                                        BOOST_CHECK(p[i] == 1.0); 
                                }
                        }
                }
        }

        void check_buffer(const std::vector<double>& buffer, const SpaceSize<ThreeDimension>& size)
        {
                const int w = size.width_;
                const int h = size.height_;
                const int d = size.depth_;
                const int a = w * h;

                for ( int k = 0, ak = a * k; k < d; ++k, ak += a ) {
                        for ( int j = 0, wj = ak + w * j; j < h; ++j, wj += w ) {
                                const double* p = &buffer[wj];
                                for ( int i = 0; i < w; ++i ) {
                                        if ( is_within_narrow_band(IntPoint3d{{i, j, k}}) ) {
                                                BOOST_CHECK(p[i] == 0.0);
                                        } else {
                                                BOOST_CHECK(p[i] == 1.0);
                                        }
                                }
                        }
                }
        }
        void test_clear_speed_within_narrow_band_2d() 
        {
                Parameters params;
                params.wband_ = 1;
                params.constant_speed_ = 1;
                params.gain_ = 2;
                
                InitialFront<TwoDimension> initial_front;
                const int left = 2;
                const int top = 3;
                const int right = 8;
                const int bottom = 7;
                initial_front.vertices_[0] = {{left, top}};
                initial_front.vertices_[1] = {{right, bottom}};
                
                SpaceSize<TwoDimension> size {11, 11};
                LevelSetMethod<TwoDimension> lsm {params, size};
                lsm.initialize_along_front(initial_front);
                lsm.initialize_over_all(initial_front);

                // narrow band
                auto& narrow_band = LevelSetMethodTester::get_narrow_band(lsm);
                set_narrow_band(narrow_band);

                // phi
                auto& dphi =  LevelSetMethodTester::get_dphi(lsm);
                std::vector<double>().swap(dphi);
                dphi.resize(size.total_, 1);

                // speed 
                auto& speed =  LevelSetMethodTester::get_speed(lsm);
                std::vector<double>().swap(speed);
                speed.resize(size.total_, 1);

                LevelSetMethodTester::clear_speed_within_narrow_band(lsm, true);
                check_buffer(speed, size);
                check_buffer(dphi, size);
        }

        void test_clear_speed_within_narrow_band_3d() 
        {
                Parameters params;
                params.wband_ = 1;
                params.constant_speed_ = 1;
                params.gain_ = 2;
                
                InitialFront<ThreeDimension> initial_front;
                const int left = 2;
                const int top = 3;
                const int right = 8;
                const int bottom = 7;
                const int front = 3;
                const int back = 7;
                initial_front.vertices_[0] = {{left, top, front}};
                initial_front.vertices_[1] = {{right, bottom, back}};
                
                SpaceSize<ThreeDimension> size {11, 11, 11};
                LevelSetMethod<ThreeDimension> lsm {params, size};
                lsm.initialize_along_front(initial_front);
                lsm.initialize_over_all(initial_front);

                // narrow band
                auto& narrow_band = LevelSetMethodTester::get_narrow_band(lsm);
                set_narrow_band(narrow_band);

                // phi
                auto& dphi =  LevelSetMethodTester::get_dphi(lsm);
                std::vector<double>().swap(dphi);
                dphi.resize(size.total_, 1);

                // speed
                auto& speed =  LevelSetMethodTester::get_speed(lsm);
                std::vector<double>().swap(speed);
                speed.resize(size.total_, 1);

                LevelSetMethodTester::clear_speed_within_narrow_band(lsm, true);
                check_buffer(speed, size);
                check_buffer(dphi, size);
       }

        
        void display_front(const Front<ThreeDimension>& front)
        {
              for ( const auto& p : front ) {
                      std::cout << "(" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";
              }
        }

        void test_front_2d()
        {
                Parameters params;
                params.wband_ = 3;
                
                InitialFront<TwoDimension> initial_front;
                initial_front.vertices_[0] = {{0, 0}};
                initial_front.vertices_[1] = {{2, 2}};
                SpaceSize<TwoDimension> size {3, 3};
                LevelSetMethod<TwoDimension> lsm {params, size};
                lsm.initialize_along_front(initial_front);
                //display_front(LevelSetMethodTester::get_front(lsm));
                BOOST_CHECK(LevelSetMethodTester::get_front(lsm).size() == 8);
        }

        void test_front_3d()
        {
                Parameters params;
                params.wband_ = 3;
                
                InitialFront<ThreeDimension> initial_front;
                initial_front.vertices_[0] = {{0, 0, 0}};
                initial_front.vertices_[1] = {{2, 2, 2}};
                SpaceSize<ThreeDimension> size {3, 3, 3};
                LevelSetMethod<ThreeDimension> lsm {params, size};
                lsm.initialize_along_front(initial_front);
                //display_front(LevelSetMethodTester::get_front(lsm));
                BOOST_CHECK(LevelSetMethodTester::get_front(lsm).size() == 26);
        }
        
        void test_set_speed_function_2d()
        {
                Parameters params;
                params.wband_ = 1;
                params.constant_speed_ = 1;
                params.gain_ = 2;
                
                InitialFront<TwoDimension> initial_front;
                const int left = 2;
                const int top = 3;
                const int right = 8;
                const int bottom = 7;
                initial_front.vertices_[0] = {{left, top}};
                initial_front.vertices_[1] = {{right, bottom}};
                
                SpaceSize<TwoDimension> size {11, 11};
                LevelSetMethod<TwoDimension> lsm {params, size};
                lsm.initialize_along_front(initial_front);
                lsm.initialize_over_all(initial_front);

                // input image
                std::vector<std::uint8_t> input = get_input_gray(size, initial_front);

                // set an input image
                auto& gray = lsm.input_object();
                gray.resize(size.total_);
                std::memcpy(&gray[0], &input[0], size.total_);
                lsm.calculate_speed_factors();
                
                
                bool resets = true;
                bool is_stopped = lsm.set_speed_function(resets);
                const auto& speed = LevelSetMethodTester::get_speed(lsm);
                BOOST_CHECK(!is_stopped);
                for ( int j = 0, wj = j * size.width_; j < size.height_; ++j, wj += size.width_ ) {
                        const double* p = &speed[wj];
                        for ( int i = 0; i < size.width_; ++i ) {
                                if ( left <= i && i <= right && j == top ) {
                                        BOOST_CHECK(p[i] != 0.0);
                                } else if ( left <= i && i <= right && j == bottom ) {
                                        BOOST_CHECK(p[i] != 0.0);
                                } else if ( i == right && top <= j && j <= bottom ) {
                                        BOOST_CHECK(p[i] != 0.0); 
                                } else if ( i == left && top <= j && j <= bottom ) {
                                        BOOST_CHECK(p[i] != 0.0); 
                                } else {
                                        BOOST_CHECK(p[i] == 0);
                                }
                        }
                }
        }

        void test_set_speed_function_3d()
        {
                Parameters params;
                params.wband_ = 1;
                params.constant_speed_ = 1;
                params.gain_ = 2;
                
                InitialFront<ThreeDimension> initial_front;
                const int left = 2;
                const int top = 3;
                const int right = 8;
                const int bottom = 7;
                const int front = 3;
                const int back = 7;

                initial_front.vertices_[0] = {{left, top, front}};
                initial_front.vertices_[1] = {{right, bottom, back}};
                
                SpaceSize<ThreeDimension> size {11, 11, 11};
                LevelSetMethod<ThreeDimension> lsm {params, size};
                lsm.initialize_along_front(initial_front);
                lsm.initialize_over_all(initial_front);

              // input image
                std::vector<std::uint8_t> input = get_input_gray(size, initial_front);

                // set an input image
                auto& gray = lsm.input_object();
                gray.resize(size.total_);
                std::memcpy(&gray[0], &input[0], size.total_);

                lsm.calculate_speed_factors();
                lsm.initialize_narrow_band();

                bool resets = true;
                bool is_stopped = lsm.set_speed_function(resets);
                BOOST_CHECK(!is_stopped);

                const auto& speed = LevelSetMethodTester::get_speed(lsm);
                const int width = size.width_;
                const int height = size.height_;
                const int area = width * height;
                const int depth = size.depth_;

                for ( int k = 0, ak = area * k; k < depth; ++k, ak += area ) { 
                        for ( int j = 0, wj = ak + width * j; j < height; ++j, j += width ) {
                                const double* p = &speed[wj];
                                for ( int i = 0; i < size.width_; ++i ) {
                                        if ( (left <= i && i <= right && top <= j && j <= bottom) 
                                                        && (k == front || k == back) ) {
                                                BOOST_CHECK(p[i] != 0.0);
                                        } else if ( (left <= i && i <= right && front <= k && k <= back) 
                                                        && (j == top || j == bottom) ) {
                                                BOOST_CHECK(p[i] != 0.0);
                                        } else if ( (top <= j && j <= bottom && front <= k && k <= back) 
                                                        && (i == left || i == right) ) { 
                                                BOOST_CHECK(p[i] != 0.0);
                                        } else {
                                                BOOST_CHECK(p[i] == 0.0);
                                        }
                                }
                        }
                }
        }

        void test_copy_nearest_speed_to_narrow_band_2d()
        {
                
                Parameters params;
                params.wband_ = 3;
                params.constant_speed_ = 1;
                params.gain_ = 2;
                params.wreset_ = 1;
                
                InitialFront<TwoDimension> initial_front;
                const int left = 2;
                const int top = 3;
                const int right = 8;
                const int bottom = 7;
                initial_front.vertices_[0] = {{left, top}};
                initial_front.vertices_[1] = {{right, bottom}};
                
                SpaceSize<TwoDimension> size {11, 11};
                LevelSetMethod<TwoDimension> lsm {params, size};
                lsm.initialize_distance_map();
                lsm.initialize_along_front(initial_front);
                lsm.initialize_over_all(initial_front);
                
                // input image
                std::vector<std::uint8_t> input = get_input_gray(size, initial_front);

                // set an input image
                auto& gray = lsm.input_object();
                gray.resize(size.total_);
                std::memcpy(&gray[0], &input[0], size.total_);
                
                lsm.calculate_speed_factors();
                
                auto& speed = LevelSetMethodTester::get_speed(lsm);
                speed.resize(size.total_);
                auto& dphi = LevelSetMethodTester::get_dphi(lsm);
                dphi.resize(size.total_);
                bool resets = true;
                LevelSetMethodTester::clear_speed_within_narrow_band(lsm, resets);
                LevelSetMethodTester::set_speed_on_front(lsm);

                const auto& statuses = lsm.get_statuses();
                LevelSetMethodTester::copy_nearest_speed_to_narrow_band(lsm, resets);
                
                const double squared_phi_answers[] = {
                        9,    10,     9,     9,     9,     9,     9,     9,     9,    10,     9,
                        8,     5,     4,     4,     4,     4,     4,     4,     4,     5,     8,
                        5,     2,     1,     1,     1,     1,     1,     1,     1,     2,     5,
                        4,     1,     0,     0,     0,     0,     0,     0,     0,     1,     4,
                        4,     1,     0,    -1,    -1,    -1,    -1,    -1,     0,     1,     4,
                        4,     1,     0,    -1,    -4,    -4,    -4,    -1,     0,     1,     4,
                        4,     1,     0,    -1,    -1,    -1,    -1,    -1,     0,     1,     4,
                        4,     1,     0,     0,     0,     0,     0,     0,     0,     1,     4,
                        5,     2,     1,     1,     1,     1,     1,     1,     1,     2,     5,
                        8,     5,     4,     4,     4,     4,     4,     4,     4,     5,     8,
                        9,    10,     9,     9,     9,     9,     9,     9,     9,    10,     9,
                };
    
                const auto& phi = lsm.get_phi();
//                display_buffer<double, double>(size, phi);
                float epsilon = 1.0e-03f;
                double phi_answer;
                for ( std::size_t i = 0; i < phi.size(); ++i ) {
                        if ( squared_phi_answers[i] < 0 ) {
                                phi_answer = -std::sqrt(-squared_phi_answers[i]);
                        } else {
                                phi_answer = std::sqrt(squared_phi_answers[i]);
                        }
                        BOOST_CHECK_CLOSE(phi[i], phi_answer, epsilon);
                }
                
                const double speed_answers[] = {
                         0,   -1.35083,   -1.35083,   0.234473,          1,          1,          1,   0.234473,   -1.35083,   -1.35083,          0,
                  -1.35083,   -1.35083,   -1.35083,   0.234473,          1,          1,          1,   0.234473,   -1.35083,   -1.35083,   -1.35083, 
                  -1.35083,   -1.35083,   -1.35083,   0.234473,          1,          1,          1,   0.234473,   -1.35083,   -1.35083,   -1.35083, 
                  -1.35083,   -1.35083,   -1.35083,   0.234473,          1,          1,          1,   0.234473,   -1.35083,   -1.35083,   -1.35083, 
                  0.234473,   0.234473,   0.234473,   0.234473,          1,          1,          1,   0.234473,   0.234473,   0.234473,   0.234473, 
                         1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1, 
                  0.234473,   0.234473,   0.234473,   0.234473,          1,          1,          1,   0.234473,   0.234473,   0.234473,   0.234473, 
                  -1.35083,   -1.35083,   -1.35083,   0.234473,          1,          1,          1,   0.234473,   -1.35083,   -1.35083,   -1.35083, 
                  -1.35083,   -1.35083,   -1.35083,   0.234473,          1,          1,          1,   0.234473,   -1.35083,   -1.35083,   -1.35083, 
                  -1.35083,   -1.35083,   -1.35083,   0.234473,          1,          1,          1,   0.234473,   -1.35083,   -1.35083,   -1.35083, 
                         0,   -1.35083,   -1.35083,   0.234473,          1,          1,          1,   0.234473,   -1.35083,   -1.35083,          0, 
                };
                
                for ( std::size_t i = 0; i < speed.size(); ++i ) {
                        BOOST_CHECK_CLOSE(speed[i], speed_answers[i], epsilon);
                }
//                display_buffer<double, double>(size, speed);

                const int status_answers[] = {
                         0,          2,          2,          2,          2,          2,          2,          2,          2,          2,          0,
                         2,          2,          1,          1,          1,          1,          1,          1,          1,          2,          2,
                         2,          1,          1,          1,          1,          1,          1,          1,          1,          1,          2,
                         1,          1,          3,          3,          3,          3,          3,          3,          3,          1,          1,
                         1,          1,          3,          1,          1,          1,          1,          1,          3,          1,          1,
                         1,          1,          3,          1,          1,          1,          1,          1,          3,          1,          1,
                         1,          1,          3,          1,          1,          1,          1,          1,          3,          1,          1,
                         1,          1,          3,          3,          3,          3,          3,          3,          3,          1,          1,
                         2,          1,          1,          1,          1,          1,          1,          1,          1,          1,          2,
                         2,          2,          1,          1,          1,          1,          1,          1,          1,          2,          2,
                         0,          2,          2,          2,          2,          2,          2,          2,          2,          2,          0,
                };
                
                for ( std::size_t i = 0; i < statuses.size(); ++i ) {
                        BOOST_CHECK(static_cast<int>(statuses[i]) == status_answers[i]);
                }
//                
        }

        template<typename T, typename U>
        void display_buffer(const SpaceSize<ThreeDimension>& size, const std::vector<T>& buffer)
        {
                const int w = size.width_;
                const int h = size.height_;
                const int d = size.depth_;
                const int a = w * h;
                std::cout << "{\n";
                for ( int k = 0, ak = k * a; k < d; ++k, ak += a ) {
                        for ( int j = 0, wj = ak + j * w; j < h; ++j, wj += w ) {
                                const T* p = &buffer[wj];
                                for ( int i = 0; i < w; ++i ) {
                                        std::cout << std::setw(10) << static_cast<U>(p[i]) << ",";
                                }
                                std::cout << std::endl;
                        }
                        std::cout << std::endl;
                }
                std::cout << "}";
        }

        void test_copy_nearest_speed_to_narrow_band_3d()
        {
                Parameters params;
                params.wband_ = 3;
                params.constant_speed_ = 1;
                params.gain_ = 2;
                params.wreset_ = 1;
                
                InitialFront<ThreeDimension> initial_front;
                const int left = 2;
                const int top = 3;
                const int right = 8;
                const int bottom = 7;
                const int front = 3;
                const int back = 7;

                initial_front.vertices_[0] = {{left, top, front}};
                initial_front.vertices_[1] = {{right, bottom, back}};
                
                SpaceSize<ThreeDimension> size {11, 11, 11};
                LevelSetMethod<ThreeDimension> lsm {params, size};
                lsm.initialize_distance_map();
                lsm.initialize_along_front(initial_front);
                lsm.initialize_over_all(initial_front);

                // input image
                std::vector<std::uint8_t> input = get_input_gray(size, initial_front);

                // set an input image
                auto& gray = lsm.input_object();
                gray.resize(size.total_);
                std::memcpy(&gray[0], &input[0], size.total_);

                lsm.calculate_speed_factors();
//                lsm.initialize_narrow_band();

                auto& speed = LevelSetMethodTester::get_speed(lsm);
                
                speed.resize(size.total_);
                auto& dphi = LevelSetMethodTester::get_dphi(lsm);
                dphi.resize(size.total_);
                bool resets = true;

                LevelSetMethodTester::clear_speed_within_narrow_band(lsm, resets);
                LevelSetMethodTester::set_speed_on_front(lsm);
                
                const auto& statuses = lsm.get_statuses();
                LevelSetMethodTester::copy_nearest_speed_to_narrow_band(lsm, resets);
                
                const int status_answers[] = {
                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,

                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,
                         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,
                         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,
                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,

                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,

                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         1,         1,         3,         3,         3,         3,         3,         3,         3,         1,         1,
                         1,         1,         3,         3,         3,         3,         3,         3,         3,         1,         1,
                         1,         1,         3,         3,         3,         3,         3,         3,         3,         1,         1,
                         1,         1,         3,         3,         3,         3,         3,         3,         3,         1,         1,
                         1,         1,         3,         3,         3,         3,         3,         3,         3,         1,         1,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,

                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         1,         1,         3,         3,         3,         3,         3,         3,         3,         1,         1,
                         1,         1,         3,         1,         1,         1,         1,         1,         3,         1,         1,
                         1,         1,         3,         1,         1,         1,         1,         1,         3,         1,         1,
                         1,         1,         3,         1,         1,         1,         1,         1,         3,         1,         1,
                         1,         1,         3,         3,         3,         3,         3,         3,         3,         1,         1,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,

                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         1,         1,         3,         3,         3,         3,         3,         3,         3,         1,         1,
                         1,         1,         3,         1,         1,         1,         1,         1,         3,         1,         1,
                         1,         1,         3,         1,         1,         1,         1,         1,         3,         1,         1,
                         1,         1,         3,         1,         1,         1,         1,         1,         3,         1,         1,
                         1,         1,         3,         3,         3,         3,         3,         3,         3,         1,         1,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,

                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         1,         1,         3,         3,         3,         3,         3,         3,         3,         1,         1,
                         1,         1,         3,         1,         1,         1,         1,         1,         3,         1,         1,
                         1,         1,         3,         1,         1,         1,         1,         1,         3,         1,         1,
                         1,         1,         3,         1,         1,         1,         1,         1,         3,         1,         1,
                         1,         1,         3,         3,         3,         3,         3,         3,         3,         1,         1,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,

                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         1,         1,         3,         3,         3,         3,         3,         3,         3,         1,         1,
                         1,         1,         3,         3,         3,         3,         3,         3,         3,         1,         1,
                         1,         1,         3,         3,         3,         3,         3,         3,         3,         1,         1,
                         1,         1,         3,         3,         3,         3,         3,         3,         3,         1,         1,
                         1,         1,         3,         3,         3,         3,         3,         3,         3,         1,         1,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,

                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         2,         1,         1,         1,         1,         1,         1,         1,         1,         1,         2,
                         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,

                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,
                         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         2,         2,         1,         1,         1,         1,         1,         1,         1,         2,         2,
                         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,
                         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,         2,
                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,

                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         0,         2,         2,         2,         2,         2,         2,         2,         2,         2,         0,
                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                };
                for ( std::size_t i = 0; i < statuses.size(); ++i ) {
                        BOOST_CHECK(static_cast<int>(statuses[i]) == status_answers[i]);
                }

                const auto& phi = lsm.get_phi();
                const double phi_answers[] = {
                         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,
                         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,
                         3,   3.31662,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.31662,         3,
                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,
                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,
                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,
                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,
                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,
                         3,   3.31662,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.31662,         3,
                         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,
                         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,

                         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,
                    3.4641,         3,   2.82843,   2.82843,   2.82843,   2.82843,   2.82843,   2.82843,   2.82843,         3,    3.4641,
                         3,   2.44949,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.44949,         3,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                         3,   2.44949,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.44949,         3,
                    3.4641,         3,   2.82843,   2.82843,   2.82843,   2.82843,   2.82843,   2.82843,   2.82843,         3,    3.4641,
                         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,

                         3,   3.31662,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.31662,         3,
                         3,   2.44949,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.44949,         3,
                   2.44949,   1.73205,   1.41421,   1.41421,   1.41421,   1.41421,   1.41421,   1.41421,   1.41421,   1.73205,   2.44949,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                   2.44949,   1.73205,   1.41421,   1.41421,   1.41421,   1.41421,   1.41421,   1.41421,   1.41421,   1.73205,   2.44949,
                         3,   2.44949,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.44949,         3,
                         3,   3.31662,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.31662,         3,

                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                         2,         1,         0,         0,         0,         0,         0,         0,         0,         1,         2,
                         2,         1,         0,         0,         0,         0,         0,         0,         0,         1,         2,
                         2,         1,         0,         0,         0,         0,         0,         0,         0,         1,         2,
                         2,         1,         0,         0,         0,         0,         0,         0,         0,         1,         2,
                         2,         1,         0,         0,         0,         0,         0,         0,         0,         1,         2,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,

                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                         2,         1,         0,         0,         0,         0,         0,         0,         0,         1,         2,
                         2,         1,         0,        -1,        -1,        -1,        -1,        -1,         0,         1,         2,
                         2,         1,         0,        -1,        -1,        -1,        -1,        -1,         0,         1,         2,
                         2,         1,         0,        -1,        -1,        -1,        -1,        -1,         0,         1,         2,
                         2,         1,         0,         0,         0,         0,         0,         0,         0,         1,         2,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,

                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                         2,         1,         0,         0,         0,         0,         0,         0,         0,         1,         2,
                         2,         1,         0,        -1,        -1,        -1,        -1,        -1,         0,         1,         2,
                         2,         1,         0,        -1,        -2,        -2,        -2,        -1,         0,         1,         2,
                         2,         1,         0,        -1,        -1,        -1,        -1,        -1,         0,         1,         2,
                         2,         1,         0,         0,         0,         0,         0,         0,         0,         1,         2,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,

                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                         2,         1,         0,         0,         0,         0,         0,         0,         0,         1,         2,
                         2,         1,         0,        -1,        -1,        -1,        -1,        -1,         0,         1,         2,
                         2,         1,         0,        -1,        -1,        -1,        -1,        -1,         0,         1,         2,
                         2,         1,         0,        -1,        -1,        -1,        -1,        -1,         0,         1,         2,
                         2,         1,         0,         0,         0,         0,         0,         0,         0,         1,         2,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,

                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                         2,         1,         0,         0,         0,         0,         0,         0,         0,         1,         2,
                         2,         1,         0,         0,         0,         0,         0,         0,         0,         1,         2,
                         2,         1,         0,         0,         0,         0,         0,         0,         0,         1,         2,
                         2,         1,         0,         0,         0,         0,         0,         0,         0,         1,         2,
                         2,         1,         0,         0,         0,         0,         0,         0,         0,         1,         2,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,

                         3,   3.31662,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.31662,         3,
                         3,   2.44949,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.44949,         3,
                   2.44949,   1.73205,   1.41421,   1.41421,   1.41421,   1.41421,   1.41421,   1.41421,   1.41421,   1.73205,   2.44949,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                   2.23607,   1.41421,         1,         1,         1,         1,         1,         1,         1,   1.41421,   2.23607,
                   2.44949,   1.73205,   1.41421,   1.41421,   1.41421,   1.41421,   1.41421,   1.41421,   1.41421,   1.73205,   2.44949,
                         3,   2.44949,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.44949,         3,
                         3,   3.31662,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.31662,         3,

                         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,
                    3.4641,         3,   2.82843,   2.82843,   2.82843,   2.82843,   2.82843,   2.82843,   2.82843,         3,    3.4641,
                         3,   2.44949,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.44949,         3,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                   2.82843,   2.23607,         2,         2,         2,         2,         2,         2,         2,   2.23607,   2.82843,
                         3,   2.44949,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.23607,   2.44949,         3,
                    3.4641,         3,   2.82843,   2.82843,   2.82843,   2.82843,   2.82843,   2.82843,   2.82843,         3,    3.4641,
                         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,

                         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,
                         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,
                         3,   3.31662,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.31662,         3,
                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,
                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,
                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,
                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,
                         3,   3.16228,         3,         3,         3,         3,         3,         3,         3,   3.16228,         3,
                         3,   3.31662,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.16228,   3.31662,         3,
                         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,
                         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,         3,
                };
                float epsilon = 1.0e-03f;
                for ( std::size_t i = 0; i < phi.size(); ++i ) {
                        BOOST_CHECK_CLOSE(phi[i], phi_answers[i], epsilon);
                }

                const double speed_answers[] = {
                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                         0,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.59222,  -3.16994,  -3.16994,         0,
                         0,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.59222,  -3.16994,  -3.16994,         0,
                         0,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.150774,  -1.59222,  -1.59222,         0,
                         0,  -1.35083,  -1.35083,  0.234473,         1,         1,         1,  0.234473,  -1.35083,  -1.35083,         0,
                         0,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.159444,  -1.60358,  -1.60358,         0,
                         0,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,         0,
                         0,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,         0,
                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,

                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.59222,  -3.16994,  -3.16994,  -3.16994,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.59222,  -3.16994,  -3.16994,  -3.16994,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.59222,  -3.16994,  -3.16994,  -3.16994,
                  -1.59222,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.150774,  -1.59222,  -1.59222,  -1.59222,
                  -1.35083,  -1.35083,  -1.35083,  0.234473,         1,         1,         1,  0.234473,  -1.35083,  -1.35083,  -1.35083,
                  -1.59222,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.159444,  -1.60358,  -1.60358,  -1.60358,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,  -3.30368,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,  -3.30368,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,  -3.30368,
                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,

                         0,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.59222,  -3.16994,  -3.16994,         0,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.59222,  -3.16994,  -3.16994,  -3.16994,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.59222,  -3.16994,  -3.16994,  -3.16994,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.59222,  -3.16994,  -3.16994,  -3.16994,
                  -1.59222,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.150774,  -1.59222,  -1.59222,  -1.59222,
                  -1.35083,  -1.35083,  -1.35083,  0.234473,         1,         1,         1,  0.234473,  -1.35083,  -1.35083,  -1.35083,
                  -1.59222,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.159444,  -1.60358,  -1.60358,  -1.60358,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,  -3.30368,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,  -3.30368,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,  -3.30368,
                         0,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,         0,

                         0,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.59222,  -3.16994,  -3.16994,         0,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.59222,  -3.16994,  -3.16994,  -3.16994,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.59222,  -3.16994,  -3.16994,  -3.16994,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.59222,  -3.16994,  -3.16994,  -3.16994,
                  -1.59222,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.150774,  -1.59222,  -1.59222,  -1.59222,
                  -1.35083,  -1.35083,  -1.35083,  0.234473,         1,         1,         1,  0.234473,  -1.35083,  -1.35083,  -1.35083,
                  -1.59222,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.159444,  -1.60358,  -1.60358,  -1.60358,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,  -3.30368,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,  -3.30368,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,  -3.30368,
                         0,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,         0,

                         0,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.150774,  -1.59222,  -1.59222,         0,
                  -1.59222,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.150774,  -1.59222,  -1.59222,  -1.59222,
                  -1.59222,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.150774,  -1.59222,  -1.59222,  -1.59222,
                  -1.59222,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.150774,  -1.59222,  -1.59222,  -1.59222,
                 -0.150774, -0.150774, -0.150774, -0.150774,  0.234473,  0.234473,  0.234473, -0.150774, -0.150774, -0.150774, -0.150774,
                  0.234473,  0.234473,  0.234473,  0.234473,         1,         1,         1,  0.234473,  0.234473,  0.234473,  0.234473,
                 -0.150774, -0.150774, -0.150774, -0.150774,  0.234473,  0.234473,  0.234473, -0.147381, -0.147381, -0.147381, -0.147381,
                  -1.59222,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.147381,  -1.59222,  -1.59222,  -1.59222,
                  -1.59222,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.147381,  -1.59222,  -1.59222,  -1.59222,
                  -1.59222,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.147381,  -1.59222,  -1.59222,  -1.59222,
                         0,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.147381,  -1.59222,  -1.59222,         0,

                         0,  -1.35083,  -1.35083,  0.234473,         1,         1,         1,  0.234473,  -1.35083,  -1.35083,         0,
                  -1.35083,  -1.35083,  -1.35083,  0.234473,         1,         1,         1,  0.234473,  -1.35083,  -1.35083,  -1.35083,
                  -1.35083,  -1.35083,  -1.35083,  0.234473,         1,         1,         1,  0.234473,  -1.35083,  -1.35083,  -1.35083,
                  -1.35083,  -1.35083,  -1.35083,  0.234473,         1,         1,         1,  0.234473,  -1.35083,  -1.35083,  -1.35083,
                  0.234473,  0.234473,  0.234473,  0.234473,         1,         1,         1,  0.234473,  0.234473,  0.234473,  0.234473,
                         1,         1,         1,         1,         1,         1,         1,         1,         1,         1,         1,
                  0.234473,  0.234473,  0.234473,  0.234473,         1,         1,         1,  0.234473,  0.234473,  0.234473,  0.234473,
                  -1.35083,  -1.35083,  -1.35083,  0.234473,         1,         1,         1,  0.234473,  -1.35083,  -1.35083,  -1.35083,
                  -1.35083,  -1.35083,  -1.35083,  0.234473,         1,         1,         1,  0.234473,  -1.35083,  -1.35083,  -1.35083,
                  -1.35083,  -1.35083,  -1.35083,  0.234473,         1,         1,         1,  0.234473,  -1.35083,  -1.35083,  -1.35083,
                         0,  -1.35083,  -1.35083,  0.234473,         1,         1,         1,  0.234473,  -1.35083,  -1.35083,         0,

                         0,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.159444,  -1.60358,  -1.60358,         0,
                  -1.59222,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.159444,  -1.60358,  -1.60358,  -1.60358,
                  -1.59222,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.159444,  -1.60358,  -1.60358,  -1.60358,
                  -1.59222,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.159444,  -1.60358,  -1.60358,  -1.60358,
                 -0.150774, -0.150774, -0.150774, -0.150774,  0.234473,  0.234473,  0.234473, -0.147381, -0.147381, -0.147381, -0.147381,
                  0.234473,  0.234473,  0.234473,  0.234473,         1,         1,         1,  0.234473,  0.234473,  0.234473,  0.234473,
                 -0.159444, -0.159444, -0.159444, -0.159444,  0.234473,  0.234473,  0.234473, -0.150546, -0.150546, -0.150546, -0.150546,
                  -1.60358,  -1.60358,  -1.60358, -0.147381,  0.234473,  0.234473,  0.234473, -0.150546,  -1.68333,  -1.68333,  -1.68333,
                  -1.60358,  -1.60358,  -1.60358, -0.147381,  0.234473,  0.234473,  0.234473, -0.150546,  -1.68333,  -1.68333,  -1.68333,
                  -1.60358,  -1.60358,  -1.60358, -0.147381,  0.234473,  0.234473,  0.234473, -0.150546,  -1.68333,  -1.68333,  -1.68333,
                         0,  -1.60358,  -1.60358, -0.147381,  0.234473,  0.234473,  0.234473, -0.150546,  -1.68333,  -1.68333,         0,

                         0,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,         0,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,  -3.30368,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,  -3.30368,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,  -3.30368,
                  -1.59222,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.147381,  -1.59222,  -1.59222,  -1.59222,
                  -1.35083,  -1.35083,  -1.35083,  0.234473,         1,         1,         1,  0.234473,  -1.35083,  -1.35083,  -1.35083,
                  -1.60358,  -1.60358,  -1.60358, -0.147381,  0.234473,  0.234473,  0.234473, -0.150546,  -1.68333,  -1.68333,  -1.68333,
                  -3.30368,  -3.30368,  -3.30368,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.68333,  -3.73411,  -3.73411,  -3.73411,
                  -3.30368,  -3.30368,  -3.30368,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.68333,  -3.73411,  -3.73411,  -3.73411,
                  -3.30368,  -3.30368,  -3.30368,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.68333,  -3.73411,  -3.73411,  -3.73411,
                         0,  -3.30368,  -3.30368,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.68333,  -3.73411,  -3.73411,         0,

                         0,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,         0,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,  -3.30368,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,  -3.30368,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,  -3.30368,
                  -1.59222,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.147381,  -1.59222,  -1.59222,  -1.59222,
                  -1.35083,  -1.35083,  -1.35083,  0.234473,         1,         1,         1,  0.234473,  -1.35083,  -1.35083,  -1.35083,
                  -1.60358,  -1.60358,  -1.60358, -0.147381,  0.234473,  0.234473,  0.234473, -0.150546,  -1.68333,  -1.68333,  -1.68333,
                  -3.30368,  -3.30368,  -3.30368,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.68333,  -3.73411,  -3.73411,  -3.73411,
                  -3.30368,  -3.30368,  -3.30368,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.68333,  -3.73411,  -3.73411,  -3.73411,
                  -3.30368,  -3.30368,  -3.30368,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.68333,  -3.73411,  -3.73411,  -3.73411,
                         0,  -3.30368,  -3.30368,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.68333,  -3.73411,  -3.73411,         0,

                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,  -3.30368,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,  -3.30368,
                  -3.16994,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,  -3.30368,
                  -1.59222,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.147381,  -1.59222,  -1.59222,  -1.59222,
                  -1.35083,  -1.35083,  -1.35083,  0.234473,         1,         1,         1,  0.234473,  -1.35083,  -1.35083,  -1.35083,
                  -1.60358,  -1.60358,  -1.60358, -0.147381,  0.234473,  0.234473,  0.234473, -0.150546,  -1.68333,  -1.68333,  -1.68333,
                  -3.30368,  -3.30368,  -3.30368,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.68333,  -3.73411,  -3.73411,  -3.73411,
                  -3.30368,  -3.30368,  -3.30368,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.68333,  -3.73411,  -3.73411,  -3.73411,
                  -3.30368,  -3.30368,  -3.30368,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.68333,  -3.73411,  -3.73411,  -3.73411,
                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,

                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                         0,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,         0,
                         0,  -3.16994,  -3.16994,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.60358,  -3.30368,  -3.30368,         0,
                         0,  -1.59222,  -1.59222, -0.150774,  0.234473,  0.234473,  0.234473, -0.147381,  -1.59222,  -1.59222,         0,
                         0,  -1.35083,  -1.35083,  0.234473,         1,         1,         1,  0.234473,  -1.35083,  -1.35083,         0,
                         0,  -1.60358,  -1.60358, -0.147381,  0.234473,  0.234473,  0.234473, -0.150546,  -1.68333,  -1.68333,         0,
                         0,  -3.30368,  -3.30368,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.68333,  -3.73411,  -3.73411,         0,
                         0,  -3.30368,  -3.30368,  -1.59222,  -1.35083,  -1.35083,  -1.35083,  -1.68333,  -3.73411,  -3.73411,         0,
                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
                         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,

                };
                
                for ( std::size_t i = 0; i < speed.size(); ++i ) {
                        BOOST_CHECK_CLOSE(speed[i], speed_answers[i], epsilon);
                }
                
//                display_buffer<double, double>(size, speed);

        }




        void test_register_to_narrow_band_2d()
        {
                Parameters params;
                params.wband_ = 3;
                params.constant_speed_ = 1;
                params.gain_ = 2;
                params.wreset_ = 1;
                
                SpaceSize<TwoDimension> size {11, 11};
                LevelSetMethod<TwoDimension> lsm {params, size};
                
                auto& statuses = LevelSetMethodTester::get_statuses(lsm);
                statuses.resize(size.total_, Status::Farway);
                
                statuses[0] = Status::Band;
                statuses[1] = Status::ResetBand;
                statuses[2] = Status::Front;
                
                LevelSetMethodTester::register_to_narrow_band(lsm);
                
                const auto& narrow_band = LevelSetMethodTester::get_narrow_band(lsm);
                BOOST_CHECK(narrow_band.size() == 3);
                BOOST_CHECK(narrow_band[0] == IntPoint2d({{0, 0}}));
                BOOST_CHECK(narrow_band[1] == IntPoint2d({{1, 0}}));
                BOOST_CHECK(narrow_band[2] == IntPoint2d({{2, 0}}));
                
        }

        void test_register_to_narrow_band_3d()
        {
                Parameters params;
                params.wband_ = 3;
                params.constant_speed_ = 1;
                params.gain_ = 2;
                params.wreset_ = 1;
                
                SpaceSize<ThreeDimension> size {11, 11, 11};
                LevelSetMethod<ThreeDimension> lsm {params, size};
                auto& statuses = LevelSetMethodTester::get_statuses(lsm);
                statuses.resize(size.total_, Status::Farway);
                
                statuses[0] = Status::Band;
                statuses[1] = Status::ResetBand;
                statuses[2] = Status::Front;
                LevelSetMethodTester::register_to_narrow_band(lsm);
                
                const auto& narrow_band = LevelSetMethodTester::get_narrow_band(lsm);
                BOOST_CHECK(narrow_band.size() == 3);
                BOOST_CHECK(narrow_band[0] == IntPoint3d({{0, 0, 0}}));
                BOOST_CHECK(narrow_band[1] == IntPoint3d({{1, 0, 0}}));
                BOOST_CHECK(narrow_band[2] == IntPoint3d({{2, 0, 0}}));
        }

        void test_propagate_front_2d()
        {
                Parameters params;
                params.wband_ = 3;
                params.constant_speed_ = 1;
                params.gain_ = 2;
                params.wreset_ = 1;
                params.time_step_ = 1;
                
                SpaceSize<TwoDimension> size {3, 3};
                LevelSetMethod<TwoDimension> lsm {params, size};

                auto& phi = LevelSetMethodTester::get_phi(lsm);
                auto& dphi = LevelSetMethodTester::get_dphi(lsm);
                auto& speed = LevelSetMethodTester::get_speed(lsm);
                auto& narrow_band = LevelSetMethodTester::get_narrow_band(lsm);

                phi.resize(size.total_, 0);
                dphi.resize(size.total_);
                speed.resize(size.total_);
                
                narrow_band.push_back({{1, 1}});
                speed[4] = 3; // positive
                const std::vector<double> sphi {
                        0, 3, 0,
                        4, 2, 6,
                        0, 5, 0,
                };
                std::copy(sphi.begin(), sphi.end(), phi.begin());
                lsm.propagate_front();
                BOOST_CHECK(phi[4] == 2.0);

                std::copy(sphi.begin(), sphi.end(), phi.begin());
                speed[4] = -3; // std::sqrt(91.0)
                lsm.propagate_front();
                BOOST_CHECK(phi[4] == 2 + 3 * std::sqrt(30.0));
                
        }
 
        void test_propagate_front_3d()
        {
                Parameters params;
                params.wband_ = 3;
                params.constant_speed_ = 1;
                params.gain_ = 2;
                params.wreset_ = 1;
                params.time_step_ = 1;
                
                SpaceSize<ThreeDimension> size {3, 3, 3};
                LevelSetMethod<ThreeDimension> lsm {params, size};

                auto& phi = LevelSetMethodTester::get_phi(lsm);
                auto& dphi = LevelSetMethodTester::get_dphi(lsm);
                auto& speed = LevelSetMethodTester::get_speed(lsm);
                auto& narrow_band = LevelSetMethodTester::get_narrow_band(lsm);

                phi.resize(size.total_, 0);
                dphi.resize(size.total_);
                speed.resize(size.total_);
                
                narrow_band.push_back({{1, 1, 1}});
                speed[13] = 3; // positive
                const std::vector<double> sphi {
                        0, 0, 0,
                        0, 7, 0,
                        0, 0, 0,

                        0, 3, 0,
                        4, 2, 6,
                        0, 5, 0,

                        0, 0, 0,
                        0, 8, 0,
                        0, 0, 0,
                };


                std::copy(sphi.begin(), sphi.end(), phi.begin());
                lsm.propagate_front();
                BOOST_CHECK(phi[13] == 2.0);

                std::copy(sphi.begin(), sphi.end(), phi.begin());
                speed[13] = -3; // std::sqrt(91.0)
                lsm.propagate_front();
                BOOST_CHECK(phi[13] == 2 + 3 * std::sqrt(91.0));
        }


        void test_2d()
        {
                std::cout << "2d\n";
                test_front_2d();
                test_initialize_along_front_2d();
                test_initialize_over_all_2d();
                test_create_speed_grid_2d();
                test_set_speed_on_front_2d();
                test_clear_speed_within_narrow_band_2d();
                test_set_speed_function_2d();
                test_copy_nearest_speed_to_narrow_band_2d();
                
                test_register_to_narrow_band_2d();
                test_propagate_front_2d();
        }

        void test_3d()
        {
                std::cout << "3d\n";
                test_front_3d();
                test_initialize_along_front_3d();
                test_initialize_over_all_3d();
                test_create_speed_grid_3d();
                test_set_speed_on_front_3d();
                test_clear_speed_within_narrow_band_3d();
                test_set_speed_function_3d();
                test_copy_nearest_speed_to_narrow_band_3d();
                test_register_to_narrow_band_3d();
                test_propagate_front_3d();
        }
}

BOOST_AUTO_TEST_CASE(TEST_LevelSetMethod)
{
        std::cout << "LevelSetMethod\n";
        test_3d();
        test_2d();
}
#endif // TEST_LevelSetMethod
