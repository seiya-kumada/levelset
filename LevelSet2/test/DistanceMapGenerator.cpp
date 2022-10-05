//
//  DistanceMapGenerator.cpp
//  LevelSetMethod2
//
//  Created by kumada on 2012/12/12.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#if(UNIT_TEST_DistanceMapGenerator)
#define BOOST_TEST_DYN_LINK

#include "DistanceMapGenerator.h"
#include <boost/test/unit_test.hpp>
#include <boost/range/algorithm/count.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <iostream>
using namespace lsm;

namespace
{
        const SpaceSize<TwoDimension> create_size(TwoDimension)
        {
                return SpaceSize<TwoDimension>{3, 3};
        }

        const SpaceSize<ThreeDimension> create_size(ThreeDimension)
        {
                return SpaceSize<ThreeDimension>{3, 3, 3};
        }
        
                template<typename D>
        struct CheckerSelector;
        
        // specilization with TwoDimension
        template<>
        struct CheckerSelector<TwoDimension>
        {
                static const int Count_ = 37;
                static const std::vector<int> Sizes_;
                static const std::vector<std::vector<DistanceMapGenerator<TwoDimension>::PointInfo>> Points_;
        };
        
        const std::vector<int> CheckerSelector<TwoDimension>::Sizes_ = {1, 4, 4, 4, 8, 4, 4, 8};
        
        const std::vector<std::vector<DistanceMapGenerator<TwoDimension>::PointInfo>> CheckerSelector<TwoDimension>::Points_ = {
                { {{{0, 0}}, 0} },
                { {{{-1,  0}}, 3}, {{{ 0, -1}}, 1}, {{{ 0,  1}}, 2}, {{{ 1, 0}}, 6} },
                { {{{-1, -1}}, 4}, {{{-1,  1}}, 5}, {{{ 1, -1}}, 7}, {{{ 1, 1}}, 8} },
                { {{{-2,  0}}, 3}, {{{ 0, -2}}, 1}, {{{ 0,  2}}, 2}, {{{ 2, 0}}, 6} },
                { {{{-2, -1}}, 4}, {{{-2,  1}}, 5}, {{{-1, -2}}, 4}, {{{-1, 2}}, 5}, {{{1, -2}}, 7}, {{{1, 2}}, 8}, {{{2, -1}}, 7}, {{{2, 1}}, 8} },
                { {{{-2, -2}}, 4}, {{{-2,  2}}, 5}, {{{ 2, -2}}, 7}, {{{ 2, 2}}, 8} },
                { {{{-3,  0}}, 3}, {{{ 0, -3}}, 1}, {{{ 0,  3}}, 2}, {{{ 3, 0}}, 6} },
                { {{{-3, -1}}, 4}, {{{-3,  1}}, 5}, {{{-1, -3}}, 4}, {{{-1, 3}}, 5}, {{{1, -3}}, 7}, {{{1, 3}}, 8}, {{{3, -1}}, 7}, {{{3, 1}}, 8} },
        };
        
        template<>
        struct CheckerSelector<ThreeDimension>
        {
                static const int Count_ = 19;
                static const std::vector<int> Sizes_;
                static const std::vector<std::vector<DistanceMapGenerator<ThreeDimension>::PointInfo>> Points_;
        };

        const std::vector<int> CheckerSelector<ThreeDimension>::Sizes_ = {1, 6, 12};
        
        const std::vector<std::vector<DistanceMapGenerator<ThreeDimension>::PointInfo>> CheckerSelector<ThreeDimension>::Points_ = {
               
                { {{{ 0,  0,  0}},  0} },
                { {{{-1,  0,  0}},  9},
                  {{{ 0, -1,  0}},  1},
                  {{{ 0,  0, -1}},  3},
                  {{{ 0,  0,  1}},  6},
                  {{{ 0,  1,  0}},  2},
                  {{{ 1,  0,  0}}, 18} },
                { {{{-1, -1,  0}}, 10},
                  {{{-1,  0, -1}}, 12},
                  {{{-1,  0,  1}}, 15},
                  {{{-1,  1,  0}}, 11},
                  {{{ 0, -1, -1}},  4},
                  {{{ 0, -1,  1}},  7},
                  {{{ 0,  1, -1}},  5},
                  {{{ 0,  1,  1}},  8},
                  {{{ 1, -1,  0}}, 19},
                  {{{ 1,  0, -1}}, 21},
                  {{{ 1,  0,  1}}, 24},
                  {{{ 1,  1,  0}}, 20} }
        };
        

        template<typename D>
        using PointInfo = typename DistanceMapGenerator<D>::PointInfo;

        inline bool operator==(const PointInfo<TwoDimension>& lhs, const PointInfo<TwoDimension>& rhs)
        {
                return lhs.label_ == rhs.label_ && lhs.point_ == rhs.point_;
        }

        inline bool operator==(const PointInfo<ThreeDimension>& lhs, const PointInfo<ThreeDimension>& rhs)
        {
                return lhs.label_ == rhs.label_ && lhs.point_ == rhs.point_;
        }
        
        template<typename D>
        void display_point(const typename DistanceMapGenerator<D>::PointInfo& p)
        {
                std::cout << "{" << p.label_ << ": ";
                for ( int i = 0; i < D::Dimension_; ++i ) {
                        std::cout << p.point_[i];
                        if ( i == D::Dimension_ - 1 ) {
                                std::cout << "}, ";
                        } else {
                                std::cout << ", ";
                        }
                }
        }

        template<typename D>
        const SpaceSize<D> get_dummy_size()
        {
                return SpaceSize<D>{3, 3};
        }
        
        template<>
        const SpaceSize<ThreeDimension> get_dummy_size()
        {
                return SpaceSize<ThreeDimension>{3, 3, 3};
        }
        

        template<typename D>
        void test_initialize_distance_map(int wband)
        {
                const SpaceSize<D> size = get_dummy_size<D>();
                std::vector<Status> statuses;
                DistanceMapGenerator<D> generator{wband, Indexer<D>{size}, statuses};
                generator.create_distance_map();
                const auto& map = generator.get_distance_map();
                
                typedef typename std::multimap<double, typename DistanceMapGenerator<D>::PointInfo>::value_type value_type;
//                std::cout << map.size() << std::endl;
                
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
                                        const PointInfo<D>& p = v.second;
                                        const PointInfo<D>& q = CheckerSelector<D>::Points_[index][k];
                                        
                                        BOOST_CHECK(q == p);
//                                        display_point<D>(p);
                                        ++k;
                                }
                        );
//                        std::cout << std::endl;
                        beg = range.second;
                        ++index;
                }
        }

        void test_select_labels_with_2d()
        {

                SpaceSize<TwoDimension> size{7, 7};
                Indexer<TwoDimension> indexer{size};
                const std::vector<Status> statuses = {
                        Status::Farway, Status::Farway, Status::Farway, Status::Farway, Status::Farway, Status::Farway, Status::Farway,
                        Status::Farway, Status::Farway, Status::Farway, Status::Farway, Status::Farway, Status::Farway, Status::Farway,
                        Status::Farway, Status::Farway, Status::Farway, Status::Front, Status::Farway, Status::Farway, Status::Farway,
                        Status::Farway, Status::Farway, Status::Farway, Status::Front, Status::Farway, Status::Farway, Status::Farway,
                        Status::Farway, Status::Farway, Status::Farway, Status::Farway, Status::Farway, Status::Farway, Status::Farway,
                        Status::Farway, Status::Farway, Status::Farway, Status::Farway, Status::Farway, Status::Farway, Status::Farway,
                        Status::Farway, Status::Farway, Status::Farway, Status::Farway, Status::Farway, Status::Farway, Status::Farway,
                };

                int wband = 3;
                DistanceMapGenerator<TwoDimension> generator{wband, indexer, statuses};
                generator.create_distance_map();
                const auto& labels = generator.select_labels({{3, 3}});
                
                BOOST_CHECK(6 == boost::count(labels, true));
                const std::vector<char> answers = {0, 2, 3, 5, 6, 8};
                for ( const auto& s : answers ) {
                        BOOST_CHECK(labels[s] == true);
                }
        }
        
        void test_select_labels_with_3d_0()
        {
                

                SpaceSize<ThreeDimension> size{3, 3, 3};
                Indexer<ThreeDimension> indexer{size};
                const std::vector<Status> statuses = {
                        Status::Farway, Status::Farway, Status::Farway,
                        Status::Farway, Status::Front, Status::Farway,
                        Status::Farway, Status::Farway, Status::Farway,

                        Status::Farway, Status::Farway, Status::Farway,
                        Status::Farway, Status::Front, Status::Farway,
                        Status::Farway, Status::Farway, Status::Farway,

                        Status::Farway, Status::Farway, Status::Farway,
                        Status::Farway, Status::Farway, Status::Farway,
                        Status::Farway, Status::Farway, Status::Farway,
                };

                int wband = 1;
                DistanceMapGenerator<ThreeDimension> generator{wband, indexer, statuses};
                generator.create_distance_map();
                
                const auto& labels = generator.select_labels({{1, 1, 1}});
                
                BOOST_CHECK(18 == boost::count(labels, true));
                const std::vector<int> answers = {0, 1, 2, 6, 7, 8, 9, 10, 11, 15, 16, 17, 18, 19, 20, 24, 25, 26};
                
                for ( std::size_t i = 0, k = 0; i < labels.size(); ++i ) {
                          if ( labels[i] ) {
                                BOOST_CHECK(i == answers[k++]);
                        }
                }
        }

        void test_select_labels_with_3d_1()
        {

                SpaceSize<ThreeDimension> size{3, 3, 3};
                Indexer<ThreeDimension> indexer{size};
                const std::vector<Status> statuses = {
                        Status::Farway, Status::Farway, Status::Farway,
                        Status::Farway, Status::Farway, Status::Farway,
                        Status::Farway, Status::Farway, Status::Farway,

                        Status::Farway, Status::Farway, Status::Farway,
                        Status::Farway, Status::Front, Status::Farway,
                        Status::Farway, Status::Farway, Status::Farway,

                        Status::Farway, Status::Farway, Status::Farway,
                        Status::Farway, Status::Farway, Status::Front,
                        Status::Farway, Status::Farway, Status::Farway,
                };

                int wband = 1;
                DistanceMapGenerator<ThreeDimension> generator{wband, indexer, statuses};
                generator.create_distance_map();


                const auto& labels = generator.select_labels({{1, 1, 1}});
                BOOST_CHECK(26 == boost::count(labels, true));
//                std::cout << labels.size() << std::endl;
                const std::vector<int> answers = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 25, 26};
                
                for ( std::size_t i = 0, k = 0; i < labels.size(); ++i ) {
                          if ( labels[i] ) {
                                BOOST_CHECK(i == answers[k++]);
                        }
                }
        }

        void test_select_labels_with_3d_2()
        {
              
                SpaceSize<ThreeDimension> size{3, 3, 3};
                Indexer<ThreeDimension> indexer{size};
                const std::vector<Status> statuses = {
                        Status::Farway, Status::Farway, Status::Farway,
                        Status::Farway, Status::Front, Status::Farway,
                        Status::Farway, Status::Farway, Status::Farway,

                        Status::Farway, Status::Farway, Status::Farway,
                        Status::Farway, Status::Front, Status::Farway,
                        Status::Farway, Status::Farway, Status::Farway,

                        Status::Farway, Status::Farway, Status::Farway,
                        Status::Farway, Status::Farway, Status::Front,
                        Status::Farway, Status::Farway, Status::Farway,
                };

                int wband = 1;
                DistanceMapGenerator<ThreeDimension> generator{wband, indexer, statuses};
                generator.create_distance_map();

                const auto& labels = generator.select_labels({{1, 1, 1}});
                BOOST_CHECK(17 == boost::count(labels, true));
//                std::cout << labels.size() << std::endl;
                const std::vector<int> answers = {0, 1, 2, 6, 7, 8, 9, 10, 11, 15, 16, 17, 18, 19, 20, 25, 26};
                
                for ( std::size_t i = 0, k = 0; i < labels.size(); ++i ) {
                          if ( labels[i] ) {
                                BOOST_CHECK(i == answers[k++]);
                        }
                }
        }
}

BOOST_AUTO_TEST_CASE(TEST_DistanceMapGenerator)
{
        std::cout << "DistanceMapGenerator\n";
        test_initialize_distance_map<TwoDimension>(3);
        test_initialize_distance_map<ThreeDimension>(1);
        test_select_labels_with_2d();
        test_select_labels_with_3d_0();
        test_select_labels_with_3d_1();
        test_select_labels_with_3d_2();
}
#endif // UNIT_TEST_DistanceMapGenerator
