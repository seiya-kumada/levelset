//
//  DistanceMapGenerator.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/12/12.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_DistanceMapGenerator_h
#define LevelSetMethod2_DistanceMapGenerator_h

#include <map>
#include "Concepts.h"
#include "NeighboringPoints.h"
#include "Status.h"
#include "Indexer.h"
#include <boost/bimap/bimap.hpp>

namespace lsm
{
        namespace detail
        {
                inline constexpr bool is_zero(int x)
                {
                        return x == 0;
                }

                inline constexpr bool is_negative(int x)
                {
                        return x < 0;
                }
                
                inline constexpr int to_symbol(int x)
                {
                        return is_zero(x) ? 0 :
                               is_negative(x) ? -1 :
                               1;
                }
                
                /// creates a table that assigns a number to each point
                /**
                 *      @tparam D the dimension type
                 *      @see the file "grid.xls"
                 */
                template<typename D>
                struct Table
                {
                        BOOST_CLASS_REQUIRE(D, lsm, DimensionConcept);
                        typedef typename boost::bimaps::bimap<IntPoint<D>, int>::value_type value_type;
                        
                        boost::bimaps::bimap<IntPoint<D>, int>           indices_;
                        std::array<IntPoint<D>, power(3, D::Dimension_)> points_;
                        
                        Table()
                        {
                                initialize();
                                for ( int i = 0; i < static_cast<int>(points_.size()); ++i ) {
                                        points_[i] = indices_.right.at(i);
                                }
                        }
                       
                        void initialize();
                        const IntPoint<D> to_symbols(const IntPoint<D>& p) const;

                        int index(const IntPoint<D>& p) const
                        {
                                return indices_.left.at(to_symbols(p));
                        }
                        
                        const IntPoint<D>& point(int index) const
                        {
                                return points_[index];
                        }
                }; // struct Table
               
                /// explicit specialization with TwoDimension
                template<>
                inline void Table<TwoDimension>::initialize()
                {
                        indices_.insert(value_type(NeighboringPoints2d( 0,  0), 0));
                        indices_.insert(value_type(NeighboringPoints2d( 0, -1), 1));
                        indices_.insert(value_type(NeighboringPoints2d( 0,  1), 2));
                        indices_.insert(value_type(NeighboringPoints2d(-1,  0), 3));
                        indices_.insert(value_type(NeighboringPoints2d(-1, -1), 4));
                        indices_.insert(value_type(NeighboringPoints2d(-1,  1), 5));
                        indices_.insert(value_type(NeighboringPoints2d( 1,  0), 6));
                        indices_.insert(value_type(NeighboringPoints2d( 1, -1), 7));
                        indices_.insert(value_type(NeighboringPoints2d( 1,  1), 8));
                }
                
                /// explicit specialization with TwoDimension
                template<>
                inline const IntPoint<TwoDimension> 
                Table<TwoDimension>::to_symbols(const IntPoint2d& p) const
                {
                        return {{to_symbol(p[0]), to_symbol(p[1])}};
                }

                /// explicit specialization with ThreeDimension
                template<>
                inline void Table<ThreeDimension>::initialize()
                {
                        indices_.insert(value_type(NeighboringPoints3d( 0,  0,  0),  0));
                        indices_.insert(value_type(NeighboringPoints3d( 0, -1,  0),  1));
                        indices_.insert(value_type(NeighboringPoints3d( 0,  1,  0),  2));
                        
                        indices_.insert(value_type(NeighboringPoints3d( 0,  0, -1),  3));
                        indices_.insert(value_type(NeighboringPoints3d( 0, -1, -1),  4));
                        indices_.insert(value_type(NeighboringPoints3d( 0,  1, -1),  5));
                        
                        indices_.insert(value_type(NeighboringPoints3d( 0,  0,  1),  6));
                        indices_.insert(value_type(NeighboringPoints3d( 0, -1,  1),  7));
                        indices_.insert(value_type(NeighboringPoints3d( 0,  1,  1),  8));
                        
                        indices_.insert(value_type(NeighboringPoints3d(-1,  0,  0),  9));
                        indices_.insert(value_type(NeighboringPoints3d(-1, -1,  0), 10));
                        indices_.insert(value_type(NeighboringPoints3d(-1,  1,  0), 11));
                        
                        indices_.insert(value_type(NeighboringPoints3d(-1,  0, -1), 12));
                        indices_.insert(value_type(NeighboringPoints3d(-1, -1, -1), 13));
                        indices_.insert(value_type(NeighboringPoints3d(-1,  1, -1), 14));
                        
                        indices_.insert(value_type(NeighboringPoints3d(-1,  0,  1), 15));
                        indices_.insert(value_type(NeighboringPoints3d(-1, -1,  1), 16));
                        indices_.insert(value_type(NeighboringPoints3d(-1,  1,  1), 17));
                        
                        indices_.insert(value_type(NeighboringPoints3d( 1,  0,  0), 18));
                        indices_.insert(value_type(NeighboringPoints3d( 1, -1,  0), 19));
                        indices_.insert(value_type(NeighboringPoints3d( 1,  1,  0), 20));
                        
                        indices_.insert(value_type(NeighboringPoints3d( 1,  0, -1), 21));
                        indices_.insert(value_type(NeighboringPoints3d( 1, -1, -1), 22));
                        indices_.insert(value_type(NeighboringPoints3d( 1,  1, -1), 23));
                        
                        indices_.insert(value_type(NeighboringPoints3d( 1,  0,  1), 24));
                        indices_.insert(value_type(NeighboringPoints3d( 1, -1,  1), 25));
                        indices_.insert(value_type(NeighboringPoints3d( 1,  1,  1), 26));
                }

                /// explicit specialization with ThreeDimension
                template<>
                inline const IntPoint<ThreeDimension> 
                Table<ThreeDimension>::to_symbols(const IntPoint<ThreeDimension>& p) const
                {
                        return {{to_symbol(p[0]), to_symbol(p[1]), to_symbol(p[2])}};
                }
        } // namespace detail

        /// calculates a reference table of distances
        /**
         *      @tparam D the dimension type
         *      @see the file "grid.xls"
         */
        template<typename D>
        class DistanceMapGenerator
        {
                BOOST_CLASS_REQUIRE(D, lsm, DimensionConcept);
                
        public:
                struct PointInfo
                {
                        IntPoint<D> point_;
                        int label_;
                };
                typedef std::multimap<double, PointInfo> DistanceMap;
                

                DistanceMapGenerator(
                        int wband,
                        const Indexer<D>&          indexer,
                        const std::vector<Status>& statuses
                )
                        : wband_{wband}
                        , squared_wband_{wband * (1 + wband)}
                        , indexer_(indexer)
                        , statuses_(statuses) 
                {
                }

                ~DistanceMapGenerator() = default;
                
                const std::multimap<double, PointInfo>& get_distance_map() const
                {
                        return distance_map_;
                }
                
                const std::vector<bool> select_labels(const IntPoint<D>& p);
                void create_distance_map();

        private:
                DistanceMapGenerator(const DistanceMapGenerator&) = delete;
                DistanceMapGenerator& operator=(const DistanceMapGenerator&) = delete;
                
                DistanceMapGenerator(DistanceMapGenerator&&) = delete;
                DistanceMapGenerator& operator=(DistanceMapGenerator&&) = delete;
               
                
                DistanceMap distance_map_;
                detail::Table<D> table_;
                int wband_;
                int squared_wband_;
                const Indexer<D>& indexer_;
                const std::vector<Status>& statuses_;
                
                bool remove(
                        const IntPoint<D>&        p,
                        int                       a,
                        const std::array<int, 3>& indices,
                        std::vector<bool>&        labels
                ) {
                        if ( statuses_[indexer_(p + table_.point(a))] == Status::Front ) {
                                labels[indices[0]] = false;
                                labels[indices[1]] = false;
                                labels[indices[2]] = false;
                                return true;
                        } else {
                                return false;
                        }
                }

                bool remove(
                        const IntPoint<D>&        p,
                        int                       a,
                        const std::array<int, 9>& indices,
                        std::vector<bool>&        labels
                ) {
                        if ( statuses_[indexer_(p + table_.point(a))] == Status::Front ) {
                                labels[indices[0]] = false;
                                labels[indices[1]] = false;
                                labels[indices[2]] = false;
                                labels[indices[3]] = false;
                                labels[indices[4]] = false;
                                labels[indices[5]] = false;
                                labels[indices[6]] = false;
                                labels[indices[7]] = false;
                                labels[indices[8]] = false;
                                return true;
                        } else {
                                return false;
                        }
                }
                
                void remove(
                        const IntPoint<D>& p,
                        int                a,
                        std::vector<bool>& labels
                ) {
                        if ( statuses_[indexer_(p + table_.point(a))] == Status::Front ) {
                                labels[a] = false;
                        }
                }
                
                void register_distance(const IntPoint<D>& p, const double& d)
                {
                        const double ds = std::sqrt(d);
                        const int index = table_.index(p);
                        distance_map_.insert(typename DistanceMap::value_type(ds, {p, index}));
                }
        }; // class DistanceMapGenerator

        /// explicit specialization with TwoDimension
        template<>
        inline const std::vector<bool>
        DistanceMapGenerator<TwoDimension>::select_labels(const IntPoint2d& p)
        {
                constexpr int size = power(3, TwoDimension::Dimension_);
                std::vector<bool> labels(size, true);
                bool flag_0 = remove(p, 1, std::array<int, 3>({{1, 4, 7}}), labels);
                bool flag_1 = remove(p, 2, std::array<int, 3>({{2, 5, 8}}), labels);
                remove(p, 3, std::array<int, 3>({{3, 4, 5}}), labels);
                remove(p, 6, std::array<int, 3>({{6, 7, 8}}), labels);
                
                if ( !flag_0 ) {
                        remove(p, 4, labels);
                        remove(p, 7, labels);
                }
                if ( !flag_1 ) {
                        remove(p, 5, labels);
                        remove(p, 8, labels);
                }
                return std::move(labels);
        }

        /// explicit specialization with ThreeDimension
        template<>
        inline const std::vector<bool>
        DistanceMapGenerator<ThreeDimension>::select_labels(const IntPoint3d& p)
        {
                constexpr int size = power(3, ThreeDimension::Dimension_);
                std::vector<bool> labels(size, true);
                
                bool flag_0 = remove(p,  1, std::array<int, 9>({{ 1,  4,  7, 10, 13, 16, 19, 22, 25}}), labels);
                bool flag_1 = remove(p,  2, std::array<int, 9>({{ 2,  5,  8, 11, 14, 17, 20, 23, 26}}), labels);
                bool flag_2 = remove(p,  3, std::array<int, 9>({{ 3,  4,  5, 12, 13, 14, 21, 22, 23}}), labels);
                
                if ( !flag_0 ) {
                        remove(p,  4, labels);
                        remove(p,  7, labels);
                        remove(p, 10, labels);
                        remove(p, 13, labels);
                        remove(p, 16, labels);
                        remove(p, 19, labels);
                        remove(p, 22, labels);
                        remove(p, 25, labels);
                }
                
                if ( !flag_1 ) {
                        remove(p,  5, labels);
                        remove(p,  8, labels);
                        remove(p, 11, labels);
                        remove(p, 14, labels);
                        remove(p, 17, labels);
                        remove(p, 20, labels);
                        remove(p, 23, labels);
                        remove(p, 26, labels);
                }
                
                if ( !flag_2 ) {
                        remove(p, 12, labels);
                        remove(p, 21, labels);
                }
                
                bool flag_3 = remove(p,  6, std::array<int, 9>({{ 6,  7,  8, 15, 16, 17, 24, 25, 26}}), labels);
                if ( !flag_3 ) {
                        remove(p, 15, labels);
                        remove(p, 24, labels);
                }
                        
                remove(p,  9, std::array<int, 9>({{ 9, 10, 11, 12, 13, 14, 15, 16, 17}}), labels);
                remove(p, 18, std::array<int, 9>({{18, 19, 20, 21, 22, 23, 24, 25, 26}}), labels);
                
                return std::move(labels);
        }

        /// explicit specialization with TwoDimension
        template<>
        inline void DistanceMapGenerator<TwoDimension>::create_distance_map()
        {
                for ( int x = -wband_; x <= wband_; ++x ) {
                        const int sx = x * x;
                        for ( int y = -wband_; y <= wband_; ++y ) {
                                const int d = sx + y * y;
                                if ( d <= squared_wband_ ) {
                                        const IntPoint2d p {{x, y}};
                                        register_distance(p, d);
                                }
                        }
                }
        }

        /// explicit specialization with ThreeDimension
        template<>
        inline void DistanceMapGenerator<ThreeDimension>::create_distance_map()
        {
                for ( int x = -wband_; x <= wband_; ++x ) {
                        const int sx = x * x;
                        for ( int y = -wband_; y <= wband_; ++y ) {
                                const int sy = y * y;
                                for ( int z = -wband_; z <= wband_; ++z ) {
                                        const int d = sx + sy + z * z;
                                        if ( d <= squared_wband_ ) {
                                                const IntPoint3d p {{x, y, z}};
                                                register_distance(p, d);
                                        }
                                }
                        }
                }
        }

} // namespace lsm

#endif // LevelSetMethod2_DistanceMapGenerator_h
