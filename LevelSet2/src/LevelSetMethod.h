//
//  LevelSetMethod.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/11/18.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_LevelSetMethod_h
#define LevelSetMethod2_LevelSetMethod_h

#include "Parameters.h"
#include "Point.h"
#include "Front.h"
#include "Indexer.h"
#include "InitialFront.h"
#include "InsideEstimator.h"
#include "Grid.h"
#include "SpeedFactor.h"
#include "UpwindScheme.h"
#include "NeighboringPoints.h"
#include "CurvatureGenerator.h"
#include "ZeroLevelSetDetector.h"
#include "DistanceMapGenerator.h"
#include "Status.h"
#include "parallel_tools.h"
#include "Debug.h"
#include "GridRange.h"
#include "StoppingCondition.h"

#include <boost/range/algorithm/fill.hpp>
#include <boost/range/algorithm/for_each.hpp>

namespace lsm
{
        /// executes the Level Set Method
        /**
         *      @tparam D the dimension type
         */
        template<typename D>
        class LevelSetMethod
        {
                friend class LevelSetMethodTester;
                BOOST_CLASS_REQUIRE(D, lsm, DimensionConcept);

        public:
                LevelSetMethod(
                        const Parameters&   parameters,
                        const SpaceSize<D>& size
                )
                        : parameters_(parameters)
                        , size_{size}
                        , indexer_{size}
                        , phi_(size.total_, 0.0)
                        , dphi_(size.total_, 0.0)
                        , speed_(size.total_, 0.0)
                        , statuses_(size.total_, Status::Farway)
                        , upwind_scheme_{indexer_, phi_}
                        , speed_factor_{indexer_, input_object_}
                        , curvature_generator_{indexer_, phi_}
                        , includes_zero_level_set_{phi_, indexer_}
                        , distance_map_generator_{parameters_.wband_, indexer_, statuses_}
                        , upper_distance_{parameters_.wband_ - parameters_.wreset_}
                        , grid_range_{size}
                        , previous_fs_{0.0}
                {
                        is_inside_space_without_edge_.set_grid(create_space_without_edge(D()));
                        is_inside_space_with_edge_.set_grid(create_space_with_edge(D()));
                }
                
                ~LevelSetMethod() = default;
                
                void initialize_narrow_band()
                {
                        narrow_band_.clear();
                        set_speed_function(true);
                }

                void calculate_speed_factors()
                {
                        speed_factor_.calculate_all(size_);
                }

                void initialize_distance_map()
                {
                        distance_map_generator_.create_distance_map();
                }
                
                const SpaceSize<D>& get_size() const
                {
                        return size_;
                }
                
                // test ok
                void initialize_along_front(const InitialFront<D>& initial_front)
                {
                        front_.clear();
                        normals_.clear();
                        create_initial_front(initial_front, D());
                        is_inside_initial_front_.set_grid(initial_front_);
                        initialize_along_front(D());
                }
                
                // test ok
                void initialize_over_all(const InitialFront<D>& initial_front)
                {
                        initialize_over_all(initial_front, D());
                }

                const std::vector<double>& get_phi() const
                {
                        return phi_;
                }
                
                const std::vector<Status>& get_statuses() const
                {
                        return statuses_;
                }
                
                const Front<D>& get_front() const
                {
                        return front_;
                }

                /// for debug
//                const Front<D>& get_non_zero_front() const
//                {
//                        return non_zero_front_;
//                }
                
                const Grid<D>& get_grid() const
                {
                        return initial_front_;
                }
                
                const Indexer<D>& get_indexer() const
                {
                        return indexer_;
                }

                const std::vector<DoublePoint<D>>& get_normals() const
                {
                        return normals_;
                }

                /// for debug
//                const std::vector<DoublePoint<D>>& get_non_zero_normals() const
//                {
//                        return non_zero_normals_;
//                }
                
                void print_verbose_description() const
                {
                        Debug::debug("[" + std::to_string(stopping_condition_.get_counter()) + "]---");
                        Debug::debug("zero speed rate = " + std::to_string(zero_count_/static_cast<double>(front_.size())));
                        Debug::debug("front size      = " + std::to_string(front_.size()));
                        Debug::debug("total speed     = " + std::to_string(total_speed_));
                }
                
                //
                bool set_speed_function(bool resets)
                {
                        clear_speed_within_narrow_band(resets);
                        total_speed_ = set_speed_on_front();
                        copy_nearest_speed_to_narrow_band(resets);
                        if ( resets ) {
                                narrow_band_.clear();
                                register_to_narrow_band(D());
                        }
                        stopping_condition_.add_total_speed(total_speed_);
                        return stopping_condition_.is_satisfied();
                }

                // test ok
                void propagate_front()
                {
                        boost::for_each(narrow_band_,
                                [&](const IntPoint<D>& p)
                                {
                                        if ( is_inside_space_without_edge_(p) ) {
                                                const int index = indexer_(p);
                                                const double& speed = speed_[index];
                                                double upwind_scheme = 0.0;
                                                if ( speed > 0.0 ) {
                                                        upwind_scheme = upwind_scheme_.template calculate<PositiveSpeed>(p, D());
                                                } else if ( speed < 0.0 ) {
                                                        upwind_scheme = upwind_scheme_.template calculate<NegativeSpeed>(p, D());
                                                } else {
                                                        upwind_scheme = 0.0;
                                                }
                                                dphi_[index] = speed * upwind_scheme * parameters_.time_step_;
                                        }
                                }
                        );

                        // a
                        boost::for_each(narrow_band_, //!!
                            [&](const IntPoint<D>& p)
                            {
                                 const int index = indexer_(p);
                                 phi_[index] -= dphi_[index];
                            }
                        );
                }
                
                void calculate_normals()
                {
                        normals_.clear();
                        boost::for_each(front_, [&](const IntPoint<D>& p)
                                {
                                        normals_.push_back(std::move(curvature_generator_.calculate_normal(p)));
                                }
                        );
                }

//                void calculate_non_zero_normals()
//                {
//                        std::vector<DoublePoint<D>>().swap(non_zero_normals_);
//                        non_zero_normals_.reserve(non_zero_front_.size());
//                        boost::for_each(non_zero_front_, [&](const IntPoint<D>& p)
//                                {
//                                        non_zero_normals_.push_back(std::move(curvature_generator_.calculate_normal(p)));
//                                }
//                        );
//                }


                bool create_labels()
                {
                        bool resets = false;
                        
                        front_.clear();
                        boost::for_each(narrow_band_,
                                [&](const IntPoint<D>& p)
                                {
                                        if ( is_inside_space_without_edge_(p) ) {
                                                const int index = indexer_(p);
                                                Status& status = statuses_[index];
                                                if ( includes_zero_level_set_(p) ) {
                                                        if ( status == Status::ResetBand ) {
                                                                resets = true;
                                                        }
                                                        status = Status::Front;
                                                        front_.push_back(p);
                                                } else {
                                                        if ( status == Status::Front ) {
                                                                status = Status::Band;
                                                        }
                                                }
                                        }
                                }
                        ); 
                        return resets;
                } 

                std::vector<std::uint8_t>& input_object()
                {
                        return input_object_;
                }

                const std::vector<std::uint8_t>& input_object() const
                {
                        return input_object_;
                }
               
        private:
                LevelSetMethod(const LevelSetMethod&) = delete;
                LevelSetMethod& operator=(const LevelSetMethod&) = delete;
                
                LevelSetMethod(LevelSetMethod&&) = delete;
                LevelSetMethod& operator=(LevelSetMethod&&) = delete;
                
                // test ok
                void register_to_narrow_band(TwoDimension)
                {
                        boost::for_each(grid_range_.y_range_, [&](int j) 
                                {
                                        boost::for_each(grid_range_.x_range_, [&](int i)
                                                {
                                                        register_to_narrow_band({{i, j}});
                                                }
                                        );
                                }
                        );
                }

                // test ok
                void register_to_narrow_band(ThreeDimension)
                {
                        boost::for_each(grid_range_.z_range_, [&](int k)
                                {
                                        boost::for_each(grid_range_.y_range_, [&](int j)
                                                {
                                                        boost::for_each(grid_range_.x_range_, [&](int i)
                                                                {
                                                                        register_to_narrow_band({{i, j, k}});
                                                                }
                                                        );
                                                }
                                        );
                                }
                        );
                }
                
                // test ok
                inline void register_to_narrow_band(const IntPoint<D>& p)
                {
                        int index = indexer_(p);
                        if ( statuses_[index] != Status::Farway ) {
                                narrow_band_.push_back(std::move(p));
                        }
                }

                const Grid<D> create_space_with_edge(TwoDimension) const
                {
                        Grid<D> grid;
                        grid.left_ = -1;
                        grid.right_ = size_.width_;
                        grid.top_ = -1;
                        grid.bottom_ = size_.height_;
                        return std::move(grid);
                }
                
                const Grid<D> create_space_with_edge(ThreeDimension) const
                {
                        Grid<D> grid;
                        grid.left_ = -1;
                        grid.right_ = size_.width_;
                        grid.top_ = -1;
                        grid.bottom_ = size_.height_;
                        grid.front_ = -1;
                        grid.back_ = size_.depth_;
                        return std::move(grid);
                }

                template<typename Pair>
                inline void copy_speed_to_narrow_band(
                        const std::vector<bool>& is_considerable,
                        const Pair&              pair,
                        const IntPoint<D>&       center,
                        bool                     resets,
                        const double&            distance,
                        const double&            center_speed
                ) {
                        typedef typename DistanceMapGenerator<D>::DistanceMap::value_type value_type;
                        
                        boost::for_each(boost::make_iterator_range(pair),
                        
                                [&](const value_type& v)
                                {
                                        const auto& info = v.second;
                                        if ( is_considerable[info.label_] ) {
                                                const IntPoint<D> p = center + info.point_;
                                                if ( is_inside_space_with_edge_(p) ) {
                                                        const int index = indexer_(p);
                                                        if ( statuses_[index] != Status::Front ) {
                                                                if ( resets ) {
                                                                        statuses_[index] = (distance > upper_distance_) ? Status::ResetBand : Status::Band;
                                                                        phi_[index] = (phi_[index] < 0) ? -distance : distance;
                                                                }
                                                                speed_[index] = center_speed;
                                                        }
                                                }
                                        }
                                }
                        
                        );
                }

                // test ok
                /**
                 *      @note is_considerable: This flag decides whether we assign the speed to a given point or not.  
                 */
                void copy_nearest_speed_to_narrow_band(bool resets)
                {
                        const auto& distance_map = distance_map_generator_.get_distance_map();
                        
                        auto beg = distance_map.rbegin();
                        const auto end = distance_map.rend();
                        std::vector<std::vector<bool>> is_considerable;
                        is_considerable.reserve(front_.size());

                        boost::for_each(front_,
                                [&](const IntPoint<D>& p)
                                {
                                        is_considerable.push_back(std::move(distance_map_generator_.select_labels(p)));
                                }
                        );
                        
                        while ( beg != end ) {
                                const double& distance = beg->first;
                                const auto range = distance_map.equal_range(distance);
                                int k = 0;
                                boost::for_each(front_, [&](const IntPoint<D>& p)
                                        {
                                                const int index = indexer_(p);
                                                if ( resets ) {
                                                        phi_[index] = 0.0;
                                                }
                                                copy_speed_to_narrow_band(is_considerable[k], range, p, resets, distance, speed_[index]);
                                                ++k;
                                        }
                                );
                                std::advance(beg, std::distance(range.first, range.second));
                        }
                }

                // test ok
                void initialize_over_all(const IntPoint<D>& p)
                {
                        const int index = indexer_(p);
                        if ( statuses_[index] != Status::Front ) {
                                phi_[index] = is_inside_initial_front_(p) ? -parameters_.wband_ : parameters_.wband_;
                        }
                }
               
                // test ok
                void initialize_over_all(const InitialFront<D>& initial_front, TwoDimension)
                {
                        boost::for_each(grid_range_.y_range_, [&](int j)
                                {
                                        boost::for_each(grid_range_.x_range_, [&](int i)
                                                {
                                                        initialize_over_all(IntPoint<D>({{i, j}}));
                                                }
                                        );
                                }
                        );
                }
               
                // test ok
                void initialize_over_all(const InitialFront<D>& initial_front, ThreeDimension)
                {
                        boost::for_each(grid_range_.z_range_, [&](int k)
                                {
                                        boost::for_each(grid_range_.y_range_, [&](int j)
                                                {
                                                        boost::for_each(grid_range_.x_range_, [&](int i)
                                                                {
                                                                        initialize_over_all(IntPoint<D>({{i, j, k}}));
                                                                }
                                                        );
                                                }
                                        );
                                }
                        );
                }

                // test ok
                inline void initialize_point_on_front(const Indexer<D>& indexer, const IntPoint<D>& p)
                {
                        const int index = indexer(p);
                        phi_[index] = 0;
                        statuses_[index] = Status::Front;
                        front_.push_back(p);
                }
                
                // test ok
                void initialize_along_front(TwoDimension)
                {
                        for ( int i = initial_front_.left_; i < initial_front_.right_; ++i ) {
                                initialize_point_on_front(indexer_, {{i, initial_front_.top_}});
                        }
                        for ( int j = initial_front_.top_; j < initial_front_.bottom_; ++j ) {
                                initialize_point_on_front(indexer_, {{initial_front_.right_, j}});
                        }
                        for ( int i = initial_front_.right_; i > initial_front_.left_; --i ) {
                                initialize_point_on_front(indexer_, {{i, initial_front_.bottom_}});
                        }
                        for ( int j = initial_front_.bottom_; j > initial_front_.top_; --j ) {
                                initialize_point_on_front(indexer_, {{initial_front_.left_, j}});
                        }
                }

                // test ok
                void initialize_along_front(ThreeDimension)
                {
                        for ( int j = initial_front_.top_; j <= initial_front_.bottom_; ++j ) {
                                for ( int i = initial_front_.left_; i <= initial_front_.right_; ++i ) {
                                        initialize_point_on_front(indexer_, {{i, j, initial_front_.front_}});
                                        initialize_point_on_front(indexer_, {{i, j, initial_front_.back_}});
                                }
                        }

                        for ( int k = initial_front_.front_ + 1; k < initial_front_.back_; ++k ) {
                                for ( int i = initial_front_.left_; i <= initial_front_.right_; ++i ) {
                                        initialize_point_on_front(indexer_, {{i, initial_front_.top_, k}});
                                        initialize_point_on_front(indexer_, {{i, initial_front_.bottom_, k}});
                                }
                        }

                        for ( int j = initial_front_.top_ + 1; j < initial_front_.bottom_; ++j ) {
                                for ( int k = initial_front_.front_ + 1; k < initial_front_.back_; ++k ) {
                                        initialize_point_on_front(indexer_, {{initial_front_.left_, j, k}});
                                        initialize_point_on_front(indexer_, {{initial_front_.right_, j, k}});
                                }
                        }
                }
                
                // test ok
                void create_initial_front(const InitialFront<D>& initial_front, TwoDimension)
                {
                        initial_front_.left_ = initial_front.vertices_[0][0];
                        initial_front_.top_ = initial_front.vertices_[0][1];
                        initial_front_.right_ = initial_front.vertices_[1][0];
                        initial_front_.bottom_ = initial_front.vertices_[1][1];
                }

                // test ok
                void create_initial_front(const InitialFront<D>& initial_front, ThreeDimension)
                {
                        create_initial_front(initial_front, TwoDimension());
                        initial_front_.front_ = initial_front.vertices_[0][2];
                        initial_front_.back_ = initial_front.vertices_[1][2];
                }
                
                // test ok
                void clear_speed_within_narrow_band(bool resets)
                {
                        // b
                        boost::for_each(narrow_band_, // ok
                                [&](const IntPoint<D>& p)
                                {
                                        const int index = indexer_(p);
                                        speed_[index] = 0.0;
                                        dphi_[index] = 0.0;
                                        if ( resets ) {
                                                if ( statuses_[index] != Status::Front ) {
                                                        statuses_[index] = Status::Farway;
                                                }
                                        }
                                }
                        );
                }
               
                // test ok
                const Grid<D> create_space_without_edge(TwoDimension) const
                {
                        Grid<D> grid;
                        grid.left_ = 0;
                        grid.right_ = size_.width_ - 1;
                        grid.top_ = 0;
                        grid.bottom_ = size_.height_ - 1;
                        return grid;
                }
               
                // test ok
                const Grid<D> create_space_without_edge(ThreeDimension) const
                {
                        Grid<D> grid;
                        grid.left_ = 0;
                        grid.right_ = size_.width_ - 1;
                        grid.top_ = 0;
                        grid.bottom_ = size_.height_ - 1;
                        grid.front_ = 0;
                        grid.back_ = size_.depth_ - 1;
                        return grid;
                }
               
                // test ok
                double set_speed_on_front()
                {
                        double fs = 0.0;
                        
                        // for debug
                        zero_count_ = 0;
//                        std::vector<IntPoint<D>>().swap(non_zero_front_);
//                        non_zero_front_.reserve(front_.size());
                        boost::for_each(front_,
                                [&](const IntPoint<D>& p) 
                                {
                                        if ( is_inside_space_without_edge_(p) ) {
                                                const int i = indexer_(p);
                                                double& speed = speed_[i];
                                                const double speed_factor = speed_factor_.value(p);
                                                speed = speed_factor;
                                                const double kappa = curvature_generator_(p);
                                                speed *= (parameters_.constant_speed_ - parameters_.gain_ * kappa);
                                                if ( std::abs(speed) < parameters_.speed_threshold_ ) {
                                                        speed = 0.0;
                                                        ++zero_count_;
                                                } else {
//                                                        non_zero_front_.push_back(p);
                                                }
                                                fs += std::abs(speed);
                                        }
                                }
                        );
                        return fs;
                }

                /// input parameters
                Parameters parameters_;
                
                /// size of the input image
                SpaceSize<D> size_;
                
                /// accessor of the array
                Indexer<D> indexer_;
                
                /// input front(zero-level set)
                Grid<D> initial_front_;
                
                /// auxiliary function
                std::vector<double> phi_;
                
                /// deviation of auxiliary function
                std::vector<double> dphi_;
                
                /// veclocity function
                std::vector<double> speed_;

                /// current statuses
                std::vector<Status> statuses_;
                
                /// front
                Front<D> front_;
                
                /// normals
                std::vector<DoublePoint<D>> normals_;
                
                /// narrow band
                std::vector<IntPoint<D>> narrow_band_;
                
                /// input image(gray image)
                std::vector<std::uint8_t> input_object_;

                UpwindScheme<D> upwind_scheme_;
                
                SpeedFactor<D> speed_factor_;
                
                CurvatureGenerator<D> curvature_generator_;
                
                InsideEstimator<D> is_inside_space_without_edge_;
                InsideEstimator<D> is_inside_space_with_edge_;
                InsideEstimator<D> is_inside_initial_front_;
                
                ZeroLevelSetDetector<D> includes_zero_level_set_;
                DistanceMapGenerator<D> distance_map_generator_;
                
                const int upper_distance_;
                
                GridRange<D> grid_range_;
                
                double previous_fs_;
                
                /// for debug
                int zero_count_ {0};
//                Front<D> non_zero_front_;
//                std::vector<DoublePoint<D>> non_zero_normals_;
                
                StoppingCondition stopping_condition_;
                
                double total_speed_ {0.0};
                
        }; // class LevelSetMethod
} // namespace lsm

#endif // LevelSetMethod2_LevelSetMethod_h
