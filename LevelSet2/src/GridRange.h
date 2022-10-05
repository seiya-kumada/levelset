//
//  GridRange.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/12/21.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_GridRange_h
#define LevelSetMethod2_GridRange_h

#include <boost/iterator.hpp>
#include <boost/range.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include "SpaceSize.h"

namespace lsm
{
        template<typename D>
        struct GridRange;

        /// explicit specialization with TwoDimension
        template<>
        struct GridRange<TwoDimension>
        {
                boost::iterator_range<boost::counting_iterator<int>> x_range_;
                boost::iterator_range<boost::counting_iterator<int>> y_range_;
                
                explicit GridRange(const SpaceSize<TwoDimension>& size)
                        : x_range_{boost::make_iterator_range(boost::counting_iterator<int>(0), boost::counting_iterator<int>(size.width_))}
                        , y_range_{boost::make_iterator_range(boost::counting_iterator<int>(0), boost::counting_iterator<int>(size.height_))}
                        {}
        }; // struct GridRange<TwoDimension>

        /// explicit specialization with ThreeDimension
        template<>
        struct GridRange<ThreeDimension>
        {
                boost::iterator_range<boost::counting_iterator<int>> x_range_;
                boost::iterator_range<boost::counting_iterator<int>> y_range_;
                boost::iterator_range<boost::counting_iterator<int>> z_range_;
                
                explicit GridRange(const SpaceSize<ThreeDimension>& size)
                        : x_range_{boost::make_iterator_range(boost::counting_iterator<int>(0), boost::counting_iterator<int>(size.width_))}
                        , y_range_{boost::make_iterator_range(boost::counting_iterator<int>(0), boost::counting_iterator<int>(size.height_))}
                        , z_range_{boost::make_iterator_range(boost::counting_iterator<int>(0), boost::counting_iterator<int>(size.depth_))}
                        {}
        }; // struct GridRange<ThreeDimension>
} // namespace lsm

#endif // LevelSetMethod2_GridRange_h
