//
//  SpeedFactor.cpp
//  LevelSetMethod2
//
//  Created by kumada on 2012/11/28.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#if(UNIT_TEST_SpeedFactor)
#define BOOST_TEST_DYN_LINK
#include "SpeedFactor.h"
#include "SpaceSize.h"
#include <boost/test/unit_test.hpp>

using namespace lsm;

namespace
{
        void test_2d()
        {
                SpaceSize<TwoDimension> size;
                size.width_ = 3;
                size.height_ = 3;
                Indexer<TwoDimension> indexer {size};
                const std::vector<std::uint8_t> gray = {50, 100, 20, 100, 0, 200, 70, 100, 30};
                SpeedFactor<TwoDimension> factor {indexer, gray};
                factor.calculate_all(size);
                const double dx = 65.0 / 4.0;
                const double dy = 15.0 / 4.0;
                const double answer = 1.0 / (1.0 + std::sqrt(dx * dx + dy * dy));
                BOOST_CHECK(answer == factor.value({{1, 1}}));
        }
        
        void test_3d()
        {
                SpaceSize<ThreeDimension> size;
                size.width_ = 3;
                size.height_ = 3;
                size.depth_ = 3;
                Indexer<ThreeDimension> indexer{size};
                const std::vector<std::uint8_t> gray = {
                        0, 100, 0, 100,   0, 100, 0, 100, 0,
                        0, 100, 0, 100,   0, 100, 0, 100, 0,
                        0, 100, 0, 100, 100, 100, 0, 100, 0,
                };
                SpeedFactor<ThreeDimension> factor {indexer, gray};
                factor.calculate_all(size);
                const double answer = 1.0 / (1.0 + std::sqrt(100 * 100 + 100 * 100 + 100 * 100));
                BOOST_CHECK(answer == factor.value({{1, 1, 1}}));
        }

}
#include <iostream>
BOOST_AUTO_TEST_CASE(TEST_SpeedFactor)
{
        std::cout << "SpeedFactor\n";
        test_2d();
//        test_3d();
}
#endif // UNIT_TEST_SpeedFactor
