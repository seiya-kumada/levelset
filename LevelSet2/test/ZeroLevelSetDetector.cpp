//
//  ZeroLevelSetDetector.cpp
//  LevelSetMethod2
//
//  Created by kumada on 2012/12/06.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#if(UNIT_TEST_ZeroLevelSetDetector)
#define BOOST_TEST_DYN_LINK

#include "ZeroLevelSetDetector.h"
#include "Parameters.h"
#include "Indexer.h"
#include <boost/test/unit_test.hpp>

using namespace lsm;

namespace
{
        void test_2d()
        {
                SpaceSize<TwoDimension> size {3, 3};
                Indexer<TwoDimension> indexer {size};
                
                std::vector<double> phi {
                        0, 3, 0,
                        4, 2, 6,
                        0, 5, 0,
                };
                
                ZeroLevelSetDetector<TwoDimension> includes_zero_level_set {phi, indexer};
                
                bool flag = includes_zero_level_set(IntPoint2d({{1, 1}}));
                BOOST_CHECK(flag == false);

                phi[1] = -10;
                flag = includes_zero_level_set(IntPoint2d({{1, 1}}));
                BOOST_CHECK(flag == true);
        }

        void test_3d()
        {
                SpaceSize<ThreeDimension> size {3, 3, 3};
                Indexer<ThreeDimension> indexer {size};

                std::vector<double> phi {
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
                
                ZeroLevelSetDetector<ThreeDimension> includes_zero_level_set {phi, indexer};

                bool flag = includes_zero_level_set(IntPoint3d({{1, 1, 1}}));
                BOOST_CHECK(flag == false);

                phi[4] = -10;
                flag = includes_zero_level_set(IntPoint3d({{1, 1, 1}}));
                BOOST_CHECK(flag == true);
        }
}
#include <iostream>

BOOST_AUTO_TEST_CASE(TEST_ZeroLevelSetDetector)
{
        std::cout << "ZeroLevelSetDetector\n";
        test_2d();
        test_3d();
        
}
#endif // UNIT_TEST_ZeroLevelSetDetector
