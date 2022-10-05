//
//  Point.cpp
//  LevelSetMethod2
//
//  Created by kumada on 2012/11/23.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#if(UNIT_TEST_Point)
#define BOOST_TEST_DYN_LINK
#include "Point.h"
#include <boost/test/unit_test.hpp>

using namespace lsm;

namespace
{
        // 2d +=
        void test_unitary_operator_0()
        {
                IntPoint2d a = {{1, 2}};
                const IntPoint2d b = {{3, 4}};
                a += b;
                BOOST_CHECK(b[0] == 3 && b[1] == 4);
                BOOST_CHECK(a[0] == 4 && a[1] == 6);
        }
        
        // 3d +=
        void test_unitary_operator_1()
        {
                IntPoint3d a = {{1, 2, 3}};
                const IntPoint3d b = {{4, 5, 6}};
                a += b;
                BOOST_CHECK(b[0] == 4 && b[1] == 5 && b[2] == 6);
                BOOST_CHECK(a[0] == 5 && a[1] == 7 && a[2] == 9);
        }
        
        // 2d -=
        void test_unitary_operator_2()
        {
                IntPoint2d a = {{1, 2}};
                const IntPoint2d b = {{3, 4}};
                a -= b;
                BOOST_CHECK(b[0] == 3 && b[1] == 4);
                BOOST_CHECK(a[0] == -2 && a[1] == -2);
        }
        
        // 3d -=
        void test_unitary_operator_3()
        {
                IntPoint3d a = {{1, 2, 3}};
                const IntPoint3d b = {{4, 5, 6}};
                a -= b;
                BOOST_CHECK(b[0] == 4 && b[1] == 5 && b[2] == 6);
                BOOST_CHECK(a[0] == -3 && a[1] == -3 && a[2] == -3);
        }

        // 2d +/-
        void test_unitary_operator_4()
        {
                const IntPoint2d a = {{1, 2}};
                const IntPoint2d b = {{3, 4}};
                const IntPoint2d c = a + b;
                BOOST_CHECK(c[0] == 4 && c[1] == 6);
                const IntPoint2d d = a - b;
                BOOST_CHECK(d[0] == -2 && d[1] == -2);
        }
        
        // 3d +/-
        void test_unitary_operator_5()
        {
                IntPoint3d a = {{1, 2, 3}};
                const IntPoint3d b = {{4, 5, 6}};
                const IntPoint3d c = a + b;
                BOOST_CHECK(c[0] == 5 && c[1] == 7 && c[2] == 9);
                const IntPoint3d d = a - b;
                BOOST_CHECK(d[0] == -3 && d[1] == -3 && d[2] == -3);
        }
        
        // 2d inner_product
        void test_inner_product_with_2d()
        {
                IntPoint2d a {{1, 2}};
                IntPoint2d b {{2, 3}};
                BOOST_CHECK(2 + 6 == inner_product(a, b));
                
                DoublePoint2d c {{1.2, 2.4}};
                DoublePoint2d d {{2.1, 3.5}};
                BOOST_CHECK(1.2 * 2.1 + 2.4 * 3.5 == inner_product(c, d));
                
        }
        
        // 3d inner_product
        void test_inner_product_with_3d()
        {
                IntPoint3d a {{1, 2, 3}};
                IntPoint3d b {{2, 3, 4}};
                BOOST_CHECK(1 * 2 + 2 * 3 + 3 * 4 == inner_product(a, b));

                DoublePoint3d c {{1.1, 2.2, 3.3}};
                DoublePoint3d d {{0.2, 0.3, 0.4}};
                BOOST_CHECK(1.1 * 0.2 + 2.2 * 0.3 + 3.3 * 0.4 == inner_product(c, d));
        }
        

        // 2d norm
        void test_norm_with_2d()
        {
                DoublePoint2d a {{1.2, 2.1}};
                double n = norm<double, TwoDimension>(a);
                BOOST_CHECK(std::sqrt(1.2 * 1.2 + 2.1 * 2.1) == n);
                
                IntPoint2d b {{1, 2}};
                double m = norm<int, TwoDimension>(b);
                BOOST_CHECK(std::sqrt(1 * 1 + 2 * 2) == m);
         }
        
        // 3d norm
        void test_norm_with_3d()
        {
                DoublePoint3d a {{1.2, 2.1, 3.1}};
                double n = norm<double, ThreeDimension>(a);
                BOOST_CHECK(std::sqrt(1.2 * 1.2 + 2.1 * 2.1 + 3.1 * 3.1) == n);
                
                IntPoint3d b {{1, 2, 3}};
                double m = norm<int, ThreeDimension>(b);
                BOOST_CHECK(std::sqrt(1 * 1 + 2 * 2 + 3 * 3) == m);
        }

        void test_outer_product_with_2d()
        {
                IntPoint2d a{{1, 0}};
                IntPoint2d b{{0, 1}};
                BOOST_CHECK(outer_product(a, b) == 1);
        }

        void test_outer_product_with_3d()
        {
                IntPoint3d a{{1, 0, 0}};
                IntPoint3d b{{0, 1, 0}};
                IntPoint3d c{{0, 0, 1}};
                BOOST_CHECK(outer_product(a, b) == c);
        }
}
#include <iostream>
BOOST_AUTO_TEST_CASE(TEST_Point)
{
        std::cout << "Point\n";
        test_unitary_operator_0();
        test_unitary_operator_1();
        test_unitary_operator_2();
        test_unitary_operator_3();
        test_unitary_operator_4();
        test_unitary_operator_5();
        test_inner_product_with_2d();
        test_inner_product_with_3d();
        test_norm_with_2d();
        test_norm_with_3d();
        test_outer_product_with_2d();
        test_outer_product_with_3d();
}
#endif // UNIT_TEST_Point
