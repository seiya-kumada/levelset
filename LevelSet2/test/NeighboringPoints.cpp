#if(UNIT_TEST_NeighboringPoints)
#define BOOST_TEST_DYN_LINK
#include "Point.h"
#include "NeighboringPoints.h"
#include <boost/test/unit_test.hpp>
#include <iostream>
using namespace lsm;

BOOST_AUTO_TEST_CASE(TEST_NeighboringPoints)
{
        std::cout << "NeighboringPoints\n";
        const NeighboringPoints<ThreeDimension> points3d;
        BOOST_CHECK(points3d(-1, -1,  0) == IntPoint3d({{-1, -1,  0}}));
        BOOST_CHECK(points3d( 0, -1,  0) == IntPoint3d({{ 0, -1,  0}}));
        BOOST_CHECK(points3d( 1, -1,  0) == IntPoint3d({{ 1, -1,  0}}));
        BOOST_CHECK(points3d(-1,  0,  0) == IntPoint3d({{-1,  0,  0}}));
        BOOST_CHECK(points3d( 0,  0,  0) == IntPoint3d({{ 0,  0,  0}}));
        BOOST_CHECK(points3d( 1,  0,  0) == IntPoint3d({{ 1,  0,  0}}));
        BOOST_CHECK(points3d(-1,  1,  0) == IntPoint3d({{-1,  1,  0}}));
        BOOST_CHECK(points3d( 0,  1,  0) == IntPoint3d({{ 0,  1,  0}}));
        BOOST_CHECK(points3d( 1,  1,  0) == IntPoint3d({{ 1,  1,  0}}));

        BOOST_CHECK(points3d(-1, -1,  -1) == IntPoint3d({{-1, -1,  -1}}));
        BOOST_CHECK(points3d( 0, -1,  -1) == IntPoint3d({{ 0, -1,  -1}}));
        BOOST_CHECK(points3d( 1, -1,  -1) == IntPoint3d({{ 1, -1,  -1}}));
        BOOST_CHECK(points3d(-1,  0,  -1) == IntPoint3d({{-1,  0,  -1}}));
        BOOST_CHECK(points3d( 0,  0,  -1) == IntPoint3d({{ 0,  0,  -1}}));
        BOOST_CHECK(points3d( 1,  0,  -1) == IntPoint3d({{ 1,  0,  -1}}));
        BOOST_CHECK(points3d(-1,  1,  -1) == IntPoint3d({{-1,  1,  -1}}));
        BOOST_CHECK(points3d( 0,  1,  -1) == IntPoint3d({{ 0,  1,  -1}}));
        BOOST_CHECK(points3d( 1,  1,  -1) == IntPoint3d({{ 1,  1,  -1}}));

        BOOST_CHECK(points3d(-1, -1,  1) == IntPoint3d({{-1, -1,  1}}));
        BOOST_CHECK(points3d( 0, -1,  1) == IntPoint3d({{ 0, -1,  1}}));
        BOOST_CHECK(points3d( 1, -1,  1) == IntPoint3d({{ 1, -1,  1}}));
        BOOST_CHECK(points3d(-1,  0,  1) == IntPoint3d({{-1,  0,  1}}));
        BOOST_CHECK(points3d( 0,  0,  1) == IntPoint3d({{ 0,  0,  1}}));
        BOOST_CHECK(points3d( 1,  0,  1) == IntPoint3d({{ 1,  0,  1}}));
        BOOST_CHECK(points3d(-1,  1,  1) == IntPoint3d({{-1,  1,  1}}));
        BOOST_CHECK(points3d( 0,  1,  1) == IntPoint3d({{ 0,  1,  1}}));
        BOOST_CHECK(points3d( 1,  1,  1) == IntPoint3d({{ 1,  1,  1}}));

        const NeighboringPoints<TwoDimension> points2d;
        BOOST_CHECK(points2d(-1, -1) == IntPoint2d({{-1, -1}}));
        BOOST_CHECK(points2d( 0, -1) == IntPoint2d({{ 0, -1}}));
        BOOST_CHECK(points2d( 1, -1) == IntPoint2d({{ 1, -1}}));
        BOOST_CHECK(points2d(-1,  0) == IntPoint2d({{-1,  0}}));
        BOOST_CHECK(points2d( 0,  0) == IntPoint2d({{ 0,  0}}));
        BOOST_CHECK(points2d( 1,  0) == IntPoint2d({{ 1,  0}}));
        BOOST_CHECK(points2d(-1,  1) == IntPoint2d({{-1,  1}}));
        BOOST_CHECK(points2d( 0,  1) == IntPoint2d({{ 0,  1}}));
        BOOST_CHECK(points2d( 1,  1) == IntPoint2d({{ 1,  1}}));

        
}
#endif // TEST_NeighboringPoints
