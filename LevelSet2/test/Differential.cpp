//
//  Differential.cpp
//  LevelSetMethod2
//
//  Created by kumada on 2012/11/25.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#if(UNIT_TEST_Differential)
#define BOOST_TEST_DYN_LINK

#include "Differential.h"
#include "SpaceSize.h"
#include "Indexer.h"
#include "DimensionTypes.h"
#include <boost/test/unit_test.hpp>

namespace lsm
{

template<typename D>
using DifferentialDouble = Differential<D, double>;

class DifferentialTester
{
public:
        static double sobel_x(DifferentialDouble<TwoDimension>& cg)
        {
                return cg.sobel_x();
        }

        static double sobel_x(DifferentialDouble<ThreeDimension>& cg)
        {
                return cg.sobel_x();
        }

        static double sobel_y(DifferentialDouble<TwoDimension>& cg)
        {
                return cg.sobel_y();
        }

        static double sobel_y(DifferentialDouble<ThreeDimension>& cg)
        {
                return cg.sobel_y();
        }

        static double sobel_z(DifferentialDouble<ThreeDimension>& cg)
        {
                return cg.sobel_z();
        }
        
        static double vx(DifferentialDouble<TwoDimension>& cg, int x, int y)
        {
                return cg.vx(x, y);
        }
        
        static double vx(DifferentialDouble<ThreeDimension>& cg, int x, int y, int z)
        {
                return cg.vx(x, y, z);
        }
        
        static double v(DifferentialDouble<TwoDimension>& cg, int x, int y)
        {
                return cg.v(x, y);
        }

        static double v(DifferentialDouble<ThreeDimension>& cg, int x, int y, int z)
        {
                return cg.v(x, y, z);
        }
        
        static double h1dx(DifferentialDouble<TwoDimension>& cg, int x, int y)
        {
                return cg.h1dx(x, y);
        }

        static double h1dx(DifferentialDouble<ThreeDimension>& cg, int x, int y, int z)
        {
                return cg.h1dx(x, y, z);
        }
        
        static double h1dx(int x, int y, int z)
        {
                return DifferentialDouble<ThreeDimension>::h1dx(x, y, z);
        }
        
        static void set_point(DifferentialDouble<TwoDimension>& cg, const IntPoint2d& p)
        {
                cg.set_point(p);
        }

        static void set_point(DifferentialDouble<ThreeDimension>& cg, const IntPoint3d& p)
        {
                cg.set_point(p);
        }
        
        static int h_total(DifferentialDouble<TwoDimension>& cg)
        {
              return cg.h0d_total_;
        }
      
        static double fx(DifferentialDouble<TwoDimension>& cg)
        {
                return cg.fx();
        }
        
        static double fy(DifferentialDouble<TwoDimension>& cg)
        {
                return cg.fy();
        }

        static double fx(DifferentialDouble<ThreeDimension>& cg)
        {
                return cg.fx();
        }
        
        static double fy(DifferentialDouble<ThreeDimension>& cg)
        {
                return cg.fy();
        }

        static double fz(DifferentialDouble<ThreeDimension>& cg)
        {
                return cg.fz();
        }
        
        static double fxy(DifferentialDouble<TwoDimension>& cg)
        {
                return cg.fxy();
        }

        static double fxy(DifferentialDouble<ThreeDimension>& cg)
        {
                return cg.fxy();
        }
        
        static double h3dxy(DifferentialDouble<ThreeDimension>& cg, int x, int y, int z)
        {
                return cg.h3dxy(x, y, z);
        }

        static double h3dxz(DifferentialDouble<ThreeDimension>& cg, int x, int y, int z)
        {
                return cg.h3dxz(x, y, z);
        }

        static double h3dyz(DifferentialDouble<ThreeDimension>& cg, int x, int y, int z)
        {
                return cg.h3dyz(x, y, z);
        }
        
        static double fxz(DifferentialDouble<ThreeDimension>& cg)
        {
                return cg.fxz();
        }
        
        static double fyz(DifferentialDouble<ThreeDimension>& cg)
        {
                return cg.fyz();
        }
        
        static constexpr double h2dx(int x, int y) 
        {
                return DifferentialDouble<TwoDimension>::h2dx(x, y);
        }

        static double fxx(DifferentialDouble<TwoDimension>& cg)
        {
                return cg.fxx();
        }

        static double h2dy(DifferentialDouble<TwoDimension>& cg, int x, int y) 
        {
                return cg.h2dy(x, y);
        }

        static double fyy(DifferentialDouble<TwoDimension>& cg)
        {
                return cg.fyy();
        }

        static constexpr double h2dx(int x, int y, int z) 
        {
                return DifferentialDouble<ThreeDimension>::h2dx(x, y, z);
        }

        static constexpr double h2dy(int x, int y, int z)
        {
                return DifferentialDouble<ThreeDimension>::h2dy(x, y, z);
        }

        static constexpr double h2dz(int x, int y, int z)
        {
                return DifferentialDouble<ThreeDimension>::h2dz(x, y, z);
        }

        static double fxx(DifferentialDouble<ThreeDimension>& cg)
        {
                return cg.fxx();
        }

        static double fyy(DifferentialDouble<ThreeDimension>& cg)
        {
                return cg.fyy();
        }

        static double fzz(DifferentialDouble<ThreeDimension>& cg)
        {
                return cg.fzz();
        }

        static double vxx(DifferentialDouble<ThreeDimension>& cg, int i, int j, int k)
        {
                return cg.vxx(i, j, k);
        }

        static double vzz(DifferentialDouble<ThreeDimension>& cg, int i, int j, int k)
        {
                return cg.vzz(i, j, k);
        }

};

} // namespace lsm

namespace
{
        using namespace lsm;
        
        typedef DifferentialDouble<TwoDimension>  DifferentialDouble2d;
        
        void test_sobel_x_2d(const std::vector<double>& input, int expected_output)
        {
                SpaceSize<TwoDimension> size {3, 3};
                Indexer<TwoDimension> indexer {size};
                DifferentialDouble2d cg {indexer, input};
                DifferentialTester::set_point(cg, IntPoint2d{{1, 1}});
                BOOST_REQUIRE_EQUAL(expected_output, DifferentialTester::sobel_x(cg));
        }

        void test_sobel_h_total_2d(const std::vector<double>& input, int expected_output)
        {
                SpaceSize<TwoDimension> size {3, 3};
                Indexer<TwoDimension> indexer {size};
                DifferentialDouble<TwoDimension> cg {indexer, input};
                DifferentialTester::set_point(cg, IntPoint2d{{1, 1}});
                BOOST_REQUIRE_EQUAL(expected_output, DifferentialTester::h_total(cg));
        }
        
        void test_fx_2d(const std::vector<double>& input, int expected_output)
        {
                SpaceSize<TwoDimension> size {3, 3};
                Indexer<TwoDimension> indexer {size};
                DifferentialDouble<TwoDimension> cg {indexer, input};
                DifferentialTester::set_point(cg, IntPoint2d{{1, 1}});
                BOOST_REQUIRE_EQUAL(expected_output, DifferentialTester::fx(cg));
        }

        void test_sobel_y_2d(const std::vector<double>& input, int expected_output)
        {
                SpaceSize<TwoDimension> size {3, 3};
                Indexer<TwoDimension> indexer {size};
                DifferentialDouble<TwoDimension> cg {indexer, input};
                DifferentialTester::set_point(cg, IntPoint2d{{1, 1}});
                BOOST_REQUIRE_EQUAL(expected_output, DifferentialTester::sobel_y(cg));
        }

        void test_fy_2d(const std::vector<double>& input, int expected_output)
        {
                SpaceSize<TwoDimension> size {3, 3};
                Indexer<TwoDimension> indexer {size};
                DifferentialDouble<TwoDimension> cg {indexer, input};
                DifferentialTester::set_point(cg, IntPoint2d{{1, 1}});
                BOOST_REQUIRE_EQUAL(expected_output, DifferentialTester::fy(cg));
        }

        void test_h1dx_3d()
        {
                constexpr int expected_output[] = {
                        -1, 0, 1,
                        -2, 0, 2,
                        -1, 0, 1,
                        -2, 0, 2,
                        -4, 0, 4,
                        -2, 0, 2,
                        -1, 0, 1,
                        -2, 0, 2,
                        -1, 0, 1,
                };
                
                BOOST_REQUIRE_EQUAL(expected_output[ 0], DifferentialTester::h1dx(-1, -1, -1));
                BOOST_REQUIRE_EQUAL(expected_output[ 1], DifferentialTester::h1dx( 0, -1, -1));
                BOOST_REQUIRE_EQUAL(expected_output[ 2], DifferentialTester::h1dx( 1, -1, -1));
                
                BOOST_REQUIRE_EQUAL(expected_output[ 3], DifferentialTester::h1dx(-1,  0, -1));
                BOOST_REQUIRE_EQUAL(expected_output[ 4], DifferentialTester::h1dx( 0,  0, -1));
                BOOST_REQUIRE_EQUAL(expected_output[ 5], DifferentialTester::h1dx( 1,  0, -1));
                
                BOOST_REQUIRE_EQUAL(expected_output[ 6], DifferentialTester::h1dx(-1,  1, -1));
                BOOST_REQUIRE_EQUAL(expected_output[ 7], DifferentialTester::h1dx( 0,  1, -1));
                BOOST_REQUIRE_EQUAL(expected_output[ 8], DifferentialTester::h1dx( 1,  1, -1));
                
                BOOST_REQUIRE_EQUAL(expected_output[ 9], DifferentialTester::h1dx(-1, -1,  0));
                BOOST_REQUIRE_EQUAL(expected_output[10], DifferentialTester::h1dx( 0, -1,  0));
                BOOST_REQUIRE_EQUAL(expected_output[11], DifferentialTester::h1dx( 1, -1,  0));
                
                BOOST_REQUIRE_EQUAL(expected_output[12], DifferentialTester::h1dx(-1,  0,  0));
                BOOST_REQUIRE_EQUAL(expected_output[13], DifferentialTester::h1dx( 0,  0,  0));
                BOOST_REQUIRE_EQUAL(expected_output[14], DifferentialTester::h1dx( 1,  0,  0));
                
                BOOST_REQUIRE_EQUAL(expected_output[15], DifferentialTester::h1dx(-1,  1,  0));
                BOOST_REQUIRE_EQUAL(expected_output[16], DifferentialTester::h1dx( 0,  1,  0));
                BOOST_REQUIRE_EQUAL(expected_output[17], DifferentialTester::h1dx( 1,  1,  0));
                
                BOOST_REQUIRE_EQUAL(expected_output[18], DifferentialTester::h1dx(-1, -1,  1));
                BOOST_REQUIRE_EQUAL(expected_output[19], DifferentialTester::h1dx( 0, -1,  1));
                BOOST_REQUIRE_EQUAL(expected_output[20], DifferentialTester::h1dx( 1, -1,  1));
                
                BOOST_REQUIRE_EQUAL(expected_output[21], DifferentialTester::h1dx(-1,  0,  1));
                BOOST_REQUIRE_EQUAL(expected_output[22], DifferentialTester::h1dx( 0,  0,  1));
                BOOST_REQUIRE_EQUAL(expected_output[23], DifferentialTester::h1dx( 1,  0,  1));
                
                BOOST_REQUIRE_EQUAL(expected_output[24], DifferentialTester::h1dx(-1,  1,  1));
                BOOST_REQUIRE_EQUAL(expected_output[25], DifferentialTester::h1dx( 0,  1,  1));
                BOOST_REQUIRE_EQUAL(expected_output[26], DifferentialTester::h1dx( 1,  1,  1));
        }
        
        void test_sobel_x_3d(const std::vector<double>& input, int expected_output)
        {
                SpaceSize<ThreeDimension> size;
                size.width_ = 3;
                size.height_ = 3;
                size.depth_ = 3;
                Indexer<ThreeDimension> indexer {size};
                DifferentialDouble<ThreeDimension> cg {indexer, input};
                DifferentialTester::set_point(cg, IntPoint3d{{1, 1, 1}});
                BOOST_REQUIRE_EQUAL(expected_output, DifferentialTester::sobel_x(cg));
        }

        void test_sobel_y_3d(const std::vector<double>& input, int expected_output)
        {
                SpaceSize<ThreeDimension> size;
                size.width_ = 3;
                size.height_ = 3;
                size.depth_ = 3;
                Indexer<ThreeDimension> indexer {size};
                DifferentialDouble<ThreeDimension> cg {indexer, input};
                DifferentialTester::set_point(cg, IntPoint3d{{1, 1, 1}});
                BOOST_REQUIRE_EQUAL(expected_output, DifferentialTester::sobel_y(cg));
        }

        void test_sobel_z_3d(const std::vector<double>& input, int expected_output)
        {
                SpaceSize<ThreeDimension> size {3, 3, 3};
                Indexer<ThreeDimension> indexer {size};
                DifferentialDouble<ThreeDimension> cg {indexer, input};
                DifferentialTester::set_point(cg, IntPoint3d{{1, 1, 1}});
                BOOST_REQUIRE_EQUAL(expected_output, DifferentialTester::sobel_z(cg));
        }
        
        void test_fx_3d(const std::vector<double>& input, double expected_output)
        {
                SpaceSize<ThreeDimension> size {3, 3, 3};
                Indexer<ThreeDimension> indexer {size};
                DifferentialDouble<ThreeDimension> cg {indexer, input};
                DifferentialTester::set_point(cg, IntPoint3d{{1, 1, 1}});
                BOOST_REQUIRE_EQUAL(expected_output, DifferentialTester::fx(cg));
        }

        void test_fy_3d(const std::vector<double>& input, double expected_output)
        {
                SpaceSize<ThreeDimension> size {3, 3, 3};
                Indexer<ThreeDimension> indexer {size};
                DifferentialDouble<ThreeDimension> cg{indexer, input};
                DifferentialTester::set_point(cg, IntPoint3d{{1, 1, 1}});
                BOOST_REQUIRE_EQUAL(expected_output, DifferentialTester::fy(cg));
        }

        void test_fz_3d(const std::vector<double>& input, double expected_output)
        {
                SpaceSize<ThreeDimension> size {3, 3, 3};
                Indexer<ThreeDimension> indexer{size};
                DifferentialDouble<ThreeDimension> cg{indexer, input};
                DifferentialTester::set_point(cg, IntPoint3d{{1, 1, 1}});
                BOOST_REQUIRE_EQUAL(expected_output, DifferentialTester::fz(cg));
        }
        
        void test_fxy_2d(const std::vector<double>& input, double expected_output)
        {
                SpaceSize<TwoDimension> size {3, 3};
                Indexer<TwoDimension> indexer {size};
                DifferentialDouble<TwoDimension> cg{indexer, input};
                DifferentialTester::set_point(cg, IntPoint2d{{1, 1}});
                BOOST_REQUIRE_EQUAL(expected_output, DifferentialTester::fxy(cg));
        }
        
        void test_fxy_3d(const std::vector<double>& input, double expected_output)
        {
                SpaceSize<ThreeDimension> size {3, 3, 3};
                Indexer<ThreeDimension> indexer {size};
                DifferentialDouble<ThreeDimension> cg{indexer, input};
                DifferentialTester::set_point(cg, IntPoint3d{{1, 1, 1}});
                BOOST_REQUIRE_EQUAL(expected_output, DifferentialTester::fxy(cg));
        }
        
        void test_fxz_3d(const std::vector<double>& input, double expected_output)
        {
                SpaceSize<ThreeDimension> size {3, 3, 3};
                Indexer<ThreeDimension> indexer {size};
                DifferentialDouble<ThreeDimension> cg{indexer, input};
                DifferentialTester::set_point(cg, IntPoint3d{{1, 1, 1}});
                BOOST_REQUIRE_EQUAL(expected_output, DifferentialTester::fxz(cg));
        }

        void test_fyz_3d(const std::vector<double>& input, double expected_output)
        {
                SpaceSize<ThreeDimension> size {3, 3, 3};
                Indexer<ThreeDimension> indexer {size};
                
                DifferentialDouble<ThreeDimension> cg{indexer, input};
                DifferentialTester::set_point(cg, IntPoint3d{{1, 1, 1}});
                BOOST_REQUIRE_EQUAL(expected_output, DifferentialTester::fyz(cg));
        }
        
        void test_fxx_2d(const std::vector<double>& input, double expected_output)
        {
                SpaceSize<TwoDimension> size {3, 3};
                Indexer<TwoDimension> indexer {size};
                DifferentialDouble<TwoDimension> cg{indexer, input};
                DifferentialTester::set_point(cg, IntPoint2d{{1, 1}});
                BOOST_REQUIRE_EQUAL(expected_output, DifferentialTester::fxx(cg));
        }
        
        
        void test_h2dx_2d()
        {
                static_assert( 1 == DifferentialTester::h2dx(-1, -1), "");
                static_assert(-2 == DifferentialTester::h2dx( 0, -1), "");
                static_assert( 1 == DifferentialTester::h2dx( 1, -1), "");
                static_assert( 2 == DifferentialTester::h2dx(-1,  0), "");
                static_assert(-4 == DifferentialTester::h2dx( 0,  0), "");
                static_assert( 2 == DifferentialTester::h2dx( 1,  0), "");
                static_assert( 1 == DifferentialTester::h2dx(-1,  1), "");
                static_assert(-2 == DifferentialTester::h2dx( 0,  1), "");
                static_assert( 1 == DifferentialTester::h2dx( 1,  1), "");
        } 
        
        void test_fyy_2d(const std::vector<double>& input, double expected_output)
        {
                SpaceSize<TwoDimension> size {3, 3};
                Indexer<TwoDimension> indexer{size};
                const std::vector<double> phi = {1, 0, 3, 2, 2, 2, 1, 0, 3};
                DifferentialDouble<TwoDimension> cg{indexer, phi};
                DifferentialTester::set_point(cg, IntPoint2d{{1, 1}});
                BOOST_REQUIRE_EQUAL(expected_output, DifferentialTester::fyy(cg));
        } 
 
        void test_fxx_3d(const std::vector<double> input, double expected_output)
        {
                SpaceSize<ThreeDimension> size {3, 3, 3};
                Indexer<ThreeDimension> indexer{size};
                DifferentialDouble<ThreeDimension> cg{indexer, input};
                DifferentialTester::set_point(cg, IntPoint3d{{1, 1, 1}});
                BOOST_REQUIRE_EQUAL(expected_output, DifferentialTester::fxx(cg));
        }
        
        void test_h2dx_3d()
        {
                static_assert( 1 == DifferentialTester::h2dx(-1, -1, -1), "");
                static_assert(-2 == DifferentialTester::h2dx( 0, -1, -1), "");
                static_assert( 1 == DifferentialTester::h2dx( 1, -1, -1), "");
                static_assert( 2 == DifferentialTester::h2dx(-1,  0, -1), "");
                static_assert(-4 == DifferentialTester::h2dx( 0,  0, -1), "");
                static_assert( 2 == DifferentialTester::h2dx( 1,  0, -1), "");
                static_assert( 1 == DifferentialTester::h2dx(-1,  1, -1), "");
                static_assert(-2 == DifferentialTester::h2dx( 0,  1, -1), "");
                static_assert( 1 == DifferentialTester::h2dx( 1,  1, -1), "");

                static_assert( 2 == DifferentialTester::h2dx(-1, -1, 0), "");
                static_assert(-4 == DifferentialTester::h2dx( 0, -1, 0), "");
                static_assert( 2 == DifferentialTester::h2dx( 1, -1, 0), "");
                static_assert( 4 == DifferentialTester::h2dx(-1,  0, 0), "");
                static_assert(-8 == DifferentialTester::h2dx( 0,  0, 0), "");
                static_assert( 4 == DifferentialTester::h2dx( 1,  0, 0), "");
                static_assert( 2 == DifferentialTester::h2dx(-1,  1, 0), "");
                static_assert(-4 == DifferentialTester::h2dx( 0,  1, 0), "");
                static_assert( 2 == DifferentialTester::h2dx( 1,  1, 0), "");

                static_assert( 1 == DifferentialTester::h2dx(-1, -1, 1), "");
                static_assert(-2 == DifferentialTester::h2dx( 0, -1, 1), "");
                static_assert( 1 == DifferentialTester::h2dx( 1, -1, 1), "");
                static_assert( 2 == DifferentialTester::h2dx(-1,  0, 1), "");
                static_assert(-4 == DifferentialTester::h2dx( 0,  0, 1), "");
                static_assert( 2 == DifferentialTester::h2dx( 1,  0, 1), "");
                static_assert( 1 == DifferentialTester::h2dx(-1,  1, 1), "");
                static_assert(-2 == DifferentialTester::h2dx( 0,  1, 1), "");
                static_assert( 1 == DifferentialTester::h2dx( 1,  1, 1), "");


        }
        
        void test_fyy_3d(const std::vector<double>& input, double expected_output)
        {
                SpaceSize<ThreeDimension> size {3, 3, 3};
                Indexer<ThreeDimension> indexer{size};
                DifferentialDouble<ThreeDimension> cg{indexer, input};
                DifferentialTester::set_point(cg, IntPoint3d{{1, 1, 1}});
                BOOST_REQUIRE_EQUAL(expected_output, DifferentialTester::fyy(cg));
        }

        void test_fzz_3d(const std::vector<double>& input, double expected_output)
        {
                SpaceSize<ThreeDimension> size {3, 3, 3};
                Indexer<ThreeDimension> indexer {size};
                DifferentialDouble<ThreeDimension> cg {indexer, input};
                DifferentialTester::set_point(cg, IntPoint3d{{1, 1, 1}});
                BOOST_CHECK(expected_output == DifferentialTester::fzz(cg));
        }

        void test_fx_and_fy_2d_with_uint8_t(
                const std::vector<std::uint8_t>& input,
                double                           expected_fx,
                double                           expected_fy
        ) {
                SpaceSize<TwoDimension> size {3, 3};
                Indexer<TwoDimension> indexer {size};
                Differential<TwoDimension, std::uint8_t> cg{indexer, input};
                cg.set_point({{1, 1}});
                BOOST_CHECK(cg.fx() == expected_fx);
                BOOST_CHECK(cg.fy() == expected_fy);
        }

        void test_fx_and_fy_3d_with_uint8_t(
                const std::vector<std::uint8_t>& gray,
                double expected_fx,
                double expected_fy,
                double expected_fz
        ) {
                SpaceSize<ThreeDimension> size {3, 3, 3};
                Indexer<ThreeDimension> indexer {size};
                Differential<ThreeDimension, std::uint8_t> cg{indexer, gray};
                cg.set_point({{1, 1, 1}});
                BOOST_CHECK(cg.fx() == expected_fx);
                BOOST_CHECK(cg.fy() == expected_fy);
                BOOST_CHECK(cg.fz() == expected_fz);
        }

        void test_2d()
        {
                test_sobel_x_2d({
                        0, 1, 2,
                        0, 1, 2,
                        0, 1, 2
                }, 8);
                
                test_sobel_h_total_2d({
                        0, 1, 2,
                        0, 1, 2,
                        0, 1, 2
                }, 4);
                
                test_fx_2d({
                        0, 1, 2,
                        0, 1, 2,
                        0, 1, 2
                }, 1);
                
                test_sobel_y_2d({
                        1, 1, 1,
                        2, 2, 2,
                        3, 3, 3
                }, 8);
                
                test_fy_2d({
                        1, 1, 1,
                        2, 2, 2,
                        3, 3, 3
                }, 1);
                
                test_fxy_2d({
                        0, 1, 2,
                        0, 1, 2,
                        0, 1, 3
                }, 0.25);
                
                test_h2dx_2d();
                
                test_fxx_2d({
                        1, 0, 3,
                        4, 0, 6,
                        7, 0, 9
                }, 40.0 / 4);
                
                test_fyy_2d({
                        1, 0, 3,
                        2, 2, 2,
                        1, 0, 3
                }, -8.0 / 4.0);
                
                test_fx_and_fy_2d_with_uint8_t({
                         50, 100,  20,
                        100,   0, 200,
                         70, 100,  30
                }, 65.0 / 4.0, 15.0 / 4.0);
        }

        void test_3d()
        {
                test_h1dx_3d();
                
                test_sobel_x_3d({
                        -1, 1, 1, /**/ -1, 1, 1, /**/ -1, 1, 1,
                        -1, 1, 1, /**/ -1, 1, 1, /**/ -1, 1, 1,
                        -1, 1, 1, /**/ -1, 1, 1, /**/ -1, 1, 1,
                }, 32);
                
                test_sobel_y_3d({
                        -1, -1, -1, /**/ 1, 1, 1, /**/ 1, 1, 1,
                        -1, -1, -1, /**/ 1, 1, 1, /**/ 1, 1, 1,
                        -1, -1, -1, /**/ 1, 1, 1, /**/ 1, 1, 1,
                }, 32);
                
                test_sobel_z_3d({
                        -1, -1, -1, /**/ -1, -1, -1, /**/ -1, -1, -1,
                         1,  1,  1, /**/  1,  1,  1, /**/  1,  1,  1,
                         1,  1,  1, /**/  1,  1,  1, /**/  1,  1,  1,
                }, 32);

                test_fx_3d({
                        0, 0, 1, /**/ 0, 0, 1, /**/ 0, 0, 1,
                        0, 0, 1, /**/ 0, 0, 1, /**/ 0, 0, 1,
                        0, 0, 1, /**/ 0, 0, 1, /**/ 0, 0, 1,
                }, 0.5);
                
                test_fy_3d({
                        0, 0, 0, /**/ 0, 0, 0, /**/ 1, 1, 1,
                        0, 0, 0, /**/ 0, 0, 0, /**/ 1, 1, 1,
                        0, 0, 0, /**/ 0, 0, 0, /**/ 1, 1, 1,
                }, 0.5);
                
                test_fz_3d({
                        0, 0, 0, /**/ 0, 0, 0, /**/ 0, 0, 0,
                        0, 0, 0, /**/ 0, 0, 0, /**/ 0, 0, 0,
                        1, 1, 1, /**/ 1, 1, 1, /**/ 1, 1, 1,
                }, 0.5);
                
                test_fxy_3d({
                        5, 0, 7, /**/ 0, 0, 0, /**/ 4, 0, 3,
                        6, 0, 8, /**/ 0, 0, 0, /**/ 5, 0, 4,
                        7, 0, 9, /**/ 0, 0, 0, /**/ 6, 0, 5,
                }, -12 / 16.0);
                
                test_fxz_3d({
                        5, 0, 7, /**/ 6, 0, 8, /**/ 7, 0, 9,
                        0, 0, 0, /**/ 0, 0, 0, /**/ 0, 0, 0,
                        4, 0, 3, /**/ 5, 0, 4, /**/ 6, 0, 5,
                }, -12 / 16.0);
                
                test_fyz_3d({
                        5, 6, 7, /**/ 0, 0, 0, /**/ 4, 5, 6,
                        0, 0, 0, /**/ 0, 0, 0, /**/ 0, 0, 0,
                        7, 8, 9, /**/ 0, 0, 0, /**/ 3, 4, 5,
                }, -12 / 16.0);
                
                test_h2dx_3d();
                
                test_fxx_3d({
                        1, -1, 1, /**/ 1, -1, 1, /**/ 1, -1, 1,
                        1, -1, 1, /**/ 1, -1, 1, /**/ 1, -1, 1,
                        1, -1, 1, /**/ 1, -1, 1, /**/ 1, -1, 1,
                }, 4.0);
                
                test_fyy_3d({
                        1, 1, 1, /**/ -1, -1, -1, /**/ 1, 1, 1,
                        1, 1, 1, /**/ -1, -1, -1, /**/ 1, 1, 1,
                        1, 1, 1, /**/ -1, -1, -1, /**/ 1, 1, 1,
                }, 4.0);
                
                test_fzz_3d({
                         1,  1,  1, /**/  1,  1,  1, /**/  1,  1,  1,
                        -1, -1, -1, /**/ -1, -1, -1, /**/ -1, -1, -1,
                         1,  1,  1, /**/  1,  1,  1, /**/  1,  1,  1,
                }, 4.0);
                
                test_fx_and_fy_3d_with_uint8_t({
                         50, 100,  20,
                        100,   0, 200,
                         70, 100,  30,
                        
                         50, 100,  20,
                        100,   0, 200,
                         70, 100,  30,
                        
                         50, 100,  20,
                        100,   0, 200,
                         70, 100,  30,
                }, 65.0 / 4.0, 15.0 / 4.0, 0);
        }
}
#include <iostream>
BOOST_AUTO_TEST_CASE(TEST_Differential)
{
        std::cout << "Differential\n";
        test_2d();
        test_3d();
}
#endif // TEST_Differential
