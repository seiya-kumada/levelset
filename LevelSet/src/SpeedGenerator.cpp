//-----------------------------------------------------------------------------------------------------------------------------
#include "SpeedGenerator.h"
#include "Parameters.h"
#include "Gaussian.h"
#include "MatrixFilter.h"
#include "MatrixForFilter.h"
#include <stdexcept>
//-----------------------------------------------------------------------------------------------------------------------------
using namespace kmd;
using namespace std;

//=============================================================================================================================
SpeedGenerator::SpeedGenerator(const Parameters* parameters)
//=============================================================================================================================
        : parameters_(parameters) {

        if ( parameters_->image_width_ == 0 || parameters_->image_height_ == 0 ) {
                throw runtime_error(Error::message("The size of the image is invalid.", __LINE__, __FILE__));
        }
}

//=============================================================================================================================
SpeedGenerator::~SpeedGenerator() {
//=============================================================================================================================
}

/*!
 *      This method should be called only once after calculating "modified_image_buffer_".
 */
 // 
//=============================================================================================================================
void SpeedGenerator::calculate_image_dependent_factors(const vector<double>& src_image) {
//=============================================================================================================================
        int w = parameters_->image_width_;
        int h = parameters_->image_height_;
        image_dependent_factors_.resize(w * h);
        vector<double> modified_image(w * h);
        modify_image(modified_image, src_image, parameters_->sigma_);
        for ( int j = 0, wj = w * j; j < h; ++j, wj += w ) {
                double* dst = &image_dependent_factors_[wj];
                const double* src = &modified_image[wj];
                for ( int i = 0; i < w; ++i ) {
                        dst[i] = 1.0 / (1.0 + src[i]);
                }
        }
}

// test ok
//=============================================================================================================================
double SpeedGenerator::calculate_curvature_dependent_speed(double curvature) const {
//=============================================================================================================================
        return -parameters_->epsilon_ * curvature;
}

// 
//=============================================================================================================================
double SpeedGenerator::calculate_curvature_dependent_speed(const double* buf, int w) const {
//=============================================================================================================================
        double curvature = calculate_curvature(buf, w);
        return calculate_curvature_dependent_speed(curvature);
}

namespace {

        // test ok
        //=====================================================================================================================
        inline double dy(const double* buf, int w) {
        //=====================================================================================================================
                return 0.5 * (buf[w] - buf[-w]);
        }

        // test ok
        //=====================================================================================================================
        inline double dx(const double* buf) {
        //=====================================================================================================================
                return 0.5 * (buf[1] - buf[-1]);
        }

        // test ok
        //=====================================================================================================================
        inline double dxdx(const double* buf) {
        //=====================================================================================================================
                return buf[1] - 2.0 * buf[0] + buf[-1];
        }

        // test ok
        //=====================================================================================================================
        inline double dydy(const double* buf, int w) {
        //=====================================================================================================================
                return buf[w] - 2.0 * buf[0] + buf[-w];
        }

        // test ok
        //=====================================================================================================================
        inline double dxdy(const double* buf, int w) {
        //=====================================================================================================================
                return 0.25 * (buf[1 + w] - buf[-1 + w] - buf[1 - w] + buf[-1 - w]);
        }
}

// test ok
//=============================================================================================================================
double SpeedGenerator::calculate_curvature(const double* buf, int w) const {
//=============================================================================================================================
        double a = dxdx(buf) * Utilities::square(dy(buf, w)) - 2.0 * dy(buf, w) * dx(buf) * dxdy(buf, w) + dydy(buf, w) * Utilities::square(dx(buf));
        double b = sqrt( Utilities::square(dx(buf)) + Utilities::square(dy(buf, w)) );
        return a / (b * b * b) / parameters_->space_step_;
}

//
//=============================================================================================================================
void SpeedGenerator::modify_image(vector<double>& dst, const vector<double>& src, double sigma) const {
//=============================================================================================================================
        int w = parameters_->image_width_;
        int h = parameters_->image_height_;
        vector<double> tmp(w * h);
        
        // src --> tmp
        if ( sigma == 0.0 ) {
                tmp = src;
        } else {
                apply_gaussian_filter(tmp, src, sigma);
        }

        // tmp --> dst
        apply_differential_filter(dst, tmp);
}

// test ok
//=============================================================================================================================
void SpeedGenerator::apply_gaussian_filter(vector<double>& dst, const vector<double>& src, double sigma) const {
//=============================================================================================================================
        int int_sigma = static_cast<int>(sigma);
        int matrix_size = 4 * int_sigma + 1;
        vector<double> matrix(Utilities::square(matrix_size));
        create_gaussian_matrix(matrix, matrix_size, sigma);

        MatrixFilter filter;
        filter.set_matrix(&matrix[0], matrix_size);
        filter.run(&dst[0], &src[0], parameters_->image_width_, parameters_->image_height_);
}

// test ok
//=============================================================================================================================
void SpeedGenerator::create_gaussian_matrix(vector<double>& matrix, int matrix_size, double sigma) const {
//=============================================================================================================================
        Gaussian gaussian(sigma);
        int range = matrix_size / 2;
        int k = 0;
        for ( int j = -range; j <= range; ++j ) {
                for ( int i = -range; i <= range; ++i ) {
                        matrix[k] = gaussian.get_value(i, j);
                        ++k;
                }
        }
}

namespace {

        // test ok
        //*********************************************************************************************************************
        struct calculate_strength : public binary_function<double, double, double> {
        //*********************************************************************************************************************
                double operator()(double x, double y) const {
                        return sqrt(Utilities::square(x) + Utilities::square(y));
                }
        };
}

// test ok
//=============================================================================================================================
void SpeedGenerator::apply_differential_filter(
        vector<double>&         dst,
        const vector<double>&   src
) const {
//=============================================================================================================================
        MatrixFilter filter;
        int w = parameters_->image_width_;
        int h = parameters_->image_height_;

        // derivative along x-direction
        vector<double> x_src(w * h);
        filter.set_matrix(MatrixForFilter::PrewittForX_, MatrixForFilter::PrewittForXSize_);
        filter.run(&x_src[0], &src[0], w, h);

        // derivative along y-direction
        vector<double> y_src(w * h);
        filter.set_matrix(MatrixForFilter::PrewittForY_, MatrixForFilter::PrewittForYSize_);
        filter.run(&y_src[0], &src[0], w, h);

        transform(x_src.begin(), x_src.end(), y_src.begin(), dst.begin(), calculate_strength());
}

#ifdef UNIT_TEST_SpeedGenerator
#define BOOST_TEST_DYN_LINK
//----------------------------------------------------------------------------------------------------------------------------
#include <kmd/without_warning/unit_test.h>
#include <kmd/without_warning/iostream.h>
#include "Front.h"
#include <fstream>
#include "TestUtilities.h"
#include "DistanceInitializer.h"
#include "TubeDomainGenerator.h"
#include "ShortestDistanceFinder.h"
#include <boost/test/floating_point_comparison.hpp>
#include "FrontDirection.h"
//----------------------------------------------------------------------------------------------------------------------------

namespace kmd {
        class SpeedGeneratorTester {
        public:
                static void test_apply_gaussian_filter(SpeedGenerator& obj, vector<double>& dst, const vector<double>& src, double sigma) {
                        return obj.apply_gaussian_filter(dst, src, sigma);
                }

                static void test_apply_differential_filter(SpeedGenerator& obj, vector<double>& dst, const vector<double>& src) {
                        return obj.apply_differential_filter(dst, src);
                }

                static void test_create_gaussian_matrix(SpeedGenerator& obj, std::vector<double>& matrix, int matrix_size, double sigma) {
                        obj.create_gaussian_matrix(matrix, matrix_size, sigma);
                }

                static double test_calculate_curvature_dependent_speed(SpeedGenerator& obj, double curvature) {
                        return obj.calculate_curvature_dependent_speed(curvature);
                }

                static double test_calculate_curvature(SpeedGenerator& obj, const double* buf, int w) {
                        return obj.calculate_curvature(buf, w);
                }

                static double test_get_image_dependent_factors(SpeedGenerator& obj, int w, int i, int j) {
                        return obj.image_dependent_factors_[w * j + i];
                }
        };

 }

//============================================================================================================================
BOOST_AUTO_TEST_CASE(TEST_SpeedGenerator) {
//============================================================================================================================

        cout << "SpeedGenerator\n";

        // dx, dy, dxdx, dydy, dxdy
        {
                double matrix[] = {
                      100, 100, 100, 100, 100,
                      100,  11,   2,   3, 100,
                      100,   4,   5,   6, 100,
                      100,   7,  48,  93, 100,
                      100, 100, 100, 100, 100,
                };

                BOOST_CHECK_EQUAL(dy(&matrix[12], 5), 23);
                BOOST_CHECK_EQUAL(dx(&matrix[12]), 1.0);
                BOOST_CHECK_EQUAL(dxdx(&matrix[12]), 0.0);
                BOOST_CHECK_EQUAL(dydy(&matrix[12], 5), 40.0);
                BOOST_CHECK_EQUAL(dxdy(&matrix[12], 5), 47.0/2.0);
        }

        // apply_gaussian_filter
        {

                Parameters parameters;
                int w = 5;
                int h = 5;
                parameters.image_width_ = w;
                parameters.image_height_ = h;
                vector<double> image_buffer(w * h, 1.0);
                parameters.epsilon_ = 0.0;

                SpeedGenerator generator(&parameters);
                vector<double> medium(w * h);
                double sigma = 2.0;
                SpeedGeneratorTester::test_apply_gaussian_filter(generator, medium, image_buffer, sigma);

                double s = 0.0;
                Gaussian gaussian(2.0);
                for ( int j = -2; j <= 2; ++j ) {
                        for ( int i = -2; i <= 2; ++i ) {
                                s += gaussian.get_value(i, j);
                        }
                }
                BOOST_CHECK_EQUAL(s, medium[12]);
        }

        // apply_differential_filter
        {


                Parameters parameters;
                int w = 5;
                int h = 5;
                parameters.image_width_ = w;
                parameters.image_height_ = h;
                vector<double> image_buffer(w * h);
                fill(image_buffer.begin(), image_buffer.begin() + 15, 1);
                fill(image_buffer.begin() + 15, image_buffer.end(), 0);
                //TestUtilities::write_image("buf", image_buffer, w, h, 1);
      
                parameters.epsilon_ = 0.0;

                SpeedGenerator generator(&parameters);
                vector<double> dst(w * h);
                SpeedGeneratorTester::test_apply_differential_filter(generator, dst, image_buffer);
                //TestUtilities::write_image("buf2", dst, w, h, 1);

                double answer[] = {
                        2.82843, 3, 3, 3, 2.82843,
                        3, 0, 0, 0, 3,
                        2.82843, 3, 3, 3, 2.82843,
                        2.23607, 3, 3, 3, 2.23607,
                        0, 0, 0, 0, 0,
                };

                for ( std::size_t i = 0, n = dst.size(); i < n; ++i ) {
                        BOOST_CHECK(fabs(answer[i] - dst[i]) < 0.00001);
                }
        }

        // calculate_strength
        {
                calculate_strength f;
                BOOST_CHECK(f(2, 4) == sqrt(20.0));

        }

        // create_gaussian_matrix
        {

                Parameters parameters;
                int w = 5;
                int h = 5;
                parameters.image_width_ = w;
                parameters.image_height_ = h;
                vector<double> image_buffer(w * h);
                parameters.epsilon_ = 0.0;

                SpeedGenerator generator(&parameters);
                double sigma = 2.0;
                int matrix_size = 5;
                vector<double> matrix(matrix_size * matrix_size);
                SpeedGeneratorTester::test_create_gaussian_matrix(generator, matrix, matrix_size, sigma);

                Gaussian gaussian(sigma);
                int range = matrix_size / 2;
                int k = 0;
                for ( int j = -range; j <= range; ++j ) {
                        for ( int i = -range; i <= range; ++i ) {
                                BOOST_CHECK_EQUAL(matrix[k], gaussian.get_value(i, j));
                                ++k;
                        }
                }
        }

        // calculate_curvature_dependent_speed
        {

                Parameters parameters;
                int w = 5;
                int h = 5;
                parameters.image_width_ = w;
                parameters.image_height_ = h;
                vector<double> image_buffer(w * h);
                parameters.epsilon_ = 2.0;

                SpeedGenerator generator(&parameters);
                BOOST_CHECK_EQUAL(-2.0, SpeedGeneratorTester::test_calculate_curvature_dependent_speed(generator, 1.0));

        }

        // calculate_curvature
        {
                // create an initial front.
                int center_x = 100;
                int center_y = 100;
                int w = 200;
                int h = 200;
                
                // create an image.
                vector<double> image_buffer(w * h, 0);
                for ( int j = 0; j < h; ++j ) {
                        for ( int i = 0; i < w; ++i ) {
                                double d = 50.0 - sqrt(Utilities::square<double>(center_x - i) + Utilities::square<double>(center_y - j));
                                if ( d > 0.0 ) {
                                        image_buffer[i + j * w] = -d;
                                }
                        }
                }
                //TestUtilities::write_image("SpeedGenerator::file0", image_buffer, w, h);

                Front front;
                double radius = 80.0;
                int step_size = 10000;
                TestUtilities::create_front(front, radius, step_size, center_x, center_y);
                Fronts fronts;
                fronts.push_back(front);
                //ofstream ofs("front");
                //TestUtilities::write_front(ofs, front);
                
                // set parameters.
                Parameters parameters;
                parameters.image_width_         = w;
                parameters.image_height_        = h;
                parameters.set_time_step(0.001);
                parameters.space_step_          = 1;
                parameters.time_step_number_    = 200;
                parameters.epsilon_             = 0.2;
                parameters.sigma_               = 2.0;
                

                //
                ShortestDistanceFinder finder(&parameters);
                int tube_width = 5;
                {
                        finder.set_tube_width(tube_width);
                        finder.run(fronts);/*
                        vector<double> dist(w * h);
                        for ( int j = 0, wj = w * j; j < h; ++j, wj += w ) {
                                double* d = &dist[wj];
                                for ( int i = 0; i < w; ++i ) {
                                        d[i] = finder.get_shortest_distance(i, j);
                                }
                        }
                        TestUtilities::write_image("dist", dist, w, h);*/
                }

                // extract tube
                TubeDomainGenerator tube_generator(&finder);
                tube_generator.run();
                const Tube& tube = tube_generator.get_tube();

                FrontDirection direction(FrontDirection::Outer);
                int signed_width = direction.attach_sign(tube_width);
                
                // create distance initializer.
                DistanceInitializer initializer(&finder);
                vector<double> distance(w * h, signed_width);
                initializer.run(&distance[0], distance.size(), fronts, tube);
                //TestUtilities::write_image("dist", distance, w, h, 2);

                vector<double> exact_distance(w * h);
                TestUtilities::calculate_exact_signed_distance(exact_distance, w, h);
                //TestUtilities::write_image("edist", exact_distance, w, h, 2);
                
                Tube::const_iterator beg = tube.begin();
                Tube::const_iterator end = tube.end();
                while ( beg != end ) {
                        int j = beg->first;
                        int wj = w * j;
                        const Segments& segments = beg->second;
                        for ( std::size_t k = 0, n = segments.size(); k < n; ++k ) {
                                const Segment& segment = segments[k];
                                int first = segment.first;
                                int second = segment.second;
                                for ( int i = first; i <= second; ++i ) {
                                       int index = wj + i;
                                       BOOST_CHECK(fabs(exact_distance[index] - distance[index]) < 1.0e-05);
                                }        
                        }
                        ++beg;        
                }

                
                // create speed generator.
                SpeedGenerator generator(&parameters);
                
      


                // see curvature except for the edge
                vector<double> curvature(w * h, 0);
                beg = tube.begin();
                while ( beg != end ) {
                        int j = beg->first;
                        int wj = w * j;
                        const Segments& segments = beg->second;
                        for ( std::size_t k = 0, n = segments.size(); k < n; ++k ) {
                                const Segment& segment = segments[k];
                                int first = segment.first;
                                int second = segment.second;
                                for ( int i = first; i <= second; ++i ) {
                                        int index = i + wj;
                                        curvature[index] = SpeedGeneratorTester::test_calculate_curvature(generator, &distance[index], w);                
                                }
                        }
                        ++beg;        
                }

                //TestUtilities::write_image("curvature", curvature, w, h);
                vector<double> exact_curvature(w * h, 0);
                TestUtilities::calculate_exact_curvature(exact_curvature, w, tube);
                //TestUtilities::write_image("exact_curvature", exact_curvature, w, h);
                
                struct Index {
                        int operator()(int x, int y) const {
                                return x + w_ * y;        
                        }
                        int     w_;
                        Index(int w)
                                : w_(w) {}
                };
                
                Index index(w);
                const int Xpos[] = { 1,  1, 1,  0, 0, -1, -1, -1};
                const int Ypos[] = {-1,  0, 1, -1, 1, -1,  0,  1};
                
                beg = tube.begin();
                while ( beg != end ) {
                        int j = beg->first;
                        const Segments& segments = beg->second;
                        for ( std::size_t k = 0, n = segments.size(); k < n; ++k ) {
                                const Segment& segment = segments[k];
                                int first = segment.first;
                                int second = segment.second;
                                for ( int i = first; i <= second; ++i ) {
                                        // consider only points inside tube
                                        if ( 
                                                (distance[index(i + Xpos[0], j + Ypos[0])] != tube_width)
                                             && (distance[index(i + Xpos[1], j + Ypos[1])] != tube_width)
                                             && (distance[index(i + Xpos[2], j + Ypos[2])] != tube_width)
                                             && (distance[index(i + Xpos[3], j + Ypos[3])] != tube_width)
                                             && (distance[index(i + Xpos[4], j + Ypos[4])] != tube_width)
                                             && (distance[index(i + Xpos[5], j + Ypos[5])] != tube_width)
                                             && (distance[index(i + Xpos[6], j + Ypos[6])] != tube_width)
                                             && (distance[index(i + Xpos[7], j + Ypos[7])] != tube_width)
                                        ) {
                                                BOOST_CHECK(fabs(curvature[index(i, j)] - exact_curvature[index(i, j)]) < 1.0e-05);
                                        }
                                }        
                        }
                        
                        ++beg;        
                }

                TestUtilities::delete_front(front);
        }
}
#endif // UNIT_TEST_SpeedGenerator


