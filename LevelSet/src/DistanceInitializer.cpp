//-----------------------------------------------------------------------------------------------------------------------------
#include "DistanceInitializer.h"
#include "ShortestDistanceFinder.h"
#include "Point.h"
#include "FrontDirection.h"
#include "Parameters.h"
#include "Pi.h"
//-----------------------------------------------------------------------------------------------------------------------------
using namespace kmd;
using namespace std;
using namespace boost;

//=============================================================================================================================
DistanceInitializer::DistanceInitializer(const ShortestDistanceFinder* distance_finder)
//=============================================================================================================================
        : distance_finder_(distance_finder) {

        assign_value_[FrontDirection::Inner] = &DistanceInitializer::assign_value_inner;
        assign_value_[FrontDirection::Outer] = &DistanceInitializer::assign_value_outer;
}

//=============================================================================================================================
DistanceInitializer::~DistanceInitializer() {
//=============================================================================================================================
}

namespace {

        const double DoublePi = 2.0 * M_PI;

        // http://kone.vis.ne.jp/diary/diaryb09.html
        // test ok
        //=====================================================================================================================
        bool is_inside2(const Point* const front[], std::size_t size, int x, int y) {
        //=====================================================================================================================
                const Point p(x, y);
                double sum = 0.0;
                for ( std::size_t i = 0; i < size; ++i ) {
                        const Point& a = *front[i];
                        std::size_t index = i + 1;
                        if ( index > size - 1 ) {
                                index = 0;
                        }
                        const Point& b = *front[index];
                        const Point ap = a - p;
                        const Point bp = b - p;
                        double s = cross_product(ap, bp); // corresponds to sin(theta)
                        double c = inner_product(ap, bp); // corresponds to cos(theta)
                        sum += atan2(s, c);
                }
                return fabs(sum - DoublePi) < 1.0e-14;
        }

        // test ok
        //=====================================================================================================================
        bool is_inside3(const Point* const front[], std::size_t size, int x, int y) {
        //=====================================================================================================================
                const Point p(x, y);
                for ( std::size_t i = 0; i < size; ++i ) {
                        const Point& a = *front[i];
                        std::size_t index = i + 1;
                        if ( index > size - 1 ) {
                                index = 0;
                        }
                        const Point& b = *front[index];

                        const Point ba = b - a;
                        const Point pa = p - a;
                        double s = cross_product(ba, pa);
                        if ( s < 0.0 ) {
                                return false;
                        }

                }
                return true;
        }

        // This function is faster than is_inside2.
        // http://kone.vis.ne.jp/diary/diaryb09.html
        // test ok
        //=====================================================================================================================
        bool is_inside(const Point* const front[], std::size_t size, int x, int y) {
        //=====================================================================================================================
                int c = 0;
                double a_x, a_y, b_x, b_y;
                double f, q;
                for ( std::size_t i = 0; i < size; ++i ) {
                        const Point* a = front[i];
                        std::size_t index = i + 1;
                        if ( index > size - 1 ) {
                                index = 0;
                        }
                        const Point* b = front[index];
                        a_y = a->y();
                        b_y = b->y();

                        if ( (a_y < y && b_y < y) || (a_y > y && b_y > y) ) {
                                continue;
                        }

                        a_x = a->x();
                        b_x = b->x();

                        if ( a_x > x && b_x > x ) {
                                continue;
                        }

                        if ( a_x < x && b_x < x ) {
                                if ( b_y == y ) {
                                        std::size_t next_index = index + 1;
                                        if ( next_index > size - 1 ) {
                                                next_index = 0;
                                        }
                                        double ny = front[next_index]->y();
                                        if ( (ny - y) * (y - a_y) > 0 ) {
                                                ++c;
                                        }
                                } else if ( a_y == y ) {
                                        // do nothing!
                                } else {
                                        ++c;
                               }
                        } else {
                                f = fabs(y - a_y) / fabs(a_y - b_y);
                                q = a_x + f * (b_x - a_x);
                                if ( q < x ) {
                                        ++c;
                               }
                        }
                }
                return c & 0x01; // Is c odd ?
        }

        bool is_inside(const Fronts& fronts, int x, int y) {
                for ( std::size_t i = 0, n = fronts.size(); i < n; ++i ) {
                        const Front& front = fronts[i];
                        if ( is_inside(&front[0], front.size(), x, y) ) {
                                return true;
                        }
                }
                return false;
        }
}

// test ok
//=============================================================================================================================
void DistanceInitializer::run(
        double*                 new_buffer,
        const double*           old_buffer,
        std::size_t             size,
        const Tube&             new_tube,
        const FrontDirection&   direction
) const {
//=============================================================================================================================
        const Parameters* parameters = distance_finder_->get_parameters();
        int w = parameters->image_width_;
        int h = parameters->image_height_;
        if ( size != static_cast<std::size_t>(w * h) ) {
                // [check] Results of __LINE__ and __FILE__ are not what you would expect
                throw runtime_error(Error::message("The size of the input buffer is invalid", __LINE__, __FILE__));
        }

        const std::vector<double>& shortest_distances = distance_finder_->get_shortest_distances();
        Tube::const_iterator beg = new_tube.begin();
        Tube::const_iterator end = new_tube.end();
        while ( beg != end ) {
                int j = beg->first;
                int wj = w * j;
                const double* dist = &shortest_distances[wj];
                double* dst = new_buffer + wj;
                const double* src = old_buffer + wj;
                const Segments& segments = beg->second;
                for ( std::size_t k = 0, n = segments.size(); k < n; ++k ) {
                        const Segment& segment = segments[k];
                        int first = segment.first;
                        int second = segment.second;
                        for ( int i = first; i <= second; ++i ) {
                                (this->*assign_value_[direction.get_value()])(dst[i], src[i], dist[i]);
                        }
                }
                ++beg;
        }
}

//=============================================================================================================================
void DistanceInitializer::assign_value_inner(double& dst, double src, double dist) const {
//=============================================================================================================================
        dst = -dist;
        if ( src > 0 ) {
                dst *= -1.0;
        }
}

//=============================================================================================================================
void DistanceInitializer::assign_value_outer(double& dst, double src, double dist) const {
//=============================================================================================================================
        dst = dist;
        if ( src < 0 ) {
                dst *= -1.0;
        }
}

//=============================================================================================================================
void DistanceInitializer::run(
        double*                 buffer,
        std::size_t             size,
        const Fronts&      fronts,
        const Tube&             tube
) const {
//=============================================================================================================================
        const Parameters* parameters = distance_finder_->get_parameters();
        int w = parameters->image_width_;
        int h = parameters->image_height_;
        if ( size != static_cast<std::size_t>(w * h) ) {
                throw runtime_error(Error::message("The size of the input buffer is invalid", __LINE__, __FILE__));
        }

        const vector<double>& shortest_distances = distance_finder_->get_shortest_distances();
        Tube::const_iterator beg = tube.begin();
        Tube::const_iterator end = tube.end();
        while ( beg != end ) {
                const int j = beg->first;
                const int wj = w * j;
                double* dst = buffer + wj;
                const double* src = &shortest_distances[wj];
                const Segments& segments = beg->second;
                for ( std::size_t k = 0, n = segments.size(); k < n; ++k ) {
                        const Segment& segment = segments[k];
                        const int first = segment.first;
                        const int second = segment.second;
                        for ( int i = first; i <= second; ++i ) {
                                dst[i] = src[i];
                                if ( is_inside(fronts, i, j) ) {
                                        dst[i] *= -1.0;
                                }
                        }
                }
                ++beg;
        }
}

#ifdef UNIT_TEST_DistanceInitializer
#define BOOST_TEST_DYN_LINK
//----------------------------------------------------------------------------------------------------------------------------
#include <kmd/without_warning/unit_test.h>
#include <kmd/without_warning/iostream.h>
#include "Parameters.h"
#include <cmath>
#include <fstream>
#include "TestUtilities.h"
#include "ZeroLevelSetDetector.h"
#include <boost/timer.hpp>
#include <set>
#include "TubeDomainGenerator.h"
#include <boost/test/floating_point_comparison.hpp>
//----------------------------------------------------------------------------------------------------------------------------
using namespace std;
using namespace boost;


//============================================================================================================================
BOOST_AUTO_TEST_CASE(TEST_DistanceInitializer) {
//============================================================================================================================

        cout << "DistanceInitializer\n";

        // is_inside/is_inside2/is_inside3
        {
                double radius = 80.0;
                int step_size = 300;
                int center_x = 100;
                int center_y = 100;
                Front front;
                TestUtilities::create_front(front, radius, step_size, center_x, center_y);

                BOOST_CHECK(is_inside(&front[0], front.size(), center_x, center_y));
                BOOST_CHECK(is_inside2(&front[0], front.size(), center_x, center_y));
                BOOST_CHECK(is_inside3(&front[0], front.size(), center_x, center_y));
                BOOST_CHECK(!is_inside(&front[0], front.size(), 0, 0));
                BOOST_CHECK(!is_inside2(&front[0], front.size(), 0, 0));
                BOOST_CHECK(!is_inside3(&front[0], front.size(), 0, 0));
                TestUtilities::delete_front(front);
        }

        // is_inside/is_inside2/is_inside3
        {
                ifstream ifs("../data/unit_test_data/front_60");
                if ( ifs.fail() ) {
                        cout << "unable to open file\n";
                        return;
                }
                const int LineSize = 100;
                char tmp[LineSize];

                vector<double> xs, ys;
                bool read = true;
                while ( ifs.getline(tmp, LineSize) ) {
                        string s(tmp);
                        string::size_type i = s.find(" ");
                        if ( i != string::npos ) {
                                if ( read ) {
                                        xs.push_back(lexical_cast<double>(s.substr(0, i)));
                                        ys.push_back(lexical_cast<double>(s.substr(i + 1, s.length())));
                                        read = false;
                                } else {
                                        read = true;
                                }
                        }
                }

                Front front;
                for ( std::size_t i = 0, n = xs.size(); i < n; ++i ) {
                        front.push_back(new Point(xs[i], ys[i]));
                }

                for ( int i = 60; i < 200; ++i ) {
                        if ( 71 <= i && i <= 167 ) {
                                BOOST_CHECK(is_inside2(&front[0], front.size(), i, 140));
                                BOOST_CHECK(is_inside(&front[0], front.size(), i, 140));
                                //BOOST_CHECK(is_inside3(&front[0], front.size(), i, 140));
                        } else {
                                BOOST_CHECK(!is_inside2(&front[0], front.size(), i, 140));
                                BOOST_CHECK(!is_inside(&front[0], front.size(), i, 140));
                                //BOOST_CHECK(!is_inside3(&front[0], front.size(), i, 140));
                        }
                }
                 TestUtilities::delete_front(front);

        }

        // run: FrontDirection::Outer
        {
                Front front;
                try {
                        double radius = 80.0;
                        int step_size = 300;
                        int center_x = 100;
                        int center_y = 100;

                        TestUtilities::create_front(front, radius, step_size, center_x, center_y);

                        Parameters parameters;
                        int w = 200;
                        int h = 200;
                        parameters.image_width_ = w;
                        parameters.image_height_ = h;
                        ShortestDistanceFinder finder(&parameters);
                        Fronts fronts;
                        fronts.push_back(front);
                        FrontDirection direction(FrontDirection::Outer);
                        int tube_width = 6;
                        int signed_tube_width = direction.attach_sign(tube_width);
                        finder.set_tube_width(tube_width);
                        finder.run(fronts);

                        //TestUtilities::write_image("nearest_infos", finder.get_shortest_distances(), w, h, 2);


                        TubeDomainGenerator tube_generator(&finder);
                        tube_generator.run();
                        const Tube& tube_domain = tube_generator.get_tube();
                        //TestUtilities::write_tube("tube", tube_domain);

                        DistanceInitializer initializer(&finder);

                        vector<double> buffer(w * h, signed_tube_width);
                        initializer.run(&buffer[0], buffer.size(), fronts, tube_domain);
                        //TestUtilities::write_image("buffer", buffer, w, h, 2);

                        BOOST_CHECK(!is_inside(&front[0], front.size(), 0, 0));
                        BOOST_CHECK(is_inside(&front[0], front.size(), 100, 100));

                        const vector<double>& shortest_distances = finder.get_shortest_distances();
                        for ( int j = 0, wj = w * j; j < h; ++j, wj += w ) {
                                const double* dis = &shortest_distances[wj];
                                const double* tmp = &buffer[wj];
                                for ( int i = 0; i < w; ++i ) {
                                        BOOST_CHECK_CLOSE(fabs(tmp[i]), dis[i], 1.0e-08);
                                }
                        }

                        vector<double> new_buffer(buffer.size(), signed_tube_width);
                        initializer.run(&new_buffer[0], &buffer[0], buffer.size(), tube_generator.get_tube(), direction);
                        for ( int j = 0, wj = w * j; j < h; ++j, wj += w ) {
                                const double* src = &buffer[wj];
                                const double* dst = &new_buffer[wj];
                                for ( int i = 0; i < w; ++i ) {
                                        BOOST_CHECK_CLOSE(src[i], dst[i], 1.0e-08);
                                }
                        }

                        //TestUtilities::write_image("new", new_buffer, w, h);

                } catch (const runtime_error& error) {
                        cout << error.what() << endl;
                }
                TestUtilities::delete_front(front);
        }

        // run: FrontDirection::Inner
        {
                Front front;
                try {
                        double radius = 80.0;
                        int step_size = 300;
                        int center_x = 100;
                        int center_y = 100;

                        TestUtilities::create_front(front, radius, step_size, center_x, center_y);

                        Parameters parameters;
                        int w = 200;
                        int h = 200;
                        parameters.image_width_ = w;
                        parameters.image_height_ = h;
                        ShortestDistanceFinder finder(&parameters);
                        Fronts fronts;
                        fronts.push_back(front);
                        FrontDirection direction(FrontDirection::Inner);
                        int tube_width = 6;
                        int signed_tube_width = direction.attach_sign(tube_width);
                        finder.set_tube_width(tube_width);
                        finder.run(fronts);

                        //TestUtilities::write_image("nearest_infos", finder.get_shortest_distances(), w, h, 2);


                        TubeDomainGenerator tube_generator(&finder);
                        tube_generator.run();
                        const Tube& tube_domain = tube_generator.get_tube();
                        //TestUtilities::write_tube("tube", tube_domain);

                        DistanceInitializer initializer(&finder);

                        vector<double> buffer(w * h, signed_tube_width);
                        initializer.run(&buffer[0], buffer.size(), fronts, tube_domain);
                        //TestUtilities::write_image("buffer", buffer, w, h, 2);

                        BOOST_CHECK(!is_inside(&front[0], front.size(), 0, 0));
                        BOOST_CHECK(is_inside(&front[0], front.size(), 100, 100));

                        const vector<double>& shortest_distances = finder.get_shortest_distances();
                        for ( int j = 0, wj = w * j; j < h; ++j, wj += w ) {
                                const double* dis = &shortest_distances[wj];
                                const double* tmp = &buffer[wj];
                                for ( int i = 0; i < w; ++i ) {
                                        BOOST_CHECK_CLOSE(fabs(tmp[i]), dis[i], 1.0e-08);
                                }
                        }

                        vector<double> new_buffer(buffer.size(), signed_tube_width);
                        initializer.run(&new_buffer[0], &buffer[0], buffer.size(), tube_generator.get_tube(), direction);
                        for ( int j = 0, wj = w * j; j < h; ++j, wj += w ) {
                                const double* src = &buffer[wj];
                                const double* dst = &new_buffer[wj];
                                for ( int i = 0; i < w; ++i ) {
                                        BOOST_CHECK_CLOSE(src[i], dst[i], 1.0e-08);
                                }
                        }

                        //TestUtilities::write_image("new", new_buffer, w, h);

                } catch (const runtime_error& error) {
                        cout << error.what() << endl;
                }
                TestUtilities::delete_front(front);
        }

}
#endif // UNIT_TEST_DistanceInitializer
