//-----------------------------------------------------------------------------------------------------------------------------
#include "ShortestDistanceFinder.h"
#include "Parameters.h"
#include <stdexcept>
#include "Point.h"
//-----------------------------------------------------------------------------------------------------------------------------
using namespace kmd;
using namespace std;

//=============================================================================================================================
ShortestDistanceFinder::ShortestDistanceFinder(const Parameters* parameters)
//=============================================================================================================================
        : parameters_(parameters),
          shortest_distances_(parameters_->image_width_ * parameters_->image_height_) {

        if ( parameters_->image_width_ == 0 || parameters_->image_height_ == 0 ) {
                throw runtime_error(Error::message("The size of the image is invalid.", __LINE__, __FILE__));
        }
}

//=============================================================================================================================
ShortestDistanceFinder::~ShortestDistanceFinder() {
//=============================================================================================================================
}

// test ok
//=============================================================================================================================
void ShortestDistanceFinder::scan_area(double* infos, int w, const Point* a, const Point* b) const {
//=============================================================================================================================
        Point c = 0.5 * ((*a) + (*b));
        int wstart = static_cast<int>(c.x()) - tube_width_;
        int wend   = static_cast<int>(c.x()) + tube_width_;
        int hstart = static_cast<int>(c.y()) - tube_width_;
        int hend   = static_cast<int>(c.y()) + tube_width_;
        double dist;
        for ( int j = hstart, wj = w * j; j <= hend; ++j, wj += w ) {
                double* tmp = &infos[wj];
                for ( int i = wstart; i <= wend; ++i ) {
                        dist = sqrt(calculate_square_distance(*a, *b, Point(i, j)));
                        if ( tmp[i] > dist ) {
                                tmp[i] = dist;
                        }
                }
        }
}

// test ok
//=============================================================================================================================
double ShortestDistanceFinder::calculate_square_distance(const Point& a, const Point& b, const Point& p) const {
//=============================================================================================================================
        const Point c = b - a;
        const Point ap = a - p;
        const double cap = inner_product(c, ap);
        const double t = -cap / c.square_norm();
        if ( 0.0 <= t && t <= 1.0 ) {
                return fabs(ap.square_norm() + t * cap);
        } else if ( t < 0.0 ) {
                return ap.square_norm();
        } else {
                const Point bp = b - p;
                return bp.square_norm();
        }
}

//=============================================================================================================================
void ShortestDistanceFinder::walk_on_front(const Point* const front[], std::size_t size, int w) {
//=============================================================================================================================
        for ( std::size_t i = 0; i < size; ++i ) {
                const Point* a = front[i];
                std::size_t index = i + 1;
                if ( index > size - 1 ) {
                        index = 0;
                }
                const Point* b = front[index];
                scan_area(&shortest_distances_[0], w, a, b);
        }
}

// test ok
//=============================================================================================================================
void ShortestDistanceFinder::run(const Fronts& fronts) {
//=============================================================================================================================
        std::fill(shortest_distances_.begin(), shortest_distances_.end(), tube_width_);
        int w = parameters_->image_width_;
        for ( std::size_t i = 0, n = fronts.size(); i < n; ++i ) {
                const Front& front = fronts[i];
                walk_on_front(&front[0], front.size(), w);
        }
}

#ifdef UNIT_TEST_ShortestDistanceFinder
#define BOOST_TEST_DYN_LINK
//----------------------------------------------------------------------------------------------------------------------------
#include <kmd/without_warning/unit_test.h>
#include <iostream>
#include "Front.h"
#include <fstream>
#include "TestUtilities.h"
#include "ZeroLevelSetDetector.h"
#include <kmd/without_warning/lexical_cast.h>
//----------------------------------------------------------------------------------------------------------------------------
using namespace boost;

namespace kmd {
        class ShortestDistanceFinderTester {
        public:
        
                static double test_calculate_square_distance(
                        ShortestDistanceFinder& obj,
                        const Point&            a,
                        const Point&            b,
                        const Point&            p
                ) {
                        return obj.calculate_square_distance(a, b, p);
                }

                static void test_scan_area(ShortestDistanceFinder& obj, double* infos, int w, const Point* a, const Point* b) {
                        obj.scan_area(infos, w, a, b);
                }
        };
}

//============================================================================================================================
BOOST_AUTO_TEST_CASE(TEST_ShortestDistanceFinder) {
//============================================================================================================================

        cout << "ShortestDistanceFinder\n";

        // calculate_square_distance
        {
                Point a(0, 0);
                Point b(2, 0);
                Parameters params;
                params.image_width_ = 100;
                params.image_height_ = 100;
                ShortestDistanceFinder finder(&params);
                BOOST_CHECK_EQUAL(1, ShortestDistanceFinderTester::test_calculate_square_distance(finder, a, b, Point(1, 1)));
                BOOST_CHECK_EQUAL(2, ShortestDistanceFinderTester::test_calculate_square_distance(finder, a, b, Point(3, 1)));
        }

        // set/get_tube_width
        {
                Parameters params;
                params.image_width_ = 100;
                params.image_height_ = 100;
                ShortestDistanceFinder finder(&params);
                finder.set_tube_width(3);
                BOOST_CHECK(3 == finder.get_tube_width());
        }

        // scan_area
        {
                int w = 10;
                int h = 10;
                vector<double> infos(w * h);
                int tube_width = 2;
                fill(infos.begin(), infos.end(), tube_width);


                Point a(4, 5);
                Point b(6, 5);

                Parameters params;
                params.image_width_ = 100;
                params.image_height_ = 100;
                ShortestDistanceFinder finder(&params);
                finder.set_tube_width(tube_width);
                ShortestDistanceFinderTester::test_scan_area(finder, &infos[0], w, &a, &b);

                for ( int j = 0, wj = 0; j < h; ++j, wj += w ) {
                        const double* tmp = &infos[wj];
                        for ( int i = 0; i < w; ++i ) {
                                if ( 3 <= i && i <= 7 && 4 <= j && j <= 6 ) {
                                        BOOST_CHECK(tmp[i] < tube_width);
                                }
                                else {
                                        BOOST_CHECK(tmp[i] == tube_width);
                                }
                        }
                }

        }

        // run
        {

                int w = 200;
                int h = 200;
                vector<double> buffer(w * h);

                int step_size = 1000;
                Front front0;
                TestUtilities::create_front(front0, 40, step_size, 50, 100);
                Front front1;
                TestUtilities::create_front(front1, 40, step_size, 150, 100);
                
                Parameters params;
                params.image_width_ = w;
                params.image_height_ = h;
                
                ShortestDistanceFinder finder(&params);
                int tube_width = 6;
                finder.set_tube_width(tube_width);
                Fronts fronts;
                fronts.push_back(front0);
                fronts.push_back(front1);
                finder.run(fronts);
                
                //TestUtilities::write_image("buf", finder.get_shortest_distances(), w, h, 2);
                const vector<double>& dist = finder.get_shortest_distances();
                for ( std::size_t i = 0, n = fronts.size(); i < n; ++i ) {
                        const Front& front = fronts[i];
                        for ( std::size_t j = 0, m = front.size(); j < m; ++j ) {
                                const Point* p = front[j];
                                int x = static_cast<int>(p->x());
                                int y = static_cast<int>(p->y());
                                BOOST_CHECK(dist[x + w * y] < tube_width);
                        }
                }
                
        }

        // run
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
                        int tube_width = 6;
                        finder.set_tube_width(tube_width);
                        Fronts fronts;
                        fronts.push_back(front);
                        finder.run(fronts);
                        const vector<double>& dist = finder.get_shortest_distances();
                        for ( std::size_t i = 0, n = fronts.size(); i < n; ++i ) {
                                const Front& front = fronts[i];
                                for ( std::size_t j = 0, m = front.size(); j < m; ++j ) {
                                        const Point* p = front[j];
                                        int x = static_cast<int>(p->x());
                                        int y = static_cast<int>(p->y());
                                        BOOST_CHECK(dist[x + w * y] < tube_width);
                                }
                        }       
                        
                        //TestUtilities::write_image("ShortestDistanceFinder::buffer", buffer, w, h);
                } catch (const runtime_error& error) {
                        cout << error.what() << endl;
                }
                TestUtilities::delete_front(front);
        }

}
#endif // UNIT_TEST

