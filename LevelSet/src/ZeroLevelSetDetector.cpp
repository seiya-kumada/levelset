//-----------------------------------------------------------------------------------------------------------------------------
#include "ZeroLevelSetDetector.h"
#include "Parameters.h"
#include "LinearInterpolator.h"
#include "Point.h"
#include "TestUtilities.h"
//-----------------------------------------------------------------------------------------------------------------------------
using namespace kmd;
using namespace std;
using namespace boost;

namespace {

        void delete_front(Front& front) {
                for ( std::size_t i = 0, n = front.size(); i < n; ++i ) {
                        delete front[i];
                }
                front.clear();
        }
}

namespace kmd {

        //*********************************************************************************************************************
        class Cell : private noncopyable {
        //*********************************************************************************************************************
                friend class CellTester;

        private:
                int     x_;
                int     y_;
                int     flag_;

                //! (x_, y_), (x_ + 1, y_), (x_ + 1, y_ + 1), (x_, y_ + 1)
                double  vertices_[4];

                LinearInterpolator      interpolator_;

                static const double     SmallSize_;



        public:

                //
                //=============================================================================================================
                Cell()
                //=============================================================================================================
                        : x_(0),
                          y_(0),
                          flag_(0) {

                        vertices_[0] = 0.0;
                        vertices_[1] = 0.0;
                        vertices_[2] = 0.0;
                        vertices_[3] = 0.0;



                }

                // test ok
                //=============================================================================================================
                int x() const {
                //=============================================================================================================
                        return x_;
                }

                // test ok
                //=============================================================================================================
                int y() const {
                //=============================================================================================================
                        return y_;
                }

                // test ok
                //=============================================================================================================
                void set(int x, int y, double v0, double v1, double v2, double v3) {
                //=============================================================================================================
                        x_ = x;
                        y_ = y;
                        vertices_[0] = v0;
                        vertices_[1] = v1;
                        vertices_[2] = v2;
                        vertices_[3] = v3;
                }

                // test ok
                //=============================================================================================================
                int get_flag() const {
                //=============================================================================================================
                        return flag_;
                }

                // test ok
                //=============================================================================================================
                void set_flag(int flag) {
                //=============================================================================================================
                        flag_ = flag;
                }

                // test ok
                //=============================================================================================================
                bool includes_zero_level_set() const {
                //=============================================================================================================
                        return *std::max_element(vertices_, vertices_ + 4) > 0 && *std::min_element(vertices_, vertices_ + 4) < 0;
                }

                /*!
                 *      structure Direction specifies the direction to go on.
                 *
                 *      "t" means a target.
                 *      2   1   0
                 *      3   t   7
                 *      4   5   6
                 */
                template<int N>
                struct Direction {
                        static const int value_ = N;
                };

                typedef Direction<0> D0;
                typedef Direction<1> D1;
                typedef Direction<2> D2;
                typedef Direction<3> D3;
                typedef Direction<4> D4;
                typedef Direction<5> D5;
                typedef Direction<6> D6;
                typedef Direction<7> D7;

                //=============================================================================================================
                template<int N>
                bool set_point(Point& p, double dx, double dy) {
                //=============================================================================================================
                        p.set_x(x_ + dx);
                        p.set_y(y_ + dy);
                        flag_ = N;
                        return true;
                }

                //=============================================================================================================
                bool set_point(D0, Point& p) {
                //=============================================================================================================
                        double dy = 1.0;
                        double dx = interpolator_.get_x_value(dy);
                        //assert(TestUtilities::print_any("dx", dx));
                        if ( is_equal_to_1(dx) ) {
                                return set_point<D0::value_>(p, dx, dy);
                        } else {
                                return false;
                        }
                }

                //=============================================================================================================
                bool set_point(D1, Point& p) {
                //=============================================================================================================
                        double dy = 1.0;
                        double dx = interpolator_.get_x_value(dy);
                        if ( is_between_0_and_1(dx) ) {
                                return set_point<D1::value_>(p, dx, dy);
                        } else {
                                return false;
                        }
                }

                //=============================================================================================================
                bool set_point(D2, Point& p) {
                //=============================================================================================================
                        double dy = 1.0;
                        double dx = interpolator_.get_x_value(dy);
                        if ( is_equal_to_0(dx) ) {
                                return set_point<D2::value_>(p, dx, dy);
                        } else {
                                return false;
                        }
                }

                //=============================================================================================================
                bool set_point(D3, Point& p) {
                //=============================================================================================================
                        double dx = 0.0;
                        double dy = interpolator_.get_y_value(dx);
                        if ( is_between_0_and_1(dy) ) {
                                return set_point<D3::value_>(p, dx, dy);
                        } else {
                                return false;
                        }
                }

                //=============================================================================================================
                bool set_point(D4, Point& p) {
                //=============================================================================================================
                        double dx = 0.0;
                        double dy = interpolator_.get_y_value(dx);
                        if ( is_equal_to_0(dy) ) {
                                return set_point<D4::value_>(p, dx, dy);
                        } else {
                                return false;
                        }
                }

                //=============================================================================================================
                bool set_point(D5, Point& p) {
                //=============================================================================================================
                        double dy = 0.0;
                        double dx = interpolator_.get_x_value(dy);
                        if ( is_between_0_and_1(dx) ) {
                                return set_point<D5::value_>(p, dx, dy);
                        } else {
                                return false;
                        }
                }

                //=============================================================================================================
                bool set_point(D6, Point& p) {
                //=============================================================================================================
                        double dy = 0.0;
                        double dx = interpolator_.get_x_value(dy);
                        if ( is_equal_to_1(dx) ) {
                                return set_point<D6::value_>(p, dx, dy);
                        } else {
                                return false;
                        }
                }

                //=============================================================================================================
                bool set_point(D7, Point& p) {
                //=============================================================================================================
                        double dx = 1.0;
                        double dy = interpolator_.get_y_value(dx);
                        if ( is_between_0_and_1(dy) ) {
                                return set_point<D7::value_>(p, dx, dy);
                        } else {
                                return false;
                        }
                }

                //=============================================================================================================
                template<int S, int T, int U>
                bool set_point(Point& p) {
                //=============================================================================================================
                        return set_point(Direction<S>(), p) || set_point(Direction<T>(), p) || set_point(Direction<U>(), p);
                }


                //=============================================================================================================
                bool trace_up(Point& p, double s) {
                //=============================================================================================================
                        if ( s > 0.0 ) {
                                return set_point<3, 2, 1>(p);
                        } else {
                                return set_point<7, 0, 1>(p);
                        }
                }

                //=============================================================================================================
                bool trace_left(Point& p, double s) {
                //=============================================================================================================
                        if ( s > 0.0 ) {
                                return set_point<1, 2, 3>(p);
                        } else {
                                return set_point<5, 4, 3>(p);
                        }
                }

                //=============================================================================================================
                bool trace_down(Point& p, double s) {
                //=============================================================================================================
                        if ( s > 0.0 ) {
                                return set_point<7, 6, 5>(p);
                        } else {
                                return set_point<3, 4, 5>(p);
                        }
                }


                //=============================================================================================================
                bool trace_right(Point& p, double s) {
                //=============================================================================================================
                        if ( s > 0.0 ) {
                                return set_point<5, 6, 7>(p);
                        } else {
                                return set_point<1, 0, 7>(p);
                        }
                }

                // test ok
                //=============================================================================================================
                bool find_next(Point& p, int flag) {
                //=============================================================================================================
                        interpolator_.set_value(vertices_[0], vertices_[1], vertices_[2], vertices_[3]);
                        interpolator_.run();
                        double s = interpolator_.get_sign();
                        //assert(TestUtilities::print_any("flag", flag));
                        //assert(TestUtilities::print_any("s", s));
                        switch ( flag ) {
                        case 0:
                                return trace_up(p, s);
                        case 1:
                                return trace_up(p, s);
                        case 2:
                                return trace_up(p, s);
                        case 3:
                                return trace_left(p, s);
                        case 4:
                                return trace_down(p, s);
                        case 5:
                                return trace_down(p, s);
                        case 6:
                                return trace_down(p, s);
                        case 7:
                                return trace_right(p, s);
                        }
                        return false;
                }

                typedef bool (Cell::*Trace)(Point& p, double s);

                //=============================================================================================================
                template<int S, int T, int U>
                bool try_starting(Front& front, double s, Trace tr) {
                //=============================================================================================================
                        Point p;
                        if ( set_point(Direction<S>(), p) || set_point(Direction<T>(), p) || set_point(Direction<U>(), p) ) {
                                front.push_back(p.clone());
                                if ( (this->*tr)(p, s) ) {
                                        front.push_back(p.clone());
                                        return true;
                                } else {
                                        return false;
                                }
                        } else {
                                return false;
                        }
                }

                //=============================================================================================================
                template<int N>
                bool try_starting(Front& front, double s, Trace tr) {
                //=============================================================================================================
                        Point p;
                        if ( set_point(Direction<N>(), p) ) {
                                front.push_back(p.clone());
                                if ( (this->*tr)(p, s) ) {
                                        front.push_back(p.clone());
                                        return true;
                                } else {
                                        return false;
                                }
                        } else {
                                return false;
                        }
                }

                // test ok
                //=============================================================================================================
                bool find_first_segment(Front& front) {
                //=============================================================================================================
                        interpolator_.set_value(vertices_[0], vertices_[1], vertices_[2], vertices_[3]);
                        interpolator_.run();
                        double s = interpolator_.get_sign();
                        if (
                                   try_starting<0, 1, 2>(front, s, &Cell::trace_down)
                                || try_starting<4, 5, 6>(front, s, &Cell::trace_up)
                                || try_starting<3>(front, s, &Cell::trace_right)
                                || try_starting<7>(front, s, &Cell::trace_left)
                        ) {
                                return true;
                        } else {
                                delete_front(front);
                                return false;
                        }
                }

                //=============================================================================================================
                bool is_between_0_and_1(double v) const {
                //=============================================================================================================
                        return (SmallSize_ < v) && (v < 1.0 - SmallSize_);
                }

                //=============================================================================================================
                bool is_equal_to_0(double v) const {
                //=============================================================================================================
                        return (0.0 <= v) && (v <= SmallSize_);
                }

                //=============================================================================================================
                bool is_equal_to_1(double v) const {
                //=============================================================================================================
                        //return (1.0 - SmallSize_ <= v) && (v <= 1.0);
                        return fabs(1.0 - v) < SmallSize_;
                }
        }; // class Cell

        const double Cell::SmallSize_ = 1.0e-08;
} // namesapce kmd



//=============================================================================================================================
ZeroLevelSetDetector::ZeroLevelSetDetector(const Parameters* parameters)
//=============================================================================================================================
        : width_(parameters->image_width_),
          height_(parameters->image_height_),
          buffer_(0),
          visit_flags_(width_ * height_),
          first_x_(-1),
          first_y_(-1) {
}

//=============================================================================================================================
ZeroLevelSetDetector::~ZeroLevelSetDetector() {
//=============================================================================================================================
        delete_fronts();
}

//=============================================================================================================================
void ZeroLevelSetDetector::eliminate_current_front() {
//=============================================================================================================================
        delete_front(fronts_.back());
        fronts_.pop_back();
}


const int ZeroLevelSetDetector::OffsetX_[] = {1, 0, -1, -1, -1,  0,  1, 1};
const int ZeroLevelSetDetector::OffsetY_[] = {1, 1,  1,  0, -1, -1, -1, 0};

const int ZeroLevelSetDetector::FlagSizes_[] = {3, 1, 3, 1, 3, 1, 3, 1};
const int ZeroLevelSetDetector::Flags_[][3] = {
        {7, 0, 1},
        {1, 0, 0},
        {1, 2, 3},
        {3, 0, 0},
        {3, 4, 5},
        {5, 0, 0},
        {5, 6, 7},
        {7, 0, 0},
};


/*!
 *
 *      [2](-1,1)  [1](0,1)  [0](1,1)
 *      [3](-1,0)  target    [7](1,0)
 *      [4](-1,-1) [5](0,-1) [6](1,-1)
 */

//=============================================================================================================================
void ZeroLevelSetDetector::trace_front(Cell* cell) {
//=============================================================================================================================

        for (;;) {

                const int cur_x = cell->x();
                const int cur_y = cell->y();
                const int flag = cell->get_flag();
                const int flag_size = FlagSizes_[flag];
                const int* flags = Flags_[flag];

                int i, j = 0;
                for ( i = 0; i < flag_size; ++i ) {
                        j = flags[i];
                        bool includes = includes_zero_level_set(cell, OffsetX_[j], OffsetY_[j], cur_x, cur_y);
                        //assert(TestUtilities::print_any("includes", includes));
                        visit_now(cell->x(), cell->y());

                        if ( includes ) {
                                break;
                        }
                }

                // unable to make closed curve
                if ( i == flag_size ) {
                        eliminate_current_front();
                        break;
                }

                // congratulations!
                if ( (cell->x() == first_x_) && (cell->y() == first_y_) ) {
                        break;
                }

                Point p;
                if ( cell->find_next(p, j) ) {
                        fronts_.back().push_back(p.clone());
                } else {
                        // unable to make closed curve
                        eliminate_current_front();
                        //assert(TestUtilities::print_any("eliminate", ""));
                        break;
                }
        }
}

//=============================================================================================================================
bool ZeroLevelSetDetector::includes_zero_level_set(Cell* cell, int offset_x, int offset_y, int cur_x, int cur_y) {
//=============================================================================================================================
        int i = cur_x + offset_x;
        int j = cur_y + offset_y;
        //assert(TestUtilities::print_point(Point(i, j)));
        cell->set(i, j, buffer(i, j), buffer(i + 1, j), buffer(i + 1, j + 1), buffer(i, j + 1));
        //assert(TestUtilities::print_any("i,j", buffer(i, j)));
        //assert(TestUtilities::print_any("i+1,j", buffer(i+1, j)));
        //assert(TestUtilities::print_any("i+1,j+1", buffer(i+1, j+1)));
        //assert(TestUtilities::print_any("i,j+1", buffer(i, j+1)));
        return cell->includes_zero_level_set();
}

//=============================================================================================================================
void ZeroLevelSetDetector::delete_fronts() {
//=============================================================================================================================
        Fronts::iterator beg = fronts_.begin();
        Fronts::iterator end = fronts_.end();
        while ( beg != end ) {
                delete_front(*beg);
                ++beg;
        }
        fronts_.clear();
}

//=============================================================================================================================
void ZeroLevelSetDetector::run(const double* buffer, const Tube& tube, int tube_width) {
//=============================================================================================================================
        delete_fronts();
        clear_visit_flags();
        Cell cell;
        buffer_ = &buffer[0];
        trace_fronts(&cell, buffer, tube, tube_width);
}

//=============================================================================================================================
void ZeroLevelSetDetector::trace_fronts(Cell* cell, const double* buffer, const Tube& tube, int tube_width) {
//=============================================================================================================================
        int w = width_;
        Tube::const_iterator beg = tube.begin();
        Tube::const_iterator end = tube.end();
        while ( beg != end ) {
                int j = beg->first;
                int wj = w * j;
                const double* buf = &buffer[wj];
                const Segments& segments = beg->second;
                for ( std::size_t k = 0, n = segments.size(); k < n; ++k ) {
                        const Segment& segment = segments[k];
                        int first = segment.first;
                        int second = segment.second;
                        for ( int i = first; i <= second; ++i ) {
                                if ( !has_already_visited(i, j) ) {
                                        visit_now(i, j);
                                        if ( (buf[i + 1] != tube_width) && (buf[i + 1 + w] != tube_width) && (buf[i + w] != tube_width) ) {
                                                cell->set(i, j, buf[i], buf[i + 1], buf[i + 1 + w], buf[i + w]);
                                                if ( cell->includes_zero_level_set() ) {
                                                        Front front;
                                                        if ( cell->find_first_segment(front) ) {
                                                                fronts_.push_back(front);
                                                                first_x_ = i;
                                                                first_y_ = j;
                                                                trace_front(cell);
                                                        }
                                                }
                                        }
                                }
                        }
                }
                ++beg;
        }
}

#ifdef UNIT_TEST_ZeroLevelSetDetector
#define BOOST_TEST_DYN_LINK
//----------------------------------------------------------------------------------------------------------------------------
#include <kmd/without_warning/unit_test.h>
#include <kmd/without_warning/iostream.h>
#include "Utilities.h"
#include <fstream>
#include <cmath>
#include "Surface.h"
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>
#include "TestUtilities.h"
#include "ShortestDistanceFinder.h"
#include "TubeDomainGenerator.h"
#include "FrontDirection.h"
#include "DistanceInitializer.h"
//----------------------------------------------------------------------------------------------------------------------------
using namespace std;
using namespace boost;

namespace {

        double calculate_distance(double x, double y) {
                return Utilities::square(x) + Utilities::square(y);
        }
        double f(int x, int y, int c) {
                return sqrt(static_cast<double>(Utilities::square(x - c) + Utilities::square(y - c))) - 30;
        }


}

namespace kmd {
        class CellTester {
        public:
                static const double* test_get_vertices(Cell& obj) {
                        return obj.vertices_;
                }
        };



        class ZeroLevelSetDetectorTester {
        public:
                static bool test_delete_fronts(ZeroLevelSetDetector& obj) {
                        obj.delete_fronts();
                        return obj.fronts_.size() == 0u;
                }



                static bool test_clear_visit_flags(ZeroLevelSetDetector& obj) {
                        obj.clear_visit_flags();
                        for ( std::size_t i = 0, n = obj.visit_flags_.size(); i < n; ++i ) {
                                if ( obj.visit_flags_[i] == true ) {
                                        return false;
                                }
                        }
                        return true;
                }

                //static void test_set_buffer(ZeroLevelSetDetector& obj, double* buffer) {
                //        obj.buffer_ = buffer;
                //}


                static void test_push_back_front(ZeroLevelSetDetector& obj, Front& front) {
                        obj.fronts_.push_back(front);
                }

                static const std::vector<bool> get_visit_flags(ZeroLevelSetDetector& obj) {
                        return obj.visit_flags_;
                }


        };

        void read_data(Fronts& fronts, const char* front_file) {

                char_separator<char> sep(" ");
                typedef tokenizer<char_separator<char> > Tokenizer;

                ifstream ifs_f(front_file);
                string str;
                Front front;
                while ( getline(ifs_f, str) ) {
                        Tokenizer tokens(str, sep);
                        Tokenizer::iterator i = tokens.begin();
                        double x = lexical_cast<double>(*i);
                        ++i;
                        double y = lexical_cast<double>(*i);
                        front.push_back(new Point(x, y));
                }
                fronts.push_back(front);
        }

        void run(const Fronts& fronts, int w, int h, std::size_t n, std::size_t s, std::size_t t) {
                Parameters params;
                params.image_width_ = w;
                params.image_height_ = h;
                ShortestDistanceFinder finder(&params);
                int tube_width = 5;
                finder.set_tube_width(tube_width);
                finder.run(fronts);
                TestUtilities::write_image("finder", finder.get_shortest_distances(), w, h, 2);

                TubeDomainGenerator tube_generator(&finder);
                tube_generator.run();
                TestUtilities::write_tube("tube", tube_generator.get_tube());

                DistanceInitializer dist_initializer(&finder);
                FrontDirection direction(FrontDirection::Outer);
                vector<double> dist(w * h, direction.attach_sign(tube_width));
                dist_initializer.run(&dist[0], dist.size(), fronts, tube_generator.get_tube());
                TestUtilities::write_image("buf", dist, w, h, 1);

                ZeroLevelSetDetector detector(&params);
                detector.run(&dist[0], tube_generator.get_tube(), direction.attach_sign(tube_width));
                const Fronts& fronts2 = detector.get_fronts();
                BOOST_CHECK(fronts2.size() == n);
                BOOST_CHECK(s == fronts2[0].size());
                //ofstream ofs2("front2-0");
                //TestUtilities::write_front(ofs2, fronts2[0]);
                if ( fronts2.size() == 2 ) {
                        BOOST_CHECK(t == fronts2[1].size());
                        //ofstream ofs2("front2-1");
                        //TestUtilities::write_front(ofs2, fronts2[1]);
                }
        }

        void test_real_data(const char* front_file, std::size_t s) {
                Fronts fronts;
                int w = 269;
                int h = 228;
                read_data(fronts, front_file);
                run(fronts, w, h, 1, s, 0);
        }

        void test_real_data(const char* front_file0, const char* front_file1, std::size_t s, std::size_t t) {
                Fronts fronts;
                int w = 269;
                int h = 228;
                read_data(fronts, front_file0);
                read_data(fronts, front_file1);
                run(fronts, w, h, 2, s, t);
        }
}

//============================================================================================================================
BOOST_AUTO_TEST_CASE(TEST_ZeroLevelSetDetector) {
//============================================================================================================================

        cout << "ZeroLevelSetDetector\n";

        // Cell::Cell/x/y/get_flag/
        {
                Cell cell;
                BOOST_CHECK(cell.x() == 0);
                BOOST_CHECK(cell.y() == 0);
                BOOST_CHECK(cell.get_flag() == 0);
                const double* vertices = CellTester::test_get_vertices(cell);
                BOOST_CHECK(vertices[0] == 0);
                BOOST_CHECK(vertices[1] == 0);
                BOOST_CHECK(vertices[2] == 0);
                BOOST_CHECK(vertices[3] == 0);
        }

        // Cell::set/includes_zero_level_set/set_flag
        {
                Cell cell;
                cell.set(3, 2, 0.02, 0.12, 0.43, 0.124);
                BOOST_CHECK(cell.x() == 3);
                BOOST_CHECK(cell.y() == 2);

                const double* vs = CellTester::test_get_vertices(cell);
                BOOST_CHECK(vs[0] == 0.02);
                BOOST_CHECK(vs[1] == 0.12);
                BOOST_CHECK(vs[2] == 0.43);
                BOOST_CHECK(vs[3] == 0.124);
                cell.set_flag(45);
                BOOST_CHECK(cell.get_flag() == 45);

                BOOST_CHECK(!cell.includes_zero_level_set());
        }

        // Cell::find_next
        {
                const double vs[] = {-1.0, 2.0, -0.5, 1.0};

                LinearInterpolator interpolation;
                interpolation.set_value(vs[0], vs[1], vs[2], vs[3]);
                interpolation.run();
                BOOST_CHECK(interpolation.get_sign() > 0.0);

                Cell cell;
                cell.set(0, 0, vs[0], vs[1], vs[2], vs[3]);

                {
                        Point p;
                        BOOST_CHECK(cell.find_next(p, 7));
                        BOOST_CHECK_EQUAL(p.x(), interpolation.get_x_value(0.0));
                        BOOST_CHECK_EQUAL(p.y(), 0.0);
                }

                {
                        Point p;
                        BOOST_CHECK(cell.find_next(p, 1));
                        BOOST_CHECK_EQUAL(p.x(), 0.0);
                        BOOST_CHECK_EQUAL(p.y(), interpolation.get_y_value(0.0));
                }

                {
                        Point p;
                        BOOST_CHECK(cell.find_next(p, 3));
                        BOOST_CHECK_EQUAL(p.x(), interpolation.get_x_value(1.0));
                        BOOST_CHECK_EQUAL(p.y(), 1.0);
                }

                {
                        Point p;
                        BOOST_CHECK(cell.find_next(p, 5));
                        BOOST_CHECK_EQUAL(p.x(), 1.0);
                        BOOST_CHECK_EQUAL(p.y(), interpolation.get_y_value(1.0));
                }
        }

        // Cell::find_first_segment
        {
                const double vs[] = {10.0, 2.0, 0.5, -1.0};

                LinearInterpolator interpolation;
                interpolation.set_value(vs[0], vs[1], vs[2], vs[3]);
                interpolation.run();
                BOOST_CHECK(interpolation.get_sign() < 0.0);

                Cell cell;
                cell.set(0, 0, vs[0], vs[1], vs[2], vs[3]);

                Front front;
                BOOST_CHECK(cell.find_first_segment(front));
                BOOST_CHECK(cell.get_flag() == 3);
                const Point* starting = front[0];
                BOOST_CHECK_EQUAL(starting->x(), interpolation.get_x_value(1.0));
                BOOST_CHECK_EQUAL(starting->y(), 1.0);

                const Point* ending = front[1];
                BOOST_CHECK_EQUAL(ending->x(), 0.0);
                BOOST_CHECK_EQUAL(ending->y(), interpolation.get_y_value(0.0));

        }

        // Cell::find_first_segment
        {
                const double vs[] = {-1.0, 2.0, -0.5, 1.0};

                LinearInterpolator interpolation;
                interpolation.set_value(vs[0], vs[1], vs[2], vs[3]);
                interpolation.run();
                BOOST_CHECK(interpolation.get_sign() > 0.0);

                Cell cell;
                cell.set(0, 0, vs[0], vs[1], vs[2], vs[3]);

                {
                        Front front;
                        BOOST_CHECK(cell.find_first_segment(front));
                        BOOST_CHECK_EQUAL(cell.get_flag(), 7);
                        const Point* starting = front[0];
                        BOOST_CHECK_EQUAL(starting->x(), interpolation.get_x_value(1.0));
                        BOOST_CHECK_EQUAL(starting->y(), 1.0);

                        const Point* ending = front[1];
                        BOOST_CHECK_EQUAL(ending->x(), 1.0);
                        BOOST_CHECK_EQUAL(ending->y(), interpolation.get_y_value(1.0));
                }
        }

        {
                const double vs[] = {0, 2.22279, -130.324, -133.221};
                /*
                0: 50 20 0
                1: 51 20 2.22279
                2: 51 21 -130.324
                3: 50 21 -133.221
                */

                LinearInterpolator interpolation;
                interpolation.set_value(vs[0], vs[1], vs[2], vs[3]);
                interpolation.run();
                BOOST_CHECK(interpolation.get_sign() < 0.0);


        }

        // Cell::find_first_segment
        {
                const double vs[] = {10.0, 2.0, 0.5, -1.0};

                LinearInterpolator interpolation;
                interpolation.set_value(vs[0], vs[1], vs[2], vs[3]);
                interpolation.run();
                BOOST_CHECK(interpolation.get_sign() < 0.0);

                Cell cell;
                cell.set(0, 0, vs[0], vs[1], vs[2], vs[3]);

                {
                        Front front;
                        BOOST_CHECK(cell.find_first_segment(front));
                        BOOST_CHECK(cell.get_flag() == 3);
                        const Point* starting = front[0];
                        BOOST_CHECK_EQUAL(starting->x(), interpolation.get_x_value(1.0));
                        BOOST_CHECK_EQUAL(starting->y(), 1.0);

                        const Point* ending = front[1];
                        BOOST_CHECK_EQUAL(ending->x(), 0.0);
                        BOOST_CHECK_EQUAL(ending->y(), interpolation.get_y_value(0.0));
                }
        }

        {
                test_real_data("../data/unit_test_data/dist_29_0", 423);
                test_real_data("../data/unit_test_data/dist_100_0", 421);
                test_real_data("../data/unit_test_data/multi_fronts_0", "../data/unit_test_data/multi_fronts_1", 236, 236);
        }
}
#endif // UNIT_TEST_ZeroLevelSetDetector

