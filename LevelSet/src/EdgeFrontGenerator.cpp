//-----------------------------------------------------------------------------------------------------------------
#include "EdgeFrontGenerator.h"
//-----------------------------------------------------------------------------------------------------------------
using namespace kmd;
using namespace std;

 //=================================================================================================================
EdgeFrontGenerator::EdgeFrontGenerator()
//=================================================================================================================
        : edge_width_(0),
          direction_() {

        is_inside_[FrontDirection::Inner] = &EdgeFrontGenerator::is_inside_in;
        is_inside_[FrontDirection::Outer] = &EdgeFrontGenerator::is_inside_out;
}

//=================================================================================================================
EdgeFrontGenerator::~EdgeFrontGenerator() {
//=================================================================================================================
}

//================================================================================================================
void EdgeFrontGenerator::set(
        const FrontDirection&   direction,
        int                     edge_width
) {
//================================================================================================================
        direction_ = direction;
        if ( direction_.is_inner() ) {
                edge_width_ = -edge_width;
        }
        else {
                edge_width_ = edge_width;
        }
}

//================================================================================================================
void EdgeFrontGenerator::run(const double* dst_buffer, const Tube& tube, int w) {
//================================================================================================================
        edge_front_.clear();
        Tube::const_iterator beg = tube.begin();
        Tube::const_iterator end = tube.end();

        while ( beg != end ) {
                int j = beg->first;
                int wj = w * j;
                const double* it = &dst_buffer[wj];
                const Segments& segments = beg->second;
                for ( std::size_t k = 0, n = segments.size(); k < n; ++k ) {
                        const Segment& segment = segments[k];
                        int first = segment.first;
                        int second = segment.second;
                        for ( int i = first; i <= second; ++i ) {
                                if ( is_inside(it[i] - edge_width_) ) {
                                        edge_front_.push_back(IntPoint(i, j));
                                }
                        }
                }
                ++beg;
        }
}

#ifdef UNIT_TEST_EdgeFrontGenerator
#define BOOST_TEST_DYN_LINK
//----------------------------------------------------------------------------------------------------------------------------
#include <kmd/without_warning/unit_test.h>
#include <iostream>
#include <fstream>
#include "Parameters.h"
#include "ShortestDistanceFinder.h"
#include "TubeDomainGenerator.h"
#include "DistanceInitializer.h"
#include "FrontDirection.h"
#include "Point.h"
#include "TestUtilities.h"
//----------------------------------------------------------------------------------------------------------------------------
using namespace kmd;

namespace kmd {

        //********************************************************************************************************************
        class EdgeFrontGeneratorTester {
        //********************************************************************************************************************
        public:
                static bool test_is_inside_in(const EdgeFrontGenerator& obj, double x) {
                        return obj.is_inside_in(x);
                }

                static bool test_is_inside_out(const EdgeFrontGenerator& obj, double x) {
                        return obj.is_inside_out(x);
                }

                static FrontDirection::Direction get_value(const EdgeFrontGenerator& obj) {
                        return obj.direction_.get_value();
                }

                static int get_off_set(const EdgeFrontGenerator& obj) {
                        return obj.edge_width_;
                }


        };
}

namespace {

        void test1() {
                const double Pi = 3.1415926535897932384;
                int r = 25;
                vector<Point*> front;
                int size = 1000;
                double step = 2 * Pi / size;
                int cx = 50;
                int cy = 50;
                for ( int i = 0; i < size; ++i ) {
                        front.push_back(new Point(cx + r * cos(i * step), cy + r * sin(i * step)));
                }
                Fronts fronts;
                fronts.push_back(front);
                TestUtilities::write_fronts("fronts", fronts);

                Parameters params;
                int w = 100;
                int h = 100;
                params.image_width_ = w;
                params.image_height_ = h;

                ShortestDistanceFinder finder(&params);
                int tube_width = 6;
                finder.set_tube_width(tube_width);
                finder.run(fronts);
                const vector<double>& shortest_distances = finder.get_shortest_distances();
                TestUtilities::write_image("shortest_distances", shortest_distances, w, h);

                TubeDomainGenerator tube_generator(&finder);
                tube_generator.run();
                const Tube& tube = tube_generator.get_tube();
                TestUtilities::write_tube("tube", tube);

                DistanceInitializer initializer(&finder);
                FrontDirection direction(FrontDirection::Outer);
                vector<double> buf(w * h, direction.attach_sign(tube_width));
                initializer.run(&buf[0], buf.size(), fronts, tube);
                TestUtilities::write_image("buf", buf, w, h);

                EdgeFrontGenerator edge_generator;
                edge_generator.set(direction, tube_width - 2);
                edge_generator.run(buf, tube, w);
                const vector<IntPoint>& edge = edge_generator.get_edge_front();
                TestUtilities::write_edge_front("edge", edge);

                for ( std::size_t i = 0, n = front.size(); i < n; ++i ) {
                        delete front[i];
                }

                front.clear();
        }

}

//============================================================================================================================
BOOST_AUTO_TEST_CASE(TEST_EdgeFrontGenerator) {
//============================================================================================================================
        cout << "EdgeFrontGenerator\n";
        FrontDirection inner(FrontDirection::Inner);
        FrontDirection outer(FrontDirection::Outer);

        // is_inside_in/out
        {
                EdgeFrontGenerator obj;
                BOOST_CHECK(EdgeFrontGeneratorTester::test_is_inside_in(obj, -1.0));
                BOOST_CHECK(EdgeFrontGeneratorTester::test_is_inside_out(obj, 1.0));
        }

        // set
        {
                EdgeFrontGenerator obj;
                obj.set(inner, 6);
                BOOST_CHECK(obj.is_inside(-1.0));
                obj.set(outer, 6);
                BOOST_CHECK(obj.is_inside(1.0));
        }

        // set
        {
                EdgeFrontGenerator obj;
                obj.set(outer, 6);
                BOOST_CHECK(FrontDirection::Outer == EdgeFrontGeneratorTester::get_value(obj));
                BOOST_CHECK(6 == EdgeFrontGeneratorTester::get_off_set(obj));
                obj.set(inner, 6);
                BOOST_CHECK(FrontDirection::Inner == EdgeFrontGeneratorTester::get_value(obj));
                BOOST_CHECK(-6 == EdgeFrontGeneratorTester::get_off_set(obj));

        }


        //test1();
}
#endif // UNIT_TEST_EdgeFrontGenerator

