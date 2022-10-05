//----------------------------------------------------------------------------------------------------------------
#include "TubeDomainGenerator.h"
#include "ShortestDistanceFinder.h"
#include "Parameters.h"
//----------------------------------------------------------------------------------------------------------------
using namespace kmd;

//================================================================================================================
TubeDomainGenerator::TubeDomainGenerator(const ShortestDistanceFinder* finder)
//================================================================================================================
        : finder_(finder) {
}

// test ok
//================================================================================================================
void TubeDomainGenerator::run() {
//================================================================================================================
        tube_.clear();
        const Parameters* params = finder_->get_parameters();
        int w = params->image_width_;
        int h = params->image_height_;
        int tube_width = finder_->get_tube_width();
        const std::vector<double>& shortest_distances = finder_->get_shortest_distances();
        for ( int j = 0, wj = w * j; j < h; ++j, wj += w ) {
                const double* it = &shortest_distances[wj];
                for ( int i = 0; i < w; ++i ) {
                        if ( it[i] != tube_width ) {
                                if ( (i == 0) || (it[i - 1] == tube_width) ) {
                                        tube_[j].push_back(Segment(i, -1));
                                }

                                if ( i == w - 1 ) {
                                        tube_[j].back().second = i;
                                }
                        }
                        else { // it[i] == tube_width
                                if ( (i != 0) && (it[i - 1] != tube_width) ) {
                                        tube_[j].back().second = i - 1;
                                }
                        }
                }
        }
}

#ifdef UNIT_TEST_TubeDomainGenerator
#define BOOST_TEST_DYN_LINK
//----------------------------------------------------------------------------------------------------------------
#include <kmd/without_warning/unit_test.h>
#include <kmd/without_warning/iostream.h>
#include "ShortestDistanceFinder.h"
#include "ZeroLevelSetDetector.h"
#include "Utilities.h"
#include "TestUtilities.h"
//----------------------------------------------------------------------------------------------------------------
using namespace std;

namespace {

        void test0() {
                Front front;
                TestUtilities::create_front(front, 40, 1000, 100, 100);
                Fronts fronts;
                fronts.push_back(front);

                Parameters params;

                int w = 200;
                int h = 200;
                params.image_width_ = w;
                params.image_height_ = h;

                ShortestDistanceFinder finder(&params);
                int tube_width = 6;
                finder.set_tube_width(tube_width);
                finder.run(fronts);

                TubeDomainGenerator generator(&finder);
                generator.run();
                const Tube& tube = generator.get_tube();

                const vector<double>& dist = finder.get_shortest_distances();
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
                                       int index = i + wj;
                                       BOOST_CHECK(dist[index] < tube_width);
                                }
                        }
                        ++beg;
                }
        }
}

//================================================================================================================
BOOST_AUTO_TEST_CASE(TEST_TubeDomainGenerator) {
//================================================================================================================
        cout << "TubeDomainGenerator\n";

        test0();
}
#endif // UNIT_TEST_TubeDomainGenerator
