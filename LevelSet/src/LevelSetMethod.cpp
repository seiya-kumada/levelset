//-----------------------------------------------------------------------------------------------------------------------------
#include "LevelSetMethod.h"
#include "SpeedGenerator.h"
#include "EvolutionEquationSolver.h"
#include "DistanceInitializer.h"
#include "ShortestDistanceFinder.h"
#include "ZeroLevelSetDetector.h"
#include "Parameters.h"
#include "EdgeFrontGenerator.h"
#include "TubeDomainGenerator.h"
#include "FrontDirection.h"
#include "TestUtilities.h"
#include <boost/progress.hpp>
#include "FrontStopper.h"
//-----------------------------------------------------------------------------------------------------------------------------
using namespace kmd;
using namespace std;
using namespace boost;
using namespace boost::gil;

//=============================================================================================================================
LevelSetMethod::LevelSetMethod(const Parameters* parameters)
//=============================================================================================================================
        : parameters_(parameters),
          distance_finder_(new ShortestDistanceFinder(parameters_)),
          speed_generator_(new SpeedGenerator(parameters)),
          distance_initializer_(new DistanceInitializer(distance_finder_.get())),
          zero_level_set_detector_(new ZeroLevelSetDetector(parameters_)),
          edge_front_generator_(new EdgeFrontGenerator),
          tube_domain_generator_(new TubeDomainGenerator(distance_finder_.get())),
          front_stopper_(new FrontStopper) {
}

//=============================================================================================================================
LevelSetMethod::~LevelSetMethod() {
//=============================================================================================================================
}

//=============================================================================================================================
void LevelSetMethod::preprocess(
        const vector<double>&   image_buffer,
        const Fronts&      fronts,
        int                     tube_width,
        double*                 curr_distance_ptr,
        int                     buffer_size,
        const FrontDirection&   direction,
        int                     w
) const {
//=============================================================================================================================

        // In advance, the source image is modified by applying Gaussian and gradient filters to it,
        // and the damping factor for the speed is calcuated over all pixels.
        speed_generator_->calculate_image_dependent_factors(image_buffer);
        // for debug
        //TestUtilities::write_image("factors", speed_generator_->get_image_dependent_factors(), parameters_->image_width_, parameters_->image_height_);

        // The shortest distances from the initial front are calcuated.
        distance_finder_->set_tube_width(tube_width);
        distance_finder_->run(fronts);
        // for debug
        //TestUtilities::write_image("finder", distance_finder_->get_shortest_distances(), parameters_->image_width_, parameters_->image_height_);

        // 2012/11/13
        // The tube domain is specified.
        // This method should be called after ShotestDistanceFinder::run.
        tube_domain_generator_->run();
        // for debug
        //TestUtilities::write_tube("tube", tube_domain_generator_->get_tube());

        // The initial signed distance function is calculated.
        distance_initializer_->run(curr_distance_ptr, buffer_size, fronts, tube_domain_generator_->get_tube());

        // The points to know the timing to reinitilize fronts are set.
        int edge_width = tube_width - 1;
        edge_front_generator_->set(direction, edge_width);
        edge_front_generator_->run(curr_distance_ptr, tube_domain_generator_->get_tube(), w);
        // for debug
        //TestUtilities::write_edge_front("edge_front", edge_front_generator_->get());
}

//=============================================================================================================================
void LevelSetMethod::run(
        const vector<double>&   image_buffer,
        const Front&            front,
        const string&           folder_path,
        const gray8c_view_t&    src_view,
        //kumada
        //bits8                   color,
        std::uint8_t            color,
        int                     interval
) {
//=============================================================================================================================
        const int w = parameters_->image_width_;
        const int h = parameters_->image_height_;
        const std::size_t buffer_size = w * h;
        const int tube_width = 6;
        FrontDirection::Direction fd = parameters_->constant_speed_ > 0 ? FrontDirection::Outer : FrontDirection::Inner;
        FrontDirection direction(fd);
        const int signed_tube_width = direction.attach_sign(tube_width);

        vector<double> distance_buffer0(buffer_size, signed_tube_width);
        vector<double> distance_buffer1(buffer_size, signed_tube_width);
        double* curr_distance_ptr = &distance_buffer0[0];
        double* next_distance_ptr = &distance_buffer1[0];

        Fronts fronts;
        fronts.push_back(front);

        // for debug
        //TestUtilities::write_fronts("initial_fronts", fronts);

        preprocess(image_buffer, fronts, tube_width, curr_distance_ptr, static_cast<int>(buffer_size), direction, w);
        // for debug
        //TestUtilities::write_image("dist", distance_buffer0, w, h);

        int count = 1;
        int upper_count = parameters_->time_step_number_;
        const Fronts* ptr_fronts = &fronts;
        progress_display pd(upper_count);

        while ( upper_count >= count ) {
                
                // develope the current buffer into next buffer
                EvolutionEquationSolver::evolve(
                        next_distance_ptr,
                        curr_distance_ptr,
                        speed_generator_.get(),
                        parameters_,
                        tube_domain_generator_->get_tube()
                );

                // using the next buffer, obtain the new fronts.
                zero_level_set_detector_->run(next_distance_ptr, tube_domain_generator_->get_tube(), signed_tube_width);
                ptr_fronts = &zero_level_set_detector_->get_fronts();

                if ( count % interval == 0 ) {
                        TestUtilities::write_intermediate_fronts(src_view, folder_path, count, *ptr_fronts, color);
                        //TestUtilities::write_next_distance(folder_path, count, next_distance_ptr, distance_buffer0, distance_buffer1, w, h, 1);
                }

                // see whether the front is inside the tube or not
                if ( !is_inside_tube(next_distance_ptr, w) ) {
                        // for debug
                        //TestUtilities::write_fronts("fronts", *ptr_fronts);
                        //TestUtilities::write_edge_front("b_edge_front", edge_front_generator_->get());
                        //TestUtilities::write_tube("b_tube", tube_domain_generator_->get_tube());


                        distance_finder_->run(*ptr_fronts);
                        tube_domain_generator_->run();
                        distance_initializer_->run(
                                next_distance_ptr,
                                curr_distance_ptr,
                                buffer_size,
                                tube_domain_generator_->get_tube(),
                                direction
                        );
                        edge_front_generator_->run(next_distance_ptr, tube_domain_generator_->get_tube(), w);

                        // for debug
                        //TestUtilities::write_edge_front("a_edge_front", edge_front_generator_->get());
                        //TestUtilities::write_tube("a_tube", tube_domain_generator_->get_tube());
                        //return;
                }

                std::swap(curr_distance_ptr, next_distance_ptr);
                ++count;
                ++pd;
        }
}

//=============================================================================================================================
bool LevelSetMethod::is_inside_tube(const double* next_distance_ptr, int w) const {
//=============================================================================================================================
        const vector<IntPoint>& edge_front = edge_front_generator_->get_edge_front();
        for ( std::size_t k = 0, n = edge_front.size(); k < n; ++k ) {
                const IntPoint& p = edge_front[k];
                int i = static_cast<int>(p.x);
                int j = static_cast<int>(p.y);
                int wj = w * j;
                if ( !edge_front_generator_->is_inside(next_distance_ptr[wj + i]) ) {
                        return false;
                }
        }
        return true;
}
