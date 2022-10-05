//-----------------------------------------------------------------------------------------------------------------------------
#ifndef LevelSetMethodH
#define LevelSetMethodH
//-----------------------------------------------------------------------------------------------------------------------------
#include <boost/noncopyable.hpp>
#include <boost/gil/typedefs.hpp>
#include "Typedef.h"
#include <list>
//-----------------------------------------------------------------------------------------------------------------------------

namespace kmd {

        struct Parameters;
        class FrontDirection;

        //*********************************************************************************************************************
        class LevelSetMethod : private boost::noncopyable {
        //*********************************************************************************************************************
        public:
                /*!
                 *      @param[in]      parameters      given parameters
                 */
                explicit LevelSetMethod(const Parameters* parameters);

                ~LevelSetMethod();

                /*!
                 *      @param[in       image_buffer
                 *      @param[in]      front           an initial front given by an user
                 *      @param[in]      folder_path     folder path to save resultants
                 *      @param[in]      src_view        source view
                 *      @param[in]      color           color to draw front
                 *      @param[in]      interval        interval to save images
                 */
                void run(
                        const std::vector<double>&              image_buffer,
                        const Front&                            front,
                        const std::string&                      folder_path,
                        const boost::gil::gray8c_view_t&        src_view,
                        //kumada
                        //boost::gil::bits8                       color,
                        std::uint8_t                            color,
                        int                                     interval
                );

        private:
                const Parameters*               parameters_;
                ShortestDistanceFinderPtr       distance_finder_;
                SpeedGeneratorPtr               speed_generator_;
                DistanceInitializerPtr          distance_initializer_;
                ZeroLevelSetDetectorPtr         zero_level_set_detector_;
                EdgeFrontGeneratorPtr           edge_front_generator_;
                TubeDomainGeneratorPtr          tube_domain_generator_;
                FrontStopperPtr                 front_stopper_;

                bool is_inside_tube(const double* next_distance_ptr, int w) const;
                void preprocess(
                        const std::vector<double>&      image_buffer,
                        const Fronts&                   front,
                        int                             tube_width,
                        double*                         curr_distance_ptr,
                        int                             buffer_size,
                        const FrontDirection&           direction,
                        int                             w
                ) const;


        };
}
#endif // LevelSetMethodH
