//------------------------------------------------------------------------------------------------------------------------
#ifndef FrontStopperH
#define FrontStopperH
//------------------------------------------------------------------------------------------------------------------------
#include <boost/noncopyable.hpp>
#include <vector>
#include <list>
#include "Typedef.h"
//------------------------------------------------------------------------------------------------------------------------
namespace kmd {

        class SpeedGenerator;

        //****************************************************************************************************************
        class FrontStopper : private boost::noncopyable {
        //****************************************************************************************************************
        public:

                /*!
                 *      @param[out]     final_fronts    a set of resultant fronts
                 *      @param[in]      moving_fronts   a set of moving fronts
                 *      @param[in]      dist            signed distance function
                 *      @param[in]      speed_generator 
                 */
                void run(
                        std::vector<Front>&             final_fronts,
                        const Fronts               moving_fronts,
                        const std::vector<double>&      dist,
                        const SpeedGenerator*           speed_generator
                ) const;
        };
}
#endif // FrontStopperH

