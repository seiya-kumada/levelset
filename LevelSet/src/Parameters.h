//-----------------------------------------------------------------------------------------------------------------------------
#ifndef PARAMETERS_H_
#define PARAMETERS_H_
//-----------------------------------------------------------------------------------------------------------------------------
#include <boost/noncopyable.hpp>
#include "Error.h"
#include <stdexcept>
//-----------------------------------------------------------------------------------------------------------------------------

namespace kmd {

        //*********************************************************************************************************************
        struct Parameters : private boost::noncopyable {
        //*********************************************************************************************************************
        private:
                double  time_step_;

        public:
                int     time_step_number_;

                //! spacial step
                double  space_step_;

                //! image size
                int     image_width_;
                int     image_height_;

                //! A Front speed is represented by "constanct_speed_ + epsilon_ * curvature".
                double  constant_speed_;

                //! value multiplied by a curvature
                double  epsilon_;

                //! deviation of gaussian
                double  sigma_;

                Parameters();

                double time_step() const {
                        if ( time_step_ == 0.0 ) {
                                throw std::runtime_error(Error::message("invalid time step", __LINE__, __FILE__));
                        }
                        return time_step_;
                }

                /*!
                 *      @param[in]      time_step
                 */
                void set_time_step(double time_step) {
                        time_step_ = time_step;
                }
        };

}
#endif /*PARAMETERS_H_*/
