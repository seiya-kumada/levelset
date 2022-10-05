//-----------------------------------------------------------------------------------------------------------------------------
#ifndef GAUSSIAN_H_
#define GAUSSIAN_H_
//-----------------------------------------------------------------------------------------------------------------------------
#include <boost/noncopyable.hpp>
#include <cmath>
#include "Utilities.h"
//-----------------------------------------------------------------------------------------------------------------------------

namespace kmd {

        /*!
         *      This class provides Gaussian on the 2 dimensional plane.
         */
        //*********************************************************************************************************************
        class Gaussian : private boost::noncopyable {
        //*********************************************************************************************************************
        public:
                /*!
                 *      @param[in]      x       x value
                 *      @param[in]      y       y value
                 *      @return                 Gaussian value
                 */
                double get_value(double x, double y) const {
                        return  coeff1_ * exp(-coeff0_ * ( Utilities::square(x) + Utilities::square(y) ));
                }

                /*!
                 *      @param[in]      sigma
                 */
                explicit Gaussian(double sigma);
                ~Gaussian() {}

        private:
                double  coeff0_;
                double  coeff1_;
        };
}
#endif /*GAUSSIAN_H_*/
