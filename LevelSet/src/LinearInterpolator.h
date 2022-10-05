//-----------------------------------------------------------------------------------------------------------------------------
#ifndef LINEARINTERPOLATOR_H_
#define LINEARINTERPOLATOR_H_
//-----------------------------------------------------------------------------------------------------------------------------
#include <boost/noncopyable.hpp>
//-----------------------------------------------------------------------------------------------------------------------------
namespace kmd {

        //*********************************************************************************************************************
        class LinearInterpolator : private boost::noncopyable {
        //*********************************************************************************************************************
        public:
                LinearInterpolator();
                ~LinearInterpolator();

                /*!
                 *      This method sets four values on the vertices.
                 * 
                 *      @param[in]      a0
                 *      @param[in]      a1
                 *      @param[in]      a2
                 *      @param[in]      a3
                 */
                void set_value(double a0, double a1, double a2, double a3);

                void run();

                /*!
                 *      @param[in]      x
                 */
                double get_y_value(double x) const {
                        return (-a0_ + x * b2_) / (b0_ + x * b1_);
                }

                /*!
                 *      @param[in]      y
                 */
                double get_x_value(double y) const {
                        return (-a0_ - y * b0_) / (-b2_ + y * b1_);
                }

                double get_sign() const {
                        return -(b0_ * b2_ / (b1_ * b1_) + a0_ / b1_);
                }

        private:
                double  a0_;
                double  a1_;
                double  a2_;
                double  a3_;

                double  b0_;
                double  b1_;
                double  b2_;
        };
}
#endif /*LINEARINTERPOLATOR_H_*/
