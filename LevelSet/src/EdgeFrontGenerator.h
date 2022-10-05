//-----------------------------------------------------------------------------------------------------------------
#ifndef EDGEFRONTGENERATOR_H_
#define EDGEFRONTGENERATOR_H_
//-----------------------------------------------------------------------------------------------------------------
#include "Typedef.h"
#include <boost/noncopyable.hpp>
#include "FrontDirection.h"
//-----------------------------------------------------------------------------------------------------------------

namespace kmd {

        //**********************************************************************************************************
        class EdgeFrontGenerator : private boost::noncopyable {
        //**********************************************************************************************************
                friend class EdgeFrontGeneratorTester;

        public:
                EdgeFrontGenerator();
                ~EdgeFrontGenerator();

                /*!
                 *      @param[in]      direction
                 *      @param[in]      edge_width
                 */
                void set(
                        const FrontDirection&   direction,
                        int                     edge_width
                );

                const std::vector<IntPoint>& get_edge_front() const {
                        return edge_front_;
                }

                /*!
                 *      @param[in]      x
                 */
                bool is_inside(double x) const {
                        return (this->*is_inside_[direction_.get_value()])(x);
                }

                /*!
                 *      @param[in]      buf
                 *      @param[in]      tube
                 *      @param[in]      w
                 */
                void run(
                        const std::vector<double>&      buf,
                        const Tube&                     tube,
                        int                             w
                ) {
                        run(&buf[0], tube, w);
                }

                /*!
                 *      @param[in]      buf
                 *      @param[in]      tube
                 *      @param[in]      w
                 */
                void run(
                        const double*   buf,
                        const Tube&     tube,
                        int             w
                );

        private:
                std::vector<IntPoint>   edge_front_;
                int                     edge_width_;
                FrontDirection          direction_;

                typedef bool (EdgeFrontGenerator::*IsInside)(double x) const;
                IsInside is_inside_[2];

                bool is_inside_in(double x) const {
                        return x < 0.0;
                }

                bool is_inside_out(double x) const {
                        return x > 0.0;
                }
        };
}
#endif /*EDGEFRONTGENERATOR_H_*/
