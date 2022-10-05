//-----------------------------------------------------------------------------------------------------------------------------
#ifndef EVOLUTIONEQUATIONSOLVER_H_
#define EVOLUTIONEQUATIONSOLVER_H_
//-----------------------------------------------------------------------------------------------------------------------------
#include <boost/noncopyable.hpp>
#include "Typedef.h"
//-----------------------------------------------------------------------------------------------------------------------------
namespace kmd {

        struct Parameters;

        //*********************************************************************************************************************
        class EvolutionEquationSolver : private boost::noncopyable {
        //*********************************************************************************************************************
        private:
                EvolutionEquationSolver();
                ~EvolutionEquationSolver();

        public:
                /*!
                 *      This method solves the differential equation.
                 *
                 *      @param[out]     dst_dist        updated distance function with t = n + 1
                 *      @param[in]      src_dist        distance function with t = n
                 *      @param[in]      speed_generator
                 *      @param[in]      parameters      paramters in relation to the solving a problem
                 *      @param[in]      tube            tube domain which includes all points to be evolved
                 */
                static void evolve(
                        double*                 dst_dist,
                        const double*           src_dist,
                        SpeedGenerator*         speed_generator,
                        const Parameters*       parameters,
                        const Tube&             tube
                );
        };
}
#endif /*EVOLUTIONEQUATIONSOLVER_H_*/
