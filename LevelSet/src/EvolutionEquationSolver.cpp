//-----------------------------------------------------------------------------------------------------------------------------
#include "EvolutionEquationSolver.h"
#include "Parameters.h"
#include "Utilities.h"
#include "SpeedGenerator.h"
//-----------------------------------------------------------------------------------------------------------------------------
using namespace kmd;
using namespace std;

namespace {

        // test ok
        //=====================================================================================================================
        inline double max_zero(double x) {
        //=====================================================================================================================
                return std::max(x, 0.0);
        }

        // test ok
        //=====================================================================================================================
        inline double min_zero(double x) {
        //=====================================================================================================================
                return std::min(x, 0.0);
        }

        // test ok
        //=====================================================================================================================
        inline double nabla_plus(int w, const double* src) {
        //=====================================================================================================================
                return std::sqrt(
                          Utilities::square( max_zero(src[0] - src[-1]) )
                        + Utilities::square( min_zero(src[1] - src[ 0]) )
                        + Utilities::square( max_zero(src[0] - src[-w]) )
                        + Utilities::square( min_zero(src[w] - src[ 0]) )
                );
        }

        // test ok
        //=====================================================================================================================
        inline double nabla_minus(int w, const double* src) {
        //=====================================================================================================================
                return std::sqrt(
                          Utilities::square( min_zero(src[0] - src[-1]) )
                        + Utilities::square( max_zero(src[1] - src[ 0]) )
                        + Utilities::square( min_zero(src[0] - src[-w]) )
                        + Utilities::square( max_zero(src[w] - src[ 0]) )
                );
        }

        // test ok
        //=====================================================================================================================
        inline double calculate_constant_part(int w, const double* src, double const_speed) {
        //=====================================================================================================================
                if ( const_speed > 0.0 ) {
                        return const_speed * nabla_plus(w, src);
                } else if ( const_speed < 0.0 ) {
                        return const_speed * nabla_minus(w, src);
                } else {
                        return 0.0;
                }
        }

        // test ok
        //=====================================================================================================================
        inline double calculate_dependent_part(int w, const double* src, double dependent_speed) {
        //=====================================================================================================================
                double d = sqrt( Utilities::square(src[1] - src[-1]) + Utilities::square(src[w] - src[-w]) ) * 0.5;
                if ( d == 0.0 ) {
                        // for getting rid of the singularity of "dependent_speed"
                        return 0.0;
                } else {
                        return dependent_speed * d;
                }
        }

        // test ok
        //=====================================================================================================================
        inline double calculate_spatial_differential_parts(int w, const double* src, double cs, double ds) {
        //=====================================================================================================================
                return calculate_constant_part(w, src, cs) + calculate_dependent_part(w, src, ds);
        }

}

//=============================================================================================================================
void EvolutionEquationSolver::evolve(
        double*                 dst_dist,
        const double*           src_dist,
        SpeedGenerator*         speed_generator,
        const Parameters*       parameters,
        const Tube&             tube
) {
//=============================================================================================================================
        int w = parameters->image_width_;
        double dt = parameters->time_step();
        double dx = parameters->space_step_;
        double d = dt / dx;

        const vector<double>& factors = speed_generator->get_image_dependent_factors();
        const double constant_speed = parameters->constant_speed_;

        Tube::const_iterator beg = tube.begin();
        Tube::const_iterator end = tube.end();
        while ( beg != end ) {
                const int wj = w * beg->first;
                double* dst_it = dst_dist + wj;
                const double* src_it = src_dist + wj;
                const double* fac = &factors[wj];
                const Segments& segments = beg->second;

                for ( std::size_t k = 0, n = segments.size(); k < n; ++k ) {
                        const Segment& segment = segments[k];
                        const int first = segment.first;
                        const int second = segment.second;
                        for ( int i = first; i <= second; ++i ) {
                                dst_it[i] = src_it[i] - d
                                        * calculate_spatial_differential_parts(
                                                w,
                                                &src_it[i],
                                                fac[i] * constant_speed,
                                                fac[i] * speed_generator->calculate_curvature_dependent_speed(&src_it[i], w)
                                        );
                        }

                }
                ++beg;
        }
}

#ifdef UNIT_TEST_EvolutionEquationSolver
#define BOOST_TEST_DYN_LINK
//----------------------------------------------------------------------------------------------------------------------------
#include <kmd/without_warning/unit_test.h>
#include <kmd/without_warning/iostream.h>
//----------------------------------------------------------------------------------------------------------------------------

//============================================================================================================================
BOOST_AUTO_TEST_CASE(TEST_EvolutionEquationSolver) {
//============================================================================================================================

        cout << "EvolutionEquationSolver\n";
        BOOST_CHECK(max_zero(3.0) == 3.0);
        BOOST_CHECK(min_zero(-3.0) == -3.0);
        BOOST_CHECK(Utilities::square(4.0) == 16.0);

        double buf[] = {
                1, 2, 3,
                4, 5, 6,
                7, 8, 9,
        };

        // nabla_plus/nabla_minus
        {
                double* src = &buf[4];
                BOOST_CHECK_EQUAL(nabla_plus(3, src), sqrt(10.0));
                BOOST_CHECK_EQUAL(nabla_minus(3, src), sqrt(10.0));
        }

        // calculate_constant_part
        {
                double* src = &buf[4];
                double const_speed = 3.0;
                BOOST_CHECK_EQUAL(calculate_constant_part(3, src, const_speed), 3.0 * sqrt(10.0));
        }

        // calculate_dependent_part(int i, int w, const double* src, double dependent_speed)
        {
                double* src = &buf[4];
                double dependent_speed = 1.0;
                BOOST_CHECK_EQUAL(calculate_dependent_part(3, src, dependent_speed), sqrt(10.0));
        }

        // calculate_spatial_differential_parts
        {
                double* src = &buf[4];
                double const_speed = 3.0;
                double dependent_speed = 1.0;
                BOOST_CHECK_EQUAL(calculate_spatial_differential_parts(3, src, const_speed, dependent_speed), 4.0 * sqrt(10.0));
        }
}
#endif // UNIT_TEST_EvolutionEquationSolver
