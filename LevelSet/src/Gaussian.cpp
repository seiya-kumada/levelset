//-----------------------------------------------------------------------------------------------------------------------------
#include "Gaussian.h"
#include "Pi.h"
//-----------------------------------------------------------------------------------------------------------------------------
using namespace kmd;

//=============================================================================================================================
Gaussian::Gaussian(double sigma)
//=============================================================================================================================
        : coeff0_(1.0 / ( 2.0 * Utilities::square(sigma) )),
          coeff1_(coeff0_ / M_PI) {
}

#ifdef UNIT_TEST_Gaussian
#define BOOST_TEST_DYN_LINK
//----------------------------------------------------------------------------------------------------------------------------
#include <kmd/without_warning/unit_test.h>
#include <kmd/without_warning/iostream.h>
//----------------------------------------------------------------------------------------------------------------------------
using namespace std;

//============================================================================================================================
BOOST_AUTO_TEST_CASE(TEST_Gaussian) {
//============================================================================================================================

        cout << "Gaussian\n";

        // get_gaussian_value
        {
                double sigma = 2.0;
                Gaussian gaussian(sigma);
                double range = 10;
                int step_number = 10;
                double step = range / step_number;
                double s = 0.0;
                for ( int j = -step_number; j <= step_number; ++j ) {
                        for ( int i = -step_number; i <= step_number; ++i ) {
                                s += gaussian.get_value(i * step, j * step) * Utilities::square(step);
                        }
                }
                BOOST_CHECK(fabs(s - 1.0) < 0.00001);
        }
}
#endif // UNIT_TEST_Gaussian


