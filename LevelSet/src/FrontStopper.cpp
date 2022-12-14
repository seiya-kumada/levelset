//------------------------------------------------------------------------------------------------------------------------
#include "FrontStopper.h"
//------------------------------------------------------------------------------------------------------------------------
using namespace kmd;
using namespace std;

//========================================================================================================================
void FrontStopper::run(
        vector<Front>&          final_fronts,
        const Fronts            moving_fronts,
        const vector<double>&   dist,
        const SpeedGenerator*   speed_generator
) const {
//========================================================================================================================
        
}

#ifdef UNIT_TEST_FrontStopper
#define BOOST_TEST_DYN_LINK
//----------------------------------------------------------------------------------------------------------------------------
#include <kmd/without_warning/unit_test.h>
#include <kmd/without_warning/iostream.h>
//----------------------------------------------------------------------------------------------------------------------------

//============================================================================================================================
BOOST_AUTO_TEST_CASE(TEST_FrontStopper) {
//============================================================================================================================

        cout << "FrontStopper\n";
}
#endif // UNIT_TEST_FrontStopper         
