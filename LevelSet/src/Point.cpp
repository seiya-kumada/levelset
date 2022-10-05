#ifdef UNIT_TEST_Point
#define BOOST_TEST_DYN_LINK
//-----------------------------------------------------------------------------------------------------------------------------
#include "Point.h"
#include <kmd/without_warning/unit_test.h>
#include <kmd/without_warning/iostream.h>
#include <numeric>
#include <cmath>
//-----------------------------------------------------------------------------------------------------------------------------
using namespace std;
using namespace kmd;

//=============================================================================================================================
BOOST_AUTO_TEST_CASE(TEST_Point) {
//=============================================================================================================================

        cout << "Point\n";

        // constructor
        {
                Point a;
                BOOST_CHECK(a.x() == 0.0 && a.y() == 0.0);
        }

        // constructor
        {
                Point a(3.43, 5.2);
                BOOST_CHECK(a.x() == 3.43 && a.y() == 5.2);
        }

        // set_x/set_y
        {
                Point a;
                a.set_x(3.2);
                BOOST_CHECK(a.x() == 3.2);
                a.set_y(5.43);
                BOOST_CHECK(a.y() == 5.43);
        }

        // x()/y()
        {
                Point a(4.2, 45.5);
                Point b(32.5, 89.2);
                a += b;
                BOOST_CHECK(a.x() == 36.7 && a.y() == 134.7);
                a -= b;
                BOOST_CHECK(fabs(a.x() - 4.2) < 1.0e-14);
                BOOST_CHECK(fabs(a.y() - 45.5) < 1.0e-13);
        }

        // operator+
        {
                Point a(1.2, 3.4);
                Point b(3.98, 12.64);
                Point c = a + b;
                BOOST_CHECK(a.x() == 1.2 && a.y() == 3.4);
                BOOST_CHECK(b.x() == 3.98 && b.y() == 12.64);
                BOOST_CHECK(c.x() == 5.18 && c.y() == 16.04);
        }

        // operator-
        {
                Point a(1.2, 3.4);
                Point b(3.98, 12.64);
                Point c = a - b;
                BOOST_CHECK(a.x() == 1.2 && a.y() == 3.4);
                BOOST_CHECK(b.x() == 3.98 && b.y() == 12.64);
                BOOST_CHECK(fabs(c.x() + 2.78) < 1.0e-15);
                BOOST_CHECK(fabs(c.y() + 9.24) < 1.0e-16);
        }

        // operator*
        {
                Point c(4.2, 3.1);
                Point d = 5 * c;
                BOOST_CHECK(d.x() == 21);
                BOOST_CHECK(d.y() == 15.5);
                BOOST_CHECK(c.x() == 4.2);
                BOOST_CHECK(c.y() == 3.1);
        }

        // operator/
        {
                Point c(4.2, 3.1);
                Point d = c / 5.0;
                BOOST_CHECK(d.x() == 4.2 / 5.0);
                BOOST_CHECK(d.y() == 3.1 / 5.0);
                BOOST_CHECK(c.x() == 4.2);
                BOOST_CHECK(c.y() == 3.1);
        }

        // cross_product/inner_product
        {
                Point a(1.0, 0.0);
                Point b(0.0, 1.0);
                BOOST_CHECK(1.0 == cross_product(a, b));
                BOOST_CHECK(0.0 == inner_product(a, b));
                BOOST_CHECK(cross_product(Point(1, 0), Point(1, 1)) > 0.0);
                BOOST_CHECK(cross_product(Point(1, 1), Point(1, 0)) < 0.0);
        }

        // operator==/operator!=
        {
                Point a(1.0, 3.3);
                Point b(1.0, 3.3);
                Point c(2.3, 5.4);
                BOOST_CHECK(a == b);
                BOOST_CHECK(a != c);
        }

}
#endif // UNIT_TEST_Point
