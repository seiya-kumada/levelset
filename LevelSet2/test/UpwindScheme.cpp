#if(UNIT_TEST_UpwindScheme)
#define BOOST_TEST_DYN_LINK

#include "UpwindScheme.h"
#include <boost/test/unit_test.hpp>
using namespace lsm;

namespace lsm
{

class UpwindSchemeTester
{
public:
        template<typename D>
        static void set_position(UpwindScheme<D>& scheme, const IntPoint<D>& p)
        {
                scheme.set_position(p);
        }

        template<typename D>
        static void calculate_with_2d_positive(UpwindScheme<D>& scheme)
        {
                scheme.calculate_with_2d(PositiveSpeed());
        }

        template<typename D>
        static void calculate_with_2d_negative(UpwindScheme<D>& scheme)
        {
                scheme.calculate_with_2d(NegativeSpeed());
        }

        template<typename D>
        static void calculate_with_3d_positive(UpwindScheme<D>& scheme)
        {
                scheme.calculate_with_3d(PositiveSpeed());
        }

        template<typename D>
        static void calculate_with_3d_negative(UpwindScheme<D>& scheme)
        {
                scheme.calculate_with_3d(NegativeSpeed());
        }



        template<typename D>
        static int left(UpwindScheme<D>& scheme)
        {
                return scheme.left_;
        }
        
        template<typename D>
        static int right(UpwindScheme<D>& scheme)
        {
                return scheme.right_;
        }

        template<typename D>
        static int top(UpwindScheme<D>& scheme)
        {
                return scheme.top_;
        }

        template<typename D>
        static int bottom(UpwindScheme<D>& scheme)
        {
                return scheme.bottom_;
        }

        template<typename D>
        static int front(UpwindScheme<D>& scheme)
        {
                return scheme.front_;
        }

        template<typename D>
        static int back(UpwindScheme<D>& scheme)
        {
                return scheme.back_;
        }

        template<typename D>
        static double fdxp(UpwindScheme<D>& scheme)
        {
                return scheme.fdxp_;
        }

        template<typename D>
        static double fdxm(UpwindScheme<D>& scheme)
        {
                return scheme.fdxm_;
        }

        template<typename D>
        static double fdyp(UpwindScheme<D>& scheme)
        {
                return scheme.fdyp_;
        }

        template<typename D>
        static double fdym(UpwindScheme<D>& scheme)
        {
                return scheme.fdym_;
        }

        template<typename D>
        static int fdzp(UpwindScheme<D>& scheme)
        {
                return scheme.fdzp_;
        }

        template<typename D>
        static int fdzm(UpwindScheme<D>& scheme)
        {
                return scheme.fdzm_;
        }
};

} // namespace lsm

namespace
{
        void test_set_position_with_2d()
        {
                const SpaceSize<TwoDimension> size {3, 3};
                const std::vector<double> phi {
                        0, 3, 0,
                        4, 2, 6,
                        0, 5, 0,
                };
                const Indexer<TwoDimension> indexer {size};
                UpwindScheme<TwoDimension> scheme {indexer, phi};
                const IntPoint2d p {{1, 1}};
                UpwindSchemeTester::set_position(scheme, p);

                BOOST_CHECK(UpwindSchemeTester::left(scheme) == 3);
                BOOST_CHECK(UpwindSchemeTester::right(scheme) == 5);
                BOOST_CHECK(UpwindSchemeTester::top(scheme) == 1);
                BOOST_CHECK(UpwindSchemeTester::bottom(scheme) == 7);
        }

        void test_calculate_with_2d_positive()
        {
                const SpaceSize<TwoDimension> size {3, 3};
                const std::vector<double> phi {
                        0, 3, 0,
                        4, 2, 6,
                        0, 5, 0,
                };
                const Indexer<TwoDimension> indexer {size};
                UpwindScheme<TwoDimension> scheme {indexer, phi};
               
                const IntPoint2d p {{1, 1}};
                UpwindSchemeTester::set_position(scheme, p);
                UpwindSchemeTester::calculate_with_2d_positive(scheme);

                BOOST_CHECK(UpwindSchemeTester::fdxm(scheme) == 0.0);
                BOOST_CHECK(UpwindSchemeTester::fdxp(scheme) == 0.0);
                BOOST_CHECK(UpwindSchemeTester::fdym(scheme) == 0.0);
                BOOST_CHECK(UpwindSchemeTester::fdyp(scheme) == 0.0);
        }

        void test_calculate_with_2d_negative()
        {
                const SpaceSize<TwoDimension> size {3, 3};
                const std::vector<double> phi {
                        0, 3, 0,
                        4, 2, 6,
                        0, 5, 0,
                };
                const Indexer<TwoDimension> indexer {size};
                UpwindScheme<TwoDimension> scheme {indexer, phi};
               
                const IntPoint2d p {{1, 1}};
                UpwindSchemeTester::set_position(scheme, p);
                UpwindSchemeTester::calculate_with_2d_negative(scheme);
                
                BOOST_CHECK(UpwindSchemeTester::fdxp(scheme) == 4.0);
                BOOST_CHECK(UpwindSchemeTester::fdxm(scheme) == -2.0);
                BOOST_CHECK(UpwindSchemeTester::fdyp(scheme) == 3.0);
                BOOST_CHECK(UpwindSchemeTester::fdym(scheme) == -1.0);
        }

        void test_calculate_2d()
        {
                const SpaceSize<TwoDimension> size {3, 3};
                const std::vector<double> phi {
                        0, 3, 0,
                        4, 2, 6,
                        0, 5, 0,
                };
                const Indexer<TwoDimension> indexer {size};
                UpwindScheme<TwoDimension> scheme {indexer, phi};
               
                const IntPoint2d p {{1, 1}};
                BOOST_CHECK(0.0 == scheme.calculate<PositiveSpeed>(p, TwoDimension()));
                BOOST_CHECK(std::sqrt(30.0) == scheme.calculate<NegativeSpeed>(p, TwoDimension()));
        }

        void test_2d()
        {
                test_set_position_with_2d();
                test_calculate_with_2d_positive();
                test_calculate_with_2d_negative();
                test_calculate_2d();
        }

        void test_3d()
        {
                const SpaceSize<ThreeDimension> size {3, 3, 3};
                const std::vector<double> phi {
                        0, 0, 0,
                        0, 7, 0,
                        0, 0, 0,

                        0, 3, 0,
                        4, 2, 6,
                        0, 5, 0,

                        0, 0, 0,
                        0, 8, 0,
                        0, 0, 0,
                };

                const Indexer<ThreeDimension> indexer {size};
                UpwindScheme<ThreeDimension> scheme {indexer, phi};

                const IntPoint3d p {{1, 1, 1}};
                BOOST_CHECK(0.0 == scheme.calculate<PositiveSpeed>(p, ThreeDimension()));
                BOOST_CHECK(std::sqrt(91.0) == scheme.calculate<NegativeSpeed>(p, ThreeDimension()));
        }
}

#include <iostream>
BOOST_AUTO_TEST_CASE(TEST_UpwindScheme)
{
        std::cout << "UpwindScheme\n";
        test_2d();
        test_3d();
}
#endif // UNIT_TEST_UpwindScheme
