//
//  CurvatureGenerator.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/11/24.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_CurvatureGenerator_h
#define LevelSetMethod2_CurvatureGenerator_h

#include "Differential.h"

namespace lsm
{
        
        /// calculates a curvature of an auxiliary function
        /**
         *      @tparam D the dimension type
         *      @see DimensionTypes
         */
        template<typename D>
        class CurvatureGenerator
        {
                friend class CurvatureGeneratorTester;
                BOOST_CLASS_REQUIRE(D, lsm, DimensionConcept);
                
        public:
                CurvatureGenerator(
                        const Indexer<D>& indexer,
                        const std::vector<double>& phi
                )
                        : differential_{indexer, phi} {}
                
                ~CurvatureGenerator() = default;

                /// calculates a curvature at a given point
                /**
                 *      @param p the given point
                 *      @return a curvature
                 */
                const double operator()(const IntPoint<D>& p);
                
                /// calculates a normal vector at a given point
                /**
                 *      @param p the given point
                 *      @return a normal vector
                 */
                const DoublePoint<D> calculate_normal(const IntPoint<D>& p);
                
        private:
                CurvatureGenerator(const CurvatureGenerator&) = delete;
                CurvatureGenerator& operator=(const CurvatureGenerator&) = delete;
                
                CurvatureGenerator(CurvatureGenerator&&) = delete;
                CurvatureGenerator& operator=(CurvatureGenerator&&) = delete;
                
                Differential<D, double> differential_;
        }; // class CurvatureGenerator

        /// explicit specialization with TwoDimension
        template<>
        inline const double CurvatureGenerator<TwoDimension>::operator()(const IntPoint2d& p)
        {
                differential_.set_point(p);
                const double dfx = differential_.fx(); // test ok
                const double dfy = differential_.fy(); // test ok
                
                const double dfxy = differential_.fxy(); // test ok
                
                const double dfxx = differential_.fxx(); // test ok
                const double dfyy = differential_.fyy(); // test ok
                
                const double dfx2 = dfx * dfx;
                const double dfy2 = dfy * dfy;
                
                const double df = std::sqrt(dfx2 + dfy2);
                return (df != 0.0) ?
                        (
                                  dfxx * dfy2
                                + dfyy * dfx2
                                - 2.0 * dfx * dfy * dfxy
                         
                        ) / (df * df * df)
                        : 0.0;
        }

        /// explicit specialization with TwoDimension
        template<>
        inline const DoublePoint2d CurvatureGenerator<TwoDimension>::calculate_normal(const IntPoint2d& p)
        {
                differential_.set_point(p);
                return {{differential_.fx(), differential_.fy()}};
        }

        /// explicit specialization with ThreeDimension
        template<>
        inline const double CurvatureGenerator<ThreeDimension>::operator()(const IntPoint3d& p)
        {
                differential_.set_point(p);
                const double dfx = differential_.fx(); // test ok
                const double dfy = differential_.fy(); // test ok
                const double dfz = differential_.fz(); // test ok
                
                const double dfxy = differential_.fxy(); // test ok
                const double dfxz = differential_.fxz(); // test ok
                const double dfyz = differential_.fyz(); // test ok
                
                const double dfxx = differential_.fxx(); // test ok
                const double dfyy = differential_.fyy(); // test ok
                const double dfzz = differential_.fzz(); // test ok

                const double dfx2 = dfx * dfx;
                const double dfy2 = dfy * dfy;
                const double dfz2 = dfz * dfz;
                
                const double df = std::sqrt(dfx2 + dfy2 + dfz2);
                return (df != 0.0) ?
                        (
                                  (dfyy + dfzz) * dfx2
                                + (dfxx + dfzz) * dfy2
                                + (dfxx + dfyy) * dfz2
                                - 2.0 * dfx * dfy * dfxy
                                - 2.0 * dfx * dfz * dfxz
                                - 2.0 * dfy * dfz * dfyz
                        ) / (df * df * df)
                        : 0.0;
        }

        /// explicit specialization with ThreeDimension
        template<>
        inline const DoublePoint3d CurvatureGenerator<ThreeDimension>::calculate_normal(const IntPoint3d& p)
        {
                differential_.set_point(p);
                return {{differential_.fx(), differential_.fy(), differential_.fz()}};
        }

} // namespace lsm

#endif // LevelSetMethod2_CurvatureGenerator_h
