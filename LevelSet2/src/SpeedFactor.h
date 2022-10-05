//
//  SpeedFactor.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/11/24.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_SpeedFactor_h
#define LevelSetMethod2_SpeedFactor_h

#include "Concepts.h"
#include "DimensionTypes.h"
#include "Indexer.h"
#include "Differential.h"
#include <boost/concept_check.hpp>

namespace lsm
{
        /// calculates factors which depends on the input image
        /**
         *      @tparam D the dimension type
         */
        template<typename D>
        class SpeedFactor
        {
                BOOST_CLASS_REQUIRE(D, lsm, DimensionConcept);
                
        public:
                SpeedFactor(const Indexer<D>& indexer, const std::vector<std::uint8_t>& gray)
                        : indexer_(indexer)
                        , differential_{indexer, gray} {}
                
                ~SpeedFactor() = default;
                
                void calculate_all(const SpaceSize<D>& size);

                const double& value(const IntPoint<D>& p) const
                {
                        return factors_[indexer_(p)];
                }

        private:
                SpeedFactor(const SpeedFactor&) = delete;
                SpeedFactor& operator=(const SpeedFactor&) = delete;
                
                SpeedFactor(SpeedFactor&&) = delete;
                SpeedFactor& operator=(SpeedFactor&&) = delete;
                
                double calculate(const IntPoint<D>& p);
                
                const Indexer<D>&               indexer_;
                Differential<D, std::uint8_t>   differential_;
                std::vector<double>             factors_;
        }; // class SpeedFactor

        /// explicit specialization with TwoDimension
        template<>
        inline double SpeedFactor<TwoDimension>::calculate(const IntPoint2d& p)
        {
                differential_.set_point(p);
                const double dx = differential_.fx();
                const double dy = differential_.fy();;
                return 1.0 / (1.0 + std::sqrt(dx * dx + dy * dy));       
        }

        /// explicit specialization with ThreeDimension
        // test ok
        template<>
        inline double SpeedFactor<ThreeDimension>::calculate(const IntPoint3d& p)
        {
                differential_.set_point(p);
                const double dx = differential_.fx();
                const double dy = differential_.fy();
                const double dz = differential_.fz();;
                return 1.0 / (1.0 + std::sqrt(dx * dx + dy * dy + dz * dz));
        }

        /// explicit specialization with TwoDimension
        template<>
        inline void SpeedFactor<TwoDimension>::calculate_all(const SpaceSize<TwoDimension>& size)
        {
                const int w = size.width_;
                const int h = size.height_;
                factors_.resize(w * h);
                for ( int j = 1, wj = w * j; j < h - 1; ++j, wj += w ) {
                        double* p = &factors_[wj];
                        for ( int i = 1; i < w - 1; ++i ) {
                                p[i] = calculate({{i, j}});
                        }
                }
        }

        /// explicit specialization with ThreeDimension
        template<>
        inline void SpeedFactor<ThreeDimension>::calculate_all(const SpaceSize<ThreeDimension>& size)
        {
                const int w = size.width_;
                const int h = size.height_;
                const int a = w * h;
                const int d = size.depth_;
                factors_.resize(a * d);
                for ( int k = 1, ak = a * k; k < d - 1; ++k, ak += a ) {
                        for ( int j = 1, wj = ak + w * j; j < h - 1; ++j, wj += w ) {
                                double* p = &factors_[wj];
                                for ( int i = 1; i < w - 1; ++i ) {
                                        p[i] = calculate({{i, j, k}});
                                }
                        }
                }
        }
} // namespace lsm

#endif // LevelSetMethod2_SpeedFactor_h
