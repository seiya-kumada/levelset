//
//  ZeroLevelSetDetector.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/12/06.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_ZeroLevelSetDetector_h
#define LevelSetMethod2_ZeroLevelSetDetector_h

#include "Concepts.h"
#include "DimensionTypes.h"
#include "Point.h"
#include "Indexer.h"
#include "NeighboringPoints.h"
#include <boost/concept_check.hpp>

namespace lsm
{
        /// finds the zero level set
        /**
         *      @tparam D the dimension type
         */
        template<typename D>
        class ZeroLevelSetDetector
        {
                BOOST_CLASS_REQUIRE(D, lsm, DimensionConcept);

        public:
                ZeroLevelSetDetector(
                        const std::vector<double>& phi,
                        const Indexer<D>& indexer
                )
                        : phi_(phi )
                        , indexer_(indexer) {}
                
                bool operator()(const IntPoint<D>& p) const;
                
        private:
                const std::vector<double>& phi_;
                const Indexer<D>& indexer_;
                
                inline bool is_negative(const double& a, const IntPoint<D>& b) const
                {
                        return a + phi_[indexer_(b)] <= 0;
                }
                
                inline bool is_positive(const double& a, const IntPoint<D>& b) const
                {
                        return a + phi_[indexer_(b)] > 0;
                }
        };

        // test ok
        /// explicit specialization with TwoDimension
        template<>
        inline bool ZeroLevelSetDetector<TwoDimension>::operator()(const IntPoint<TwoDimension>& p) const
        {
                const double phi_p = phi_[indexer_(p)];
                if ( phi_p >= 0 ) {
                        if ( is_negative(phi_p, p + NeighboringPoints2d(-1,  0)) ) {
                                return true;
                        }
                        if ( is_negative(phi_p, p + NeighboringPoints2d( 1,  0)) ) {
                                return true;
                        }
                        if ( is_negative(phi_p, p + NeighboringPoints2d( 0, -1)) ) {
                                return true;
                        }
                        if ( is_negative(phi_p, p + NeighboringPoints2d( 0,  1)) ) {
                                return true;
                        }
                        return false;
                } else {
                        if ( is_positive(phi_p, p + NeighboringPoints2d(-1,  0)) ) {
                                return true;
                        }
                        if ( is_positive(phi_p, p + NeighboringPoints2d( 1,  0)) ) {
                                return true;
                        }
                        if ( is_positive(phi_p, p + NeighboringPoints2d( 0, -1)) ) {
                                return true;
                        }
                        if ( is_positive(phi_p, p + NeighboringPoints2d( 0,  1)) ) {
                                return true;
                        }
                        return false;
                }
        }

        // test ok
        /// explicit specialization with ThreeDimension
        template<>
        inline bool ZeroLevelSetDetector<ThreeDimension>::operator()(const IntPoint<ThreeDimension>& p) const
        {
                const double phi_p = phi_[indexer_(p)];
                if ( phi_p >= 0 ) {
                        if ( is_negative(phi_p, p + NeighboringPoints3d(-1,  0,  0)) ) {
                                return true;
                        }
                        if ( is_negative(phi_p, p + NeighboringPoints3d( 1,  0,  0)) ) {
                                return true;
                        }
                        if ( is_negative(phi_p, p + NeighboringPoints3d( 0, -1,  0)) ) {
                                return true;
                        }
                        if ( is_negative(phi_p, p + NeighboringPoints3d( 0,  1,  0)) ) {
                                return true;
                        }
                        if ( is_negative(phi_p, p + NeighboringPoints3d( 0,  0,  1)) ) {
                                return true;
                        }
                        if ( is_negative(phi_p, p + NeighboringPoints3d( 0,  0, -1)) ) {
                                return true;
                        }
                        return false;
                        
                } else {
                        if ( is_positive(phi_p, p + NeighboringPoints3d(-1,  0,  0)) ) {
                                return true;
                        }
                        if ( is_positive(phi_p, p + NeighboringPoints3d( 1,  0,  0)) ) {
                                return true;
                        }
                        if ( is_positive(phi_p, p + NeighboringPoints3d( 0, -1,  0)) ) {
                                return true;
                        }
                        if ( is_positive(phi_p, p + NeighboringPoints3d( 0,  1,  0)) ) {
                                return true;
                        }
                        if ( is_positive(phi_p, p + NeighboringPoints3d( 0,  0,  1)) ) {
                                return true;
                        }
                        if ( is_positive(phi_p, p + NeighboringPoints3d( 0,  0, -1)) ) {
                                return true;
                        }
                        return false;
                }
        }
} // namespace lsm
#endif // LevelSetMethod2_ZeroLevelSetDetector_h
