//
//  InitialFront.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/11/23.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_InitialFront_h
#define LevelSetMethod2_InitialFront_h

#include "Point.h"
#include "Concepts.h"
#include <boost/concept_check.hpp>

namespace lsm
{
        /// represents an initial front
        /**
         *      @tparam D the dimension type
         */
        template<typename D>
        struct InitialFront
        {
                BOOST_CLASS_REQUIRE(D, lsm, DimensionConcept);
                std::array<IntPoint<D>, 2> vertices_;
        }; // struct InitialFront
} // namespace lsm

#endif // LevelSetMethod2_InitialFront_h
