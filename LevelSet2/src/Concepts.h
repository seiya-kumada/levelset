//
//  Concepts.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/11/18.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_Concepts_h
#define LevelSetMethod2_Concepts_h

namespace lsm
{

/// Concept for the dimension
/**
 *      @tparam T the dimension type
 *      @see DimensionTypes.h
 */
template<typename T>
struct DimensionConcept
{
        void constraints()
        {
                (void)T::Dimension_;
                static_assert(T::Dimension_ == 2 || T::Dimension_ == 3, "Invalid Dimension Type");
        }
}; // struct DimensionConcept

} // namespace lsm
#endif // LevelSetMethod2_Concepts_h
