//
//  GeometryDescriptor.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/12/22.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_GeometryDescriptor_h
#define LevelSetMethod2_GeometryDescriptor_h

namespace lsm
{
        /// helper class
        class GeometryDescriptor
        {
        public:
                explicit GeometryDescriptor(GLenum enm)
                {
                        glBegin(enm);
                }
                
                ~GeometryDescriptor()
                {
                        glEnd();
                }
        }; // class GeometryDescriptor
} // namespace lsm

#endif // LevelSetMethod2_GeometryDescriptor_h
