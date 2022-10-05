//
//  MatrixStack.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/12/22.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_MatrixStack_h
#define LevelSetMethod2_MatrixStack_h
#ifdef _WIN32
#include <Windows.h>
#endif
#include <gl/gl.h>

namespace lsm
{
        /// helper class
        class MatrixStack
        {
        public:
                MatrixStack()
                {
                        glPushMatrix();
                }
                
                ~MatrixStack()
                {
                        glPopMatrix();
                }
        }; // class MatrixStack
} // namespace lsm
#endif
