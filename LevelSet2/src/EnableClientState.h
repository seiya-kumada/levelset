//
//  Enable.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/12/26.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_EnableClientState_h
#define LevelSetMethod2_EnableClientState_h
#ifdef _WIN32
#include <windows.h>
#endif
#include <gl/GL.h>

namespace lsm
{

   
/// helper class
class EnableClientState
{
        GLenum enm_;
public:
        explicit EnableClientState(GLenum enm)
                : enm_{enm}
        {
                glEnableClientState(enm_);
        }
        
        ~EnableClientState()
        {
                glDisableClientState(enm_);
        }
}; // class EnableClientState

} // namespace lsm

#endif // LevelSetMethod2_EnableClientState_h
