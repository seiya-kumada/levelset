//
//  Cube.cpp
//  LevelSetMethod2
//
//  Created by kumada on 2013/01/25.
//  Copyright (c) 2013年 kumada. All rights reserved.
//

#include "Cube.h"
#include "MatrixStack.h"
#ifdef _WIN32
#include <windows.h>
#endif
#include <gl/glut.h>
#include <gl/GL.h>
using namespace lsm;

namespace
{
        constexpr int Left = 0;
        constexpr int Top = 1;
        constexpr int Front = 2;
        
        constexpr int Right = 0;
        constexpr int Bottom = 1;
        constexpr int Back = 2;
}

Cube::Cube(
        bool              is_solid,
        const IntPoint3d& center,
        int               size
)
        : Geometry{is_solid}
        , center_(center)
        , size_(size)
        , left_top_front_({{center[0] - size / 2, center[1] - size / 2, center[2] - size / 2}})
        , right_bottom_back_({{center[0] + size / 2, center[1] + size / 2, center[2] + size / 2}}) {}


void Cube::assign_value(const IntPoint3d& p, std::uint8_t& v) const
{
        if ( is_solid_ ) {
                if ( is_inside(p) ) {
                        v = 128;
                }
        } else {
                if ( !is_inside(p) ) {
                        v = 128;
                }
        }
}

bool Cube::is_inside(const IntPoint3d& p) const
{
        const int x = p[0];
        const int y = p[1];
        const int z = p[2];
        return (left_top_front_[Left]  < x && x < right_bottom_back_[Right]) &&
               (left_top_front_[Top]   < y && y < right_bottom_back_[Bottom]) &&
               (left_top_front_[Front] < z && z < right_bottom_back_[Back]);
}

void Cube::draw_geometry(const double matrix[16]) const
{
        MatrixStack stack;
        glColor3d(0, 0, 1);
        glMultMatrixd(matrix);
        glTranslated(center_[0], center_[1], center_[2]);
        if ( is_solid_ ) {
                glutSolidCube(size_);
        } else {
                glutWireCube(size_);
        }
}