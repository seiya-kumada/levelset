//
//  Sphere.cpp
//  LevelSetMethod2
//
//  Created by kumada on 2013/01/24.
//  Copyright (c) 2013年 kumada. All rights reserved.
//

#include "Sphere.h"
#include "MatrixStack.h"
#ifdef _WIN32
#include <Windows.h>
#endif
#include <gl/glut.h>

using namespace lsm;

Sphere::Sphere(bool is_solid, const IntPoint3d& center, int radius)
        : Geometry(is_solid)
        , center_(center)
        , radius_(radius) {}

void Sphere::assign_value(const IntPoint3d& p, std::uint8_t& v) const 
{
        const int dx = p[0] - center_[0];
        const int dy = p[1] - center_[1];
        const int dz = p[2] - center_[2];
        const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
        if ( is_solid_ ) {
                if ( r < radius_ ) {
                        v = 128;
                }
        } else {
                if ( r > radius_ ) {
                        v = 128;
                }
        }
}

void Sphere::draw_geometry(const double matrix[16]) const
{
        MatrixStack stack;
        glColor3d(0, 0, 1);
        glMultMatrixd(matrix);
        glTranslated(center_[0], center_[1], center_[2]);
        if ( is_solid_ ) {
                glutSolidSphere(radius_, 32, 32);
        } else {
                glutWireSphere(radius_, 16, 16);
        }
}