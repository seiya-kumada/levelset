//
//  LoadedGeometry.cpp
//  LevelSetMethod2
//
//  Created by kumada on 2013/02/10.
//  Copyright (c) 2013年 kumada. All rights reserved.
//

#include "LoadedGeometry.h"
#include "MatrixStack.h"
#include "EnableClientState.h"
#ifdef _WIN32
#include <Windows.h>
#endif
#include <gl/glut.h>

using namespace lsm;

LoadedGeometry::LoadedGeometry(
       bool                             is_solid,
       const std::vector<FloatPoint3d>& mixed_vertices
)
        : Geometry{is_solid}
        , mixed_vertices_(mixed_vertices)
        , vertex_size_(mixed_vertices.size() / 2) {}

void LoadedGeometry::draw_geometry(const double matrix[16]) const
{
        MatrixStack stack;
        EnableClientState enable_vertex {GL_VERTEX_ARRAY};
        glColor3d(0, 0, 1);
        glMultMatrixd(matrix);
        
        glInterleavedArrays(GL_N3F_V3F, 0, &mixed_vertices_[0]);
        glDrawArrays(GL_POINTS, 0, static_cast<int>(vertex_size_));
}