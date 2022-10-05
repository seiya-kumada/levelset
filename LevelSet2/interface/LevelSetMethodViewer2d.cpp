//
//  LevelSetMethodViewer2d.cpp
//  LevelSetMethod2
//
//  Created by kumada on 2013/01/12.
//  Copyright (c) 2013年 kumada. All rights reserved.
//

#include "LevelSetMethodViewer2d.h"
#include "MatrixStack.h"
#include "EnableClientState.h"
#include "GeometryDescriptor.h"

using namespace lsm;

/**
 *      @param[in] dummy not-used argument
 *      @param[in] front front drawn on the viewer
 */
void LevelSetMethodViewer2d::display_front(
        const std::vector<DoublePoint2d>& dummy,
        const Front2d&                    front
) const {
//glClear(GL_COLOR_BUFFER_BIT);
        glPointSize(1);
        {
                glColor3d(0, 0, 0);
                MatrixStack stack;
                glMultMatrixd(translation_matrix_);
                EnableClientState enable_vertex {GL_VERTEX_ARRAY};
                glVertexPointer(2, GL_INT, sizeof(IntPoint2d), &front[0]);
                glDrawArrays(GL_POINTS, 0, static_cast<int>(front.size()));
        }
}

void LevelSetMethodViewer2d::display_front()
{
        const auto& front = level_set_method_->get_front();
        display_front_if_enables(previous_normals_, front); // previous_normals_ is not used.
        previous_front_ = std::move(front);
        future_ = std::async(std::launch::async, std::bind(&LevelSetMethodViewer2d::evolve_front, this));
}

void LevelSetMethodViewer2d::preprocess_to_display_front(bool enables_objects)
{
        glClear(GL_COLOR_BUFFER_BIT);
        if ( enables_objects ) {
                MatrixStack stack;
                glMultMatrixd(translation_matrix_);
                
                glColor3d(1.0, 1.0, 1.0);
                GeometryDescriptor gd(GL_POLYGON);
                
                glTexCoord2d(0, 0);
                glVertex3d(0, 0, 0);
                
                glTexCoord2f(1, 0);
                glVertex3d(space_size_.width_ - 1, 0, 0);
                
                glTexCoord2f(1, 1);
                glVertex3d(space_size_.width_ - 1, space_size_.height_ - 1, 0);

                glTexCoord2f(0, 1);
                glVertex3d(0, space_size_.height_ - 1, 0);
        }
}

void LevelSetMethodViewer2d::keyboard(unsigned char key, int x, int y)
{
        Base::keyboard_base(key, x, y);
}

/**
 *      @param[in] argc         a well-known argument
 *      @param[in] argv         a well-known argument
 *      @param[in] space_size      
 *      @param[in] title        viewer title
 *      @param[in] display_fun  callback for OpenGL
 *      @param[in] keyboard_fun callback for OpenGL
 */
void LevelSetMethodViewer2d::initialize_viewer_derived(
        int                argc,
        char*              argv[],
        const SpaceSize2d& space_size,
        const std::string& title,
        void (*display_fun)(),
        void (*keyboard_fun)(unsigned char key, int x, int y)
) {
        space_size_ = space_size;
        set_translation_matrix();
        
        glutInit(&argc, argv);
        glutInitWindowPosition(100, 50);
        glutInitWindowSize(space_size.width_, space_size.height_);
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
        glutCreateWindow(title.c_str());
        
        //_/_/_/
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(0, space_size.width_ - 1, -(space_size.height_ - 1), 0,  2, 4);
        
        glutIdleFunc(glutPostRedisplay);
        glutDisplayFunc(display_fun);
        glutKeyboardFunc(keyboard_fun);

        glClearColor(1, 1, 1, 1);

        // set a texture
        glEnable(GL_TEXTURE_2D);
        glGenTextures(1, &texture_name_);
        glBindTexture(GL_TEXTURE_2D , texture_name_);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        
        glTexImage2D(
                GL_TEXTURE_2D,    // target
                0,                // level
                GL_LUMINANCE,     // internal format
                space_size_.width_,
                space_size_.height_,
                0,                // border
                GL_LUMINANCE,     // format
                GL_UNSIGNED_BYTE, // type
                texture_          // data
        );
        Debug::debug(">> method: initialize_viewer done");
}

void LevelSetMethodViewer2d::set_translation_matrix()
{
        /*
                (x,y,z,w) is converted to (X,Y,Z,W) by the following matrix:
        
                X       1,  0,  0,  0    x
                Y       0, -1,  0,  0    y
                Z       0,  0,  0, -2    z
                W       0,  0,  0,  1    1
                
                It should be noticed that the alignment of elements is different from one we expect.
         */

        set_translation_row(0, 1,  0,  0, 0);
        set_translation_row(1, 0, -1,  0, 0);
        set_translation_row(2, 0,  0,  0, 0);
        set_translation_row(3, 0,  0, -2, 1);
}
