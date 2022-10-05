//
//  LevelSetMethodViewer3d.cpp
//  LevelSetMethod2
//
//  Created by kumada on 2013/01/12.
//  Copyright (c) 2013年 kumada. All rights reserved.
//

#include "LevelSetMethodViewer3d.h"
#include "MatrixStack.h"
#include "EnableClientState.h"
#include "GeometryGenerator.h"
#include "Geometry.h"

using namespace lsm;

void LevelSetMethodViewer3d::set_translation_matrix()
{
        const double a = static_cast<double>(-space_size_.width_ / 2);
        const double b = static_cast<double>(space_size_.height_ / 2);
        const double c = incremental_point_of_view_;
        /*
                (x,y,z,w) is converted to (X,Y,Z,W) by the following matrix:
        
                X       1,  0,  0, a    x
                Y       0, -1,  0, b    y
                Z       0,  0, -1, c    z
                W       0,  0,  0, 1    1
                
                It should be noticed that the alignment of elements is different from one we expect.
         */
        
        set_translation_row(0, 1,  0,  0, 0);
        set_translation_row(1, 0, -1,  0, 0);
        set_translation_row(2, 0,  0, -1, 0);
        set_translation_row(3, a,  b,  c, 1);
}

void LevelSetMethodViewer3d::keyboard(unsigned char key, int x, int y)
{
        Base::keyboard_base(key, x, y);
        switch ( key ) {
        case 60: // <
                ++incremental_point_of_view_;
                break;
        case 62: // >
                --incremental_point_of_view_;
                break;
        case 117: // U
                x_angle_ += 1.0;
                break;
        case 110: // N
                x_angle_ -= 1.0;
                break;
        case 104: // H
                y_angle_ += 1.0;
                break;
        case 106: // J
                y_angle_ -= 1.0;
                break;
        default:
                break;
        };
}

void LevelSetMethodViewer3d::set_geometries(const std::vector<std::shared_ptr<Geometry>>* geometries)
{
        geometries_ = geometries;
}

/**
 *      @param[in] argc         a well-known argument
 *      @param[in] argv         a well-known argument
 *      @param[in] space_size      
 *      @param[in] title        viewer title
 *      @param[in] display_fun  callback for OpenGL
 *      @param[in] keyboard_fun callback for OpenGL
 */
void LevelSetMethodViewer3d::initialize_viewer_derived(
        int                argc,
        char*              argv[],
        const SpaceSize3d& space_size,
        const std::string& title,
        void (*display_fun)(),
        void (*keyboard_fun)(unsigned char key, int x, int y))
{
        space_size_ = space_size;
        set_translation_matrix();
        
        glutInit(&argc, argv);
        glutInitWindowPosition(100, 50);
        glutInitWindowSize(space_size.width_, space_size.height_);
        glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
        glutCreateWindow(title.c_str());
        glutIdleFunc(glutPostRedisplay);
        glutDisplayFunc(display_fun);
        glutKeyboardFunc(keyboard_fun);

        //_/_/_/
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        double ar = static_cast<double>(space_size.width_) / space_size.height_;
        gluPerspective(fovy_, ar, 1, 3 * space_size.depth_);
        gluLookAt(
                0.0, 0.0, 2 * space_size.depth_ + incremental_point_of_view_, // point of view
                0.0, 0.0, 0.0, // point of center
                0.0, 1.0, 0.0  // up direction
        );
        
        //_/_/_/
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        glClearColor(1, 1, 1, 1);
        
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_LIGHT0);
        glEnable(GL_NORMALIZE);
        glEnable(GL_COLOR_MATERIAL);
        glEnable(GL_LIGHTING);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        
        constexpr GLfloat light_ambient[]  = { 0.0f, 0.0f, 0.0f, 1.0f };
        constexpr GLfloat light_diffuse[]  = { 1.0f, 1.0f, 1.0f, 1.0f };
        constexpr GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
        const GLfloat light_position[]     = { 0.0f, 0.0f, static_cast<GLfloat>(space_size.depth_), 0.0f };
        
        glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
        glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
        glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
        glLightfv(GL_LIGHT0, GL_POSITION, light_position);

        constexpr GLfloat mat_ambient[]    = { 0.7f, 0.7f, 0.7f, 1.0f };
        constexpr GLfloat mat_diffuse[]    = { 0.8f, 0.8f, 0.8f, 1.0f };
        constexpr GLfloat mat_specular[]   = { 1.0f, 1.0f, 1.0f, 1.0f };
        constexpr GLfloat high_shininess[] = { 100.0f };

        glMaterialfv(GL_FRONT, GL_AMBIENT,   mat_ambient);
        glMaterialfv(GL_FRONT, GL_DIFFUSE,   mat_diffuse);
        glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
        glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess);
        
        Debug::debug(">> method: initialize_viewer done");

}

void LevelSetMethodViewer3d::display_front()
{
        level_set_method_->calculate_normals();
        level_set_method_->calculate_normals();
        const auto& normals = level_set_method_->get_normals();
        const auto& front = level_set_method_->get_front();
        display_front_if_enables(normals, front);
        previous_front_ = std::move(front);
        previous_normals_ = std::move(normals);
        future_ = std::async(std::launch::async, std::bind(&LevelSetMethodViewer3d::evolve_front, this));
}

void LevelSetMethodViewer3d::preprocess_to_display_front(bool enables_objects)
{
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glPushMatrix();
        set_translation_matrix();
        glRotated(x_angle_, 1.0, 0.0, 0.0);
        glRotated(y_angle_, 0.0, 1.0, 0.0);
        
        if ( enables_objects ) {
                for ( const auto& g : *geometries_ ) {
                        g->draw(translation_matrix_);
                }
        }
}

void LevelSetMethodViewer3d::postprocess_to_display_front() const
{
        glPopMatrix();
}

/**
 *      @param[in] normals normal vectors for the points that constitute the front
 *      @param[in] front   the front drawn on the viewer
 */
void LevelSetMethodViewer3d::display_front(
        const std::vector<DoublePoint3d>& normals,
        const Front3d&                    front
) const {
        glPointSize(1);
        {
                glColor4d(1, 0, 0, 0.5);
                MatrixStack stack;
                glMultMatrixd(translation_matrix_);
                
                EnableClientState enable_vertex {GL_VERTEX_ARRAY};
                EnableClientState enable_normal {GL_NORMAL_ARRAY};

                glVertexPointer(3 , GL_INT , sizeof(IntPoint3d), &front[0]);
                glNormalPointer(GL_DOUBLE , sizeof(DoublePoint3d), &normals[0]);
                glDrawArrays(GL_POINTS , 0 , static_cast<int>(front.size()));
        }
}