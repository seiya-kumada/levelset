//
//  LevelSetMethodViewer3d.h
//  LevelSetMethod2
//
//  Created by kumada on 2013/01/12.
//  Copyright (c) 2013年 kumada. All rights reserved.
//

#ifndef __LevelSetMethod2__LevelSetMethodViewer3d__
#define __LevelSetMethod2__LevelSetMethodViewer3d__

#include "LevelSetMethodViewer.h"

namespace lsm
{
        class Sphere;
        class Geometry;
        
        class LevelSetMethodViewer3d
                : public LevelSetMethodViewer<
                                ThreeDimension,
                                LevelSetMethodViewer3d
                         >
        {
                typedef LevelSetMethodViewer<ThreeDimension, LevelSetMethodViewer3d> Base;
                friend class LevelSetMethodViewer<ThreeDimension, LevelSetMethodViewer3d>;
        
        public:
                LevelSetMethodViewer3d() = default;
                ~LevelSetMethodViewer3d() = default;
                
                void set_geometries(const std::vector<std::shared_ptr<Geometry>>* geometries);
                void keyboard(unsigned char key, int x, int y);
                
        private:
                LevelSetMethodViewer3d(const LevelSetMethodViewer3d&) = delete;
                LevelSetMethodViewer3d& operator=(const LevelSetMethodViewer3d&) = delete;
                LevelSetMethodViewer3d(LevelSetMethodViewer3d&&) = delete;
                LevelSetMethodViewer3d& operator=(LevelSetMethodViewer3d&&) = delete;
                
                /// used to change a view point
                double incremental_point_of_view_ {0.0};

                /// angle around x axis
                double x_angle_ {0.0};
                
                /// angle around y axis
                double y_angle_ {0.0};
                
                /// field of view along y axis
                double fovy_ {46.0};

                /// objects drawn in the viewer
                const std::vector<std::shared_ptr<Geometry>>* geometries_;
                
                /**
                 *      @param[in] argc         a well-known argument
                 *      @param[in] argv         a well-known argument
                 *      @param[in] space_size      
                 *      @param[in] title        viewer title
                 *      @param[in] display_fun  callback for OpenGL
                 *      @param[in] keyboard_fun callback for OpenGL
                 */
                void initialize_viewer_derived(
                        int                argc,
                        char*              argv[],
                        const SpaceSize3d& space_size,
                        const std::string& title,
                        void (*display_fun)(),
                        void (*keyboard_fun)(unsigned char key, int x, int y)
                );

                void display_front();
                
                void preprocess_to_display_front(bool enables_objects);
                
                void postprocess_to_display_front() const;
                
                /**
                 *      @param[in] normals vectors for the points that constitute the front
                 *      @param[in] front   the front drawn on the viewer
                 */
                void display_front(
                        const std::vector<DoublePoint3d>& normals,
                        const Front3d&                    front
                ) const;
                
                void set_translation_matrix();
        }; // class LevelSetMethodViewer3d
} // namespace lsm
#endif /* defined(__LevelSetMethod2__LevelSetMethodViewer3d__) */
