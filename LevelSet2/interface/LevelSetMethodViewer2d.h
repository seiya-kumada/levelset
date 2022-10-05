//
//  LevelSetMethodViewer2d.h
//  LevelSetMethod2
//
//  Created by kumada on 2013/01/12.
//  Copyright (c) 2013年 kumada. All rights reserved.
//

#ifndef __LevelSetMethod2__LevelSetMethodViewer2d__
#define __LevelSetMethod2__LevelSetMethodViewer2d__

#include "LevelSetMethodViewer.h"

namespace lsm
{
        class LevelSetMethodViewer2d
                : public LevelSetMethodViewer<
                                TwoDimension,
                                LevelSetMethodViewer2d
                         >
        {
                typedef LevelSetMethodViewer<TwoDimension, LevelSetMethodViewer2d> Base;
                friend class LevelSetMethodViewer<TwoDimension, LevelSetMethodViewer2d>;
                
        public:
                LevelSetMethodViewer2d() = default;
                ~LevelSetMethodViewer2d() = default;
                
                void set_texture(const std::uint8_t* texture)
                {
                        texture_ = texture;
                }
                void keyboard(unsigned char key, int x, int y);
                
        private:
                LevelSetMethodViewer2d(const LevelSetMethodViewer2d&) = delete;
                LevelSetMethodViewer2d& operator=(const LevelSetMethodViewer2d&) = delete;
                LevelSetMethodViewer2d(LevelSetMethodViewer2d&&) = delete;
                LevelSetMethodViewer2d& operator=(LevelSetMethodViewer2d&&) = delete;

                const std::uint8_t* texture_;
                
                std::uint32_t texture_name_;
                
                void display_front(
                        const std::vector<DoublePoint2d>& dummy,
                        const Front2d&                    front
                ) const;
                
                void display_front();
                void preprocess_to_display_front(bool enables_objects);
                void postprocess_to_display_front() const {}
                
                void set_translation_matrix();
                
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
                        const SpaceSize2d& size,
                        const std::string& title,
                        void (*display_fun)(),
                        void (*keyboard_fun)(unsigned char key, int x, int y)
                );

        }; // class LevelSetMethodViewer2d
        
} // namespace lsm
#endif /* defined(__LevelSetMethod2__LevelSetMethodViewer2d__) */
