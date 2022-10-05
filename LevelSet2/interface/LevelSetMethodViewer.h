//
//  LevelSetMethodViewer.h
//  LevelSetMethod2
//
//  Created by kumada on 2013/01/12.
//  Copyright (c) 2013年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_LevelSetMethodViewer_h
#define LevelSetMethod2_LevelSetMethodViewer_h

#include <boost/concept_check.hpp>
#include <memory>
#include <thread>

#include "DimensionTypes.h"
#include "LevelSetMethod.h"
#include "SpaceSize.h"
#ifdef _WIN32
#include <Windows.h>
#endif
#include <gl/gl.h>
#ifdef LINUX
        #include <GL/glut.h>
#else
        #include <gl/glut.h>
#endif

namespace lsm
{
        /// Viewer for displaying a sequence of fronts calculated by the Level Set Method
        /**
         *      @tparam Dimension the dimension type
         *      @tparam Derived   the derived type
         *      @see CRTP(Curiously Reccursive Template Pattern)
         */
        template<typename Dimension, typename Derived>
        class LevelSetMethodViewer
        {
                BOOST_CLASS_REQUIRE(Dimension, lsm, DimensionConcept);
                typedef LevelSetMethodViewer Self;
                
        public:
                LevelSetMethodViewer() = default;
                ~LevelSetMethodViewer() = default;
                
                /**
                 *      @param[in] argc         a well-known argument
                 *      @param[in] argv         a well-known argument
                 *      @param[in] space_size      
                 *      @param[in] title        viewer title
                 *      @param[in] display_fun  callback for OpenGL
                 *      @param[in] keyboard_fun callback for OpenGL
                 */
                void initialize_viewer(
                        int                         argc,
                        char*                       argv[],
                        const SpaceSize<Dimension>& space_size,
                        const std::string&          title,
                        void (*display_fun)(),
                        void (*keyboard_fun)(unsigned char key, int x, int y))
                {
                        derived().initialize_viewer_derived(argc, argv, space_size, title, display_fun, keyboard_fun);
                }
                
                void set_verbose(bool is_verbose)
                {
                        is_verbose_ = is_verbose;
                }
                
                void start()
                {
                        start_time_ = std::chrono::steady_clock::now();
                        glutMainLoop();
                }

                /**
                 *      @param[in] params        command-line paramters
                 *      @param[in] space_size
                 *      @param[in] object        input object
                 *      @param[in] initial_front The initial front must embrace an input object.
                 */
                void initialize_level_set_method(
                        const Parameters&              params,
                        const SpaceSize<Dimension>&    space_size,
                        const std::uint8_t*            object,
                        const InitialFront<Dimension>& initial_front
                ) {
                        level_set_method_.reset(new LevelSetMethod<Dimension>{params, space_size});
                        
                        // set an input object
                        auto& input_object = level_set_method_->input_object();
                        input_object.resize(space_size.total_);
                        std::memcpy(&input_object[0], object, space_size.total_);
                        
                        level_set_method_->initialize_distance_map();
                        level_set_method_->initialize_along_front(initial_front);
                        level_set_method_->initialize_over_all(initial_front);
                        level_set_method_->calculate_speed_factors();
                        level_set_method_->initialize_narrow_band();
                        
                        future_ = std::async(std::launch::async, std::bind(&Self::evolve_front, this));
                        Debug::debug(">> method: initialize_level_set_method done");
                }
                
                /// used in a callback function for OpenGL
                void display()
                {
                        derived().preprocess_to_display_front(enables_objects_);
                        display_core();
                        glutSwapBuffers();
                        const_derived().postprocess_to_display_front();
                }

                /// set (v0,v1,v2,v3) to a index-th row
                void set_translation_row(int index, double v0, double v1, double v2, double v3)
                {
                        const int i = index << 2;
                        translation_matrix_[i]     = v0;
                        translation_matrix_[i + 1] = v1;
                        translation_matrix_[i + 2] = v2;
                        translation_matrix_[i + 3] = v3;
                }
                
        protected:
                std::unique_ptr<LevelSetMethod<Dimension>> level_set_method_;
                
                Front<Dimension> previous_front_;
                
                std::vector<DoublePoint<Dimension>> previous_normals_;
                
                std::future<bool> future_;
                
                SpaceSize<Dimension> space_size_;
                
                /// a 4x4 matrix to move objects to the center of the viewer
                double translation_matrix_[16];

                /// used in a callback function for OpenGL
                void keyboard_base(unsigned char key, int x, int y)
                {
//                        Debug::debug(std::to_string(static_cast<int>(key)));
                        switch ( key ) {
                        case 27: // ESC key
                                quits_viewer_ = true;
                                Debug::debug(">> After finishing current process, the viewer will force-quit.");
                                break;
                                
                        case 112: // P
                                pauses_level_set_method_ = !pauses_level_set_method_;
                                if ( pauses_level_set_method_ ) {
                                        Debug::debug(">> The Level Set Method is paused");
                                } else {
                                        Debug::debug(">> The Level Set Method restarts");
                                }
                                break;
                                
                        case 102: // F
                                enables_front_ = !enables_front_;
                                break;
                                
                        case 111: // O
                                enables_objects_ = !enables_objects_;
                                break;
                                
                        default:
                                break;
                        }
                }
        
                /**
                 *      @param[in]      normals normal vectors for each point that constitutes the front
                 *      @param[in]      front   the front drawn on the viewer
                 */
                void display_front_if_enables(
                        const std::vector<DoublePoint<Dimension>>& normals,
                        const Front<Dimension>&                    front
                ) const {
                        if ( enables_front_ ) {
                                const_derived().display_front(normals, front);
                        }
                }
//        public:
                bool evolve_front()
                {
                        if ( quits_viewer_ ) {
                                return false;
                        }
                        
                        completes_level_set_method_ = level_set_method_->set_speed_function(initializes_phi_);
                        if ( completes_level_set_method_ ) {
                                auto end = std::chrono::steady_clock::now();
                                std::cout << ">> computation time: " << std::chrono::duration<double, std::milli>(end - start_time_).count() / 1000 << "[sec]\n";
                                return false;
                        }
                        level_set_method_->propagate_front();
                        initializes_phi_ = level_set_method_->create_labels();
                        return true;
                }

        private:        
                /// a flag that the user sets to pause the viewer
                bool pauses_level_set_method_ {true};

                /// a flag to display a front
                bool enables_front_ {true};
 
                /// a flag to display objects
                bool enables_objects_ {true};
                
                /// a flag that the user sets to force-quit the viewer
                bool quits_viewer_ {false};
                
                /// starting time of the Level Set Method
                std::chrono::time_point<std::chrono::steady_clock> start_time_;
                
                /// a flag representing completion of the Level Set Method
                bool completes_level_set_method_ {false};
                
                /// a flag used for initializing the function "phi"
                bool initializes_phi_ {true};
                
                /// a flag to decide whether a verbose description enables or not
                bool is_verbose_ {false};
                
                LevelSetMethodViewer(const LevelSetMethodViewer&) = delete;
                LevelSetMethodViewer& operator=(const LevelSetMethodViewer&) = delete;
                LevelSetMethodViewer(LevelSetMethodViewer&&) = delete;
                LevelSetMethodViewer& operator=(LevelSetMethodViewer&&) = delete;

                Derived& derived()
                {
                        return static_cast<Derived&>(*this);
                }
                
                const Derived& const_derived() const
                {
                        return static_cast<const Derived&>(*this);
                }
                
                void display_core()
                {
                        if ( pauses_level_set_method_ ) {
                                display_front_if_enables(previous_normals_, previous_front_);
                                return;
                        }
                        
                        if ( !future_.valid() ) {
                                display_front_with_exit(">> force-quit(The Level Set Method was over.)", "");
                                return;
                        }
                        
                        if ( future_.wait_for(std::chrono::seconds(0)) == std::future_status::ready ) {
                                if ( future_.get() ) {
                                        if ( is_verbose_ ) {
                                                level_set_method_->print_verbose_description();
                                        }
                                        derived().display_front();
                                        
                                } else {
                                        display_front_with_exit(">> force-quit", ">> The Level Set Method has just ended.");
                                }
                                return;
                        }
                        display_front_if_enables(previous_normals_, previous_front_);
                }

                void display_front_with_exit(
                        const std::string& message0,
                        const std::string& message1
                ) const {
                        if ( quits_viewer_ ) {
                                if ( !message0.empty() ) {
                                        Debug::debug(message0);
                                }
                                std::exit(0);
                                return;
                        }
                        
                        if ( !message1.empty() ) {
                                Debug::debug(message1);
                        }
                        display_front_if_enables(previous_normals_, previous_front_);
                
                }
        }; // class LevelSetMethodViewer
} // namespace lsm

#endif // LevelSetMethod2_LevelSetMethodViewer_h
