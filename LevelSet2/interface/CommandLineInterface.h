//
//  CommandLineInterface.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/12/18.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_CommandLineInterface_h
#define LevelSetMethod2_CommandLineInterface_h

#include <memory>
#include <boost/program_options.hpp>

namespace lsm
{
        class LevelSetMethodViewer2d;
        class LevelSetMethodViewer3d;
        struct Parameters;
        class GeometryGenerator;
        class GeometryLoader;
        
        /// class that implements command-line interfaces
        class CommandLineInterface
        {
        private:
                /// used for 2d model
                static std::unique_ptr<LevelSetMethodViewer2d> viewer2d_;
                static void display_2d();
                static void keyboard_2d(unsigned char key, int x, int y);
                static void execute_level_set_method_in_2d(
                        int                                          argc,
                        char*                                        argv[],
                        const boost::program_options::variables_map& vm,
                        const Parameters&                            params
                );
        
                /// used for 3d model
                static std::unique_ptr<LevelSetMethodViewer3d> viewer3d_;
                static std::unique_ptr<GeometryGenerator>      geometry_generator_;
                static std::unique_ptr<GeometryLoader>         geometry_loader_;
                static void display_3d();
                static void keyboard_3d(unsigned char key, int x, int y);
                static void execute_level_set_method_in_3d(
                        int                                          argc,
                        char*                                        argv[],
                        const boost::program_options::variables_map& vm,
                        const Parameters&                            params
                );
                
        public:
                static void execute_level_set_method(int argc, char* argv[]);
        }; // class CommandLineInterface

} // namespace lsm

#endif // LevelSetMethod2_CommandLineInterface_h
