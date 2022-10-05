//
//  CommandLineInterface.cpp
//  LevelSetMethod2
//
//  Created by kumada on 2012/12/18.
//  Copyright (c) 2012年 kumada. All rights reserved.
//
#ifdef _WIN32
#define NOMINMAX
#endif
#include "CommandLineInterface.h"
#include "LevelSetMethodViewer2d.h"
#include "LevelSetMethodViewer3d.h"
#include "GeometryGenerator.h"
#include "GeometryLoader.h"
#include "Sphere.h"
#include "Cube.h"

#include <memory>
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>

#include <boost/program_options.hpp>
#include <opencv2/opencv.hpp>

namespace program_options = boost::program_options;
using namespace lsm;

namespace
{
        /// extract a parameter from "vm"
        /**
         *      @param[in] name          parameter name to extract
         *      @param[in] error_message error message in case of failure
         *      @param[in] vm            argument container
         */
        template<typename T>
        inline const T extract_parameter(
                const std::string&                    name,
                const std::string&                    error_message,
                const program_options::variables_map& vm
        ) {
                if ( vm.count(name) ) {
                        return vm[name].as<T>();
                } else {
                        throw std::runtime_error(error_message);
                }
        }
        
        /// load an input image used for the Level Set Method in 2D
        /**
         *      @param[in]  filename   input file name
         *      @param[out] image      input image
         *      @param[out] space_size image size
         */
        void load_input_image(
                const std::string&         filename,
                std::vector<std::uint8_t>& image,
                SpaceSize<TwoDimension>&   space_size
        ) {
                cv::Mat gray = cv::imread(filename, CV_8UC1);
                if ( gray.empty() ) {
                        throw std::runtime_error("unable to load an input image");
                }
                space_size.width_ = gray.cols;
                space_size.height_ = gray.rows;
                space_size.total_ = gray.cols * gray.rows;
                std::copy(gray.data, gray.data + space_size.total_, std::back_inserter(image));
        }
        
        void set_common_arguments(program_options::options_description& desc)
        {
                desc.add_options()
                        ("help",                                                   "produce help message")
                        ("dim",             program_options::value<int>(),         "set either 2 or 3")
                        ("verbose",                                                "print verbose description")
                        ("input",           program_options::value<std::string>(), "set GRAY image path in 2D/set pattern in 3D")
                        ("wband",           program_options::value<int>(),         "set width of band")
                        ("wreset",          program_options::value<int>(),         "set width to reset")
                        ("time_step",       program_options::value<double>(),      "set time step")
                        ("gain",            program_options::value<double>(),      "set gain")
                        ("constant_speed",  program_options::value<double>(),      "set constant speed")
                        ("speed_threshold", program_options::value<double>(),      "set speed threshold")
                        ("left",            program_options::value<int>(),         "set left of initial rectangle")
                        ("top",             program_options::value<int>(),         "set top of initial rectangle")
                        ("right",           program_options::value<int>(),         "set right of initial rectangle")
                        ("bottom",          program_options::value<int>(),         "set bottom of initial rectangle")
                        ("front",           program_options::value<int>(),         "set front of initial rectangle (not used in 2D)")
                        ("back",            program_options::value<int>(),         "set back of initial rectangle (not used in 2D)")
                ;
        }
        
        void set_parameters(const program_options::variables_map& vm, Parameters& params)
        {
                params.wband_           = extract_parameter<int>(   "wband",           "invalid wband",           vm);
                params.wreset_          = extract_parameter<int>(   "wreset",          "invalid wreset",          vm);
                params.time_step_       = extract_parameter<double>("time_step",       "invalid time step",       vm);
                params.gain_            = extract_parameter<double>("gain",            "invalid gain",            vm);
                params.constant_speed_  = extract_parameter<double>("constant_speed",  "invalid constant speed",  vm);
                params.speed_threshold_ = extract_parameter<double>("speed_threshold", "invalid speed threshold", vm);
        }
}

std::unique_ptr<LevelSetMethodViewer2d> CommandLineInterface::viewer2d_;

void CommandLineInterface::execute_level_set_method(int argc, char* argv[])
{
        program_options::options_description desc("allowed options");
        set_common_arguments(desc);

        program_options::variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);

        if ( vm.empty() || vm.count("help") ) {
                std::cout << desc << "\n";
                return;
        }

        Parameters params;
        set_parameters(vm, params);
        
        int dim = extract_parameter<int>("dim", "invalid dimension", vm);
        if ( dim == 2 ) {
                execute_level_set_method_in_2d(argc, argv, vm, params);
                return;
        }
        if ( dim == 3 ) {
                execute_level_set_method_in_3d(argc, argv, vm, params);
                return;
        }
        Debug::debug("invalid dimension");
}

void CommandLineInterface::display_2d()
{
        viewer2d_->display();
}

void CommandLineInterface::keyboard_2d(unsigned char key, int x, int y)
{
        viewer2d_->keyboard(key, x, y);
}

/// execute the Level Set Method in 2D
void CommandLineInterface::execute_level_set_method_in_2d(
        int                                   argc,
        char*                                 argv[],
        const program_options::variables_map& vm,
        const Parameters&                     params
) {
        // set an initial front
        const int left   = extract_parameter<int>("left",   "invalid left",   vm);
        const int top    = extract_parameter<int>("top",    "invalid top",    vm);
        const int right  = extract_parameter<int>("right",  "invalid left",   vm);
        const int bottom = extract_parameter<int>("bottom", "invalid bottom", vm);
        InitialFront<TwoDimension> initial_front;
        initial_front.vertices_[0] = {{left, top}};
        initial_front.vertices_[1] = {{right, bottom}};

        // load an input image
        const std::string filename = extract_parameter<std::string>("input", "invalid input", vm);
        SpaceSize2d space_size;
        std::vector<std::uint8_t> pixels;
        load_input_image(filename, pixels, space_size);
        
        viewer2d_.reset(new LevelSetMethodViewer2d);
        viewer2d_->set_texture(&pixels[0]);

        if ( vm.count("verbose") ) {
                viewer2d_->set_verbose(true);
        }

        viewer2d_->initialize_viewer(argc, argv, space_size, "Level Set Method in 2D", display_2d, keyboard_2d);
        viewer2d_->initialize_level_set_method(params, space_size, &pixels[0], initial_front);
        viewer2d_->start();
}

namespace
{
        enum class Pattern
        {
                Pattern0,
                Pattern1,
                Pattern2,
                Pattern3,
                Pattern4,
                Pattern5,
                NonPattern,
        };
        
        inline void add_sphere(GeometryGenerator& generator, bool is_solid, const IntPoint3d& center, int radius)
        {
                generator.add_geometry(std::make_shared<Sphere>(is_solid, center, radius));
        }

        inline void add_cube(GeometryGenerator& generator, bool is_solid, const IntPoint3d& center, int size)
        {
                generator.add_geometry(std::make_shared<Cube>(is_solid, center, size));
        }
        
        void create_pattern0(
                GeometryGenerator&         geometry_generator,
                int                        w,
                int                        h,
                int                        d,
                int                        factor,
                std::vector<std::uint8_t>& src
        ) {
                const int radius {15 * factor};
                bool is_solid {true};
                add_sphere(geometry_generator, is_solid, {{    w / 4,     h / 4,     d / 4}}, radius);
                add_sphere(geometry_generator, is_solid, {{3 * w / 4, 3 * h / 4, 3 * d / 4}}, radius);
                geometry_generator.generate(src);
        }

        void create_pattern1(
                GeometryGenerator&         geometry_generator,
                int                        w,
                int                        h,
                int                        d,
                int                        factor,
                int                        wband,
                std::vector<std::uint8_t>& src
        ) {
                const int radius {15 * factor};
                bool is_solid {true};
                add_sphere(geometry_generator, is_solid,  {{w / 2,     h / 2, d / 4}}, radius);
                add_sphere(geometry_generator, is_solid,  {{w / 2, 3 * h / 4, d / 2}}, radius);
                add_sphere(geometry_generator, !is_solid, {{w / 2,     h / 2, d / 2}}, std::min(w, std::min(h, d)) / 2 - wband);
                geometry_generator.generate(src);
        }
        
        void create_pattern2(
                GeometryGenerator&         geometry_generator,
                int                        w,
                int                        h,
                int                        d,
                int                        factor,
                std::vector<std::uint8_t>& src
        ) {
                const int radius {15 * factor};
                bool is_solid {true};
                add_sphere(geometry_generator, is_solid, {{w / 4,     h / 4,     d / 4}}, radius);
                add_sphere(geometry_generator, is_solid, {{w / 4,     h / 4, 3 * d / 4}}, radius);
                add_sphere(geometry_generator, is_solid, {{w / 4, 3 * h / 4, 3 * d / 4}}, radius);
                add_sphere(geometry_generator, is_solid, {{w / 4, 3 * h / 4,     d / 4}}, radius);
                add_sphere(geometry_generator, is_solid, {{3 * w / 4, 3 * h / 4, 3 * d / 4}}, radius);
                add_sphere(geometry_generator, is_solid, {{3 * w / 4, 3 * h / 4,     d / 4}}, radius);
                add_sphere(geometry_generator, is_solid, {{3 * w / 4,     h / 4, 3 * d / 4}}, radius);
                add_sphere(geometry_generator, is_solid, {{3 * w / 4,     h / 4,     d / 4}}, radius);

                geometry_generator.generate(src);
        }
        
        void create_pattern3(
                GeometryGenerator&         geometry_generator,
                int                        w,
                int                        h,
                int                        d,
                int                        factor,
                int                        wband,
                std::vector<std::uint8_t>& src
        ) {
                const int radius {15 * factor};
                bool is_solid {true};
                
                add_sphere(geometry_generator, is_solid, {{w / 4,     h / 4,     d / 4}}, radius);
                add_sphere(geometry_generator, is_solid, {{w / 4,     h / 4, 3 * d / 4}}, radius);
                add_sphere(geometry_generator, is_solid, {{w / 4, 3 * h / 4, 3 * d / 4}}, radius);
                add_sphere(geometry_generator, is_solid, {{w / 4, 3 * h / 4,     d / 4}}, radius);
                
                add_sphere(geometry_generator, is_solid, {{3 * w / 4, 3 * h / 4, 3 * d / 4}}, radius);
                add_sphere(geometry_generator, is_solid, {{3 * w / 4, 3 * h / 4,     d / 4}}, radius);
                add_sphere(geometry_generator, is_solid, {{3 * w / 4,     h / 4, 3 * d / 4}}, radius);
                add_sphere(geometry_generator, is_solid, {{3 * w / 4,     h / 4,     d / 4}}, radius);

                add_sphere(geometry_generator, !is_solid, {{w / 2, h / 2, d / 2}}, std::min(w, std::min(h, d)) / 2 - wband);
                geometry_generator.generate(src);
        }

        void create_pattern4(
                GeometryGenerator&         geometry_generator,
                int                        w,
                int                        h,
                int                        d,
                int                        factor,
                std::vector<std::uint8_t>& src
        ) {
                const int radius {15 * factor};
                bool is_solid {true};
                add_cube(geometry_generator, is_solid, {{    w / 4,     h / 4,     d / 4}}, 2 * radius);
                add_cube(geometry_generator, is_solid, {{3 * w / 4, 3 * h / 4, 3 * d / 4}}, 2 * radius);
                geometry_generator.generate(src);
        }

        void create_pattern5(
                GeometryGenerator&         geometry_generator,
                int                        w,
                int                        h,
                int                        d,
                int                        factor,
                int                        wband,
                std::vector<std::uint8_t>& src
        ) {
                const int radius {15 * factor};
                bool is_solid {true};
                add_cube(geometry_generator, is_solid,  {{w / 2,     h / 2, d / 4}}, 2 * radius);
                add_cube(geometry_generator, is_solid,  {{w / 2, 3 * h / 4, d / 2}}, 2 * radius);
                add_cube(geometry_generator, !is_solid, {{w / 2,     h / 2, d / 2}}, 2* (std::min(w, std::min(h, d)) / 2 - wband));
                geometry_generator.generate(src);
        }

        /// create 3D objects
        void create_objects(
                Pattern                    pattern,
                const Parameters&          params,
                GeometryGenerator&         geometry_generator,
                int                        factor,
                std::vector<std::uint8_t>& src
        ) {
                const auto& space_size = geometry_generator.get_space_size();
                const int w = space_size.width_;
                const int h = space_size.height_;
                const int d = space_size.depth_;

                switch ( pattern ) {
                case Pattern::Pattern0:
                        create_pattern0(geometry_generator, w, h, d, factor, src);
                        break;
                
                case Pattern::Pattern1:
                        create_pattern1(geometry_generator, w, h, d, factor, params.wband_, src);
                        break;

                case Pattern::Pattern2:
                        create_pattern2(geometry_generator, w, h, d, factor, src);
                        break;
                        
                case Pattern::Pattern3:
                        create_pattern3(geometry_generator, w, h, d, factor, params.wband_, src);
                        break;
                        
                case Pattern::Pattern4:
                        create_pattern4(geometry_generator, w, h, d, factor, src);
                        break;

                case Pattern::Pattern5:
                        create_pattern5(geometry_generator, w, h, d, factor, params.wband_, src);
                        break;

                default:
                        break;
                }
        }
        
        inline const Pattern extract_pattern(const std::string& input)
        {
                if ( input == "pattern0" ) {
                        return Pattern::Pattern0;
                }
                if ( input == "pattern1" ) {
                        return  Pattern::Pattern1;
                }
                if ( input == "pattern2" ) {
                        return  Pattern::Pattern2;
                }
                if ( input == "pattern3" ) {
                        return Pattern::Pattern3;
                }
                if ( input == "pattern4" ) {
                        return Pattern::Pattern4;
                }
                if ( input == "pattern5" ) {
                        return Pattern::Pattern5;
                }
                
                return Pattern::NonPattern;
        }
}

std::unique_ptr<LevelSetMethodViewer3d> CommandLineInterface::viewer3d_;
std::unique_ptr<GeometryLoader>         CommandLineInterface::geometry_loader_;
std::unique_ptr<GeometryGenerator>      CommandLineInterface::geometry_generator_;

void CommandLineInterface::display_3d()
{
        viewer3d_->display();
}

void CommandLineInterface::keyboard_3d(unsigned char key, int x, int y)
{
        viewer3d_->keyboard(key, x, y);
}

/// execute the Level Set Method in 3D
void CommandLineInterface::execute_level_set_method_in_3d(
        int                                   argc,
        char*                                 argv[],
        const program_options::variables_map& vm,
        const Parameters&                     params
) {
        // set an initial front
        const int left   = extract_parameter<int>("left",   "invalid left",   vm);
        const int top    = extract_parameter<int>("top",    "invalid top",    vm);
        const int right  = extract_parameter<int>("right",  "invalid left",   vm);
        const int bottom = extract_parameter<int>("bottom", "invalid bottom", vm);
        const int front  = extract_parameter<int>("front",  "invalid front",  vm);
        const int back   = extract_parameter<int>("back",   "invalid back",   vm);
        InitialFront<ThreeDimension> initial_front;
        initial_front.vertices_[0] = {{left, top, front}};
        initial_front.vertices_[1] = {{right, bottom, back}};
        
        // set a 3D-space size
        constexpr int factor = 2;
        constexpr int width  = factor * 100;
        constexpr int height = factor * 100;
        constexpr int depth  = factor * 100;
        SpaceSize<ThreeDimension> space_size {width, height, depth};

        const std::string input = extract_parameter<std::string>("input", "invalid input", vm);
        Pattern pattern = extract_pattern(input);
        
        viewer3d_.reset(new LevelSetMethodViewer3d);
        std::vector<std::uint8_t> pixels(space_size.total_, 255);
        
        // create 3D objects
        if ( pattern != Pattern::NonPattern ) {
                geometry_generator_.reset(new GeometryGenerator(space_size));
                create_objects(pattern, params, *geometry_generator_, factor, pixels);
                viewer3d_->set_geometries(&geometry_generator_->get_geometries());
        // load 3D objects
        } else {
                geometry_loader_.reset(new GeometryLoader());
                geometry_loader_->read_file(input);
                geometry_loader_->scale_scene(width * 0.85);
                const FloatPoint3d center {{width/2, height/2, depth/2}};
                geometry_loader_->move_to_center(center);
                geometry_loader_->convert(space_size, pixels);
                viewer3d_->set_geometries(&geometry_loader_->get_geometries());
        }

        if ( vm.count("verbose") ) {
                viewer3d_->set_verbose(true);
        }
        
        viewer3d_->initialize_viewer(argc, argv, space_size, "Level Set Method in 3D", display_3d, keyboard_3d);
        viewer3d_->initialize_level_set_method(params, space_size, &pixels[0], initial_front);
        viewer3d_->start();
}