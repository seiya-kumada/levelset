#ifdef UNIT_TEST_main
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <kmd/without_warning/unit_test.h>
#else
//----------------------------------------------------------------------------------------------------------------------------------------------------
#include <boost/program_options.hpp>
#include <string>
#include <iostream>
//kumada
#include <boost/gil/extension/io/jpeg.hpp>
#include <iomanip>
#include "TestUtilities.h"
#include "LevelSetMethod.h"
#include "Parameters.h"
#include <fstream>
#include "Pi.h"
#include <boost/chrono.hpp>
//----------------------------------------------------------------------------------------------------------------------------------------------------
using namespace boost::program_options;
using namespace std;
using namespace boost::gil;
using namespace boost;
using namespace kmd;

namespace {

        //============================================================================================================================================
        void run_with_mode_info(variables_map& vm) {
        //============================================================================================================================================
                const string file_path = vm.count("input") ? vm["input"].as<string>() : "";
                gray8_image_t img;
                read_image(file_path, img, jpeg_tag{});
                
                cout << setw(10) << "width : " << setw(5) << img.width() << endl;
                cout << setw(10) << "height : " << setw(5) << img.height() << endl;
        }

        //============================================================================================================================================
        void write_loop_to_view(const gray8_view_t& view, int cx, int cy, int r, int color) {
        //============================================================================================================================================

                const int size = 1000;
                double step = 2.0 * M_PI / size;
                for ( int i = 0; i < size; ++i ) {
                        int x = static_cast<int>(r * cos(step * i) + cx + 0.5);
                        int y = static_cast<int>(r * sin(step * i) + cy + 0.5);
                        view(x, y) = color;
                }
        }

        //============================================================================================================================================
        void get_loop_parameters(int* x, int* y, int* r, variables_map& vm) {
        //============================================================================================================================================
                *x = vm.count("x") ? vm["x"].as<int>() : -1;
                *y = vm.count("y") ? vm["y"].as<int>() : -1;
                *r = vm.count("radius") ? vm["radius"].as<int>() : -1;
                if ( *x == -1 || *y == -1 || *r == -1 ) {
                        throw runtime_error("invalid parameters for the initial loop");
                }
        }

        //============================================================================================================================================
        void write_circle_parameters(const string& path, int x, int y, int r, int color) {
        //============================================================================================================================================
                ofstream ofs(path.c_str());
                ofs << setw(10) << "x : " << setw(10) << x << endl;
                ofs << setw(10) << "y : " << setw(10) << y << endl;
                ofs << setw(10) << "r : " << setw(10) << r << endl;
                ofs << setw(10) << "color : " << setw(10) << color << endl;
        }

        //============================================================================================================================================
        void run_with_mode_initial(variables_map& vm) {
        //============================================================================================================================================
                int x, y, r;
                get_loop_parameters(&x, &y, &r, vm);
                int color = vm.count("color") ? vm["color"].as<int>() : 255;

                const string file_path = vm.count("input") ? vm["input"].as<string>() : "";
                gray8_image_t img;
                read_image(file_path, img, jpeg_tag{});
                write_loop_to_view(view(img), x, y, r, color);

                const string folder_path = vm.count("output") ? vm["output"].as<string>() : "./";
                write_view(folder_path + "/initial_loop.jpg", const_view(img), jpeg_tag{});
                write_circle_parameters(folder_path + "/circle_parameters.txt", x, y, r, color);
        }

        //=====================================================================================================================
        void show_parameters(const Parameters& parameters) {
        //=====================================================================================================================
                cout << setw(20) << "time step:"        << setw(10) << parameters.time_step()           << endl;
                cout << setw(20) << "time step number:" << setw(10) << parameters.time_step_number_     << endl;
                cout << setw(20) << "space step:"       << setw(10) << parameters.space_step_           << endl;
                cout << setw(20) << "width:"            << setw(10) << parameters.image_width_          << endl;
                cout << setw(20) << "height:"           << setw(10) << parameters.image_height_         << endl;
                cout << setw(20) << "constant speed:"   << setw(10) << parameters.constant_speed_       << endl;
                cout << setw(20) << "epsilon:"          << setw(10) << parameters.epsilon_              << endl;
                cout << setw(20) << "sigma:"            << setw(10) << parameters.sigma_                << endl;
        }

        //============================================================================================================================================
        void run_with_mode_final(variables_map& vm) {
        //============================================================================================================================================
                int x, y, r;
                get_loop_parameters(&x, &y, &r, vm);
                
                const string file_path = vm.count("input") ? vm["input"].as<string>() : "";
                gray8_image_t img;
                read_image(file_path, img, jpeg_tag{});
                
                vector<double> buffer;
                copy(const_view(img).begin(), const_view(img).end(), back_inserter(buffer));

                const string config_path = vm.count("config") ? vm["config"].as<string>() : "";
                ifstream ifs(config_path.c_str());
                if ( ifs.fail() ) {
                        throw runtime_error("invalid config file path");
                }

                Parameters parameters;
                TestUtilities::read_parameters(ifs, &parameters, static_cast<int>(const_view(img).width()), static_cast<int>(const_view(img).height()));
                show_parameters(parameters);

                vector<Point*> front;
                int step_size = 1000;

                LevelSetMethod method(&parameters);
                TestUtilities::create_front(front, r, step_size, x, y);

                const string folder_path = vm.count("output") ? vm["output"].as<string>() : "./";
                int color = vm.count("color") ? vm["color"].as<int>() : 255;
                int interval = vm.count("interval") ? vm["interval"].as<int>() : 20;
                const boost::chrono::steady_clock::time_point start = boost::chrono::steady_clock::now();
                method.run(buffer, front, folder_path, const_view(img), std::uint8_t(color), interval);
                const boost::chrono::steady_clock::time_point end = boost::chrono::steady_clock::now();
                std::cout << boost::chrono::duration<double, boost::milli>(end - start).count() << "[ms]\n";
                TestUtilities::delete_front(front);
        }
}

//====================================================================================================================================================
int main(int argc, char* argv[]) {
//====================================================================================================================================================
        try {
                options_description desc("allowed options");
                desc.add_options()
                        ("help,h", "produce help message")
                        ("input,I", value<string>(), "set input GRAY image path")
                        ("output,O", value<string>(), "set output folder path")
                        ("mode,M", value<string>(), "set mode to \"info\", \"initial\", or \"final\"")
                        ("x,X", value<int>(), "set x of circle center")
                        ("y,Y", value<int>(), "set y of circle center")
                        ("radius,R", value<int>(), "set radius of circle")
                        ("color,C", value<int>(), "set gray level within [0,255] to draw circle")
                        ("config,F", value<string>(), "set config file path")
                        ("interval,T", value<int>(), "set interval to save images")
                ;

                variables_map vm;
                store(parse_command_line(argc, argv, desc), vm);
                notify(vm);

                if ( vm.empty() || vm.count("help") ) {
                        cout << desc << "\n";
                        return 1;
                }

                string mode = vm.count("mode") ? vm["mode"].as<string>() : "";
                if ( mode == "info" ) {
                        run_with_mode_info(vm);
                        return 1;
                }
                if ( mode == "initial" ) {
                        run_with_mode_initial(vm);
                        return 1;
                }
                if ( mode == "final" ) {
                        run_with_mode_final(vm);
                        return 1;
                }

        }
        catch (const unknown_option& error) {
                cout << "<<unknown_option>> " << error.what() << endl;
        }
        catch (const invalid_option_value& error) {
                cout << "<<invalid_option_value>> " << error.what() << endl;
        }
        catch (const ios_base::failure& error) {
                cout << "<<ios_base::failure>> " << error.what() << endl;
        }
        catch (const runtime_error& error) {
                cout << "<<runtime_error>> " << error.what() << endl;
        }
        return 0;
}
#endif
