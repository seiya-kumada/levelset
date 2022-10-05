//------------------------------------------------------------------------------------------------------------------------------------
#ifndef TESTUTILITIES_H_
#define TESTUTILITIES_H_
//------------------------------------------------------------------------------------------------------------------------------------
#include <vector>
#include "Point.h"
#include <string>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/gil/utilities.hpp>
#include <boost/gil/typedefs.hpp>
#include "Typedef.h"
#include <list>
//------------------------------------------------------------------------------------------------------------------------------------

namespace kmd {

        struct Parameters;

        class TestUtilities {
        public:
                static void create_front(Front& closed_points, double radius, int step_size, int center_x, int center_y);
                static void create_front(Front& closed_points, std::ifstream& ifs);
                static bool write_image(const std::string& filename, const std::vector<double>& image, int w, int h, int space = 2);

                static bool write_images(const std::string& filename2d, const std::string& filename3d, const std::vector<double>& image, int w, int h);
                static bool write_tube(const std::string& filename, const Tube& tube);
                static void calculate_exact_signed_distance(std::vector<double>& distance, int w, int h);
                static void calculate_exact_distance(std::vector<double>& distance, int w, int h);
                static void calculate_exact_distance2(std::vector<double>& distance, int w, int h);
                static void calculate_exact_curvature(std::vector<double>& curvature, int w, const Tube& tube);
                static const std::string str(const char* msg, int count);

                static void create_gray_image(
                        std::vector<double>& gray_image, 
                        std::vector<Point>& closed_points, 
                        int w, 
                        int h, 
                        int major_radius, 
                        int minor_radius, 
                        int step_size, 
                        int center_x, 
                        int center_y
                );
                static void create_closed_points(std::vector<Point>& closed_points, int major_radius, int minor_radius, int step_size, int center_x, int center_y);
                //static void print_ctr(const char* message);
                //static void print_dtr(const char* message);
                static void read_buffer(int* w, int* h, std::vector<double>& image, const char* file_name);
                static void read_image(int* w, int* h, std::vector<double>& image, const char* file_name);
                static void delete_front(Front& front);
                template<typename T>
                static bool print_any(const std::string& message, T t) {
                        std::cout << message << " : " << boost::lexical_cast<std::string>(t) << std::endl;
                        return true;
                }

                static bool print_point(const Point& p) {
                        std::cout << p.x() << " " << p.y() << std::endl;
                        return true;
                }
                static bool write_front(std::ofstream& ofs, const Front& front);
                static bool write_front(const char* filename, const Front& front);
                static bool write_edge_front(const char* filename, const std::vector<IntPoint>& edge_front);
                static bool write_fronts(std::ofstream& ofs, const Fronts& fronts);
                static bool write_fronts(const char* filename, const Fronts& fronts);
                static bool write_buffer(const std::string& path, const std::vector<double>& buffer, int w, int h);
                static bool write_intermediate_image(
                        const std::string&              folder_path,
                        int                             count,
                        const Fronts&         fronts,
                        const double*                   ptr,
                        const std::vector<double>&      buf0,
                        const std::vector<double>&      buf1,
                        int                             w,
                        int                             h,
                        int                             space
                );

                static bool write_intermediate_fronts(
                        const boost::gil::gray8c_view_t&        view,
                        const std::string&                      folder_path,
                        int                                     count,
                        const Fronts&                 fronts,
                        //kumada
                        //boost::gil::bits8                       color
                        std::uint8_t                       color
                );

                static bool write_next_distance(
                        const std::string&              folder_path,
                        int                             count,
                        const double*                   ptr,
                        const std::vector<double>&      buf0,
                        const std::vector<double>&      buf1,
                        int                             w,
                        int                             h,
                        int                             space
                );

                static void read_parameters(std::ifstream& ifs, Parameters* parameters, int w, int h);


        };
}

#endif /*TESTUTILITIES_H_*/
