//------------------------------------------------------------------------------------------------------------------------------------
#include "TestUtilities.h"
#include "Front.h"
#include <string>
#include <cmath>
#include <fstream>
#include "Error.h"
#include "Ellipse.h"
#include <stdexcept>
#include <boost/gil/image_view_factory.hpp>
#include <boost/gil/extension/io/jpeg.hpp>
#include <boost/gil/typedefs.hpp>
#include <boost/algorithm/minmax_element.hpp>
#include <boost/tokenizer.hpp>
#include "Parameters.h"
#include "Pi.h"
//------------------------------------------------------------------------------------------------------------------------------------
using namespace std;
using namespace kmd;
using namespace boost::gil;
using namespace boost;

namespace {

        void show_debug(const string& str) {
                cout << "<<DEBUG>> : " << str << endl;
        }

        std::size_t get_file_size(ifstream& ifs) {
                ifs.seekg(0, ios::end);
                std::size_t data_size = ifs.tellg();
                ifs.seekg(0, ios::beg);
                return data_size;
        }
}


//====================================================================================================================================
const string TestUtilities::str(const char* msg, int count) {
//====================================================================================================================================
        string s(msg);
        return msg + lexical_cast<string>(count);
}

//====================================================================================================================================
void TestUtilities::create_front(Front& closed_points, double radius, int step_size, int center_x, int center_y) {
//====================================================================================================================================
//        show_debug("TestUtilities::create_front");
        test::Front front;
        front.set_radius(radius);
        double step = 2 * M_PI / step_size;
        closed_points.reserve(step_size);
        front.set_center(center_x, center_y);
        double radian;
        for ( int i = 0; i < step_size; ++i ) {
                radian = i * step;
                closed_points.push_back(new Point(front.x(radian), front.y(radian)));
        }
}

//====================================================================================================================================
void TestUtilities::create_front(Front& closed_points, ifstream& ifs) {
//====================================================================================================================================
        std::size_t file_size = get_file_size(ifs);

        vector<char> data(file_size);
        ifs.read(&data[0], file_size);

        char_separator<char> sep(" \t\n");
        string file_str(&data[0]);
        tokenizer<char_separator<char> > tokens(file_str, sep);

        for ( tokenizer<char_separator<char> >::iterator i = tokens.begin(), e = tokens.end(); i != e; ++i ) {
                double x = lexical_cast<double>(*i);
                ++i;
                double y = lexical_cast<double>(*i);
                closed_points.push_back(new Point(x, y));
        }

}

//====================================================================================================================================
bool TestUtilities::write_image(const string& filename, const vector<double>& image, int w, int h, int space) {
//====================================================================================================================================
        show_debug("TestUtilities::write_image");
        ofstream ofs(filename.c_str());
        for ( int j = 0; j < h; j += space ) {
                for ( int i = 0; i < w; i += space ) {
                        ofs << i << " " << j << " " << image[i + w * j] << endl;
                }
        }
        return true;
}

namespace {

        //****************************************************************************************************************************
        struct reduce_pixel : public unary_function<double, unsigned char> {
        //****************************************************************************************************************************
                reduce_pixel(double min_v, double max_v)
                : min_v_(min_v), max_v_(max_v) {}

                unsigned char operator()(double x) const {
                        return clip(255 * (x - min_v_) / (max_v_ - min_v_));
                }

                unsigned char clip(double x) const {
                        if ( x < 0 ) {
                                return 0;
                        }
                        if ( x > 255 ) {
                                return 255;
                        }
                        return static_cast<unsigned char>(x);
                }

        private:
                double min_v_;
                double max_v_;
        };
}

namespace {
        struct inverse : unary_function<double, double> {
                void operator()(double& x) const {
                        x = 1.0 / (1.0 + x);
                }
        };
}


//====================================================================================================================================
bool TestUtilities::write_buffer(const string& folder_path, const vector<double>& buffer, int w, int h) {
//====================================================================================================================================
        show_debug("TestUtilities::write_buffer");

        vector<double> tmp0(buffer);
        for_each(tmp0.begin(), tmp0.end(), inverse());

        write_image(folder_path + "/modified", tmp0, w, h, 2);

        typedef vector<double>::const_iterator Diter;
        pair<Diter, Diter> r = boost::minmax_element(tmp0.begin(), tmp0.end());

        cout << "min : " << *r.first << ", max : " << *r.second << endl;

        vector<unsigned char> tmp1(tmp0.size());
        transform(tmp0.begin(), tmp0.end(), tmp1.begin(), reduce_pixel(*r.first, *r.second));
        gray8c_view_t v = interleaved_view(w, h, (gray8c_pixel_t*)&tmp1[0], w);
        gil::write_view(folder_path + "/modified.jpg", v, jpeg_tag{});
        return true;
}



//============================================================================================================================
bool TestUtilities::write_next_distance(
        const string&           folder_path,
        int                     count,
        const double*           ptr,
        const vector<double>&   buf0,
        const vector<double>&   buf1,
        int                     w,
        int                     h,
        int                     space
) {
//============================================================================================================================
        //show_debug("TestUtilities::write_next_distance");
        if ( ptr == &buf0[0] ) {
                write_image(folder_path + (string("/dist_") + boost::lexical_cast<string>(count)).c_str(), buf0, w, h, space);
        } else {
                write_image(folder_path + (string("/dist_") + boost::lexical_cast<string>(count)).c_str(), buf1, w, h, space);
        }
        return true;
}



//====================================================================================================================================
bool TestUtilities::write_intermediate_image(
        const string&           folder_path,
        int                     count,
        const Fronts&      fronts,
        const double*           ptr,
        const vector<double>&   buf0,
        const vector<double>&   buf1,
        int                     w,
        int                     h,
        int                     space
) {
//====================================================================================================================================
        show_debug("TestUtilities::write_intermediate_image");
        string str(folder_path + "fronts_");
        str += boost::lexical_cast<string>(count);
        ofstream ofs(str.c_str());
        TestUtilities::write_fronts(ofs, fronts);
        write_next_distance(folder_path, count, ptr, buf0, buf1, w, h, space);
        return true;
}

//====================================================================================================================================
bool TestUtilities::write_fronts(const char* filename, const Fronts& fronts) {
//====================================================================================================================================
        ofstream ofs(filename);
        return write_fronts(ofs, fronts);
}

//====================================================================================================================================
bool TestUtilities::write_fronts(ofstream& ofs, const Fronts& fronts) {
//====================================================================================================================================
        show_debug("TestUtilities::write_fronts");
        Fronts::const_iterator beg = fronts.begin();
        Fronts::const_iterator end = fronts.end();
        while ( beg != end ) {
                write_front(ofs, *beg);
                ofs << endl;
                ofs << endl;
                ++beg;
        }
        return true;
}

namespace {

        //============================================================================================================================
        //kumada
        //void write_intermediate_front(const gray8_view_t& v, const Front& front, bits8 color) {
        void write_intermediate_front(const gray8_view_t& v, const Front& front, std::uint8_t color) {
        //============================================================================================================================
                for ( std::size_t i = 0, n = front.size(); i < n; ++i ) {
                        v(static_cast<int>(front[i]->x() + 0.5), static_cast<int>(front[i]->y() + 0.5)) = color;
                }
        }
}

//====================================================================================================================================
bool TestUtilities::write_intermediate_fronts(
        const gray8c_view_t&    v,
        const string&           folder_path,
        int                     count,
        //const vector<Front>&    fronts,
        const Fronts&      fronts,
        //kumada
        //bits8                   color
        std::uint8_t                   color

) {
//====================================================================================================================================
        show_debug("TestUtilities::write_intermediate_fronts");
        gray8_image_t img(v.dimensions());
        copy_pixels(v, view(img));
        Fronts::const_iterator beg = fronts.begin();
        Fronts::const_iterator end = fronts.end();
        while ( beg != end ) {
                write_intermediate_front(view(img), *beg, color);
                ++beg;
        }
        write_view(folder_path + "/front_" + lexical_cast<string>(count) + ".jpg", const_view(img), jpeg_tag{});
        return true;
}


//====================================================================================================================================
bool TestUtilities::write_images(const string& filename2d, const string& filename3d, const vector<double>& image, int w, int h) {
//====================================================================================================================================
        show_debug("TestUtilities::write_images");
        ofstream ofs2(filename2d.c_str());
        ofstream ofs3(filename3d.c_str());

        ofs2.precision(20);
        ofs2.setf(ios::fixed, ios::floatfield);
        ofs3.precision(20);
        ofs3.setf(ios::fixed, ios::floatfield);

        for ( int j = 0; j < h; ++j ) {
                for ( int i = 0; i < w; ++i ) {
                        double v = image[i + w * j];
                        if ( v == 0.0 ) {
                                ofs2 << i << " " << j << endl;
                        }
                        ofs3 << i << " " << j << " " << v << endl;
                }
        }
        return true;
}

//====================================================================================================================================
void TestUtilities::calculate_exact_signed_distance(vector<double>& distance, int w, int h) {
//====================================================================================================================================
        show_debug("TestUtilities::calculate_exact_signed_distance");
        for ( int j = 0, wj = w * j; j < h; ++j, wj += w ) {
                double* tmp = &distance[wj];
                for ( int i = 0; i < w; ++i ) {
                        double x = i - 100.0;
                        double y = j - 100.0;
                        tmp[i] = -80.0 + sqrt(x*x + y*y);
                }
        }
}

//====================================================================================================================================
void TestUtilities::calculate_exact_distance(vector<double>& distance, int w, int h) {
//====================================================================================================================================
        show_debug("TestUtilities::calculate_exact_distance");
        for ( int j = 0, wj = w * j; j < h; ++j, wj += w ) {
                double* tmp = &distance[wj];
                for ( int i = 0; i < w; ++i ) {
                        double x = i - 100.0;
                        double y = j - 100.0;
                        tmp[i] = -80.0 + sqrt(x*x + y*y);
                        if ( tmp[i] < 0.0 ) {
                                tmp[i] *= -1.0;
                        }
                }
        }
}

//====================================================================================================================================
void TestUtilities::calculate_exact_distance2(vector<double>& distance, int w, int h) {
//====================================================================================================================================
        show_debug("TestUtilities::calculate_exact_distance2");
        for ( int j = 0, wj = w * j; j < h; ++j, wj += w ) {
                double* tmp = &distance[wj];
                for ( int i = 0; i < w; ++i ) {
                        double x = i - 50;
                        double y = j - 50;
                        double d0 = fabs(sqrt(Utilities::square(x) + Utilities::square(y)) - 30);

                        x = i - 150;
                        y = j - 150;
                        double d1 = fabs(sqrt(Utilities::square(x) + Utilities::square(y)) - 30);

                        tmp[i] = d0 < d1 ? d0 : d1;
                }
        }
}

//====================================================================================================================================
void TestUtilities::calculate_exact_curvature(std::vector<double>& curvature, int w, const Tube& tube) {
//====================================================================================================================================
        show_debug("TestUtilities::calculate_exact_curvature");
        Tube::const_iterator beg = tube.begin();
        Tube::const_iterator end = tube.end();
        while ( beg != end ) {
                int j = beg->first;
                int wj = w * j;
                const Segments& segments = beg->second;
                for ( std::size_t k = 0, n = segments.size(); k < n; ++k ) {
                        const Segment& segment = segments[k];
                        int first = segment.first;
                        int second = segment.second;
                        for ( int i = first; i <= second; ++i ) {
                                int index = i + wj;
                                double x = i - 100.0;
                                double y = j - 100.0;
                                curvature[index] = 1.0 / sqrt(x * x + y * y);
                        }
                }
                ++beg;
        }
 }

//====================================================================================================================================
void TestUtilities::create_gray_image(vector<double>& gray_image, vector<Point>& closed_points, int w, int h, int major_radius, int minor_radius, int step_size, int center_x, int center_y) {
//====================================================================================================================================
        show_debug("TestUtilities::create_gray_image");
        Ellipse ellipse;
        ellipse.set_radius(major_radius, minor_radius);
        ellipse.set_center(center_x, center_y);
        gray_image.resize(w * h, 255);
        double step = 2 * M_PI / step_size;
        double radian;
        for ( int i = 0; i < step_size; ++i ) {
                radian = i * step;
                int x = static_cast<int>(ellipse.x(radian) + 0.5);
                int y = static_cast<int>(ellipse.y(radian) + 0.5);
                closed_points.push_back(Point(x, y));
                gray_image[x + y * w] = 0;
        }
}

//====================================================================================================================================
void TestUtilities::create_closed_points(vector<Point>& closed_points, int major_radius, int minor_radius, int step_size, int center_x, int center_y) {
//====================================================================================================================================
        show_debug("TestUtilities::create_closed_points");
        Ellipse ellipse;
        ellipse.set_radius(major_radius, minor_radius);
        ellipse.set_center(center_x, center_y);
        double step = 2 * M_PI / step_size;
        double radian;
        for ( int i = 0; i < step_size; ++i ) {
                radian = i * step;
                closed_points.push_back(Point(ellipse.x(radian), ellipse.y(radian)));
        }
}

//====================================================================================================================================
void TestUtilities::delete_front(Front& front) {
//====================================================================================================================================
//        show_debug("TestUtilities::delete_front");
        for ( std::size_t i = 0, n = front.size(); i < n; ++i ) {
                delete front[i];
        }
}

//====================================================================================================================================
void TestUtilities::read_buffer(int* w, int* h, vector<double>& image, const char* file_name) {
//====================================================================================================================================
        show_debug("TestUtilities::read_buffer");
        ifstream ifs(file_name);
        if ( ifs.fail() ) {
                throw runtime_error(Error::message("TestUtilities::read_buffer", __LINE__, __FILE__));
        }
        ifs >> *w;
        ifs >> *h;
        int v;
        int size = (*w) * (*h);
        image.reserve(size);
        for ( int i = 0; i < size; ++i ) {
                ifs >> v;
                image.push_back(v);
        }
}



//====================================================================================================================================
void TestUtilities::read_parameters(ifstream& ifs, Parameters* parameters, int w, int h) {
//====================================================================================================================================
//        show_debug("TestUtilities::read_parameters");
        std::size_t data_size = get_file_size(ifs);
        vector<char> data(data_size);
        ifs.read(&data[0], data_size);

        char_separator<char> sep(" \t\n");
        string data_str(&data[0]);
        tokenizer<char_separator<char> > tokens(data_str, sep);

        typedef tokenizer<char_separator<char> >::iterator Iter;
        parameters->image_width_ = w;
        parameters->image_height_ = h;
        for ( Iter i = tokens.begin(), end = tokens.end(); i != end; ++i ) {

                if ( *i == "time_step" ) {
                        ++i;
                        parameters->set_time_step(lexical_cast<double>(*i));
                }

                if ( *i == "time_step_number" ) {
                        ++i;
                        parameters->time_step_number_ = lexical_cast<int>(*i);
                }

                if ( *i == "space_step" ) {
                        ++i;
                        parameters->space_step_ = lexical_cast<double>(*i);
                }

                if ( *i == "constant_speed" ) {
                        ++i;
                        parameters->constant_speed_ = lexical_cast<double>(*i);
                }

                if ( *i == "epsilon" ) {
                        ++i;
                        parameters->epsilon_ = lexical_cast<double>(*i);
                }

                if ( *i == "sigma" ) {
                        ++i;
                        parameters->sigma_ = lexical_cast<double>(*i);
                }
        }
}

//====================================================================================================================================
void TestUtilities::read_image(int* w, int* h, vector<double>& image, const char* file_name) {
//====================================================================================================================================
        show_debug("TestUtilities::read_image");
        gray8_image_t img;
        gil::read_image(file_name, img, jpeg_tag{});
        copy(const_view(img).begin(), const_view(img).end(), back_inserter(image));
        *w = static_cast<int>(img.width());
        *h = static_cast<int>(img.height());
}

#if 0
//====================================================================================================================================
void TestUtilities::print_ctr(const char*) {
//====================================================================================================================================
        //cout << message << ": constructor\n";
}

//====================================================================================================================================
void TestUtilities::print_dtr(const char*) {
//====================================================================================================================================
        //cout << message << ": destructor\n";
}
#endif

//====================================================================================================================================
bool TestUtilities::write_front(const char* filename, const Front& front) {
//====================================================================================================================================
        ofstream ofs(filename);
        return write_front(ofs, front);
}

//====================================================================================================================================
bool TestUtilities::write_front(ofstream& ofs, const Front& front) {
//====================================================================================================================================
        show_debug("TestUtilities::write_front");
        for ( std::size_t i = 0, n = front.size(); i < n; ++i ) {
                const Point* p = front[i];
                ofs << p->x() << " " << p->y() << endl;

        }
        return true;
}

//====================================================================================================================================
bool TestUtilities::write_edge_front(const char* filename, const std::vector<IntPoint>& edge_front) {
//====================================================================================================================================
        show_debug("TestUtilities::write_edge_front");
        ofstream ofs(filename);
        for ( std::size_t i = 0, n = edge_front.size(); i < n; ++i ) {
                const IntPoint& p = edge_front[i];
                ofs << p.x << " " << p.y << endl;
        }
        return true;
}

//====================================================================================================================================
bool TestUtilities::write_tube(const std::string& filename, const Tube& tube) {
//====================================================================================================================================
        show_debug("TestUtilities::write_tube");
        ofstream ofs(filename.c_str());
        Tube::const_iterator beg = tube.begin();
        Tube::const_iterator end = tube.end();
        while ( beg != end ) {
                int j = beg->first;
                const Segments& segments = beg->second;
                for ( std::size_t k = 0, n = segments.size(); k < n; ++k ) {
                        const Segment& segment = segments[k];
                        int first = segment.first;
                        int second = segment.second;
                        for ( int i = first; i <= second; ++i ) {
                                ofs << i << " " << j << endl;
                        }
                }
                ++beg;
        }
        return true;
}
