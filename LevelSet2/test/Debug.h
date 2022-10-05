#ifndef LevelSetMethod2_Debug_h
#define LevelSetMethod2_Debug_h

#include "Point.h"
#include "SpaceSize.h"
#include "Status.h"
#include "Front.h"
#include "Grid.h"
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>
#include <boost/range/algorithm/for_each.hpp>
#include <opencv2/opencv.hpp>

namespace lsm
{

/// static functions for debugging
class Debug
{
public:
        static void debug(const std::string& message)
        {
                std::cout << message << std::endl;
        }

        template<typename D>
        static void display_distance_map(const std::multimap<double, IntPoint<D>>& map)
        {
                typedef typename std::multimap<double, IntPoint<D>>::value_type value_type;
                std::cout << map.size() << std::endl;
                auto beg = map.begin();
                auto end = map.end();
                int index = 0;
                while ( beg != end ) {
                        double key = beg->first;
                        std::cout << key << " : ";
                        auto range = map.equal_range(key);
                        const std::ptrdiff_t size = std::distance(range.first, range.second);
                        std::cout << "[" << size << "]: ";
                        int k = 0;
                        boost::for_each(range,
                                [&](const value_type& v)
                                {
                                        const auto& p = v.second;
                                        display_point<D>(p);
                                        ++k;
                                }
                        );
                        std::cout << std::endl;
                        beg = range.second;
                        ++index;
                }
        }

        template<typename T, typename U>
        static void display_buffer(const SpaceSize<TwoDimension>& size, const std::vector<T>& buffer)
        {
                std::cout << "{\n";
                for ( int j = 0, wj = j * size.width_; j < size.height_; ++j, wj += size.width_ ) {
                        const T* p = &buffer[wj];
                        for ( int i = 0; i < size.width_; ++i ) {
                                std::cout << std::setw(10) << static_cast<U>(p[i]) << ", ";
                        }
                        std::cout << std::endl;
                }
                std::cout << "}\n";
        }
        
        static void draw_phi(const std::vector<double>& phi, int w, int h)
        {
                cv::Mat dst(h, w, CV_8UC1);
                for ( int j = 0, wj = w * j; j < h; ++j, wj += w ) {
                        std::uint8_t* p = dst.ptr<std::uint8_t>(j);
                        const double* s = &phi[wj];
                        for ( int i = 0; i < w; ++i ) {
                                const double& si = s[i];
                                if ( si < 0.0 ) {
                                        p[i] = 100;
                                } else {
                                        p[i] = 250;
                                }
                        }
                }
                
                cv::imwrite("/Users/kumada/Data/level-set-method/phi.jpg", dst);
        }

        static void display_grid(const Grid<TwoDimension>& grid)
        {
                std::cout << "(left, top, right, bottom) = ("
                        << grid.left_ << ", " << grid.top_ << ", "
                        << grid.right_ << ", " << grid.bottom_ << ")\n";
        }
        static void display_grid(const Grid<ThreeDimension>& grid)
        {
        }
        
        static void draw_statuses(const std::vector<Status>& statuses, const SpaceSize<TwoDimension>& size)
        {
                const int w = size.width_;
                const int h = size.height_;
                cv::Mat dst(h, w, CV_8UC1);
                for ( int j = 0, wj = w * j; j < h; ++j, wj += w ) {
                        std::uint8_t* p = dst.ptr<std::uint8_t>(j);
                        const Status* s = &statuses[wj];
                        for ( int i = 0; i < w; ++i ) {
                                const Status& si = s[i];
                                if ( si == Status::Farway ) {
                                        p[i] = 0;
                                } else if ( si == Status::Band ) {
                                        p[i] = 100;
                                } else if ( si == Status::ResetBand ) {
                                        p[i] = 180;
                                } else { // Status::Front
                                        p[i] = 255;
                                }
                        }
                }
                cv::imwrite("/Users/kumada/Data/level-set-method/statuses/status_" + std::to_string(index_) + ".jpg", dst);
                ++index_;
        }
        
        static int index_;
        static void draw_contour(const Front<TwoDimension>& front, const SpaceSize<TwoDimension>& size)
        {
                cv::Mat img(size.height_, size.width_, CV_8UC1);
                std::fill(img.data, img.data + size.width_ * size.height_, 0);
                for ( const auto& p : front ) {
                        cv::Point p0(static_cast<int>(p[0]), static_cast<int>(p[1]));
                        cv::Point p1(static_cast<int>(p[0]) + 1, static_cast<int>(p[1]) + 1);
                        cv::rectangle(img, p0, p1, cv::Scalar(255));
                }
                std::cout << "index_ = " << index_ << std::endl;
                cv::imwrite("/Users/kumada/Data/level-set-method/fronts/front_" + std::to_string(index_) + ".jpg", img);
                ++index_;
        }

        static void draw_contour(const Front<ThreeDimension>& front, const SpaceSize<ThreeDimension>& size) {}
        

private:
        template<typename D>
        static void display_point(const IntPoint<D>& p)
        {
                std::cout << "{";
                for ( int i = 0; i < D::Dimension_; ++i ) {
                        std::cout << p[i];
                        if ( i == D::Dimension_ - 1 ) {
                                std::cout << "}, ";
                        } else {
                                std::cout << ", ";
                        }
                }
        }
};

} // namespace lsm

#endif // LevelSetMethod2_Debug_h
