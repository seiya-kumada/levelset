//-----------------------------------------------------------------------------------------------------------------------------
#ifndef POINT_H_
#define POINT_H_
//-----------------------------------------------------------------------------------------------------------------------------
#include <cmath>
#include "Utilities.h"
//-----------------------------------------------------------------------------------------------------------------------------

namespace kmd {

        //*********************************************************************************************************************
        class Point {
        //*********************************************************************************************************************
        private:
                double     x_;
                double     y_;

        public:
                // test ok
                //=============================================================================================================
                Point(double x, double y)
                //=============================================================================================================
                        : x_(x),
                          y_(y) {
                }

                //=============================================================================================================
                Point()
                //=============================================================================================================
                        : x_(0.0),
                          y_(0.0) {
                }

                // test ok
                //=============================================================================================================
                void set_x(double x) {
                //=============================================================================================================
                        x_ = x;
                }

                // test ok
                //=============================================================================================================
                void set_y(double y) {
                //=============================================================================================================
                        y_ = y;
                }

                // test ok
                //=============================================================================================================
                double x() const {
                //=============================================================================================================
                        return x_;
                }

                // test ok
                //=============================================================================================================
                double y() const {
                //=============================================================================================================
                        return y_;
                }

                // test ok
                //=============================================================================================================
                double& x() {
                //=============================================================================================================
                        return x_;
                }

                // test ok
                //=============================================================================================================
                double& y() {
                //=============================================================================================================
                        return y_;
                }


                // test ok
                //=============================================================================================================
                Point& operator+=(const Point& other) {
                //=============================================================================================================
                        x_ += other.x_;
                        y_ += other.y_;
                        return *this;
                }

                // test ok
                //=============================================================================================================
                Point& operator-=(const Point& other) {
                //=============================================================================================================
                        x_ -= other.x_;
                        y_ -= other.y_;
                        return *this;
                }

                //=============================================================================================================
                double square_norm() const {
                //=============================================================================================================
                        return Utilities::square(x_) + Utilities::square(y_);
                }

                Point* clone() const {
                        return new Point(x_, y_);
                }



        };

        struct LessPoint {
                bool operator()(const Point& a, const Point& b) const {
                        if ( a.x() == b.x() ) {
                                return a.y() < b.y();
                        } else {
                                return a.x() < b.x();
                        }
                }
        };

        // test ok
        //=====================================================================================================================
        inline const Point operator+(Point a, const Point& b) {
        //=====================================================================================================================
                return a += b;
        }

        // test ok
        //=====================================================================================================================
        inline const Point operator-(Point a, const Point& b) {
        //=====================================================================================================================
                return a -= b;
        }

        // test ok
        //=====================================================================================================================
        inline const Point operator*(double a, const Point& b) {
        //=====================================================================================================================
                return Point(a * b.x(), a * b.y());
        }

        // test ok
        //=====================================================================================================================
        inline const Point operator/(const Point& a, double b) {
        //=====================================================================================================================
                return Point(a.x() / b, a.y() / b);
        }

        // test ok
        //=====================================================================================================================
        inline double cross_product(const Point& a, const Point& b) {
        //=====================================================================================================================
                return a.x() * b.y() - a.y() * b.x();
        }

        // test ok
        //=====================================================================================================================
        inline double inner_product(const Point& a, const Point& b) {
        //=====================================================================================================================
                return a.x() * b.x() + a.y() * b.y();
        }



        // test ok
        //=====================================================================================================================
        inline bool operator==(const Point& a, const Point& b) {
        //=====================================================================================================================
                return (fabs(a.x() - b.x()) < std::numeric_limits<double>::epsilon())
                    && (fabs(a.y() - b.y()) < std::numeric_limits<double>::epsilon());
        }

        //=====================================================================================================================
        inline bool operator!=(const Point& a, const Point& b) {
        //=====================================================================================================================
                return !operator==(a, b);
        }



}
#endif /*POINT_H_*/
