//------------------------------------------------------------------------------------------------------------------------------------
#ifndef FRONT_H_
#define FRONT_H_
//------------------------------------------------------------------------------------------------------------------------------------
#include <cmath>
//------------------------------------------------------------------------------------------------------------------------------------

namespace kmd {

        namespace test {
                class Front {
                public:
                        void set_radius(double radius) {
                                radius_ = radius;
                        }

                        double x(double radian) const {
                                return center_x_ + radius_ * cos(radian);
                        }

                        double y(double radian) const {
                                return center_y_ + radius_ * sin(radian);
                        }

                        void set_center(double center_x, double center_y) {
                                center_x_ = center_x;
                                center_y_ = center_y;
                        }

                private:
                        double  radius_;
                        double  center_x_;
                        double  center_y_;
                };
        }
}
#endif /*FRONT_H_*/
