//------------------------------------------------------------------------------------------------------------------------------------
#ifndef ELLIPSE_H_
#define ELLIPSE_H_
//------------------------------------------------------------------------------------------------------------------------------------
namespace kmd {
        class Ellipse {
        public:
                void set_radius(int major_radius, int minor_radius) {
                        major_radius_ = major_radius;
                        minor_radius_ = minor_radius;
                }

                double x(double radian) const {
                        return center_x_ + major_radius_ * cos(radian);
                }

                double y(double radian) const {
                        return center_y_ + minor_radius_ * sin(radian);
                }

                void set_center(double center_x, double center_y) {
                        center_x_ = center_x;
                        center_y_ = center_y;
                }

        private:
                double  major_radius_;
                double  minor_radius_;
                double  center_x_;
                double  center_y_;
        };
}
#endif /*ELLIPSE_H_*/
