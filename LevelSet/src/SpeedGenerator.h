//-----------------------------------------------------------------------------------------------------------------------------
#ifndef SPEEDGENERATOR_H_
#define SPEEDGENERATOR_H_
//-----------------------------------------------------------------------------------------------------------------------------
#include <vector>
#include <boost/noncopyable.hpp>
#include "Typedef.h"
#//-----------------------------------------------------------------------------------------------------------------------------

namespace kmd {

        struct Parameters;

        //*********************************************************************************************************************
        class SpeedGenerator : private boost::noncopyable {
        //*********************************************************************************************************************
                friend class SpeedGeneratorTester;

        public:
                /*!
                 *      @param[in]      parameters
                 */
                explicit SpeedGenerator(const Parameters* parameters);

                ~SpeedGenerator();

                /*!
                 *      @param[in]      src     source image buffer
                 */
                void calculate_image_dependent_factors(const std::vector<double>& src);
                
                /*!
                 *      @return pointer to a head of "image_dependent_factors_"
                 */
                const std::vector<double>& get_image_dependent_factors() const {
                        return image_dependent_factors_;        
                }
                
                /*!
                 *      @param[in]      buf     address of a target pixel
                 *      @param[in]      w       width of image
                 *      @return                 curvature-dependent speed
                 */
                double calculate_curvature_dependent_speed(const double* buf, int w) const;
                
        private:
                const Parameters*       parameters_;
                std::vector<double>     image_dependent_factors_;
        
                double calculate_curvature(const double* buf, int w) const;
                void create_gaussian_matrix(std::vector<double>& matrix, int matrix_size, double sigma) const;
                void apply_gaussian_filter(std::vector<double>& medium, const std::vector<double>& src, double sigma) const;
                void apply_differential_filter(std::vector<double>& dst, const std::vector<double>& medium) const;
                double calculate_curvature_dependent_speed(double curvature) const;
                void modify_image(std::vector<double>& dst, const std::vector<double>& src, double sigma) const;
        };
}
#endif /*SPEEDGENERATOR_H_*/
