//-----------------------------------------------------------------------------------------------------------------------------
#ifndef MATRIXFILTER_H_
#define MATRIXFILTER_H_
//-----------------------------------------------------------------------------------------------------------------------------
#include <boost/noncopyable.hpp>
#include <cstddef>
//-----------------------------------------------------------------------------------------------------------------------------
namespace kmd {

        /*!
         *      This class applies a filter to a given image.
         */
        //*********************************************************************************************************************
        class MatrixFilter : private boost::noncopyable {
        //*********************************************************************************************************************
                friend class MatrixFilterTester;

        public:
                MatrixFilter();
                ~MatrixFilter();

                /*!
                 *      This method sets the filter which is applied to an image.
                 *      This may throw std::runtime_error.
                 *
                 *      @param[in]      matrix
                 *      @param[in]      matrix_size
                 */
                void set_matrix(const double* matrix, std::size_t matrix_size);

                /*!
                 *      This method applies a filter to a given image.
                 *      This may throw std::runtime_error.
                 *
                 *      @param[out]     dst_buf output image
                 *      @param[in]      src_buf input image
                 *      @param[in]      w       image width
                 *      @param[in]      h       image height
                 */
                void run(double* dst_buf, const double* src_buf, int w, int h) const;

        private:
                const double*   matrix_;
                std::size_t     matrix_size_;
                int             range_;

                double apply_matrix_to_local_area(const double* src, int w, int h, int center_x, int center_y) const;
        };
}
#endif /*MATRIXFILTER_H_*/

