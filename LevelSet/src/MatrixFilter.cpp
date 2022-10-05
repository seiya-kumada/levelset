//-----------------------------------------------------------------------------------------------------------------------------
#include "Error.h"
#include "MatrixFilter.h"
#include <stdexcept>
//-----------------------------------------------------------------------------------------------------------------------------
using namespace std;
using namespace kmd;

//=============================================================================================================================
MatrixFilter::MatrixFilter()
//=============================================================================================================================
        : matrix_(0),
          matrix_size_(0),
          range_(0) {
 }

//=============================================================================================================================
MatrixFilter::~MatrixFilter() {}
//=============================================================================================================================

// test ok
//=============================================================================================================================
void MatrixFilter::set_matrix(const double* matrix, std::size_t matrix_size) {
//=============================================================================================================================
        if ( matrix_size % 2 == 0 ) {
                throw runtime_error(Error::message("The size of matrix is invalid.", __LINE__, __FILE__));
        }
        matrix_ = matrix;
        matrix_size_ = matrix_size;
        range_ = static_cast<int>(matrix_size / 2);
}

// test ok
//=============================================================================================================================
double MatrixFilter::apply_matrix_to_local_area(const double* src, int w, int h, int center_x, int center_y) const {
//=============================================================================================================================
        double s = 0.0;
        int k = 0;
        for ( int j = -range_, wj = w * j; j <= range_; ++j, wj += w ) {
                const double* tmp = &src[wj];
                for ( int i = -range_; i <= range_; ++i ) {
                        int x_pos = i + center_x;
                        int y_pos = j + center_y;
                        if ( 0 <= x_pos && x_pos < w && 0 <= y_pos && y_pos < h ) {
                                // if (i, j) is inside the image
                                s += tmp[i] * matrix_[k];
                        }
                        ++k;
                }
        }
        return s;
}

// test ok
//=============================================================================================================================
void MatrixFilter::run(double* dst_buf, const double* src_buf, int w, int h) const {
//=============================================================================================================================
        for ( int j = 0, wj = w * j; j < h; ++j, wj += w ) {
                double* dst = &dst_buf[wj];
                const double* src = &src_buf[wj];
                for ( int i = 0; i < w; ++i ) {
                        dst[i] = apply_matrix_to_local_area(&src[i], w, h, i, j);
                }
        }
}



#ifdef UNIT_TEST_MatrixFilter
#define BOOST_TEST_DYN_LINK
//-----------------------------------------------------------------------------------------------------------------------------
#include <kmd/without_warning/unit_test.h>
#include <kmd/without_warning/iostream.h>
#include <cmath>
#include "MatrixForFilter.h"
//-----------------------------------------------------------------------------------------------------------------------------
using namespace std;

namespace kmd {

        //*********************************************************************************************************************
        class MatrixFilterTester {
        //*********************************************************************************************************************
        public:


                //=============================================================================================================
                static const double* test_matrix(MatrixFilter& filter) {
                //=============================================================================================================
                        return filter.matrix_;
                }

                //=============================================================================================================
                static std::size_t test_matrix_size(MatrixFilter& filter) {
                //=============================================================================================================
                        return filter.matrix_size_;
                }

                //=============================================================================================================
                static int test_range(MatrixFilter& filter) {
                //=============================================================================================================
                        return filter.range_;
                }


        };
}

//=============================================================================================================================
BOOST_AUTO_TEST_CASE(TEST_MatrixFilter) {
//=============================================================================================================================

        cout << "MatrixFilter\n";

        // exception behaviour
        {

                MatrixFilter filter;
                bool flag = false;
                try {
                        filter.set_matrix(0, 4);
                } catch (const runtime_error&) {
                        flag = true;
                }
                BOOST_CHECK(flag);
        }

        // matrix_size_/matrix_/set_matrix/run
        {
                MatrixFilter filter;
                int w = 3;
                vector<double> matrix(w * w, 5.0 / 45.0);
                filter.set_matrix(&matrix[0], w);
                BOOST_CHECK(MatrixFilterTester::test_matrix_size(filter) == static_cast<std::size_t>(w));
                for ( std::size_t i = 0, n = matrix.size(); i < n; ++i ) {
                        BOOST_CHECK(matrix[i] == 5 / 45.0);
                }

                int image_w = 10;
                int image_h = 10;
                vector<double> src_buf(image_w * image_h, 1.0);
                vector<double> dst_buf(image_w * image_h, 0.0);

                //
                filter.run(&dst_buf[0], &src_buf[0], image_w, image_h);

                double epsilon = 0.000000001;
                for ( int j = 0; j < image_h; ++j ) {
                        for ( int i = 0; i < image_w; ++i ) {

                                if (
                                           (i == 0 && j == 0)
                                        || (i == 0 && j == image_h - 1)
                                        || (i == image_w - 1 && j == 0)
                                        || (i == image_w - 1 && j == image_h - 1)
                                ) {
                                        BOOST_CHECK(fabs(dst_buf[i + image_w * j] - 4.0 / 9.0) < epsilon);
                                } else if ( i == 0 || i == image_w - 1 || j == 0 || j == image_h - 1 ) {
                                        BOOST_CHECK(fabs(dst_buf[i + image_w * j] - 2.0 / 3.0) < epsilon);
                                } else {
                                        BOOST_CHECK(fabs(dst_buf[i + image_w * j] - 1.0) < epsilon);
                                }
                        }
                }
        }

}
#endif // UNIT_TEST_MatrixFilter
