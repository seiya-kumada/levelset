//-----------------------------------------------------------------------------------------------------------------------------
#include "MatrixForFilter.h"
#include <cstddef>
//-----------------------------------------------------------------------------------------------------------------------------
using namespace kmd;

const double MatrixForFilter::PrewittForX_[] = {
        1, 0, -1,
        1, 0, -1,
        1, 0, -1,
};

const double MatrixForFilter::PrewittForY_[] = {
        1,  1,  1,
        0,  0,  0,
       -1, -1, -1,
};

const double MatrixForFilter::SobelForX_[] = {
        1, 0, -1,
        2, 0, -2,
        1, 0, -1,
};

const double MatrixForFilter::SobelForY_[] = {
        1,  2,  1,
        0,  0,  0,
       -1, -2, -1,
};

const double MatrixForFilter::EdgeDetectForX_[] = {
        -1, 1, -1,
        -1, 1, -1,
        -1, 1, -1,
};

const double MatrixForFilter::EdgeDetectForY_[] = {
        -1, -1, -1,
         1,  1,  1,
        -1, -1, -1,
};

const double MatrixForFilter::Laplacian_[] = {
        0,  1, 0,
        1, -4, 1,
        0,  1, 0,
};

const int MatrixForFilter::PrewittForXSize_    = 3;
const int MatrixForFilter::PrewittForYSize_    = 3;
const int MatrixForFilter::SobelForXSize_      = 3;
const int MatrixForFilter::SobelForYSize_      = 3;
const int MatrixForFilter::EdgeDetectForXSize_ = 3;
const int MatrixForFilter::EdgeDetectForYSize_ = 3;
const int MatrixForFilter::LaplacianSize_      = 3;



