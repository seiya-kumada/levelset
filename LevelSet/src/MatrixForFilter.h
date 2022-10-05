//-----------------------------------------------------------------------------------------------------------------------------
#ifndef MATRIXFORFILTER_H_
#define MATRIXFORFILTER_H_
//-----------------------------------------------------------------------------------------------------------------------------
#include <boost/noncopyable.hpp>
//-----------------------------------------------------------------------------------------------------------------------------
namespace kmd {



        //*********************************************************************************************************************
        class MatrixForFilter : private boost::noncopyable {
        //*********************************************************************************************************************
        public:
                static const double PrewittForX_[];
                static const double PrewittForY_[];

                static const double SobelForX_[];
                static const double SobelForY_[];

                static const double EdgeDetectForX_[];
                static const double EdgeDetectForY_[];

                static const double Laplacian_[];

                static const int PrewittForXSize_;
                static const int PrewittForYSize_;

                static const int SobelForXSize_;
                static const int SobelForYSize_;

                static const int EdgeDetectForXSize_;
                static const int EdgeDetectForYSize_;

                static const int LaplacianSize_;

        private:
                MatrixForFilter();
                ~MatrixForFilter();
        };
}
#endif /*MATRIXFORFILTER_H_*/
