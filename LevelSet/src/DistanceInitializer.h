//-----------------------------------------------------------------------------------------------------------------------------
#ifndef DISTANCEINITIALIZER_H_
#define DISTANCEINITIALIZER_H_
//-----------------------------------------------------------------------------------------------------------------------------
#include <boost/noncopyable.hpp>
#include "Typedef.h"
#include <list>
//-----------------------------------------------------------------------------------------------------------------------------

namespace kmd {
        class FrontDirection;
        
        /*!
         *      This class initializes a signed distance function.
         */
        //*********************************************************************************************************************
        class DistanceInitializer : private boost::noncopyable {
        //*********************************************************************************************************************
                friend class DistanceInitializerTester;

        public:
                /*!
                 *      @param[in]      distance_finder
                 */
                explicit DistanceInitializer(
                        const ShortestDistanceFinder* distance_finder
                );
                
                ~DistanceInitializer();
                                
                /*!
                 *      This method may throw std::runtime_error.
                 *
                 *      @param[out]     buffer
                 *      @param[in]      size    size of buffer
                 *      @param[in]      fronts
                 *      @param[in]      tube    a set of points which are inside the tube domain
                 */
                void run(
                        double*         buffer,
                        std::size_t     size,
                        const Fronts&   fronts,
                        const Tube&     tube
                ) const;
                
                /*!
                 *      @param[out]     new_buffer
                 *      @param[in]      old_buffer      
                 *      @param[in]      size            size of buffer
                 *      @param[in]      new_tube        a set of points which are inside the tube domain
                 *      @param[in]      direction       direction for the front to develop
                 */
                void run(
                        double*                 new_buffer,
                        const double*           old_buffer,
                        std::size_t             size,
                        const Tube&             new_tube,
                        const FrontDirection&   direction
                ) const;

        private:
                typedef void (DistanceInitializer::*AssignValue)(double& dst, double src, double dist) const;
                
                const ShortestDistanceFinder*   distance_finder_;
                AssignValue                     assign_value_[2];
                
                
                void assign_value_inner(double& dst, double src, double dist) const;
                void assign_value_outer(double& dst, double src, double dist) const;
        };
}
#endif /*DISTANCEINITIALIZER_H_*/
