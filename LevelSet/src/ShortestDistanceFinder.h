//-----------------------------------------------------------------------------------------------------------------------------
#ifndef SHORTESTDISTANCEFINDER_H_
#define SHORTESTDISTANCEFINDER_H_
//-----------------------------------------------------------------------------------------------------------------------------
#include <vector>
#include <boost/noncopyable.hpp>
#include "Typedef.h"
#include <list>
//-----------------------------------------------------------------------------------------------------------------------------
namespace kmd {

        struct Parameters;

        /*!
         *      This class finds out a point which lies on a front and is the closest to a given location.
         */
        //*********************************************************************************************************************
        class ShortestDistanceFinder : private boost::noncopyable {
        //*********************************************************************************************************************
                friend class ShortestDistanceFinderTester;

        public:
                /*!
                 *      @param[in]      parameters
                 */
                explicit ShortestDistanceFinder(const Parameters* parameters);

                ~ShortestDistanceFinder();

                // test ok
                /*!
                 *      @param[in]      tube_width
                 */
                void set_tube_width(int tube_width) {
                        tube_width_ = tube_width;
                }

                // test ok
                int get_tube_width() const {
                        return tube_width_;
                }

                const std::vector<double>& get_shortest_distances() const {
                        return shortest_distances_;
                }

                /*!
                 *      @param[in]      fronts
                 */
                void run(const Fronts& fronts);

                const Parameters* get_parameters() const {
                        return parameters_;
                }

        private:
                const Parameters*       parameters_;
                std::vector<double>     shortest_distances_;
                int                     tube_width_;

                double calculate_square_distance(const Point& a, const Point& b, const Point& p) const;
                void scan_area(double* infos, int w, const Point* a, const Point* b) const;
                void walk_on_front(const Point* const front[], std::size_t size, int w);
        };
}
#endif /*SHORTESTDISTANCEFINDER_H_*/

