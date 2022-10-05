//
//  Differential.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/11/24.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_Differential_h
#define LevelSetMethod2_Differential_h

#include "Concepts.h"
#include "DimensionTypes.h"
#include "NeighboringPoints.h"
#include "Indexer.h"
#include <boost/concept_check.hpp>

namespace lsm
{

        namespace detail
        {
                /// These constants are used to extend the 2D soble filter to the 3D one.
                /**
                 *      @see "Extension to other dimensions" 
                 *      in http://en.wikipedia.org/wiki/Sobel_operator
                 */
                constexpr double h0d_[] = { 1,  2,  1};
                constexpr double h1d_[] = {-1,  0,  1};
                constexpr double h2d_[] = { 1, -2,  1};
                constexpr double h3d_[] = { 1,  0, -1};
                
                /// This class is used to implement the Differential class.
                /**
                 *      @tparam D the dimension type
                 *      @tparam T the value type (std::uint8_t or double)
                 *
                 *      @see http://en.wikipedia.org/wiki/Sobel_operator
                 */
                template<typename D, typename T>
                class DifferentialTool
                {
                        BOOST_CLASS_REQUIRE(D, lsm, DimensionConcept);
                        
                protected:
                        // test ok
                        static constexpr int h0d_total_ = 1 + 2 + 1;
                
                        static constexpr int index1(int x)
                        {
                                return x + 1;
                        };
                        
                        static constexpr double h(int x)
                        {
                                return h0d_[index1(x)];
                        }

                        static constexpr double h1d(int x)
                        {
                                return h1d_[index1(x)];
                        }

                        static constexpr double h2d(int x)
                        {
                                return h2d_[index1(x)];
                        }

                        static constexpr double h3d(int x)
                        {
                                return h3d_[index1(x)];
                        }

                        DifferentialTool(
                                const Indexer<D>& indexer,
                                const std::vector<T>& buffer
                        )
                                : indexer_(indexer)
                                , buffer_(buffer) {}
                        
                        ~DifferentialTool() = default;

                        const Indexer<D>& indexer_;
                        const std::vector<T>& buffer_;
                        std::array<T, power(3, D::Dimension_)> values_;

                        const T& value(const IntPoint<D>& p) const
                        {
                                return buffer_[indexer_(p)];
                        }
                        
                private:
                        DifferentialTool(const DifferentialTool&) = delete;
                        DifferentialTool& operator=(const DifferentialTool&) = delete;
                        
                        DifferentialTool(DifferentialTool&&) = delete;
                        DifferentialTool& operator=(DifferentialTool&&) = delete;
                }; // class DifferentialTool
        } // namespace detail

        /// calculates differential
        /**
         *      @tparam D the dimension type
         *      @tparam T the value type (std::uint8_t or double)
         */
        template<typename D, typename T>
        class Differential;

        /// explicit specialization with TwoDimension
        /**
         *      @tparam T the value type (std::uint8_t or double)
         */
        template<typename T>
        class Differential<TwoDimension, T> : private detail::DifferentialTool<TwoDimension, T>
        {
                typedef TwoDimension D;
                typedef detail::DifferentialTool<D, T> Base;

                using Base::value;
                using Base::h0d_total_;
                
                friend class DifferentialTester;
                
        public:
                Differential(
                        const Indexer<D>&     indexer,
                        const std::vector<T>& buffer
                )
                        : Base{indexer, buffer} {}
                
                ~Differential() = default;
         
                // test ok
                const double sobel_x() const
                {
                        return    vx(-1, -1) + vx(0, -1) + vx(1, -1)
                                + vx(-1,  0) + vx(0,  0) + vx(1,  0)
                                + vx(-1,  1) + vx(0,  1) + vx(1,  1);
                }

                // test ok
                const double sobel_y() const
                {
                        return    vy(-1, -1) + vy(0, -1) + vy(1, -1)
                                + vy(-1,  0) + vy(0,  0) + vy(1,  0)
                                + vy(-1,  1) + vy(0,  1) + vy(1,  1);
                }

                // test ok
                const double fx() const
                {
                        return sobel_x() / (2 * h0d_total_);
                }
                
                // test ok
                const double fy() const
                {
                        return sobel_y() / (2 * h0d_total_);
                }

                // test ok
                const double fxx() const
                {
                        return   (vxx(-1, -1) + vxx(0, -1) + vxx(1, -1)
                                + vxx(-1,  0) + vxx(0,  0) + vxx(1,  0)
                                + vxx(-1,  1) + vxx(0,  1) + vxx(1,  1)) / h0d_total_;
                }

                // test ok
                const double fyy() const
                {
                        return   (vyy(-1, -1) + vyy(0, -1) + vyy(1, -1)
                                + vyy(-1,  0) + vyy(0,  0) + vyy(1,  0)
                                + vyy(-1,  1) + vyy(0,  1) + vyy(1,  1)) / h0d_total_ ;
                }

                // test ok
                const double fxy() const
                {
                        return   (vxy(-1, -1) + vxy(0, -1) + vxy(1, -1)
                                + vxy(-1,  0) + vxy(0,  0) + vxy(1,  0)
                                + vxy(-1,  1) + vxy(0,  1) + vxy(1,  1)) / 4.0;
                }

                /// set points around a given point
                void set_point(const IntPoint<D>& p)
                {
                        v(-1, -1) = value(p + NeighboringPoints2d(-1, -1));
                        v( 0, -1) = value(p + NeighboringPoints2d( 0, -1));
                        v( 1, -1) = value(p + NeighboringPoints2d( 1, -1));

                        v(-1,  0) = value(p + NeighboringPoints2d(-1,  0));
                        v( 0,  0) = value(p + NeighboringPoints2d( 0,  0));
                        v( 1,  0) = value(p + NeighboringPoints2d( 1,  0));
                        
                        v(-1,  1) = value(p + NeighboringPoints2d(-1,  1));
                        v( 0,  1) = value(p + NeighboringPoints2d( 0,  1));
                        v( 1,  1) = value(p + NeighboringPoints2d( 1,  1));
                }

        private:
                Differential(const Differential&) = delete;
                Differential& operator=(const Differential&) = delete;
                
                Differential(Differential&&) = delete;
                Differential& operator=(Differential&&) = delete;
               
                static constexpr int index2(int i, int j)
                {
                        return (i + 1) + 3 * (j + 1);
                }

                static constexpr double h1dx(int x, int y)
                {
                        return Base::h1d(x) * Base::h(y);
                }

                static constexpr double h1dy(int x, int y)
                {
                        return Base::h(x) * Base::h1d(y);
                }

                // test ok
                static constexpr double h2dx(int x, int y)
                {
                        return Base::h2d(x) * Base::h(y);
                }

                static constexpr double h2dy(int x, int y)
                {
                        return Base::h(x) * Base::h2d(y);
                } 

                static constexpr double h3dxy(int x, int y)
                {
                        return Base::h3d(x) * Base::h3d(y);
                }

                T& v(int x, int y)
                {
                        return Base::values_[index2(x, y)];
                }

                const T& v(int x, int y) const
                {
                        return Base::values_[index2(x, y)];
                }

                const double vx(int x, int y) const
                {
                        return v(x, y) * h1dx(x, y);
                }

                const double vy(int x, int y) const
                {
                        return v(x, y) * h1dy(x, y);
                }

                const double vxx(int x, int y) const
                {
                        return v(x, y) * h2dx(x, y);
                }

                const double vyy(int x, int y) const
                {
                        return v(x, y) * h2dy(x, y);
                }

                const double vxy(int x, int y) const
                {
                        return v(x, y) * h3dxy(x, y);
                }
        }; // class Differential<TwoDimension, T>

        /// explicit specialization with ThreeDimension
        /**
         *      @tparam T the value type (std::uint8_t or double)
         */
        template<typename T>
        class Differential<ThreeDimension, T> : private detail::DifferentialTool<ThreeDimension, T>
        {
                typedef ThreeDimension D;
                typedef detail::DifferentialTool<D, T> Base;
                
                using Base::value;
                using Base::h0d_total_;
                
                friend class DifferentialTester;

        public:
                Differential(
                        const Indexer<D>&     indexer,
                        const std::vector<T>& buffer
                )
                        : Base{indexer, buffer} {}
         
                ~Differential() = default;
         
                // test ok
                const double sobel_x() const
                {
                        return    vx(-1, -1, -1)
                                + vx( 0, -1, -1)
                                + vx( 1, -1, -1)
                                + vx(-1,  0, -1)
                                + vx( 0,  0, -1)
                                + vx( 1,  0, -1)
                                + vx(-1,  1, -1)
                                + vx( 0,  1, -1)
                                + vx( 1,  1, -1)
                        
                                + vx(-1, -1,  0)
                                + vx( 0, -1,  0)
                                + vx( 1, -1,  0)
                                + vx(-1,  0,  0)
                                + vx( 0,  0,  0)
                                + vx( 1,  0,  0)
                                + vx(-1,  1,  0)
                                + vx( 0,  1,  0)
                                + vx( 1,  1,  0)
                        
                                + vx(-1, -1,  1)
                                + vx( 0, -1,  1)
                                + vx( 1, -1,  1)
                                + vx(-1,  0,  1)
                                + vx( 0,  0,  1)
                                + vx( 1,  0,  1)
                                + vx(-1,  1,  1)
                                + vx( 0,  1,  1)
                                + vx( 1,  1,  1);
                }

                // test ok
                const double sobel_y() const
                {
                        return    vy(-1, -1, -1)
                                + vy( 0, -1, -1)
                                + vy( 1, -1, -1)
                                + vy(-1,  0, -1)
                                + vy( 0,  0, -1)
                                + vy( 1,  0, -1)
                                + vy(-1,  1, -1)
                                + vy( 0,  1, -1)
                                + vy( 1,  1, -1)

                                + vy(-1, -1,  0)
                                + vy( 0, -1,  0)
                                + vy( 1, -1,  0)
                                + vy(-1,  0,  0)
                                + vy( 0,  0,  0)
                                + vy( 1,  0,  0)
                                + vy(-1,  1,  0)
                                + vy( 0,  1,  0)
                                + vy( 1,  1,  0)
                        
                                + vy(-1, -1,  1)
                                + vy( 0, -1,  1)
                                + vy( 1, -1,  1)
                                + vy(-1,  0,  1)
                                + vy( 0,  0,  1)
                                + vy( 1,  0,  1)
                                + vy(-1,  1,  1)
                                + vy( 0,  1,  1)
                                + vy( 1,  1,  1);
                }

                // test ok
                const double sobel_z() const
                {
                        return    vz(-1, -1, -1)
                                + vz( 0, -1, -1)
                                + vz( 1, -1, -1)
                                + vz(-1,  0, -1)
                                + vz( 0,  0, -1)
                                + vz( 1,  0, -1)
                                + vz(-1,  1, -1)
                                + vz( 0,  1, -1)
                                + vz( 1,  1, -1)
                        
                                + vz(-1, -1,  0)
                                + vz( 0, -1,  0)
                                + vz( 1, -1,  0)
                                + vz(-1,  0,  0)
                                + vz( 0,  0,  0)
                                + vz( 1,  0,  0)
                                + vz(-1,  1,  0)
                                + vz( 0,  1,  0)
                                + vz( 1,  1,  0)
                        
                                + vz(-1, -1,  1)
                                + vz( 0, -1,  1)
                                + vz( 1, -1,  1)
                                + vz(-1,  0,  1)
                                + vz( 0,  0,  1)
                                + vz( 1,  0,  1)
                                + vz(-1,  1,  1)
                                + vz( 0,  1,  1)
                                + vz( 1,  1,  1);
                }

                // test ok
                const double fx() const
                {
                        return sobel_x() / 32.0;
                }

                // test ok
                const double fy() const
                {
                        return sobel_y() / 32.0;
                }

                // test ok
                const double fz() const
                {
                        return sobel_z() / 32.0;
                }
                
                // test ok
                const double fxx() const
                {
                        return   (vxx(-1, -1, -1)
                                + vxx( 0, -1, -1)
                                + vxx( 1, -1, -1)
                                + vxx(-1,  0, -1)
                                + vxx( 0,  0, -1)
                                + vxx( 1,  0, -1)
                                + vxx(-1,  1, -1)
                                + vxx( 0,  1, -1)
                                + vxx( 1,  1, -1)
                        
                                + vxx(-1, -1,  0)
                                + vxx( 0, -1,  0)
                                + vxx( 1, -1,  0)
                                + vxx(-1,  0,  0)
                                + vxx( 0,  0,  0)
                                + vxx( 1,  0,  0)
                                + vxx(-1,  1,  0)
                                + vxx( 0,  1,  0)
                                + vxx( 1,  1,  0)
                        
                                + vxx(-1, -1,  1)
                                + vxx( 0, -1,  1)
                                + vxx( 1, -1,  1)
                                + vxx(-1,  0,  1)
                                + vxx( 0,  0,  1)
                                + vxx( 1,  0,  1)
                                + vxx(-1,  1,  1)
                                + vxx( 0,  1,  1)
                                + vxx( 1,  1,  1)) / 16.0;
                }

                // test ok
                const double fyy() const
                {
                        return   (vyy(-1, -1, -1)
                                + vyy( 0, -1, -1)
                                + vyy( 1, -1, -1)
                                + vyy(-1,  0, -1)
                                + vyy( 0,  0, -1)
                                + vyy( 1,  0, -1)
                                + vyy(-1,  1, -1)
                                + vyy( 0,  1, -1)
                                + vyy( 1,  1, -1)
                        
                                + vyy(-1, -1,  0)
                                + vyy( 0, -1,  0)
                                + vyy( 1, -1,  0)
                                + vyy(-1,  0,  0)
                                + vyy( 0,  0,  0)
                                + vyy( 1,  0,  0)
                                + vyy(-1,  1,  0)
                                + vyy( 0,  1,  0)
                                + vyy( 1,  1,  0)
                        
                                + vyy(-1, -1,  1)
                                + vyy( 0, -1,  1)
                                + vyy( 1, -1,  1)
                                + vyy(-1,  0,  1)
                                + vyy( 0,  0,  1)
                                + vyy( 1,  0,  1)
                                + vyy(-1,  1,  1)
                                + vyy( 0,  1,  1)
                                + vyy( 1,  1,  1)) / 16.0;
                }

                // test ok
                const double fzz() const
                {
                        return   (vzz(-1, -1, -1)
                                + vzz( 0, -1, -1)
                                + vzz( 1, -1, -1)
                                + vzz(-1,  0, -1)
                                + vzz( 0,  0, -1)
                                + vzz( 1,  0, -1)
                                + vzz(-1,  1, -1)
                                + vzz( 0,  1, -1)
                                + vzz( 1,  1, -1)
                        
                                + vzz(-1, -1,  0)
                                + vzz( 0, -1,  0)
                                + vzz( 1, -1,  0)
                                + vzz(-1,  0,  0)
                                + vzz( 0,  0,  0)
                                + vzz( 1,  0,  0)
                                + vzz(-1,  1,  0)
                                + vzz( 0,  1,  0)
                                + vzz( 1,  1,  0)
                        
                                + vzz(-1, -1,  1)
                                + vzz( 0, -1,  1)
                                + vzz( 1, -1,  1)
                                + vzz(-1,  0,  1)
                                + vzz( 0,  0,  1)
                                + vzz( 1,  0,  1)
                                + vzz(-1,  1,  1)
                                + vzz( 0,  1,  1)
                                + vzz( 1,  1,  1)) / 16.0;
                }

                // test ok
                const double fxy() const
                {
                        return   (vxy(-1, -1, -1)
                                + vxy( 0, -1, -1)
                                + vxy( 1, -1, -1)
                                + vxy(-1,  0, -1)
                                + vxy( 0,  0, -1)
                                + vxy( 1,  0, -1)
                                + vxy(-1,  1, -1)
                                + vxy( 0,  1, -1)
                                + vxy( 1,  1, -1)
                        
                                + vxy(-1, -1,  0)
                                + vxy( 0, -1,  0)
                                + vxy( 1, -1,  0)
                                + vxy(-1,  0,  0)
                                + vxy( 0,  0,  0)
                                + vxy( 1,  0,  0)
                                + vxy(-1,  1,  0)
                                + vxy( 0,  1,  0)
                                + vxy( 1,  1,  0)
                        
                                + vxy(-1, -1,  1)
                                + vxy( 0, -1,  1)
                                + vxy( 1, -1,  1)
                                + vxy(-1,  0,  1)
                                + vxy( 0,  0,  1)
                                + vxy( 1,  0,  1)
                                + vxy(-1,  1,  1)
                                + vxy( 0,  1,  1)
                                + vxy( 1,  1,  1)) / (4.0 * h0d_total_);
                }

                // test ok
                const double fxz() const
                {
                        return   (vxz(-1, -1, -1)
                                + vxz( 0, -1, -1)
                                + vxz( 1, -1, -1)
                                + vxz(-1,  0, -1)
                                + vxz( 0,  0, -1)
                                + vxz( 1,  0, -1)
                                + vxz(-1,  1, -1)
                                + vxz( 0,  1, -1)
                                + vxz( 1,  1, -1)
                        
                                + vxz(-1, -1,  0)
                                + vxz( 0, -1,  0)
                                + vxz( 1, -1,  0)
                                + vxz(-1,  0,  0)
                                + vxz( 0,  0,  0)
                                + vxz( 1,  0,  0)
                                + vxz(-1,  1,  0)
                                + vxz( 0,  1,  0)
                                + vxz( 1,  1,  0)
                                  
                                + vxz(-1, -1,  1)
                                + vxz( 0, -1,  1)
                                + vxz( 1, -1,  1)
                                + vxz(-1,  0,  1)
                                + vxz( 0,  0,  1)
                                + vxz( 1,  0,  1)
                                + vxz(-1,  1,  1)
                                + vxz( 0,  1,  1)
                                + vxz( 1,  1,  1)) / (4.0 * h0d_total_);
                }

                // test ok
                const double fyz() const
                {
                        return   (vyz(-1, -1, -1)
                                + vyz( 0, -1, -1)
                                + vyz( 1, -1, -1)
                                + vyz(-1,  0, -1)
                                + vyz( 0,  0, -1)
                                + vyz( 1,  0, -1)
                                + vyz(-1,  1, -1)
                                + vyz( 0,  1, -1)
                                + vyz( 1,  1, -1)
                        
                                + vyz(-1, -1,  0)
                                + vyz( 0, -1,  0)
                                + vyz( 1, -1,  0)
                                + vyz(-1,  0,  0)
                                + vyz( 0,  0,  0)
                                + vyz( 1,  0,  0)
                                + vyz(-1,  1,  0)
                                + vyz( 0,  1,  0)
                                + vyz( 1,  1,  0)
                        
                                + vyz(-1, -1,  1)
                                + vyz( 0, -1,  1)
                                + vyz( 1, -1,  1)
                                + vyz(-1,  0,  1)
                                + vyz( 0,  0,  1)
                                + vyz( 1,  0,  1)
                                + vyz(-1,  1,  1)
                                + vyz( 0,  1,  1)
                                + vyz( 1,  1,  1)) / (4.0 * h0d_total_);
                }

                /// set points around a given point
                void set_point(const IntPoint<D>& p)
                {
                        v(-1, -1, -1) = value(p + NeighboringPoints3d(-1, -1, -1));
                        v( 0, -1, -1) = value(p + NeighboringPoints3d( 0, -1, -1));
                        v( 1, -1, -1) = value(p + NeighboringPoints3d( 1, -1, -1));
                        v(-1,  0, -1) = value(p + NeighboringPoints3d(-1,  0, -1));
                        v( 0,  0, -1) = value(p + NeighboringPoints3d( 0,  0, -1));
                        v( 1,  0, -1) = value(p + NeighboringPoints3d( 1,  0, -1));
                        v(-1,  1, -1) = value(p + NeighboringPoints3d(-1,  1, -1));
                        v( 0,  1, -1) = value(p + NeighboringPoints3d( 0,  1, -1));
                        v( 1,  1, -1) = value(p + NeighboringPoints3d( 1,  1, -1));

                        v(-1, -1,  0) = value(p + NeighboringPoints3d(-1, -1,  0));
                        v( 0, -1,  0) = value(p + NeighboringPoints3d( 0, -1,  0));
                        v( 1, -1,  0) = value(p + NeighboringPoints3d( 1, -1,  0));
                        v(-1,  0,  0) = value(p + NeighboringPoints3d(-1,  0,  0));
                        v( 0,  0,  0) = value(p + NeighboringPoints3d( 0,  0,  0));
                        v( 1,  0,  0) = value(p + NeighboringPoints3d( 1,  0,  0));
                        v(-1,  1,  0) = value(p + NeighboringPoints3d(-1,  1,  0));
                        v( 0,  1,  0) = value(p + NeighboringPoints3d( 0,  1,  0));
                        v( 1,  1,  0) = value(p + NeighboringPoints3d( 1,  1,  0));

                        v(-1, -1,  1) = value(p + NeighboringPoints3d(-1, -1,  1));
                        v( 0, -1,  1) = value(p + NeighboringPoints3d( 0, -1,  1));
                        v( 1, -1,  1) = value(p + NeighboringPoints3d( 1, -1,  1));
                        v(-1,  0,  1) = value(p + NeighboringPoints3d(-1,  0,  1));
                        v( 0,  0,  1) = value(p + NeighboringPoints3d( 0,  0,  1));
                        v( 1,  0,  1) = value(p + NeighboringPoints3d( 1,  0,  1));
                        v(-1,  1,  1) = value(p + NeighboringPoints3d(-1,  1,  1));
                        v( 0,  1,  1) = value(p + NeighboringPoints3d( 0,  1,  1));
                        v( 1,  1,  1) = value(p + NeighboringPoints3d( 1,  1,  1));
                }

        private:
                Differential(const Differential&) = delete;
                Differential& operator=(const Differential&) = delete;
                
                Differential(Differential&&) = delete;
                Differential& operator=(Differential&&) = delete;
               
                static constexpr int index3(int i, int j, int k)
                {
                        return (i + 1) + 3 * (j + 1) + 9 * (k + 1);
                }

                static constexpr double h1dx(int x, int y, int z)
                {
                        return Base::h1d(x) * Base::h(y) * Base::h(z);
                }

                static constexpr double h1dy(int x, int y, int z)
                {
                        return Base::h(x) * Base::h1d(y) * Base::h(z);
                }

                // test ok
                static constexpr double h2dx(int x, int y, int z)
                {
                        return Base::h2d(x) * Base::h(y) * Base::h(z);
                }

                static constexpr double h2dy(int x, int y, int z)
                {
                        return Base::h(x) * Base::h2d(y) * Base::h(z);
                }

                static constexpr double h3dxy(int x, int y, int z)
                {
                        return Base::h3d(x) * Base::h3d(y) * Base::h(z);
                }

                static constexpr double hdz(int x, int y, int z)
                {
                        return Base::h(x) * Base::h(y) * Base::h1d(z);
                }

                static constexpr double h2dz(int x, int y, int z)
                {
                        return Base::h(x) * Base::h(y) * Base::h2d(z);
                }

                static constexpr double h3dxz(int x, int y, int z)
                {
                        return Base::h3d(x) * Base::h(y) * Base::h3d(z);
                }

                static constexpr double h3dyz(int x, int y, int z)
                {
                        return Base::h(x) * Base::h3d(y) * Base::h3d(z);
                }

                T& v(int x, int y, int z)
                {
                        return Base::values_[index3(x, y, z)];
                }

                const T& v(int x, int y, int z) const
                {
                        return Base::values_[index3(x, y, z)];
                }

                const double vx(int x, int y, int z) const
                {
                        return v(x, y, z) * h1dx(x, y, z);
                }

                const double vy(int x, int y, int z) const
                {
                        return v(x, y, z) * h1dy(x, y, z);
                }
         
                const double vz(int x, int y, int z) const
                {
                        return v(x, y, z) * hdz(x, y, z);
                }

                const double vxx(int x, int y, int z) const
                {
                        return v(x, y, z) * h2dx(x, y, z);
                }

                const double vyy(int x, int y, int z) const
                {
                        return v(x, y, z) * h2dy(x, y, z);
                }

                const double vzz(int x, int y, int z) const
                {
                        return v(x, y, z) * h2dz(x, y, z);
                }

                const double vxy(int x, int y, int z) const
                {
                        return v(x, y, z) * h3dxy(x, y, z);
                }
                        
                const double vxz(int x, int y, int z) const
                {
                        return v(x, y, z) * h3dxz(x, y, z);
                }
                
                const double vyz(int x, int y, int z) const
                {
                        return v(x, y, z) * h3dyz(x, y, z);
                }
        }; // class Differential<ThreeDimension, T>

} // namespace lsm

#endif // LevelSetMethod2_Differential_h
