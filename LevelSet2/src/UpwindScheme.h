#ifndef LevelSetMethod2_UpwindSheme_h
#define LevelSetMethod2_UpwindSheme_h

#include "NeighboringPoints.h"
#include "Indexer.h"

namespace lsm
{

        struct PositiveSpeed {};
        struct NegativeSpeed {};

        /// implements the upwind scheme
        /**
         *      @tparam D the dimension type
         *      @see Level Set Methods and Fast Marching Methods, J. A. Sethian, 
         *      Cambridge University Press, ISBN 0-521-64557-3, p65
         */
        template<typename D>
        class UpwindScheme
        {
                BOOST_CLASS_REQUIRE(D, lsm, DimensionConcept);
                friend class UpwindSchemeTester;

        private:
                const Indexer<D>& indexer_;
                const std::vector<double>& phi_;

        public:
                UpwindScheme(
                        const Indexer<D>&               indexer,
                        const std::vector<double>&      phi)
                        : indexer_(indexer)
                        , phi_(phi) {}

                ~UpwindScheme() = default;

                // test ok
                template<typename U>
                const double calculate(const IntPoint<D>& p, TwoDimension)
                {
                        set_position(p);
                        calculate_with_2d(U());
                        return std::sqrt(
                                          fdxm_ * fdxm_
                                        + fdxp_ * fdxp_
                                        + fdym_ * fdym_
                                        + fdyp_ * fdyp_
                                );
                }

                // test ok
                template<typename U>
                const double calculate(const IntPoint<D>& p, ThreeDimension)
                {
                        set_position(p);
                        calculate_with_3d(U());
                        return std::sqrt(
                                          fdxm_ * fdxm_
                                        + fdxp_ * fdxp_
                                        + fdym_ * fdym_
                                        + fdyp_ * fdyp_
                                        + fdzp_ * fdzp_
                                        + fdzm_ * fdzm_
                                );
                }

                template<typename U>
                void calculate_(const IntPoint<D>& p, TwoDimension)
                {
                        set_position(p);
                        calculate_with_2d(U());
                        upwind_scheme_difference_ =
                                std::sqrt(
                                          fdxm_ * fdxm_
                                        + fdxp_ * fdxp_
                                        + fdym_ * fdym_
                                        + fdyp_ * fdyp_
                                );
                        
                        calculate_central_difference_with_2d(); 
                        central_difference_ = std::sqrt(dx_ * dx_ + dy_ * dy_);
                }

                
                void calculate_only_central_difference(const IntPoint<D>& p, TwoDimension)
                {
                        set_position(p);
                        upwind_scheme_difference_ = 0;
                        calculate_central_difference_with_2d(); 
                        central_difference_ = std::sqrt(dx_ * dx_ + dy_ * dy_);
                }

                template<typename U>
                void calculate_(const IntPoint<D>& p, ThreeDimension)
                {
                        set_position(p);
                        calculate_with_3d(U());
                        upwind_scheme_difference_ =
                                std::sqrt(
                                          fdxm_ * fdxm_
                                        + fdxp_ * fdxp_
                                        + fdym_ * fdym_
                                        + fdyp_ * fdyp_
                                        + fdzp_ * fdzp_
                                        + fdzm_ * fdzm_
                                );
                        
                        calculate_central_difference_with_3d(); 
                        central_difference_ = std::sqrt(dx_ * dx_ + dy_ * dy_ + dz_ * dz_);
                }

                
                void calculate_only_central_difference(const IntPoint<D>& p, ThreeDimension)
                {
                        set_position(p);
                        upwind_scheme_difference_ = 0;
                        calculate_central_difference_with_3d(); 
                        central_difference_ = std::sqrt(dx_ * dx_ + dy_ * dy_ + dz_ * dz_);
                }

                const double& get_central_difference() const
                {
                        return central_difference_;
                }

                const double& get_upwind_scheme_difference() const
                {
                        return upwind_scheme_difference_;
                }

        private:
                int left_;
                int right_;
                int top_;
                int self_;
                int bottom_;
                int front_;
                int back_;

                /// upwind scheme
                double fdxm_;
                double fdxp_;
                double fdym_;
                double fdyp_;
                double fdzm_;
                double fdzp_;

                /// central difference
                double dx_;
                double dy_;
                double dz_;

                double central_difference_;
                double upwind_scheme_difference_;

                // test ok
                void set_position(const IntPoint<TwoDimension>& p)
                {
                        left_   = indexer_(p + NeighboringPoints2d(-1,  0));
                        right_  = indexer_(p + NeighboringPoints2d( 1,  0));
                        self_   = indexer_(p);
                        top_    = indexer_(p + NeighboringPoints2d( 0, -1));
                        bottom_ = indexer_(p + NeighboringPoints2d( 0,  1));
                }
                
                void set_position(const IntPoint<ThreeDimension>& p)
                {
                        left_   = indexer_(p + NeighboringPoints3d(-1,  0,  0));
                        right_  = indexer_(p + NeighboringPoints3d( 1,  0,  0));
                        self_   = indexer_(p);
                        top_    = indexer_(p + NeighboringPoints3d( 0, -1,  0));
                        bottom_ = indexer_(p + NeighboringPoints3d( 0,  1,  0));
                        front_  = indexer_(p + NeighboringPoints3d( 0,  0, -1));
                        back_   = indexer_(p + NeighboringPoints3d( 0,  0,  1));
                }

                // test ok
                void calculate_with_2d(PositiveSpeed)
                {
                        fdxm_ = std::max(phi_[self_]   - phi_[left_],   0.0);
                        fdxp_ = std::min(phi_[right_]  - phi_[self_],   0.0);
                        fdym_ = std::max(phi_[self_]   - phi_[top_],    0.0);
                        fdyp_ = std::min(phi_[bottom_] - phi_[self_],   0.0);
                }
               
                // test ok
                void calculate_with_2d(NegativeSpeed)
                {
                        fdxp_ = std::max(phi_[right_]  - phi_[self_],   0.0);
                        fdxm_ = std::min(phi_[self_]   - phi_[left_],   0.0);
                        fdyp_ = std::max(phi_[bottom_] - phi_[self_],   0.0);
                        fdym_ = std::min(phi_[self_]   - phi_[top_],    0.0);
                }

                void calculate_central_difference_with_2d()
                {
                        dx_ = 0.5 * (phi_[right_]  - phi_[left_]);
                        dy_ = 0.5 * (phi_[bottom_] - phi_[top_]);
                }
                
                void calculate_with_3d(PositiveSpeed)
                {
                        calculate_with_2d(PositiveSpeed()); 
                        fdzm_ = std::max(phi_[self_] - phi_[front_],    0.0);
                        fdzp_ = std::min(phi_[back_] - phi_[self_],     0.0);
                }

                void calculate_with_3d(NegativeSpeed)
                {
                        calculate_with_2d(NegativeSpeed());
                        fdzp_ = std::max(phi_[back_] - phi_[self_],     0.0);
                        fdzm_ = std::min(phi_[self_] - phi_[front_],    0.0);
                }
         
                void calculate_central_difference_with_3d()
                {
                        calculate_central_difference_with_2d();
                        dz_ = 0.5 * (phi_[back_] - phi_[front_]);
                }
                
                UpwindScheme(const UpwindScheme&) = delete;
                UpwindScheme& operator=(const UpwindScheme&) = delete;
                 
                UpwindScheme(UpwindScheme&&) = delete;
                UpwindScheme& operator=(UpwindScheme&&) = delete;
        }; // class UpwindScheme
} // namespace lsm
#endif // LevelSetMethod2_UpwindSheme_h 
