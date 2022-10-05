//--------------------------------------------------------------------------------------------------------------------------
#ifndef FRONTDIRECTION_H_
#define FRONTDIRECTION_H_
//--------------------------------------------------------------------------------------------------------------------------

namespace kmd {
        //******************************************************************************************************************
        class FrontDirection {
        //******************************************************************************************************************
        public:
                enum Direction {
                        Inner,
                        Outer,
                };
                
                /*!
                 *      @param[in]      kind
                 */
                FrontDirection(Direction kind = Inner)
                        : direction_(kind) {}
                
                Direction get_value() const {
                        return direction_;
                }      
                
                bool is_inner() const {
                        return direction_ == Inner;
                }
                
                bool is_outer() const {
                        return direction_ == Outer;        
                }
                
                /*!
                 *      @param[in]      x
                 */
                int attach_sign(int x) const {
                        if ( direction_ == Inner ) {
                                return -x;
                        }
                        else {
                                return x;
                        }
                }
                
        private:
                Direction       direction_;
                FrontDirection(const FrontDirection&);
        };
}

#endif /*FRONTDIRECTION_H_*/
