//----------------------------------------------------------------------------------------------------------------
#ifndef TUBEDOMAINGENERATOR_H_
#define TUBEDOMAINGENERATOR_H_
//----------------------------------------------------------------------------------------------------------------
#include <boost/noncopyable.hpp>
#include "Typedef.h"
//----------------------------------------------------------------------------------------------------------------
namespace kmd {

        class ShortestDistanceFinder;

        //********************************************************************************************************
        class TubeDomainGenerator : private boost::noncopyable {
        //********************************************************************************************************
        public:
                /*!
                 *      @param[in]      finder
                 */
                explicit TubeDomainGenerator(const ShortestDistanceFinder* finder);

                const Tube& get_tube() const {
                        return tube_;
                }

                void run();

        private:
                Tube                            tube_;
                const ShortestDistanceFinder*   finder_;
        };
}
#endif /*TUBEDOMAINGENERATOR_H_*/
