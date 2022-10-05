//-----------------------------------------------------------------------------------------------------------------------------
#ifndef TypedefH
#define TypedefH
//-----------------------------------------------------------------------------------------------------------------------------
#include <boost/shared_ptr.hpp>
#include <vector>
//#include <boost/gil/utilities.hpp>
#include <boost/gil/point.hpp>
#include <map>
//-----------------------------------------------------------------------------------------------------------------------------
namespace kmd {
        
        class Point;
        class ShortestDistanceFinder;
        class SpeedGenerator;
        class DistanceInitializer;
        class ZeroLevelSetDetector;
        class EdgeFrontGenerator;
        class TubeDomainGenerator;
        class FrontStopper;
        
        typedef std::vector<Point*>                             Front;
        typedef std::vector<Front>                              Fronts;
        typedef boost::shared_ptr<ShortestDistanceFinder>       ShortestDistanceFinderPtr;
        typedef boost::shared_ptr<SpeedGenerator>               SpeedGeneratorPtr;
        typedef boost::shared_ptr<DistanceInitializer>          DistanceInitializerPtr;
        typedef boost::shared_ptr<ZeroLevelSetDetector>         ZeroLevelSetDetectorPtr;
        typedef boost::gil::point<std::ptrdiff_t>              IntPoint;
        //kumada
        //typedef boost::gil::point2<std::ptrdiff_t>              IntPoint;
        typedef boost::shared_ptr<EdgeFrontGenerator>           EdgeFrontGeneratorPtr;
        typedef boost::shared_ptr<TubeDomainGenerator>          TubeDomainGeneratorPtr;
        typedef boost::shared_ptr<FrontStopper>                 FrontStopperPtr;
        typedef std::pair<int, int>                             Segment;
        typedef std::vector<Segment>                            Segments;
        typedef std::map<int, Segments>                         Tube;
                
}
#endif // TypedefH
