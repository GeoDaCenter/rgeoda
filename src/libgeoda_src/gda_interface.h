#ifndef __GDA_INTERFACE_H__
#define __GDA_INTERFACE_H__

#include <vector>
//#include "geofeature.h"

// forward declaration
namespace gda {
    struct PointContents;
    struct MainMap;
}

class AbstractGeoDa {
public:
    AbstractGeoDa() {};
    virtual ~AbstractGeoDa() {};

    virtual int GetNumObs() const = 0;

    virtual const std::vector<gda::PointContents*>& GetCentroids() = 0;

    virtual int GetMapType() = 0;

    virtual gda::MainMap& GetMainMap() = 0;

};
#endif
