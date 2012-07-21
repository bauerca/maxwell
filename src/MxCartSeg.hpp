#ifndef MX_CART_SEG
#define MX_CART_SEG

#include "MxDimVector.hpp"
#include "MxDynMatrix.hpp"
#include "MxSegment.h"

template<size_t DIM>
class MxCartSeg : public MxSegment<DIM> {
  public:

    MxCartSeg(char aType, double aLen) : MxSegment<DIM>(aLen, MxDimVector<double, DIM>(1)) {
      this->dir = 0;
      switch (aType) {
        case 'x':
          this->dir[0] = 1;
          break;
        case 'y':
          this->dir[1] = 1;
          break;
        case 'z':
          this->dir[2] = 1;
          break;
        default:
          std::cout << "MxCartSeg: invalid type given. Should be 'x', 'y', or 'z'.";
          throw 1;
      }
    }

    ~MxCartSeg() {}
};

#endif
