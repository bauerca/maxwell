#ifndef MX_CRAB_CAV
#define MX_CRAB_CAV

#include <string>
#include <cmath>

#include "MxShape.hpp"
#include "MxShapeSubtract.hpp"
#include "MxShapeUnion.hpp"
#include "MxShapeIntersection.hpp"
#include "MxTorus.hpp"
#include "MxCylinder.hpp"
#include "MxSlab.hpp"
#include "MxCone.hpp"
#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"

class MxCrabCav : public MxShape<3> {
  public:
    MxCrabCav(MxDimVector<double, 3> midPt, int myNumCells, double myCellLength, double myIrisRadius, double myCavRadius, double myIrisRho, double myCavRho);

    virtual ~MxCrabCav();

    virtual double func(const MxDimVector<double, 3> & p) const;

    virtual MxDimVector<double, 3> gradFunc(const MxDimVector<double, 3> & p) const;

    virtual bool hasGradFunc() const {return true;}

    virtual MxShape<3>::DimVecDblPair boundingBox() const {return bbox;}

  private:

    MxShape<3>::DimVecDblPair bbox;

    // inputs
    MxDimVector<double, 3> crabCavMidPt;
    int numCells;
    double cellLen, halfCellLen, irisRad, cavRad, irisRho, cavRho;

    // calculated vars
    MxDimVector<double, 3> zhat;
    double irisTorMajRad, irisTorMinRad, cavTorMajRad, cavTorMinRad;
    double cosTheta, theta, sinTheta, tanTheta, cotTheta;
    double coneVertOffset;

    void setVars();

    MxShapeUnion<3> * halfCell;

    MxShape<3> * slab;

    std::vector<MxShape<3> *> shapes;

    MxShape<3> * create(MxShape<3> * aShape);
};


MxShape<3> * MxCrabCav::create(MxShape<3> * aShape) {
  shapes.push_back(aShape);
  return aShape;
}
 

void MxCrabCav::setVars() {
  this->name = "crab_cavity";

  zhat = 0; zhat[2] = 1;

  irisTorMajRad = irisRad + irisRho;
  irisTorMinRad = irisRho;

  cavTorMajRad = cavRad - cavRho;
  cavTorMinRad = cavRho;

  // find angle!
  double rhoSum  = cavRho + irisRho;
  double rhoSum2 = rhoSum * rhoSum;
  double radDiff = cavRad - irisRad;
  double halfCellLen2 = 0.25 * cellLen * cellLen;
  double diff2 = pow(radDiff - rhoSum, 2);

  cosTheta = (rhoSum - radDiff) * rhoSum;
  cosTheta += sqrt(halfCellLen2 * (diff2 - rhoSum2 + halfCellLen2));
  cosTheta /= halfCellLen2 + diff2;

  theta = acos(cosTheta);
  sinTheta = sqrt(1 - cosTheta * cosTheta);
  tanTheta = sinTheta / cosTheta;
  cotTheta = 1.0 / tanTheta;
  

  //std::cout << "sinTheta = " << sinTheta << "\n";
  //std::cout << "cosTheta = " << cosTheta << "\n";
  //std::cout << "tanTheta = " << tanTheta << "\n";
  //std::cout << "cotTheta = " << cotTheta << "\n";

  coneVertOffset = 0.5 * cellLen - irisRho * sinTheta + (irisRad + irisRho * (1.0 - cosTheta)) * cotTheta;
  //std::cout << "coneVertOffset = " << coneVertOffset << "\n";

  bbox.first[2] = crabCavMidPt[2] - 0.5 * double(numCells) * cellLen;
  bbox.first[1] = crabCavMidPt[1] - cavRad;
  bbox.first[0] = crabCavMidPt[0] - cavRad;

  bbox.second[2] = crabCavMidPt[2] + 0.5 * double(numCells) * cellLen;
  bbox.second[1] = crabCavMidPt[1] + cavRad;
  bbox.second[0] = crabCavMidPt[0] + cavRad;
}

MxCrabCav::~MxCrabCav() {
  for (size_t i = 0; i < shapes.size(); ++i)
    delete shapes[i];
}

MxCrabCav::MxCrabCav(MxDimVector<double, 3> midPt, int myNumCells, double myCellLength, 
double myIrisRadius, double myCavRadius, double myIrisRho, double myCavRho) : 
crabCavMidPt(midPt),
numCells(myNumCells), 
cellLen(myCellLength), 
halfCellLen(0.5 * myCellLength), 
irisRad(myIrisRadius),
cavRad(myCavRadius), 
irisRho(myIrisRho), 
cavRho(myCavRho) {

  setVars();


  MxShape<3> * irisTube = create(new MxCylinder(0, zhat, irisRad + irisRho * (1.0 - cosTheta)));
  MxShape<3> * irisTorus = create(new MxTorus<3>(0.5 * cellLen * zhat, zhat, irisTorMajRad, irisTorMinRad));

  MxShapeSubtract<3> * corrugatedIrisTube = new MxShapeSubtract<3>();
  corrugatedIrisTube->setBaseShape(irisTube);
  corrugatedIrisTube->subtractShape(irisTorus);
  create(corrugatedIrisTube);


  MxShape<3> * cavTube = create(new MxCylinder(0, zhat, cavRad - cavRho * (1.0 - cosTheta)));
  MxShape<3> * cavCone = create(new MxCone(coneVertOffset * zhat, zhat, theta));
  MxShape<3> * cavTorus = create(new MxTorus<3>(0, zhat, cavTorMajRad, cavTorMinRad));
  
  MxShapeIntersection<3> * preCav = new MxShapeIntersection<3>();
  preCav->add(cavCone);
  preCav->add(cavTube);
  create(preCav);

  halfCell = new MxShapeUnion<3>();
  halfCell->add(preCav);
  halfCell->add(cavTorus);
  halfCell->add(corrugatedIrisTube);
  create(halfCell);

  slab = create(new MxSlab<3>(crabCavMidPt, zhat, double(numCells) * cellLen));

}

  
double MxCrabCav::func(const MxDimVector<double, 3> & p) const {

  MxDimVector<double, 3> pcopy(p - crabCavMidPt);
  pcopy[2] = halfCellLen - fabs(fmod(fabs(p[2] - crabCavMidPt[2]), cellLen) - halfCellLen);

  double fcell = halfCell->func(pcopy);
  double fslab = slab->func(p);

  // return intersection of fcell and fslab
  return fcell < fslab ? fcell : fslab;
}


MxDimVector<double, 3> MxCrabCav::gradFunc(const MxDimVector<double, 3> & p) const {

  MxDimVector<double, 3> pcopy(p - crabCavMidPt);
  pcopy[2] = halfCellLen - fabs(fmod(p[2] - crabCavMidPt[2], cellLen) - halfCellLen);

  double fcell = halfCell->func(pcopy);
  double fslab = slab->func(p);

  MxDimVector<double, 3> gradfCell(halfCell->gradFunc(pcopy));
  MxDimVector<double, 3> gradfSlab(slab->gradFunc(p));

  return fcell < fslab ? gradfCell : gradfSlab;
}



#if 0

  MxShapeSubtract<3> * corrugatedIrisTube = new MxShapeSubtract<3>();
  MxShape<3> * irisTube = create(new MxCylinder(crabCavMidPt, zhat, irisRad + irisRho * (1.0 - cosTheta)));
  corrugatedIrisTube->setBaseShape(irisTube);

  // get all the tori that describe the irises

  MxDimVector<double, 3> irisMidPt = crabCavMidPt;
  irisMidPt[2] -= 0.5 * (double(numCells) - 1.0) * cellLen;
  for (int cell = 0; cell < numCells; ++cell) {
    corrugatedIrisTube->subtractShape(create(new MxTorus(irisMidPt, zhat, irisTorMajRad, irisTorMinRad)));
    irisMidPt[2] += cellLen;
  }
  create(corrugatedIrisTube);

  crabCav = corrugatedIrisTube;


#if 1

  // Now we form the cavity pieces

  MxDimVector<double, 3> cavMidPt = crabCavMidPt;
  cavMidPt[2] -= 0.5 * double(numCells) * cellLen;

  MxShape<3> * cavTube = create(new MxCylinder(crabCavMidPt, zhat, cavRad - cavRho * (1.0 - cosTheta)));
  MxShape<3> * cavTorus, * cavCone1, * cavCone2, * slab;
  MxShapeIntersection<3> * preCav;
  MxShapeUnion<3> * cav;

  MxShapeUnion<3> * allCavs = new MxShapeUnion<3>();

  for (int cell = 0; cell <= numCells; ++cell) {
    // first make a torus
    cavTorus = create(new MxTorus(cavMidPt, zhat, cavTorMajRad, cavTorMinRad));
    //std::cout << "cell " << cell << ", cavTorus.func = " << cavTorus->func(crabCavMidPt) << "\n";
    
    // make two cones
    cavCone1 = create(new MxCone(cavMidPt - coneVertOffset * zhat, zhat, theta));
    cavCone2 = create(new MxCone(cavMidPt + coneVertOffset * zhat, zhat, theta));
    //std::cout << "cell " << cell << ", cone1.func = " << cavCone1->func(crabCavMidPt) << "\n";
    //std::cout << "cell " << cell << ", cone2.func = " << cavCone2->func(crabCavMidPt) << "\n";

    // bound cones at vertices
    slab = create(new MxSlab<3>(cavMidPt, zhat, 2.0 * coneVertOffset));
    //std::cout << "cell " << cell << ", slab.func = " << slab->func(crabCavMidPt) << "\n";

    // need to know about derived member functions
    preCav = new MxShapeIntersection<3>();

    preCav->add(cavCone1);
    preCav->add(cavCone2);
    preCav->add(slab);
    preCav->add(cavTube);
    //std::cout << "cell " << cell << ", preCav.func = " << preCav->func(crabCavMidPt) << "\n";

    create(preCav);

    // need to know about derived member functions
    cav = new MxShapeUnion<3>();

    cav->add(preCav);
    cav->add(cavTorus);
    //std::cout << "cell " << cell << ", cav.func = " << cav->func(crabCavMidPt) << "\n";

    create(cav);

    allCavs->add(cav);

    cavMidPt[2] += cellLen;
  }

  create(allCavs);
  //std::cout << "allCavs.func = " << allCavs->func(crabCavMidPt) << "\n";


  // merge allCavs with corrugatedIrisTube
  MxShapeUnion<3> * openCrabCav = new MxShapeUnion<3>();
  openCrabCav->add(allCavs);
  openCrabCav->add(corrugatedIrisTube);
  create(openCrabCav);
  
  // cap ends of the thing
  MxShapeIntersection<3> * crabCavPre = new MxShapeIntersection<3>();
  crabCavPre->add(openCrabCav);
  crabCavPre->add(create(new MxSlab<3>(crabCavMidPt, zhat, double(numCells) * cellLen)));

  crabCav = create(crabCavPre);
#endif

  this->name = "crab_cavity";

}
#endif


#endif
