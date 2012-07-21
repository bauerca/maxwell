
#include <map>

#include "MxTypes.h"


class A {
  private:
    int i;
};


class AContainer {
  public:
    AContainer() : count(0) {}

    ~AContainer() {
      typename std::map<int, RCP<A> >::const_iterator iter;
      for (iter = mMap.begin(); iter != mMap.end(); iter++)
        std::cout << iter->first << "\n"
            << "  strong: " << iter->second.strong_count() << "\n"
            << "  weak: " << iter->second.weak_count() << "\n"
            << "  has own: " << iter->second.has_ownership() << "\n"
            << "  address: " << iter->second.get() << "\n";
    }

    void addA(bool bad) {
      if (bad)
        mMap.insert(std::make_pair(count, new A()));
      else
        mMap.insert(std::make_pair(count, rcp(new A())));
      count++;
      //return a;
    }

    RCP<A> getA(int i) const {
      std::map<int, RCP<A> >::const_iterator it;
      it = mMap.find(i);
      return it->second;
    }


  private:
    std::map<int, RCP<A> > mMap;
    
    int count;
};


class AFusser {
  public:
    AFusser() : mACont(rcp(new AContainer())), mA(Teuchos::null) {
      mACont->addA(false);
      mACont->addA(true);
      mA = mACont->getA(0);
    }

  private:
    RCP<A> mA;
    RCP<AContainer> mACont;
};




int main() {
  
  //RCP<A> a;

  //RCP<AContainer> cont = rcp(new AContainer());
  //a = cont->addA();

  ////RCP<A> a = rcp(new A());
  ////RCP<A> a2 = a;
  //RCP<A> a2;


  //a2 = cont->getA(0);

  //std::cout << "a  has own: " << a.has_ownership() << "\n";
  //std::cout << "a2 has own: " << a2.has_ownership() << "\n";

  //cont = Teuchos::null;

  RCP<AFusser> fusser = rcp(new AFusser());

  //RCP<MxEMOps<3, double> > emOps = rcp(new MxEMOps<3, double>(

  return 0;
}
