#include "dihed.h"
#include <cassert>
#include <cmath>

using namespace std;

namespace MD
{
  //CONSTRUCTORS
  dihed::dihed(const partk* atom1, const partk* atom2, const partk* atom3,
               const partk* atom4, const dtype* sType) {
    assert(!(atom1 == atom2 || atom1 == atom3 || atom1 == atom4 ||
             atom2 == atom3 || atom2 == atom4 || atom3 == atom4) &&
           (*atom1).obRsid() == (*atom2).obRsid() &&
           (*atom2).obRsid() == (*atom3).obRsid() &&
           (*atom3).obRsid() == (*atom4).obRsid());
    uppr = atom1;
    left = atom2;
    righ = atom3;
    lowr = atom4;
    bLef = obBond(uppr, left);
    bMid = obBond(left, righ);
    bRig = obBond(righ, lowr);
    aLef = obAngl(uppr, left, righ);
    aRig = obAngl(left, righ, lowr);
    mole = (*left).obRsid();
    type = sType;
  }

  dihed::dihed(istream& iput, const world* syst) {
    string name1, name2, name3, name4;
    iput >> name1 >> name2 >> name3 >> name4;
    uppr = (*syst).obPart(name1);
    left = (*syst).obPart(name2);
    righ = (*syst).obPart(name3);
    lowr = (*syst).obPart(name4);
    if(uppr != NULL && left != NULL && righ != NULL && lowr != NULL ) {
      assert((*uppr).obRsid() == (*left).obRsid() &&
             (*left).obRsid() == (*righ).obRsid() &&
             (*righ).obRsid() == (*lowr).obRsid());
      mole = (*left).obRsid();
      bLef = obBond(uppr, left);
      bMid = obBond(left, righ);
      bRig = obBond(righ, lowr);
      aLef = obAngl(uppr, left, righ);
      aRig = obAngl(left, righ, lowr);
      const string type1 = (*uppr).obMark();
      const string type2 = (*left).obMark();
      const string type3 = (*righ).obMark();
      const string type4 = (*lowr).obMark();
      type = (*((*mole).obRtyp())).obDtyp((*uppr).obMark(), (*left).obMark(),
                                          (*righ).obMark(), (*lowr).obMark());
    }
  }

  //MODIFIERS
  const dihed& dihed::operator =(const dihed& dhed) {
    if(this != &dhed) {
      type = dhed.obType();
    }
    return *this;
  }

  //OBSERVERS
  //BUG: need energy formula
  myCoor dihed::obEner() const { 
  }

  bool dihed::operator ==(const dihed& dhed) const {
    return (uppr == dhed.obUppr() && left == dhed.obLeft() &&
            righ == dhed.obRigh() && lowr == dhed.obLowr());
  }
}
