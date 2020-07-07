/****************************************************************************
 * Container class for the simulation
 ***************************************************************************/

#ifndef WORLD_H
#define WORLD_H

#include <vector>
#include "const.h"
#include "resid.h"
#include "partk.h"
#include "memor.h"
#include "bound.h"
#include "angle.h"
#include "dihed.h"
#include "tlist.h"

namespace MD
{
  class tlist;
  class resid;
  class partk;
  class bound;
  class angle;
  class dihed;
  class pmemo;
  class auxMm;
  class world
  {
    public:
    //CONSTRUCTORS
      explicit world();
      world(world& syst);
      ~world();

    //MODIFIERS
      void adRsid(istream& iput, const tlist* types);
      void adPart(istream& iput);
      void adBond(istream& iput);
      void adAngl(istream& iput);
      void adDhed(istream& iput);
      void fitSze();

      void mkBcon(const int sBcon);
      void mkEcon(const int sEcon);
      void mkCcon(const int sCcon);
      void mkLcon(const int sLcon);
      void mkDcon(const int sDcon);
      void mkVcon(const int sVcon);
      void mkIcon(const int sIcon);

      void mkPmOf(const int sPmOf);
      void mkOmOf(const int sOmOf);

      void mkOnLj(const myCoor sOnLj);
      void mkOfLj(const myCoor sOfLj);
      void mkOfVl(const myCoor sOfVl);
      void mkRtVl(const myCoor sRtVl);
      void mkOfDp(const myCoor sOfDp);
      void mkRtUp(const myCoor sRtUp);
      void mkTemp(const myCoor sTemp);
      void mkTime(const myCoor sTime);

      void mkEdge(const myVect* sEdge)
      { myVectCpy(edge, sEdge); } 

      void mkCond(ostream& oput, istream& iput);

      void iniCon(istream& iput);
      void iniSys(istream& iput, const tlist* types);
      void iniTrj(istream& iput);
      void skpTrj(istream& iput);

      void adMemo();
      void iniMem(istream& iput);

      void iniRnd();
      void rndCrd();
      void rndVel();

      void nwTraj();
      void nwMemo();
      void nwVlst();
      void nwAccl();
      void impuls();

      void apBcon(myVect* uVec);
      void apVcon();
      void apAMem();
      void apVMem();

      //NOTE: for clarity, the calculation is done in a visualizable route
      //      it may be changed to a more efficient approach in the future
      void ljonSz(const partk* atom1, const partk* atom2,
                  myCoor& dist, myVect* len1);
      void bondSz(bound* bond, myVect* len1);
      void anglSz(angle* angl, myVect* len1, myVect* len2);
      void dhedSz(dihed* dhed, myVect* len1, myVect* len2,
                  myVect* len3, myVect* nrm1, myVect* nrm2);
     
      void nbndFc(const partk* atom1, const partk* atom2, 
                  const int ctrl, myVect* fVec);
      void bondFc(bound* bond, const int ctrl, myVect* for1);
      void anglFc(angle* angl, const int ctrl, myVect* for1, myVect* for2);
      void dhedFc(dihed* dhed, const int ctrl, myVect* for1, myVect* for2, 
                                               myVect* for3);

      void ljonFc(const partk* atom1, const partk* atom2, const int ctrl,
                  const myCoor dist, const myVect* uVec, myVect* fVec);
      void dispFc(const partk* atom1, const partk* atom2,
                  const myCoor dist, const myVect* uVec, myVect* fVec);
      void randFc(const partk* atom1, const partk* atom2,
                  const myCoor dist, const myVect* uVec, myVect* fVec);

      myCoor ljonSt(const partk* atom1, const partk* atom2, const int ctrl,
                    const myCoor dist);
      myCoor bondSt(const bound* bond, const int ctrl);
      myCoor anglSt(const angle* angl, const int ctrl);
      myCoor dhedSt(const dihed* dhed, const int ctrl);
      myCoor pDpdSt(const partk* atom1, const partk* atom2,
                    const myCoor dist);
      myCoor oDpdSt(const partk* atom1, const partk* atom2,
                    const myCoor dist);

      void verlet();
      void lpFrgA();
      void lpFrgB();
      void vVerlA();
      void vVerlB();

    //OBSERVERS
      void wrCond(ostream& oput) const;
      void wrTraj(ostream& oput) const;
      void wrCoor(ostream& oput) const;

      const vector<resid>& obRsidLs() const { return rsids; }
      const vector<partk>& obPartLs() const { return parts; }
      const vector<bound>& obBondLs() const { return bonds; }
      const vector<angle>& obAnglLs() const { return angls; }
      const vector<dihed>& obDhedLs() const { return dheds; }

      const resid* obRsid(const string name) const;
      const partk* obPart(const string name) const;
      const bound* obBond(const partk* atom1, const partk* atom2)
                          const;
      const angle* obAngl(const partk* atom1, const partk* atom2,
                          const partk* atom3) const;
      const dihed* obDhed(const partk* atom1, const partk* atom2,
                          const partk* atom3, const partk* atom4) const;

      const resid* obRsid(const int indx) const;
      const partk* obPart(const int indx) const;
      const bound* obBond(const int indx) const;
      const angle* obAngl(const int indx) const;
      const dihed* obDhed(const int indx) const;

      int obNresid() const { return rsids.size(); }
      int obNpartk() const { return parts.size(); }
      int obNbound() const { return bonds.size(); }
      int obNangle() const { return angls.size(); }
      int obNdihed() const { return dheds.size(); }

      int obBcon() const { return boundCon; }
      int obEcon() const { return ensemCon; }
      int obCcon() const { return colmbCon; }
      int obLcon() const { return ljPotCon; }
      int obDcon() const { return disspCon; }
      int obVcon() const { return verltCon; }
      int obIcon() const { return intgrCon; }

      int obPmOf() const { return pmOf; }
      int obOmOf() const { return omOf; }

      myCoor obOnLj() const { return onLj; }
      myCoor obOfLj() const { return ofLj; }
      myCoor obOfVl() const { return ofVl; }
      myCoor obRtVl() const { return rtVl; }
      myCoor obOfDp() const { return ofDp; }
      myCoor obRtUp() const { return rtUp; }
      myCoor obTime() const { return time; }
      myCoor obTemp() const { return temp; }
      myCoor obSysT() const
      { return (this->obKine() / (parts.size() * DIM / 2)); }

      const myVect* obEdge() const { return edge; }

      const myVect* obMomt() const;
      myCoor obKine() const;
      myCoor obPote() const;
      myCoor obEner() const;
      myCoor obPFVC(const vector<myVect*> *oVel) const;
      myCoor obPVVC(const vector<myVect*> *oVel) const;
      myCoor obtVAC(const vector<myVect*> *oVel) const;
      myCoor obtMSD(const vector<myVect*> *oPos) const;
      vector<int> obLjDs(const ptype* ptyp1, const ptype* ptyp2,
                         const int bSze) const;
      vector<int> ob12Ds(const btype* btyp, const int bSze) const;
      vector<int> ob13Ds(const atype* atyp, const int bSze) const;
    private:
      vector<resid> rsids;
      vector<partk> parts;
      vector<bound> bonds;
      vector<angle> angls;
      vector<dihed> dheds;
      //boundary control: 1 for periodic, 2 for rigd, 3 for repulsive
      //Ensemble control: 1 for NVE, 2 for NVT, 3 for NPT 
      //Coulomb force control: 1 for on
      //LJ force control: 1 for shifted potential, 2 for smooth potential, 
      //                  3 for truncated force, 4 for shifted force,
      //                  5 for smooth force
      //disspative force control: 1 for uncorrelated, 2 for parallel,
      //                          3 for orthogonal, 4 for isotropic,
      //                          5 for parallel MZ, 6 for orthogonal MZ,
      //                          7 for full MZ, 8 for parallel auxiliary MZ,
      //                          9 for orthogonal auxiliary MZ, 
      //                          10 for full auxiliary MZ
      //verlet list control: 1 for constant update rate, 
      //                     2 for automatic rate, 
      //                     3 for using with cell list,
      //                     4 for using automatic rate AND cell list 
      //integrator control: 1 for Verlet, 
      //                    2 for leapfrog, 21 for DPD-LF,
      //                    22 for impulsive DPD-LF,
      //                    3 for velocity verlet, 31 for DPD-VV,
      //                    32 for consistent DPD-VV, 33 for Lowe method
      //                    34 for impulsive DPD-VV
      int boundCon, ensemCon, colmbCon, ljPotCon, disspCon, 
          verltCon, intgrCon;
      //Memory kernel factor
      vector<myCoor> pDCor,    //parallel friction factor
                     oDCor,    //orthogonal friction factor
                     pRCor,    //parallel random factor
                     oRCor;    //orthogonal random factor
      //Auxiliary variables factor
      //A_sp, A_ss and B_sr in parallel and orthogonal direction
      //B_sr is diagonal, so it's sufficient to use vector
      myVect *pAsp, *oAsp, *pDsr, *oDsr;
      myMatr *oAss, *oBsr;
      vector<myMatr*> pAss, pBsr;
      int pmOf, omOf;          //size of the memory kernels or
                               //number of auxiliary variables
      myCoor onLj, ofLj,       //cuton and cutoff for LJ potential        
             ofVl, rtVl,       //cutoff and rate for Verlet list
             ofDp,             //cutoff for DPD
             rtUp,             //rate for integrator
             temp,             //temperature and inverse
             time;
      myVect *edge;            //lattice vector

    //SUPPORT VARIABLES (no access to outside)
      //force control : 2 for conservative force only,
      //                3 for friction force only, 
      //                5 for random force only,
      //                6 for conservative and friction force,
      //                10 for conservative and random force,
      //                15 for friction and random force,
      //                30 for all 3 forces,
      //homogeneity control: 1 for on
      int forceCon, homogCon;
      //memory kernel for MZ (as circular-linked lists)
      pmemo xMemo, uMemo, vMemo, pRMem, oRMem;
      //auxiliary variables for MZ
      vector<auxMm> pAxMm, oAxMm;
      //supports for Verlet list
      //BUG: maybe it's better to store the position in the parts list
      //     instead
      //BUG: maybe it's structurally compatible with MZ by storing a
      //     boolean list instead
      vector<partk*> vList;
      vector<int> vLens;
      vector<myVect*> vPost;
      myCoor invT, sqrT;              //inverse and square root of temperature
      myVect *fTot, *forc, *unit;     //dummy vectors for force calculation
      myVect *velo, *post;            //dummy vectors for MZ
      myVect *pAcc, *oAcc, 
             *pFor, *oFor;            //dummy vectors for auxiliary
      myVect *pbc1, *pbc2;            //dummy vectors for PBC
  };
  
  //NON-MEMBER FUNCTIONS
  template <class type> int obPInd(const vector<type>* parts,
                                   const type* atom1, const type* atom2);
}
#endif
