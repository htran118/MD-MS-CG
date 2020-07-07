#include <cassert>
#include <cmath>
#include <random>
#include <stdlib.h>
#include "world.h"

using namespace std;

namespace MD {
  //CONSTRUCTORS
  world::world() {
    boundCon = ensemCon = colmbCon = ljPotCon = verltCon = intgrCon = 0;
    forceCon = 0b0000000000000000;
    omOf = pmOf = 0;
    onLj = ofLj = ofVl = rtVl = ofDp = rtUp = time = temp = 0;
    edge = myVectMalloc(DIM);
    fTot = myVectMalloc(DIM);
    forc = myVectMalloc(DIM);
    unit = myVectMalloc(DIM);
    pbc1 = myVectMalloc(DIM);
    pbc2 = myVectMalloc(DIM);
  }

  world::world(world& systm) {
    rsids = systm.obRsidLs();
    parts = systm.obPartLs();
    bonds = systm.obBondLs();
    angls = systm.obAnglLs();
    boundCon = systm.obBcon();
    ensemCon = systm.obEcon();
    colmbCon = systm.obCcon();
    ljPotCon = systm.obLcon();
    disspCon = systm.obDcon();
    verltCon = systm.obVcon();
    intgrCon = systm.obIcon();
    forceCon = 0b0000000000000000;
    if(colmbCon != 0 || ljPotCon != 0)
      forceCon |= 0b0000000000000001;
    if(intgrCon != 24 && intgrCon != 33 && intgrCon == 34)
      if(disspCon != 0)
        forceCon |= 0b0000000000001110; 
    pmOf = systm.obPmOf();
    omOf = systm.obOmOf();
    onLj = systm.obOnLj();
    ofLj = systm.obOfLj();
    ofVl = systm.obOfVl();
    rtVl = systm.obRtVl();
    ofDp = systm.obOfDp();
    rtUp = systm.obRtUp();
    time = systm.obTime();
    temp = systm.obTemp();
    invT = 1 / BOLTZM / temp;
    sqrT = sqrt(BOLTZM * temp);
    edge = myVectMalloc(DIM);
    myVectCpy(edge, systm.obEdge());
    fTot = myVectMalloc(DIM);
    forc = myVectMalloc(DIM);
    unit = myVectMalloc(DIM);
    pbc1 = myVectMalloc(DIM);
    pbc2 = myVectMalloc(DIM);
    if(verltCon == 2)
      for(int i = 0; i < parts.size(); ++i) {
        vLens.push_back(0);
        vPost.push_back(myVectMalloc(DIM));
      }
  }

  world::~world() {
    myVectMemdel(edge);
    myVectMemdel(fTot);
    myVectMemdel(forc);
    myVectMemdel(unit);
    myVectMemdel(pbc1);
    myVectMemdel(pbc2);
    if(disspCon != 0 && disspCon != 1) {
      myVectMemdel(post);
      myVectMemdel(velo); 
      if(disspCon >= 8 && disspCon < 11) {
        myVectMemdel(pAcc);
        myVectMemdel(oAcc);
        myVectMemdel(pFor);
        myVectMemdel(oFor);
        for(int i = 0; i < 100; ++i) {
          myMatrMemdel(pAss.at(i));
          myMatrMemdel(pBsr.at(i));
        }
        myMatrMemdel(oAss);
        myMatrMemdel(oBsr);
      }
    }
    if(verltCon == 2)
      for(int i = 0; i < parts.size(); ++i)
        myVectMemdel(vPost.at(i));
  } 

  //MODIFIERS
  //add resid class
  void world::adRsid(istream& iput, const tlist* types) {
    rsids.push_back(resid(iput, types));
  }

  //add partk class
  void world::adPart(istream& iput) {
    parts.push_back(partk(iput, this));
    if(verltCon != 0) {
      vLens.push_back(0);
      if(verltCon == 2) {
        vPost.push_back(myVectMalloc(DIM));
        myVectCpy(vPost.at(vPost.size() - 1),
                                (parts.at(parts.size() - 1)).obPost());
      }
    }
  }

  //add bound class
  void world::adBond(istream& iput) {
    bonds.push_back(bound(iput, this));
  }

  //add angle class
  void world::adAngl(istream& iput) {
    angls.push_back(angle(iput, this));
  }

  inline void world::fitSze() {
    rsids.shrink_to_fit();
    parts.shrink_to_fit();
    bonds.shrink_to_fit();
    angls.shrink_to_fit();
    vPost.shrink_to_fit();
    vLens.shrink_to_fit();
    pAxMm.shrink_to_fit();
    xMemo.fitSze();
    uMemo.fitSze();
    vMemo.fitSze();
    pRMem.fitSze();
    oRMem.fitSze();
  }

  //set boundary control
  void world::mkBcon(const int sBcon) {
    assert(sBcon >= 1 && sBcon < 4);
    boundCon = sBcon;
  }

  //set ensemble control
  void world::mkEcon(const int sEcon) {
    assert(sEcon >= 1 && sEcon < 4);
    ensemCon = sEcon;
  }
  
  //set coulomb force control
  void world::mkCcon(const int sCcon) {
    assert(sCcon >= 0 && sCcon < 2);
    colmbCon = sCcon;
    if(colmbCon != 0)
      forceCon |= 0b0000000000000001;
  }

  //set LJ force control
  void world::mkLcon(const int sLcon) {
    assert(LJONPART && sLcon >= 0 && sLcon < 7);
    ljPotCon = sLcon;
    if(ljPotCon != 0)
      forceCon |= 0b0000000000000001;
  }

  //set dissipative force control
  void world::mkDcon(const int sDcon) {
    assert(DPDPARTK && sDcon >= 0 && sDcon < 11);
    disspCon = sDcon;
    if(disspCon != 0)
      forceCon |= 0b0000000000001110;
  }

  //set Verlet list control
  void world::mkVcon(const int sVcon) {
    assert(sVcon >= 0 && sVcon < 5);
    verltCon = sVcon;
  }

  //set integrator control
  void world::mkIcon(const int sIcon) {
    assert((sIcon >= 1 && sIcon < 4) || (sIcon >= 21 && sIcon < 25)
                                     || (sIcon >= 31 && sIcon < 35));
    intgrCon = sIcon;
    if(intgrCon == 24 || intgrCon == 33 || intgrCon == 34)
      if(disspCon != 0)
        forceCon &= 0b1111111111111001;
  }

  //set cutoff for memory kernel
  void world::mkPmOf(const int sPmOf) {
    assert (DPDPARTK && (disspCon == 5 || disspCon == 7 ||
                         disspCon == 8 || disspCon == 10) && sPmOf > 0);
    pmOf = sPmOf;
  }

  void world::mkOmOf(const int sOmOf) {
    assert (DPDPARTK && (disspCon == 6 || disspCon == 7 ||
                         disspCon == 9 || disspCon == 10) && sOmOf > 0);
    omOf = sOmOf;
  }

  //set cuton for LJ force
  void world::mkOnLj(const myCoor sOnLj) {
    assert(LJONPART && (ljPotCon == 3 || ljPotCon == 6) && sOnLj > 0);
    onLj = sOnLj;
  }

  //set cutoff for LJ force
  void world::mkOfLj(const myCoor sOfLj) {
    assert(LJONPART && (ljPotCon != 0 || ljPotCon != 1) && sOfLj > 0);
    ofLj = sOfLj;
  }

  //set cutoff for Verlet list
  void world::mkOfVl(const myCoor sOfVl) {
    assert(verltCon != 0 && sOfVl > 0);
    ofVl = sOfVl;
  }

  //set update rate for Verlet list
  void world::mkRtVl(const myCoor sRtVl) {
    assert(verltCon != 0 && sRtVl > 0);
    rtVl = sRtVl;
  }

  //set cutoff for dissipative force
  void world::mkOfDp(const myCoor sOfDp) {
    assert(DPDPARTK && (disspCon != 0 || disspCon != 1) && sOfDp > 0);
    ofDp = sOfDp;
  }

  //set update rate for heat bath
  void world::mkRtUp(const myCoor sRtUp) {
    assert((intgrCon == 32 || intgrCon == 33) && sRtUp > 0);
    rtUp = sRtUp;
  }

  //set temperature
  void world::mkTemp(const myCoor sTemp) {
    assert(sTemp > 0);
    temp = sTemp;
    invT = 1 / BOLTZM * temp;
    sqrT = sqrt(BOLTZM * temp);
  }

  //set time
  void world::mkTime(const myCoor sTime) {
    assert(sTime >= 0);
    time = sTime;
  }

  //initialize world
  //set up the topology and structure of the world
  void world::iniSys(istream& iput, const tlist* types) {
    string ctrl, prev;
    while(iput >> ctrl) {
      if(ctrl == ctrlIgnr) {
        iput.ignore(256, '\n');
        ctrl = prev;
      }
      else if(ctrl == ctrlRsid)
        rsids.push_back(resid(iput, types));
      else if(ctrl == ctrlAtom) 
        parts.push_back(partk(iput, this));
      else if(ctrl == ctrlBond && BONDFORC > 0)
        bonds.push_back(bound(iput, this));
      else if(ctrl == ctrlAngl && BONDFORC > 1)
        angls.push_back(angle(iput, this));
      else if(ctrl == ctrlDhed && BONDFORC > 2)
        dheds.push_back(dihed(iput, this));
      else if(ctrl == ctrlStop) {
        if(prev == ctrlAtom)
          for(int i = 0; i < parts.size(); ++i)
            (*(const_cast<resid*> ((parts.at(i)).obRsid()))).
            adPart(&(parts.at(i)));
        else if(prev == ctrlBond && BONDFORC > 0)
          for(int i = 0; i < bonds.size(); ++i) {
            bound *bond = &(bonds.at(i));
            (*(const_cast<resid*> ((*bond).obRsid()))).adBond(bond);
            (*(const_cast<partk*> ((*bond).obLeft()))).adBond(bond);
            (*(const_cast<partk*> ((*bond).obRigh()))).adBond(bond);
          }
        else if(prev == ctrlAngl && BONDFORC > 1)
          for(int i = 0; i < angls.size(); ++i) {
            angle *angl = &(angls.at(i));
            (*(const_cast<resid*> ((*angl).obRsid()))).adAngl(angl);
            (*(const_cast<partk*> ((*angl).obLeft()))).adAngl(angl);
            (*(const_cast<partk*> ((*angl).obMidl()))).adAngl(angl);
            (*(const_cast<partk*> ((*angl).obRigh()))).adAngl(angl);
          }
        else if(prev == ctrlDhed && BONDFORC > 2)
          for(int i = 0; i < dheds.size(); ++i) {
            dihed *dhed = &(dheds.at(i));
            (*(const_cast<resid*> ((*dhed).obRsid()))).adDhed(dhed);
            (*(const_cast<partk*> ((*dhed).obUppr()))).adDhed(dhed);
            (*(const_cast<partk*> ((*dhed).obLeft()))).adDhed(dhed);
            (*(const_cast<partk*> ((*dhed).obRigh()))).adDhed(dhed);
            (*(const_cast<partk*> ((*dhed).obLowr()))).adDhed(dhed);
          }
      }
      else if(ctrl == ctrlEnd)
        break;
      prev = ctrl;
    }
  }

  //BUG: check this line-by-line for sync
  //initialize simulation conditions (with output prompt)
  void world::mkCond(ostream& oput, istream& iput) {
    myCoor data;
    oput << "Choose boundary control: 1 for periodic, 2 for rigd, "
         << "3 for repulsive" << endl;
    iput >> boundCon;
    oput << boundCon << endl;
    assert(boundCon >= 1 && boundCon < 4); 
    oput << "Choose the box vector" << endl;
    for(int k = 0; k < DIM; ++k){
      iput >> data;
      oput << data << endl;
      assert(data > 0);
      myVectSet(edge, k, data);
    }
    oput << "Choose ensemble control: 1 for NVE, 2 for NVT, 3 for NPT" 
         << endl;
    iput >> ensemCon;
    oput << ensemCon << endl;
    assert(ensemCon >= 1 && ensemCon < 4);
    oput << "Choose Coulomb force control: 1 for on" << endl;
    iput >> colmbCon;
    oput << colmbCon << endl;
    assert(colmbCon >= 0 && colmbCon < 2);
    if(colmbCon != 0)
      forceCon |= 0b0000000000000001;
    oput << "Choose LJ force control: 1 for shifted potential, "
         << "2 for smooth potential, 3 for truncated force, "
         << "4 for shifted force, 5 for smooth force" << endl;
    iput >> ljPotCon;
    oput << ljPotCon << endl;
    assert(LJONPART || ljPotCon == 0);
    assert(ljPotCon >= 0 && ljPotCon < 6);
    if(ljPotCon != 0) {
      forceCon |= 0b0000000000000001;
      if(ljPotCon != 1) {
        oput << "Choose cutoff distance" << endl;
        iput >> ofLj;
        oput << ofLj << endl;
        assert(ofLj > 0);
        if(ljPotCon == 3 || ljPotCon == 6) {
          oput << "Choose cuton distance" << endl;
          iput >> onLj;
          oput << onLj << endl;
          assert(onLj > 0);
        }
      }
    }
    oput << "Choose dissipative force control: 1 for uncorrelated, "
         << "2 for parallel, 3 for orthogonal, 4 for isotropic, "
         << "5 for parallel MZ, 6 for orthogonal MZ, "
         << "7 for full MZ, 8 for parallel auxiliary MZ, "
         << "9 for orthogonal auxiliary MZ, 10 for full auxiliary MZ"
         << endl;
    iput >> disspCon;
    oput << disspCon << endl;
    assert(DPDPARTK || disspCon == 0);
    assert(disspCon >= 0 && disspCon < 11);
    if(disspCon != 0) {
      forceCon |= 0b0000000000001110;
      if(disspCon != 1) {
        oput << "Choose cutoff distance for the force correlation"
             << endl;
        iput >> ofDp;
        oput << ofDp << endl;
        assert(ofDp > 0);
      }
      if(disspCon >= 5 && disspCon < 8) {
        if(disspCon != 6) {
          oput << "Choose cutoff size for the parallel memory kernel" << endl;
          iput >> pmOf;
          oput << pmOf << endl;
          assert(pmOf > 0);
        }
        if(disspCon != 5) {
          oput << "Choose cutoff size for the orthogonal memory kernel" << endl;
          iput >> omOf;
          oput << omOf << endl;
          assert(omOf > 0);
        }
      }
      else if(disspCon >= 8 && disspCon < 11) {
        oput << "WARNING: Auxiliary MZ requires thermalized system"
             << endl;
         if(disspCon != 9) {
          oput << "Choose number of auxiliary variables in parallel " 
               << "direction" << endl;
          iput >> pmOf;
          oput << pmOf << endl;
          assert(pmOf > 0);
        }
        if(disspCon != 8) {
          oput << "Choose number of auxiliary variables in orthogonal "
               << "direction" << endl;
          iput >> omOf;
          oput << omOf << endl;
          assert(omOf > 0);
        }
      }
    }
    oput << "Choose Verlet list control: 1 for constant update rate, "
         << "2 for automatic rate, 3 for using with cell list, "
         << "4 for using automatic rate AND cell list" << endl;
    iput >> verltCon;
    oput << verltCon << endl;
    assert(verltCon >= 0 && verltCon < 5);
    if(verltCon != 0) {
      oput << "Choose cutoff distance for Verlet list" << endl;
      iput >> ofVl;
      oput << ofVl << endl;
      assert(ofVl > 0);
      oput << "Choose update rate (for automatic update, it is the ratio "
           << "between greatest displacement and Verlet 'shell' thickness)"
           << endl;
      iput >> rtVl;
      oput << rtVl << endl;
      assert(rtVl > 0);
    }
    oput << "Choose integrator control: 1 for Verlet, "
         << "2 for leapfrog, 21 for DPD-LF, 24 for impulse-LF,"
         << "3 for velocity verlet, 31 for DPD-VV, "
         << "32 for consistent DPD-VV, 33 for Lowe method, "
         << "34 for impulsive DPD-VV " 
         << endl;
    iput >> intgrCon;
    oput << intgrCon << endl;
    assert((intgrCon >= 1 && intgrCon < 4) ||
           (intgrCon >= 21 && intgrCon < 25) || 
           (intgrCon >= 31 && intgrCon < 35));
    if(intgrCon == 32) {
      assert(disspCon != 0);
      oput << "WARNING: Self-consistent DPD-VV requires thermalized system"
           << endl;
      oput << "Choose limiting ratio" << endl;
      iput >> data;
      assert(data > 0);
      rtUp = data;
    }
    if(intgrCon == 33) {
      oput << "Choose coupling rate with the heat bath "
           << "(the collision probability per unit of time)"
           << endl;
      iput >> rtUp;
      assert(rtUp > 0);
      rtUp *= STPSZE;
    }
    if(intgrCon == 24 || intgrCon == 33 || intgrCon == 34)
      if(disspCon != 0)
        forceCon &= 0b1111111111111001;
    oput << "Choose temperature" << endl;
    iput >> temp;
    oput << temp << endl;
    assert(temp > 0);
    invT = 1 / BOLTZM / temp;
    sqrT = sqrt(BOLTZM * temp);
  }

  //initialize simulation conditions (from input stream)
  void world::iniCon(istream& iput) {
    myCoor dist;
    iput >> boundCon;
    for(int k = 0; k < DIM; ++k) {
      iput >> dist;
      myVectSet(edge, k, dist);
    }
    iput >> ensemCon;
    iput >> colmbCon;
    iput >> ljPotCon;
    if(LJONPART && (ljPotCon != 0 || ljPotCon != 1)) {
      iput >> ofLj;
      if(ljPotCon == 3 || ljPotCon == 6)
        iput >> onLj;
    }
    iput >> disspCon;
    if(DPDPARTK && (disspCon != 0 || disspCon != 1)) {
      iput >> ofDp;
      if(disspCon >= 5 && disspCon < 8) {
        if(disspCon != 6)
          iput >> pmOf;
        if(disspCon != 5)
          iput >> omOf;
      }
      else if(disspCon >= 8 && disspCon < 11) {
        if(disspCon != 9)
          iput >> pmOf;
        if(disspCon != 8)
          iput >> omOf;
      }
    }
    iput >> verltCon;
    if(verltCon != 0)
      iput >> ofVl >> rtVl;
    iput >> intgrCon;
    if(intgrCon == 32 || intgrCon == 33)
      iput >> rtUp;
    if(colmbCon != 0 || ljPotCon != 0)
      forceCon |= 0b0000000000000001;
    if(disspCon != 0)
      if(intgrCon != 24 && intgrCon != 33 && intgrCon != 34)
        forceCon |= 0b0000000000001110;
    iput >> temp;
    invT = 1 / BOLTZM / temp;
    sqrT = sqrt(BOLTZM * temp);
  }

  //initialize trajectory
  void world::iniTrj(istream& iput) {
    iput >> time;
    for(int i = 0; i < parts.size(); ++i)
      (parts.at(i)).iniTrj(iput);
    //create Verlet list for the first time
    for(int i = 0; i < parts.size(); ++i) {
      vLens.push_back(0);
      vPost.push_back(myVectMalloc(DIM));
    }
    this->apVcon();
    //BUG: is normal comparison a better choice here?
    if(BONDFORC > 0)
      for(int i = 0; i < bonds.size(); ++i)
        this->bondSz(&(bonds.at(i)), LEN1);
    if(BONDFORC > 1)
      for(int i = 0; i < angls.size(); ++i)
        this->anglSz(&(angls.at(i)), LEN1, LEN2);
    if(BONDFORC > 2)
      for(int i = 0; i < dheds.size(); ++i)
        this->dhedSz(&(dheds.at(i)), LEN1, LEN2, LEN3, NRM1, NRM2);
  }

  void world::skpTrj(istream& iput) {
    //BUG: if using ignore() for the line of time, it misses one line. why?
    iput >> time;
    for(int i = 0; i <= parts.size(); ++i)
      iput.ignore(256, '\n');
  }

  //add memory kernel
  void world::adMemo() {
    partk *atom1, *atom2;
    int memO;
    //memory for DPD-VV
    if(intgrCon == 31 || intgrCon == 32)
      for(int i = 0; i < parts.size(); ++i)
        (parts.at(i)).adAMem(1);
    //memory kernel for MZ
    if(disspCon != 0 && disspCon != 1) {
      myVectMemdel(post);
      myVectMemdel(velo);
      post = myVectMalloc(DIM);
      velo = myVectMalloc(DIM);
      if(disspCon >= 5 && disspCon < 8) {
        memO = fmax(pmOf, omOf);
        //clear memory
        pDCor.clear();
        pRCor.clear();
        oDCor.clear();
        oRCor.clear();
        xMemo.adPmem(0, &parts);
        vMemo.adPmem(0, &parts);
        uMemo.adPmem(0, &parts);
        pRMem.adPmem(0, &parts);
        oRMem.adPmem(0, &parts);
        //allocate memory
        xMemo.adPmem(memO, &parts);
        vMemo.adPmem(memO, &parts);
        uMemo.adPmem(memO, &parts);
        if(disspCon != 6) {
          pRMem.adPmem(memO, &parts);
          for(int i = 0; i < pmOf; ++i) {
            pDCor.push_back(0.0);
            pRCor.push_back(0.0);
          }
        }
        if(disspCon != 5) {
          oRMem.adPmem(omOf, &parts);
          for(int i = 0; i < omOf; ++i) {
            oDCor.push_back(0.0);
            oRCor.push_back(0.0);
          }
        }
      }
      //auxiliary variables for MZ
      else if(disspCon >= 8 && disspCon < 11) {
        //clear memory
        myVectMemdel(pAcc);
        myVectMemdel(oAcc);
        myVectMemdel(pFor);
        myVectMemdel(oFor);
        for(int i = 0; i < pAss.size(); ++i) {
          myMatrMemdel(pAss.at(i));
          myMatrMemdel(pBsr.at(i));
        }
        myMatrMemdel(oAss);
        myMatrMemdel(oBsr);
        pAxMm.clear();
        oAxMm.clear();
        //allocate memory
        post = myVectMalloc(DIM);
        velo = myVectMalloc(DIM);
        if(disspCon != 9) {
          pAcc = myVectMalloc(pmOf);
          pFor = myVectMalloc(pmOf);
          for(int i = 0; i < 100; ++i) {
            pAss.push_back(myMatrMalloc(pmOf, pmOf));
            pBsr.push_back(myMatrMalloc(pmOf, pmOf));
          }
          for(int i = 0; i < (parts.size() * (parts.size() - 1) / 2); ++i)
            pAxMm.push_back(auxMm(pmOf));
        }
        if(disspCon != 8) {
          oAcc = myVectMalloc(omOf);
          oFor = myVectMalloc(omOf);
          oAss = myMatrMalloc(omOf, omOf);
          oBsr = myMatrMalloc(omOf, omOf);
          for(int i = 0; i < (parts.size() * (parts.size() - 1) / 2); ++i)
            oAxMm.push_back(auxMm(omOf));
        }
      }
    }
  }

  //initialize memory kernel
  void world::iniMem(istream& iput) {
    assert(disspCon >= 5 && disspCon < 11);
    myCoor dumy, memO;
    partk *atom1, *atom2;
    //BUG: need support for case where the trajectory is restarted
    if(disspCon >= 5 && disspCon < 8) {
      memO = fmax(pmOf, omOf);
      for(int k = 0; k < DIM; ++k) {
        myVectSet(POST, k, ofDp);
        myVectZro(VELO);
      }
      for(int i = 0; i < parts.size(); ++i) {
        atom1 = &(parts.at(i));
        for(int j = i + 1; j < parts.size(); ++j) {
          atom2 = &(parts.at(j));
          obXdif(atom2, atom1, UNIT);
          for(int k = 0; k < memO; ++k) {
            xMemo.mkPmem(UNIT, &parts, atom1, atom2, k);
            this->apBcon(UNIT);
            myVectUnt(UNIT, dumy);
            uMemo.mkPmem(UNIT, &parts, atom1, atom2, k);
            vMemo.mkPmem(VELO, &parts, atom1, atom2, k);
            if(disspCon != 6)
              pRMem.mkPmem(VELO, &parts, atom1, atom2, k);
            if(disspCon != 5)
              oRMem.mkPmem(VELO, &parts, atom1, atom2, k);
          }
        }
      }
      if(disspCon != 6) {
        iput >> dumy;
        for(int k = 0; k < pmOf; ++k) {
          pRCor.at(k) = pow(dumy, k);
          pDCor.at(k) = 0;
          pDCor.at(k) = -pow(dumy, k) / (1 - dumy * dumy);
        }
        //for(int k = 0; k < pmOf; ++k)
        // for(int l = 0; l < (pmOf - k); ++l)
        //    pDCor.at(k) += pRCor.at(k + l) * pRCor.at(l);
        //for(int k = 0; k < pmOf; ++k) {
        //  iput >> dumy;
        //  pDCor.at(k) = dumy;
        //  iput >> dumy;
        //  pRCor.at(k) = dumy;
        //} 
      }
      if(disspCon != 5) {
        iput >> dumy;
        for(int k = 0; k < omOf; ++k) {
          oRCor.at(k) = pow(dumy, k);
          oDCor.at(k) = 0;
          oDCor.at(k) = -pow(dumy, k) / (1 - dumy * dumy);
        }
        //for(int k = 0; k < omOf; ++k)
        //  for(int l = 0; l < (omOf - k); ++l)
        //    oDCor.at(k) = oRCor.at(k + l) * oRCor.at(l);
        //for(int k = 0; k < omOf; ++k) {
        //  iput >> dumy;
        //  oDCor.at(k) = dumy;
        //  iput >> dumy;
        //  oRCor.at(k) = dumy;
        //} 
      }
      for(int k = 0; k < (!(pmOf * omOf) ? fmax(pmOf, omOf)
                                         : min(pmOf, omOf)) ; ++k) {
        if(disspCon != 6)
          COUTPT << pRCor.at(k) << " " << pDCor.at(k) << " ";
        if(disspCon != 5)
          COUTPT << oRCor.at(k) << " " << oDCor.at(k) << " ";
        COUTPT << endl;
      }
    }
    else if(disspCon >= 8 && disspCon < 11) {
      if(disspCon != 9) {
        for(int i = 0; i < 100; ++i) {
          COUTPT << endl << i * ofDp / 100 << endl;
          for(int j = 0; j < pmOf; ++j) {
            COUTPT << endl << "pAss= ";
            myVectInp(pAcc, iput);
            myVectOut(pAcc);
            myMatrSetrow(pAss.at(i), j, pAcc);
          }
          for(int j = 0; j < pmOf; ++j) {
            COUTPT << endl << "pBsr= ";
            myVectInp(pAcc, iput);
            myVectOut(pAcc);
            myMatrSetrow(pBsr.at(i), j, pAcc);
          }
        }
      }
      if(disspCon != 8) {
        COUTPT << endl << "oAss= ";
        for(int i = 0; i < omOf; ++i) {
          myVectInp(oAcc, iput);
          myVectOut(oAcc);
          myMatrSetrow(oAss, i, oAcc);
        }
        COUTPT << endl << "oBsr= ";
        for(int i = 0; i < omOf; ++i) {
          myVectInp(oAcc, iput);
          myVectOut(oAcc);
          myMatrSetrow(oBsr, i, oAcc);
        }
      }
    }
  }

  //initialize trajectory randomly
  void world::iniRnd() {
    cerr << "Not supporting random initialization" << endl;
    terminate();
  }
 
  //initialize coordinate randomly
  inline void world::rndCrd() {
    myCoor boxV = 1, rsdV;
    resid *mole;
    mole = &(rsids.at(0));
    for(int i = 1; i < rsids.size(); ++i)
      assert(&(rsids.at(i)) == mole);
    for(int k = 0; k < DIM; ++k)
      boxV *= myVectGet(edge, k);
    boxV /= parts.size();
    //constructing the relaxed residue
    cerr << "Not supporting randomization of position" << endl;
    terminate();
  }

  //initialize velocity randomly
  inline void world::rndVel() {
    myCoor tmpMs = 0, rVel;
    resid *mole;
    mole = &(rsids.at(0));
    for(int i = 1; i < rsids.size(); ++i)
      assert(&(rsids.at(i)) == mole);
    assert(temp > 0);
    for(int i = 0; i < (*mole).obNpartk(); ++i)
      tmpMs += (*(*mole).obPart(i)).obMass();
    tmpMs *= invT;
    //assign velocity s.t. the COM velocity i 0
    cerr << "Not supporting randomization of velocity" << endl;
    terminate();
  } 
 
  //update trajectory
  void world::nwTraj() {
    //NOTE: be careful of using dummy vectors. Some of them may not be
    //dummy in the super-routine. Need to check all dummy cases.
    switch(intgrCon) {
      case 1:
        this->nwMemo();
        this->nwVlst();
        this->nwAccl();
        this->verlet();
        break;
      case 2:
        this->nwMemo();
        this->nwVlst();
        this->nwAccl();
        this->lpFrgA();
        this->lpFrgB();
        break;
      case 21:
        cerr << "Not supporting DPD-LF" << endl;
        terminate();
        break;
      case 24:
        this->nwMemo();
        this->nwVlst();
        this->nwAccl();
        this->lpFrgA();
        this->impuls();
        this->lpFrgB();
        break; 
      case 3:
        this->vVerlA();
        this->nwMemo();
        this->nwVlst();
        this->nwAccl();
        this->vVerlB();
        break;
      case 31:
        this->vVerlA();
        this->nwMemo();
        this->nwVlst();
        forceCon = (colmbCon != 0 || ljPotCon != 0) ? 0b0000000000000001
                                                    : 0b0000000000000000;
        forceCon = (disspCon != 0) ? (forceCon | 0b0000000000001100)
                                   : forceCon;
        this->nwAccl();
        this->apAMem();
        forceCon &= 0b1111111111111010;
        forceCon |= 0b0000000000000010;
        this->nwAccl();
        this->apAMem();
        this->vVerlB();
        this->apVMem();
        this->nwAccl();
        this->apAMem();
        break;
      case 32:
        this->vVerlA();
        this->nwMemo();
        this->nwVlst();
        forceCon = (colmbCon != 0 || ljPotCon != 0) ? 0b0000000000000001
                                                    : 0b0000000000000000;
        forceCon = (disspCon != 0) ? (forceCon | 0b0000000000001100)
                                   : forceCon;
        this->nwAccl();
        this->apAMem();
        forceCon &= 0b1111111111111010;
        forceCon |= 0b0000000000000010;
        do {
          forceCon += 0b0000010000000000;
          this->nwAccl();
          this->apAMem();
          this->vVerlB();
          this->apVMem();
          //if(forceCon & 0b1000000000000000) {
          //  cerr << "Fail to thermalize at T=" << this->obSysT()
          //       << ". Terminating..." << endl;
          //  terminate();
          //}
        //} while(forceCon < (0b0000010000000000 * rtUp )); 
        } while(abs((this->obSysT() / temp) - 1) > rtUp);
        this->nwAccl();
        this->apAMem();
        break;
      case 33:
      case 34:
        this->nwMemo();
        this->vVerlA();
        this->nwVlst();
        this->nwAccl();
        this->vVerlB();
        this->impuls();
        break;
      }
    time += STPSZE;
  }

  //update memory kernel
  inline void world::nwMemo() {
    if(!(disspCon >= 5 && disspCon < 11)) 
      return;
    partk *atom1, *atom2;
    myCoor dist;
    //BUG: Verlet list does not support MZ
    //BUG: distance calculation is repeated later. Maybe it's better to have
    //     distance memory
    if(disspCon >= 5 && disspCon < 8) {
      for(int i = 0; i < parts.size(); ++i) {
        atom1 = &(parts.at(i));
        for(int j = i + 1; j < parts.size(); ++j) {
          atom2 = &(parts.at(j));
          obXdif(atom2, atom1, UNIT);
          xMemo.nwPmem(UNIT, &parts, atom1, atom2);
          this->apBcon(UNIT);
          myVectUnt(UNIT, dist);
          uMemo.nwPmem(UNIT, &parts, atom1, atom2);
          obVdif(atom2, atom1, VELO);
          vMemo.nwPmem(VELO, &parts, atom1, atom2);
          if(disspCon != 6) {
            myVectCpy(POST, UNIT);
            myVectScl(POST, NRMDIS(ENGINE));
            pRMem.nwPmem(POST, &parts, atom1, atom2);
          }
          if(disspCon != 5) {
            myVectOrd(POST, UNIT);
            oRMem.nwPmem(POST, &parts, atom1, atom2);
          }
        }
      }
      xMemo.nwMpos();
      uMemo.nwMpos();
      vMemo.nwMpos();
      if(disspCon != 6)
        pRMem.nwMpos();
      if(disspCon != 5)
        oRMem.nwMpos();
    }
  }

  //update Verlet list
  inline void world::nwVlst() {
    if(verltCon == 0)
      return;
    myCoor dist, maxD, ctOf = fmax(ofDp, ofLj);
    bool updt = false;
    //check if the list needs update
    switch(verltCon) {
      case 1:
        updt = fmod(time, STPSZE * rtVl) < STPSZE;
        break;
      case 2:
        maxD = 0;
        for(int i = 0; i < parts.size(); ++i) {
          myVectCpy(UNIT, vPost.at(i));
          myVectSub(UNIT, (parts.at(i)).obPost());
          dist = myVectNrm(UNIT);
          if(dist > maxD)
            maxD = dist;
        }
        updt = maxD > (rtVl * (ofVl - ctOf));
        break;
    }
    //update the list
    if(updt)
      this->apVcon();
  }

  //update acceleration (using normal method)
  inline void world::nwAccl() {
    partk *atom1, *atom2, *atom3, *atom4;
    bound* bond;
    angle* angl;
    dihed* dhed;
    int cPtr = 0, cMax, cMin;
    //wipe acceleration
    myVectZro(FORC);
    for(int i = 0; i < parts.size(); ++i)
      (parts.at(i)).mkAccl(FORC);
    if(disspCon == 8 || disspCon == 9 || disspCon == 10) {
      if(disspCon != 9) {
        myVectZro(pAcc);
        for(int i = 0; i < (parts.size() * (parts.size() - 1) / 2); ++i)
          (pAxMm.at(i)).mkVelo(pAcc);
      }
      if(disspCon != 8) {
        myVectZro(oAcc);
        for(int i = 0; i < (parts.size() * (parts.size() - 1) / 2); ++i)
          (oAxMm.at(i)).mkVelo(oAcc);
      }
    }
    //nonbonded interactions
    //BUG: If there are only bonded interactions, this will fail
    //if(!(forceCon & 0b0000000000000111))
    //  return;
    for(int i = 0; i < parts.size(); ++i) {
      atom1 = &(parts.at(i));
      //uncorrelated disspative force
      if(disspCon == 1) {
        myVectZro(FTOT);
        this->dispFc(atom1, NULL, 0.0, NULL, FORC);
        (*atom1).adAccl(FORC);
        this->randFc(atom1, NULL, 0.0, NULL, FORC);
        (*atom1).adAccl(FORC);
      }
      //inter-particle force
      if((ljPotCon != 0  && ljPotCon != 1) ||
         (disspCon != 0  && disspCon != 1  &&
          intgrCon != 24 && intgrCon != 33 && intgrCon != 34)) {
        //get the next atom from Verlet list or in the sequence
        cMax = (verltCon == 0) ? parts.size() : vLens.at(i);
        cMin = (verltCon == 0) ? (i + 1) : 0;
        for(int j = cMin; j < cMax; ++j) {
          if(verltCon != 0) {
            atom2 = vList.at(cPtr);
            ++cPtr;
          }
          else {
            atom2 = &(parts.at(j));
            if(samTop(atom1, atom2))
              continue;
          }
          this->nbndFc(atom1, atom2, 0, FORC);
          myVectInv(FORC, (*atom1).obMass());
          (*atom1).adAccl(FORC);
          myVectScl(FORC, -((*atom1).obMass() / (*atom2).obMass()));
          (*atom2).adAccl(FORC);
        }
      }
    }
    //bonded interactions
    //BUG: is normal comparison a better choice here?
    if(BONDFORC > 0)
      for(int i = 0; i < bonds.size(); ++i) {
        bond = &(bonds.at(i));
        atom1 = const_cast<partk*> ((*bond).obLeft());
        atom2 = const_cast<partk*> ((*bond).obRigh());
        this->bondFc(bond, 1, FORC);
        myVectInv(FORC, (*atom1).obMass());
        (*atom1).adAccl(FORC);
        myVectScl(FORC, -((*atom1).obMass() / (*atom2).obMass()));
        (*atom2).adAccl(FORC);
      } 
    if(BONDFORC > 1)
      for(int i = 0; i < angls.size(); ++i) {
        angl = &(angls.at(i)); 
        atom1 = const_cast<partk*> ((*angl).obLeft());
        atom2 = const_cast<partk*> ((*angl).obMidl());
        atom3 = const_cast<partk*> ((*angl).obRigh());
        this->anglFc(angl, 1, FOR1, FOR2);
        myVectInv(FOR1, (*atom1).obMass());
        (*atom1).adAccl(FOR1);
        myVectScl(FOR1, -((*atom1).obMass() / (*atom2).obMass()));
        (*atom2).adAccl(FOR1);
        myVectInv(FOR2, (*atom3).obMass());
        (*atom3).adAccl(FOR2);
        myVectScl(FOR2, -((*atom3).obMass() / (*atom2).obMass()));
        (*atom2).adAccl(FOR2);
      }
    if(BONDFORC > 2)
      for(int i = 0; i < dheds.size(); ++i) {
        cerr << "Not supporting dihedral angle force" << endl;
        terminate();
    }
  }

  //update acceleration (using impulse method)
  inline void world::impuls() {
    if(disspCon == 0)
      return;
    myCoor dist, mass, sdot, pVel, oVel,
           pGma = 0, oGma = 0, pSig = 0, oSig = 0;
    int cPtr = 0, cMax, cMin, pos, rPtr;
    partk *atom1, *atom2;
    bool updt = true;
    for(int i = 0; i < parts.size(); ++i) {
        atom1 = &(parts.at(i));
      switch(disspCon) {
        case 1:
          //impulse frequency
          switch(intgrCon) {
            case 24:
            case 34:
              pGma = exp(-0.5 * (*atom1).obSigm() * (*atom1).obSigm()
                              * invT * STPSZE);
              break;
            case 33:
              pGma = 0;
              break;
          }
          pSig = sqrt(BOLTZM * temp * (1 - pGma * pGma) / (*atom1).obMass());
          pGma -= 1;
          //impulse
          myVectCpy(FORC, (*atom1).obVelo());
          myVectScl(FORC, pGma);
          (*atom1).adVelo(FORC);
          for(int k = 0; k < DIM; ++k)
            myVectSet(FORC, k, pSig * NRMDIS(ENGINE)); 
          (*atom1).adVelo(FORC);
          break;
        case 2:
        case 3:
        case 4:
          cMax = (verltCon == 0) ? parts.size() : vLens.at(i);
          cMin = (verltCon == 0) ? (i + 1) : 0;
          for(int j = cMin; j < cMax; ++j) {
            if(verltCon != 0) {
              atom2 = vList.at(cPtr);
              ++cPtr;
            }
            else {
              atom2 = &(parts.at(j));
              if(samTop(atom1, atom2))
                continue;
            }
            //unit vector pointing from 1 to 2
            obXdif(atom2, atom1, UNIT);
            this->apBcon(UNIT);
            myVectUnt(UNIT, dist);
            switch(intgrCon) {
              case 24:
              case 34:
                updt = dist < ofDp && 
                       UNIDIS(ENGINE) <= this->pDpdSt(atom1, atom2, dist);
                break;
              case 33:
                updt = dist < ofDp && UNIDIS(ENGINE) <= rtUp;
                break;
            }
            if(updt) {
              //parallel projection
              obVdif(atom2, atom1, FORC);
              myVectDot(FORC, UNIT, &pVel);
              //impulse frequency and reduced mass 
              mass = (*atom1).obMass() * (*atom2).obMass()
                     / ((*atom1).obMass() + (*atom2).obMass()); 
              switch(intgrCon) {
                case 24:
                case 34:
                  pGma = exp(-0.333 * (*atom1).obSigm() * (*atom2).obSigm()
                                    * invT * STPSZE);
                  break;
                case 33:
                  pGma = 0;
                  break;
              }
              pSig = sqrt(BOLTZM * temp * (1 - pGma * pGma) / mass);
              pGma -= 1;
              //impulse
              if(disspCon != 2) {
                cerr << "Not supporting orthogonal damping" << endl;
                terminate();
              }
              if(disspCon != 3) {
                pVel *= pGma;
                pSig *= NRMDIS(ENGINE);
                myVectScl(UNIT, (pSig + pVel));
                myVectScl(UNIT, -mass / (*atom2).obMass());
                (*atom2).adVelo(UNIT);
                myVectScl(UNIT, -(*atom2).obMass() / (*atom1).obMass());
                (*atom1).adVelo(UNIT);
              }
            }
          }
          break;
        case 5:
        case 6:
        case 7:
          cerr << "Not supporting MZ" << endl;
          terminate();
          break;
        case 8:
        case 9:
        case 10:
          cMax = (verltCon == 0) ? parts.size() : vLens.at(i);
          cMin = (verltCon == 0) ? (i + 1) : 0;
          for(int j = cMin; j < cMax; ++j) {
            if(verltCon != 0) {
              atom2 = vList.at(cPtr);
              ++cPtr;
            }
            else
              atom2 = &(parts.at(j));
            pos = obPInd(&parts, atom1, atom2);
            //BUG: may need velocity reduction factor
            //impulse correlation and reduced mass
            mass = (*atom1).obMass() * (*atom2).obMass()
                   / ((*atom1).obMass() + (*atom2).obMass());
            pGma = -this->pDpdSt(atom1, atom2, dist);
            oGma = -this->oDpdSt(atom1, atom2, dist);
            //unit vector pointing from 1 to 2
            ljonSz(atom1, atom2, dist, UNIT);
            if(dist < ofDp) {
              rPtr = floor(dist / ofDp * 100);
              //parallel projection
              obVdif(atom2, atom1, VELO);
              myVectDot(VELO, UNIT, &pVel);
              if(disspCon != 8) {
                cerr << "Not supporting orthogonal auxiliary MZ" << endl;
                terminate();
              }
              if(disspCon != 9) {
                //S(0) = v_ij
                myVectCpy(pFor, (pAxMm.at(pos)).obPost());
                myVectSet(pFor, 0, pVel);
                //S = A_sp*S + B_sr*xi
                myMatrVecmul(CblasNoTrans, 1.0, pAss.at(rPtr), pFor, 0.0, pAcc);
                pVel -= myVectGet(pAcc, 0);
                for(int k = 0; k < pmOf; ++k)
                  myVectSet(pFor, k, NRMDIS(ENGINE));
                myMatrVecmul(CblasNoTrans, sqrT, pBsr.at(rPtr), pFor, 
                                           1.0, pAcc);
                (pAxMm.at(pos)).mkPost(pAcc);
                //impulse
                myVectCpy(VELO, UNIT);
                myVectScl(VELO, -0.5 * myVectGet((pAxMm.at(pos)).obPost(), 0));
                myVectScl(UNIT, -mass / (*atom2).obMass());
                (*atom2).adVelo(UNIT);
                myVectScl(UNIT, -(*atom2).obMass() / (*atom1).obMass());
                (*atom1).adVelo(UNIT);
              }
            }
          }
          break;
      }
    }
  }

  //apply boundary condition
  //POSTCON: uVec is a vector satisfying the boundary condition
  inline void world::apBcon(myVect* uVec) {
    switch(boundCon) {
      case 1:
        myVectCpy(PBC1, uVec);
        myVectDiv(PBC1, edge);
        for(int k = 0; k < DIM; ++k)
          myVectSet(PBC1, k, myVectGet(edge, k) * round(myVectGet(PBC1, k)));
        myVectSub(uVec, PBC1);
    }
  }

  //apply Verlet list condition
  inline void world::apVcon() {
    partk *atom1, *atom2;
    int cPtr;
    myCoor dist;
    vList.clear();
    for(int i = 0; i < parts.size(); ++i) {
      atom1 = &(parts.at(i));
      cPtr = 0;
      if(verltCon == 2)
        myVectCpy(vPost.at(i), (*atom1).obPost());
      for(int j = i + 1; j < parts.size(); ++j) {
        atom2 = &(parts.at(j));
        if(!samTop(atom1, atom2)) {
          obXdif(atom2, atom1, UNIT);
          this->apBcon(UNIT);
          dist = myVectNrm(UNIT);
          if(dist < ofVl) {
            vList.push_back(atom2);
            ++cPtr;
          }
        }
      }
      vLens.at(i) = cPtr;
    }
  }

  //apply acceleration memory for DPD-VV
  inline void world::apAMem() {
    if(forceCon & 0b0000000000000101)
      for(int i = 0; i < parts.size(); ++i)
        (parts.at(i)).nwAMem((parts.at(i)).obAccl());
    else if(forceCon & 0b0000000000000010)
      for(int i = 0; i < parts.size(); ++i)
        (parts.at(i)).adAccl((parts.at(i)).obAMem(0));
    else {
      myVectZro(FTOT);
      for(int i = 0; i < parts.size(); ++i)
        (parts.at(i)).nwAMem(FTOT);
    }
  }

  //apply velocity memory for DPD-VV
  inline void world::apVMem() {
    if(disspCon >= 5 && disspCon < 8) 
      for(int i = 0; i < parts.size(); ++i)
        for(int j = i + 1; j < parts.size(); ++j) {
          obVdif(&(parts.at(j)), &(parts.at(i)), VELO);
          vMemo.mkPmem(VELO, &parts, &(parts.at(i)), &(parts.at(j)), 0);
        }
  }

  //calculate the distance between two particles
  //POSTCON: len1 is unit vector pointing from atom1 to atom2
  void world::ljonSz(const partk* atom1, const partk* atom2,
                     myCoor& dist, myVect* len1) {
    obXdif(atom2, atom1, len1);
    this->apBcon(len1);
    myVectUnt(len1, dist);
  }

  //update bound size
  //BUG: this calculation is also done in nonbonded force, so this is
  //redundant, but is kept for now to make things clear
  //POSTCON: len1 is unit vector pointing from left to righ
  inline void world::bondSz(bound* bond, myVect* len1) {
    myCoor leng;
    obXdif((*bond).obLeft(), (*bond).obRigh(), len1);
    leng = myVectNrm(len1);
    myVectInv(len1, leng);
    (*bond).mkSize(leng);
  }

  //update angle size
  //POSTCON: len1, len2 are unit vector pointing from midl to left, righ
  inline void world::anglSz(angle* angl, myVect* len1, myVect* len2) {
    myCoor sdot;
    obXdif((*angl).obMidl(), (*angl).obLeft(), len1);
    obXdif((*angl).obMidl(), (*angl).obRigh(), len2);
    myVectInv(len1, (*((*angl).obBlef())).obSize());
    myVectInv(len2, (*((*angl).obBrig())).obSize());
    myVectDot(len1, len2, &sdot);
    if(abs(sdot) > (1 - MINZRO))
      (*angl).mkSize(acos(copysign(1.0, sdot)));
    else
      (*angl).mkSize(acos(sdot));
  }    

  //update dihed size
  //POSTCON: len1, len2, len3 are unit vector pointing from uppr to left,
  //      left to righ, righ to lowr
  //      nrm1, nrm2 are normal unit vector for the uppr, lowr planes
  inline void world::dhedSz(dihed* dhed, myVect* len1, myVect* len2,
                            myVect* len3, myVect* nrm1, myVect* nrm2) {
    myCoor sdot;
    obXdif((*dhed).obUppr(), (*dhed).obLeft(), len1);
    obXdif((*dhed).obLeft(), (*dhed).obRigh(), len2);
    obXdif((*dhed).obRigh(), (*dhed).obLowr(), len3);
    myVectInv(len1, (*((*dhed).obBlef())).obSize());
    myVectInv(len2, (*((*dhed).obBmid())).obSize());
    myVectInv(len3, (*((*dhed).obBrig())).obSize());
    myVectVecCrs(len2, len1, nrm1);
    myVectVecCrs(len2, len3, nrm2);
    myVectInv(nrm1, sin((*((*dhed).obAlef())).obSize()));
    myVectInv(nrm2, sin((*((*dhed).obArig())).obSize()));
    myVectDot(nrm1, nrm2, &sdot);
    (*dhed).mkSize(acos(fmin(fmax(sdot, -1.0), 1.0)));
  }

  //update non-bonded force
  void world::nbndFc(const partk* atom1, const partk* atom2,
                     const int ctrl, myVect* fVec) {
    myCoor dist;
    myVectZro(fVec);
    //unit vector pointing from 1 to 2
    if(disspCon != 5 && disspCon != 6 && disspCon != 7)
      this->ljonSz(atom1, atom2, dist, UNIT);
    else {
      myVectCpy(UNIT, xMemo.obPmem(&parts, atom1, atom2, 0));
      this->apBcon(UNIT);
      dist = myVectNrm(UNIT);
      myVectCpy(UNIT, uMemo.obPmem(&parts, atom1, atom2, 0));
    }
    //conservative force
    if(ljPotCon != 0 && (ljPotCon == 1 || dist < ofLj)) {
      this->ljonFc(atom1, atom2, ctrl, dist, UNIT, LEN1);
      myVectAdd(fVec, LEN1);
    }
    //DPD force
    //BUG: MZ uses memory kernel, which means some pairs
    //may not be included because of Verlet list
    if(disspCon != 0 && disspCon != 1) {
    //if(disspCon != 0 && disspCon != 1 && dist < ofDp) {
      this->dispFc(atom1, atom2, dist, UNIT, LEN1);
      myVectAdd(fVec, LEN1);
      this->randFc(atom1, atom2, dist, UNIT, LEN1);
      myVectAdd(fVec, LEN1);
    }
  }

  //update bond force
  void world::bondFc(bound* bond, const int ctrl, myVect* fVec) {
    this->bondSz(bond, fVec);
    myVectScl(fVec, this->bondSt(bond, ctrl));
  }

  //update angle force
  //PRECON: size of bLef, bRig is updated
  //POSTCON: for1, for2 are the forces acting on left and righ
  void world::anglFc(angle* angl, const int ctrl, myVect* for1,
                                                  myVect* for2) {
    myCoor aFor;
    this->anglSz(angl, LEN1, LEN2);
    aFor = this->anglSt(angl, ctrl);
    myVectVecCrs(LEN1, LEN2, NRM1);
    myVectVecCrs(LEN1, NRM1, for1);
    myVectScl(for1, aFor / myVectNrm(for1)
                         / (*((*angl).obBrig())).obSize());
    myVectVecCrs(NRM1, LEN2, for2);
    myVectScl(for2, aFor / myVectNrm(for2)
                         / (*((*angl).obBlef())).obSize());
  }

  //update dihed force
  //PRECON: size of bLef, bMid, bRig is updated
  //        size of aLef, aRig is updated
  //BUG: A mess
  void world::dhedFc(dihed* dhed, const int ctrl, 
                     myVect* for1, myVect* for2, myVect* for3) {
  }

  //update LJ force between 2 partk
  //POSTCON: fVec stores LJ force vector from atom2 to atom1
  void world::ljonFc(const partk* atom1, const partk* atom2,
                     const int ctrl, const myCoor dist,
                     const myVect* uVec, myVect* fVec) {
    myVectZro(fVec);
    if(!(forceCon & 0b0000000000000001))
      return;
    else {
      myVectCpy(fVec, uVec);
      myVectScl(fVec, this->ljonSt(atom1, atom2, ctrl, dist));
    }
  }

  //update dissipative force between 2 partk
  //POSTCON: fVec stores dissipative force vector from atom2 to atom1
  void world::dispFc(const partk* atom1, const partk* atom2,
                     const myCoor dist, const myVect* uVec,
                     myVect* fVec) {
    myVectZro(fVec);
    if(!(forceCon & 0b0000000000000010) || !(forceCon & 0b0000000000001000))
      return;
    else {
      myCoor pVel, oVel, sdot, pGma = 0, oGma = 0;
      int rPtr;
      int pos;
      switch(disspCon) {
        case 1:
          pGma = -0.5 * (*atom1).obSigm() * (*atom1).obSigm() * invT;
          myVectCpy(fVec, (*atom1).obVelo());
          myVectScl(fVec, pGma);
          break;
        case 2:
        case 3:
        case 4:
          //parallel projection
          myVectCpy(POST, uVec);
          obVdif(atom2, atom1, VELO);
          myVectDot(VELO, POST, &pVel);
          //orthogonal force
          if(disspCon != 2) {
            //orthogonal projection
            myVectScl(POST, pVel);
            myVectSub(VELO, POST);
            oVel = myVectNrm(VELO);
            oGma = -0.5 * (*atom1).obSigm() * (*atom2).obSigm() * invT *
                          this->oDpdSt(atom1, atom2, dist) *
                          this->oDpdSt(atom1, atom2, dist);
            oGma *= oVel;
            //orthogonal unit vector
            while(oVel < MINZRO) {
              myVectOrd(VELO, uVec);
              oVel = myVectNrm(VELO);
            }
            myVectInv(VELO, oVel);
            myVectScl(VELO, oGma);
            myVectAdd(fVec, VELO);
          }
          //parallel force
          if(disspCon != 3) {
            pGma = -0.5 * (*atom1).obSigm() * (*atom2).obSigm() * invT *
                          this->pDpdSt(atom1, atom2, dist) *
                          this->pDpdSt(atom1, atom2, dist);
            pGma *= pVel;
            myVectCpy(VELO, uVec);
            myVectScl(VELO, pGma);
            myVectAdd(fVec, VELO);
          }
          break;
        case 5:
        case 6:
        case 7:
          pGma = -this->pDpdSt(atom1, atom2, dist) *
                  this->pDpdSt(atom1, atom2, dist) * invT;
          oGma = -this->oDpdSt(atom1, atom2, dist) *
                  this->oDpdSt(atom1, atom2, dist) * invT;
          pos = fmax(pmOf, omOf);
          for(int k = 0; k < pos; ++k) {
            myVectCpy(VELO, vMemo.obPmem(&parts, atom1, atom2, k));
            myVectCpy(POST, uMemo.obPmem(&parts, atom1, atom2, k));
            //parallel projection
            myVectDot(VELO, POST, &pVel);
            //orthogonal force
            if(disspCon != 5 && k < omOf) {
              myVectScl(POST, pVel);
              myVectSub(VELO, POST);
              myVectScl(VELO, oGma * oDCor.at(k));
              myVectAdd(fVec, VELO);
            }
            //parallel force
            if(disspCon != 6 && k < pmOf) {
              myVectCpy(VELO, POST);
              myVectScl(VELO, pGma * pDCor.at(k));
              myVectAdd(fVec, VELO);
            }
          }
          break;
        case 8:
        case 9:
        case 10:
          if(dist < ofDp) {
            rPtr = floor(dist / ofDp * 100);
            pGma = -this->pDpdSt(atom1, atom2, dist);
            oGma = -this->oDpdSt(atom1, atom2, dist);
            pos = obPInd(&parts, atom1, atom2);
            //parallel projection
            obVdif(atom2, atom1, VELO);
            myVectCpy(POST, uVec);
            myVectDot(VELO, POST, &pVel);
            //orthogonal force
            if(disspCon != 8) {
              //orthogonal projection
               myVectScl(POST, pVel);
              myVectSub(VELO, POST);
              oVel = myVectNrm(VELO);
              //S(0) = v_ij
              myVectCpy(oFor, (oAxMm.at(pos)).obPost());
              myVectSet(oFor, 0, oVel);
              myVectCpy(oAcc, oFor);
              //dS = (A_ss-1)*S + B_sr*xi
              myMatrVecmul(CblasNoTrans, 1.0, oAss, oFor, -1.0, oAcc);
              for(int k = 0; k < omOf; ++k)
                myVectSet(oFor, k, NRMDIS(ENGINE));
              myMatrVecmul(CblasNoTrans, sqrT, oBsr, oFor, 1.0, oAcc);
              //orthogonal unit vector
              while(oVel < MINZRO) {
                myVectOrd(VELO, uVec);
                oVel = myVectNrm(VELO);
              } 
              myVectInv(VELO, oVel);
              (oAxMm.at(pos)).adPost(oAcc);
              myVectScl(VELO, myVectGet(oAcc, 0) * STPINV);
              myVectAdd(fVec, VELO);
            }
            //parallel force
            if(disspCon != 9) {
              //S(0) = v_ij
              myVectCpy(pFor, (pAxMm.at(pos)).obPost());
              myVectSet(pFor, 0, pVel);
              myVectCpy(pAcc, pFor);
              //dS = (A_ss-1)*S + B_sr*xi
              myMatrVecmul(CblasNoTrans, 1.0, pAss.at(rPtr), pFor, -1.0, pAcc);
              for(int k = 0; k < pmOf; ++k)
                myVectSet(pFor, k, NRMDIS(ENGINE));
              myMatrVecmul(CblasNoTrans, sqrT, pBsr.at(rPtr), pFor, 1.0, pAcc);
              (pAxMm.at(pos)).adPost(pAcc);
              //parallel unit vector
              myVectCpy(VELO, uVec);
              myVectScl(VELO, myVectGet(pAcc, 0) * STPINV);
              myVectAdd(fVec, VELO);
            }
          }
          break;
      }
    }
  }

  //update random force between 2 partk
  //POSTCON: fVec stores random force vector from atom2 to atom1
  void world::randFc(const partk* atom1, const partk* atom2,
                     const myCoor dist, const myVect* uVec,
                     myVect* fVec) {
    myVectZro(fVec);
    if(!(forceCon & 0b0000000000000100) || !(forceCon & 0b0000000000001000))
      return;
    else {
      myCoor pSig = 0, oSig = 0;
      int pos;
      switch(disspCon) {
        case 1:
          pSig = STPIRT * (*atom1).obSigm();
          pSig *= NRMDIS(ENGINE);
          switch(DIM) {
            case 1:
              pSig *= (ENGINE() % 2 == 0) ? 1.0 : -1.0;
              myVectSet(fVec, 0, pSig);
              break;
            case 2:
              break;
            case 3:
              for(int k = 0; k < DIM; ++k)
                myVectSet(fVec, k, pSig * NRMDIS(ENGINE));
              break;
          }
          break;
        case 2:
        case 3:
        case 4:
          if(disspCon != 3) {
            pSig = STPIRT * sqrt((*atom1).obSigm() * (*atom2).obSigm())
                          * this->pDpdSt(atom1, atom2, dist); 
            pSig *= NRMDIS(ENGINE);
            myVectCpy(POST, uVec);
            myVectScl(POST, pSig);
            myVectAdd(fVec, POST);
          }
          if(disspCon != 2) {
            oSig = STPIRT * sqrt((*atom1).obSigm() * (*atom2).obSigm())
                          * this->oDpdSt(atom1, atom2, dist);
            myVectOrd(POST, uVec);
            myVectScl(POST, oSig);
            myVectAdd(fVec, POST);
          }
          break;
        case 5:
        case 6:
        case 7:
          //parallel force
          if(disspCon != 6) {
            pSig = STPIRT * this->pDpdSt(atom1, atom2, dist);
            for(int k = 0; k < pmOf; ++k) {
              myVectCpy(POST, pRMem.obPmem(&parts, atom1, atom2, k));
              myVectScl(POST, pSig * pRCor.at(k));
              myVectAdd(fVec, POST);
            }
          }
          //orthogonal force
          if(disspCon != 5) {
            oSig = STPIRT * this->oDpdSt(atom1, atom2, dist);
            for(int k = 0; k < omOf; ++k) {
              myVectCpy(POST, oRMem.obPmem(&parts, atom1, atom2, k));
              myVectScl(POST, oSig * oRCor.at(k));
              myVectAdd(fVec, POST);
            }
          }
          break;
        case 8:
        case 9:
        case 10:
          break;
      }
    }
  }

  inline myCoor world::ljonSt(const partk* atom1, const partk* atom2,
                              const int ctrl, const myCoor dist) {
    myCoor rmin, epsl, ljDist, ctOfLj;
    rmin = ((*atom1).obRmin() + (*atom2).obRmin());
    epsl = sqrt((*atom1).obEpsl() * (*atom2).obEpsl());
    ljDist = pow(rmin / dist, 6);
    switch(ljPotCon) {
      case 2:
        cerr << "Not supporting shifted lj potential" << endl;
        terminate();
        break;
      case 3:
        cerr << "Not supporting smooth lj potential" << endl;
        terminate();
        break;
      case 4:
        switch(ctrl) {
          case 0:
            return (12 * epsl * (ljDist - 1) * ljDist / dist);
            break;
          case 1:
            return (15.7282 * dist * exp(-4.7395 * dist * dist));
            break;
          case 2:
            return -(568.791 * pow(dist, 5) - 3156.54 * pow(dist, 4) +
                     6680.22 * pow(dist, 3) - 6777.74 * pow(dist, 2) +
                     3304.80 * dist - 619.931);
            break;
          case 3:
            return (759.69 * (1 + 4 * dist / ofLj)
                           * pow((1 - dist / ofLj), 4));
            break;
        }
        break;
      case 5:
        ctOfLj = pow(rmin / ofLj, 6);
        return (12 * epsl * (ljDist - 1) * ljDist / dist -
                12 * epsl * (ctOfLj - 1) * ctOfLj / ofLj);
        break;
      case 6:
        cerr << "Not supporting smooth LJ force" << endl;
        terminate();
        break;
    }
  }

  //update DPD force strength
  //BUG: Maybe we need different cutoff for 2 direction?
  inline myCoor world::pDpdSt(const partk* atom1, const partk* atom2,
                              const myCoor dist) {
    if(dist >= ofDp)
      return 0.0;
    else if(ofDp >= MINZRO)
      return sqrt(4064.42 * exp(-0.1436 * dist)); 
      //return sqrt(1448.37 * (1 + 1.78 * dist/ofDp)
      //                    * pow((1 - dist/ofDp), 1.78));
      //return (1 - dist / ofDp);
    else
      switch(intgrCon) {
        case 24:
        case 34:
          return 0.5;
          break;
        default:
          return 1.0;
          break;
      }
  }

  inline myCoor world::oDpdSt(const partk* atom1, const partk* atom2,
                              const myCoor dist) {
    if(dist >= ofDp)
      return 0.0;
    else if(ofDp >= MINZRO)
      //return sqrt(934.76 * (1 + 2.62 * dist/ofDp)
      //                   * pow((1 - dist/ofDp), 2.62));
      return (1 - dist / ofDp);
    else
      switch(intgrCon) {
        case 24:
        case 34:
          return 0.5;
          break;
        default:
          return 1.0;
          break;
      }
 }

  inline myCoor world::bondSt(const bound* bond, const int ctrl) {
    myCoor dist;
    switch(ctrl) {
      case 0:
        return ((*bond).obRigd() * ((*bond).obSize() - (*bond).obRelx()));
        break;
      case 1: 
        dist = (*bond).obSize();
        return (-2022.1 + 9526.6 * dist - 14912. * pow(dist, 2)
                                        + 7670.8 * pow(dist, 3));
        break;
      case 2: 
        dist = (*bond).obSize();
        return (71.1900 * pow(dist, 5) - 440.175 * pow(dist, 4) +
                1083.32 * pow(dist, 3) - 1314.21 * pow(dist, 2) +
                782.440 * dist - 184.030);
        break;
      case 3: 
        dist = (*bond).obSize();
        return (46.8684 * pow(dist, 5) - 297.015 * pow(dist, 4) +
                754.680 * pow(dist, 3) - 946.080 * pow(dist, 2) +
                579.340 * dist - 139.300);
        break;
    }
  }

  inline myCoor world::anglSt(const angle* angl, const int ctrl) {
    myCoor dist = (*angl).obSize();
    switch(ctrl) {
      case 0:
        return (-(*angl).obRigd() * ((*angl).obSize() - (*angl).obRelx()));
        break;
      case 1:
        dist = (*angl).obSize();
        return -(53.010 * pow(dist, 5) - 515.00 * pow(dist, 4) +
                1931.2 * pow(dist, 3) - 3474.0 * pow(dist, 2) +
                2978.0 * pow(dist, 1) - 974.5); 
        return (840.889 * pow(dist, 7) + 1104.74 * pow(dist, 6) -
                269.228 * pow(dist, 5) - 701.450 * pow(dist, 4) -
                76.2428 * pow(dist, 3) + 94.9230 * pow(dist, 2) +
                16.8579 * dist + 7.03344) * sin((*angl).obSize());
        break;
      case 2:
        dist = cos((*angl).obSize());
        return ((exp(-3.2726 + 3.0931 * dist - 2.0000 * pow(dist, 2) -
                      2.7968 * pow(dist, 3) + 1.7396 * pow(dist, 4) +
                      3.4985 * pow(dist, 5)) *
                 (17.492 * pow(dist, 4) + 6.9584 * pow(dist, 3) -
                  8.3904 * pow(dist, 2) - 4.0000 * dist + 3.0931) +
                 0.91375) * sin((*angl).obSize()));
        break;
    }
  }

  inline myCoor world::dhedSt(const dihed* dhed, const int ctrl) {
    switch(ctrl) {
      case 0:
        return (0.5 * (2 * (*dhed).obRig2() * sin(2 * (*dhed).obSize()) +
                       3 * (*dhed).obRig2() * sin(3 * (*dhed).obSize()) -
                       (*dhed).obRig1() * sin((*dhed).obSize())));
        break;
    }
  }

  //update using Verlet integrator
  inline void world::verlet() {
    for(int i = 0; i < parts.size(); ++i)
      (parts.at(i)).verlet();
  }

  //update using leap frog integrator
  inline void world::lpFrgA() {
    for(int i = 0; i < parts.size(); ++i)
      (parts.at(i)).lpFrgA();
  }

  inline void world::lpFrgB() {
    for(int i = 0; i < parts.size(); ++i)
      (parts.at(i)).lpFrgB();
  }

  //update using velocity Verlet integrator
  inline void world::vVerlA() {
    for(int i = 0; i < parts.size(); ++i)
      (parts.at(i)).vVerlA();
  }

  inline void world::vVerlB() {
    for(int i = 0; i < parts.size(); ++i) 
      (parts.at(i)).vVerlB();
  }

  //OBSERVERS
  //write simulation conditions
  void world::wrCond(ostream& oput) const {
    oput << boundCon;
    for(int k = 0; k < DIM; ++k)
      oput << '\t' << myVectGet(edge, k);
    oput << '\n' << ensemCon;
    oput << '\n' << colmbCon;
    oput << '\n' << ljPotCon;
    if(LJONPART && (ljPotCon != 0 || ljPotCon != 1)) {
      oput << '\t' << ofLj;
      if(ljPotCon == 3 || ljPotCon == 6)
        oput << '\t' << onLj;
    }
    oput << '\n' << disspCon;
    if(DPDPARTK && (disspCon != 0 || disspCon != 1)) {
      oput << '\t' << ofDp;
      if(disspCon >= 5 && disspCon < 8) {
        if(disspCon != 6)
          oput << '\t' << pmOf;
        if(disspCon != 5)
          oput << '\t' << omOf;
      }
      else if(disspCon >= 8 && disspCon < 11) {
        if(disspCon != 9)
          oput << '\t' << pmOf;
        if(disspCon != 8)
          oput << '\t' << omOf;
      }
    }
    oput << '\n' << verltCon;
    if(verltCon != 0)
      oput << '\t' << ofVl << '\t' << rtVl;
    oput << '\n' << intgrCon;
    if(intgrCon == 32 || intgrCon == 33)
      oput << '\t' << rtUp;
    oput << '\n' << temp;
    oput << '\n';
  }

  //write trajectory
  void world::wrTraj(ostream& oput) const {
    oput << time << endl;
    for(int i = 0; i < parts.size(); ++i)
      (parts.at(i)).wrTraj(oput);
    oput << endl;
  }

  //write coordinate
  void world::wrCoor(ostream& oput) const {
    oput << time << endl;
    for(int i = 0; i < parts.size(); ++i)
      (parts.at(i)).wrCoor(oput);
  }

  const resid* world::obRsid(const string name) const {
    for(int pos = 0; pos < rsids.size(); ++pos)
      if((rsids.at(pos)).obName() == name)
        return &rsids.at(pos);
    return NULL;
  }

  const partk* world::obPart(const string name) const {
    for(int pos = 0; pos < parts.size(); ++pos)
      if((parts.at(pos)).obName() == name)
        return &parts.at(pos);
    return NULL;
  }

  const bound* world::obBond(const partk* atom1, const partk* atom2)
        const {
    for(int pos = 0; pos < bonds.size(); ++pos)
      if((bonds.at(pos)).obLeft() == atom1 &&
         (bonds.at(pos)).obRigh() == atom2)
        return &bonds.at(pos);
    return NULL;
  }

  const angle* world::obAngl(const partk* atom1, const partk* atom2,
                             const partk* atom3) const  {
    for(int pos = 0; pos < angls.size(); ++pos)
      if((angls.at(pos)).obLeft() == atom1 &&
         (angls.at(pos)).obMidl() == atom2 &&
         (angls.at(pos)).obRigh() == atom3)
        return &angls.at(pos);
    return NULL;
  }

  const resid* world::obRsid(const int indx) const {
    assert(indx >= 0 && indx < rsids.size());
    return &rsids.at(indx);
  }

  const partk* world::obPart(const int indx) const {
    assert(indx >= 0 && indx < parts.size());
    return &parts.at(indx);
  }

  const bound* world::obBond(const int indx) const {
    assert(indx >= 0 && indx < bonds.size());
    return &bonds.at(indx);
  }

  const angle* world::obAngl(const int indx) const {
    assert(indx >= 0 && indx < angls.size());
    return &angls.at(indx);
  }

  const dihed* world::obDhed(const int indx) const {
    assert(indx >= 0 && indx < dheds.size());
    return &dheds.at(indx);
  }

  const myVect* world::obMomt() const {
    myVectScl(FORC, 0);
    for(int i = 0; i < parts.size(); ++i) {
      myVectCpy(UNIT, (parts.at(i)).obVelo());
      myVectScl(UNIT, (parts.at(i)).obMass());
      myVectAdd(FORC, UNIT);
    }
    return FORC;
  }

  myCoor world::obKine() const {
    myCoor ener = 0, vSqr;
    switch(intgrCon) {
      case 1:
        for(int i = 0; i < parts.size(); ++i) {
          myVectCpy(UNIT, (parts.at(i)).obPost());
          myVectSub(UNIT, (parts.at(i)).obVelo());
          myVectInv(UNIT, STPSZE); 
          myVectDot(UNIT, UNIT, &vSqr);
          ener += 0.5 * vSqr * (parts.at(i)).obMass();
        }
        break;
      case 2:
      case 21:
      case 24:
      case 3:
      case 31:
      case 32:
      case 33:
      case 34:
        for(int i = 0; i < parts.size(); ++i) {
          myVectCpy(UNIT, (parts.at(i)).obVelo());
          myVectDot(UNIT, UNIT, &vSqr);
          ener += 0.5 * vSqr * (parts.at(i)).obMass();
        }
        break;
    }
    return ener;
  }

  myCoor world::obPote() const {
    myCoor ener = 0, dist, rmin, epsl, LjDist, ctOfLj, delt;
    const partk *atom1, *atom2;
    const bound *bond;
    for(int i = 0; i < parts.size(); ++i) {
      atom1 = &(parts.at(i));
      for(int j = i + 1; j < parts.size(); ++j) {
        //calculate distance
        atom2 = &(parts.at(j));
        obXdif(atom1, atom2, UNIT);
        (const_cast<world*> (this))->apBcon(UNIT);
        dist = myVectNrm(UNIT); 
        //calculate LJ Potential
        if(ljPotCon != 0 && (ljPotCon == 1 || dist < ofLj)) {
          rmin = ((*atom1).obRmin() + (*atom2).obRmin());
          epsl = sqrt((*atom1).obEpsl() * (*atom2).obEpsl());
          LjDist = pow(rmin / dist, 6);
          switch(ljPotCon) {
            case 1:
            case 4:
              ener += epsl * (2 - LjDist) * LjDist;
              break;
            case 5:
              ctOfLj = pow(rmin / ofLj, 6);
              ener += epsl * (2 - LjDist) * LjDist -
                      epsl * (2 - ctOfLj) * ctOfLj;
              break;
          }
        }
        //calculate Intramolecular potential
        if(samRsd(atom1, atom2)) {
          bond = obBond(atom1, atom2);
          if(bond != NULL) {
            delt = dist - (*bond).obRelx();
            ener += 0.5 * (*bond).obRigd() * delt * delt;
          }
        }
      }
    }
    return ener;
  }
  
  myCoor world::obPFVC(const vector<myVect*> *oVel) const {
    const partk *atom1, *atom2;
    int count, cPtr;
    myCoor aveg = 0, fDot;
    for(int i = 0; i < parts.size(); ++i) {
      atom1 = &(parts.at(i));
      for(int j = i + 1; j < parts.size(); ++j) {
        atom2 = &(parts.at(j));
        //BUG: no idea how to handle pair force-velocity correlation
        //     for bonded interaction
        if(samTop(atom1, atom2))
          continue;
        else {
          ++count;
          cPtr = obPInd(&parts, atom1, atom2);
          (const_cast<world*> (this))->nbndFc(atom1, atom2, 0, FTOT);
          myVectDot(FTOT, (*oVel).at(cPtr), &fDot);
          aveg += fDot;
        }
      }
    }
    return (aveg / count);    
  }

  //calculate pair velocity-velocity correlation
  myCoor world::obPVVC(const vector<myVect*> *oVel) const {
    const partk *atom1, *atom2;
    int count, cPtr;
    myCoor aveg = 0, vDot;
    for(int i = 0; i < parts.size(); ++i) {
      atom1 = &(parts.at(i));
      for(int j = i + 1; j < parts.size(); ++j) {
        atom2 = &(parts.at(j));
        if(samBnd(atom1, atom2)) {
          ++count;
          cPtr = obPInd(&parts, atom1, atom2);
          obVdif(atom2, atom1, VELO);
          myVectDot(VELO, (*oVel).at(cPtr), &vDot);
          aveg += vDot;
        }
        else {
          continue;
        }
      }
    }
    return (aveg / count);    
  }

  //calculate velocity autocorrelation
  myCoor world::obtVAC(const vector<myVect*> *oVel) const {
    myCoor aveg = 0, vDot;
    for(int i = 0; i < parts.size(); ++i) {
      myVectDot((parts.at(i)).obVelo(), (*oVel).at(i), &vDot);
      aveg += vDot;
    }
    return (aveg / parts.size());    
  }

  //calculate mean-squared displacement
  myCoor world::obtMSD(const vector<myVect*> *oPos) const {
    myCoor aveg = 0, xDot;
    for(int i = 0; i < parts.size(); ++i) {
      myVectCpy(FORC, (parts.at(i)).obPost());
      myVectSub(FORC, (*oPos).at(i));
      myVectDot(FORC, FORC, &xDot);
      aveg += xDot;
    }
    return (aveg / parts.size());
  }

  //calculate non-bonded distribution
  vector<int> world::obLjDs(const ptype* ptyp1, const ptype* ptyp2,
                            const int bSze) const {
    myCoor ctOf, dist;
    vector<int> bins;
    ctOf = fmax(ofLj, ofDp);
    ctOf = (ctOf < MINZRO) ? (myVectMin(edge) / 2) : (ctOf * 2);
    for(int k = 0; k < bSze; ++k)
      bins.push_back(0);
    for(int i = 0; i < parts.size(); ++i)
      for(int j = i + 1; j < parts.size(); ++j)
        if((parts.at(i)).obType() == ptyp1 &&
           (parts.at(j)).obType() == ptyp2 &&
           !samTop(&(parts.at(i)), &(parts.at(j)))) {
          obXdif(&(parts.at(i)), &(parts.at(j)), UNIT);
          (const_cast<world*> (this))->apBcon(UNIT);
          dist = myVectNrm(UNIT);
          if(dist < ctOf)
            ++bins.at((int) floor(dist / ctOf * bSze));
        }
    return bins;
  }

  //calculate bond distribution
  vector<int> world::ob12Ds(const btype* btyp, const int bSze) const {
    const bound *bond;
    vector<int> bins;
    for(int k = 0; k < bSze; ++k)
      bins.push_back(0);
    for(int i = 0; i < bonds.size(); ++i) {
      bond = &(bonds.at(i));
      if((*bond).obType() == btyp &&
         (*bond).obSize() < 2 * (*bond).obRelx())
        ++bins.at((int) floor((*bond).obSize() / (2 * (*bond).obRelx())
                                               * bSze));
    }
    return bins;
  }

  //calculate angle distribution
  vector<int> world::ob13Ds(const atype* atyp, const int bSze) const {
    const angle *angl;
    vector<int> bins;
    for(int k = 0; k < bSze; ++k)
      bins.push_back(0);
    for(int i = 0; i < angls.size(); ++i) {
      angl = &(angls.at(i));
      if((*angl).obType() == atyp &&
         (*angl).obSize() < ONEPIE) {
        ++bins.at((int) floor((cos((*angl).obSize()) / 2 + 0.5) * bSze));
      }
      else if((*angl).obSize() >= ONEPIE)
        cerr << "Angle size " << (*angl).obSize() / ONEPIE * 180.0
             << " is greater than 180.0 degrees" << endl;
    }
    return bins;
  }

  //NON-MEMBER FUNCTIONS
  template <class type> int obPInd(const vector<type>* parts,
                                   const type* atom1, const type* atom2) {
    int pos1 = atom1 - &((*parts).front()),
        pos2 = atom2 - &((*parts).front());
    if(pos1 > pos2) {
      pos1 = pos1 + pos2;
      pos2 = pos1 - pos2;
      pos1 = pos1 - pos2;
    }
    //NOTE: this is slower, but more elegant
    //int pos1 = min((atom1 - &((*parts).front())),
    //               (atom2 - &((*parts).front()))),
    //    pos2 = max((atom1 - &((*parts).front())),
    //               (atom2 - &((*parts).front())));
    return (pos1 * ((*parts).size() - 1) + pos2 - 1 - pos1 * (pos1 + 1) / 2);
  }
  template int obPInd(const vector<partk>*, const partk*, const partk*);

}
