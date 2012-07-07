#include "TVector3.h"
#include "TLorentzVector.h"

class MyGenParticle { 
 public:
  int pdg;
  TLorentzVector p4;
  unsigned label;
};

class MyRecoElectron { 
 public:
  unsigned label;
  TLorentzVector p4;
  double SCenergy, SCeta, SCphi;
  unsigned trkDriven, ecalDriven;
  // ID variables
  double sigieie;
  double deta;
  double dphi;
  double HoverE;
  double eSCoverP;
  double Itrk;
  double Iecal;
  double Ihcal;
  double convDist;
  double convDcot;
  double misHits;
  double fBrem;
  double convR;
  double mva_epi; // for PF
  unsigned isEB, isEE;
};

class MyHLTObj { 
 public:
  TVector3 p3;
  unsigned label;
  std::string triggername;
  std::string filtername;
};

class MyHLT { 
 public: 
  unsigned long long bits; 
      // 0x1   HLT_Jet30U
      // 0x2   HLT_Jet50U
      // 0x4   HLT_Jet70U
      // 0x8   HLT_Jet100U
};

class MyMET {
 public : 
  TVector2 pf;
  TVector2 calo;
};


class MyGsfTrack { 
 public:
  unsigned lostHits;
  unsigned npixHits;
  unsigned pixHitMask;
};

class MyPFElectron { 
 public:
  unsigned gsfIndex;
  double pt;
  unsigned ecalDriven;
  unsigned ecalDrivenSeed;
  unsigned tkDrivenSeed;
  unsigned passMvaPreselection;
  unsigned passCutBasedPreselection;
  double mva_e_pi;
  unsigned lostHits;
  unsigned npixHits;
  unsigned pixHitMask;
  double convR;
};


class MyPFJet { 
 public:
  unsigned index;
  TLorentzVector p4;
  TLorentzVector p4cor;
  float et;
  float area;
  float etaeta;
  float phiphi;
  float etaphi;
  double dR;
  double dRminpT;
  double Itrk, Ical;
  double Et_0_2;
  unsigned ntrks;
  unsigned nele;
  unsigned nphoton;
  float emf;
  float cr;
  float crtrans;
  float had;
  float em;
  unsigned gsfeleindices;
  unsigned nGsfTracks;
  int pfele_min_index;
  int pfele_max_index;
  float electron_energy, photon_energy;
  float electron_energyFrac, photon_energyFrac;
  float neutralHadFrac, chargedHadFrac;
  float neutralEMFrac, chargedEMFrac;
  float ecalConstitEnergy, hcalConstitEnergy;
};


class MyCaloJet { 
 public:
  unsigned index;
  TLorentzVector p4;
  TLorentzVector p4cor;
  float energy;
  float et;
  float area;
  float etaeta;
  float phiphi;
  float etaphi;
  double dR;
  double Itrk, Ical;
  double Et_0_2;
  unsigned ntrks;
  float emf;
  float cr;
  float crtrans;
  float had;
  float em;
  float btag_secvtx;
  float btag_softele;
  float btag_jetprob;
  unsigned ntrkAbove1, ntrkAbove2_5, ntrkAbove5;
  float n90Hits, fHPD;


};

class MyEvent { 
 public : 
  unsigned run;
  unsigned event;
  unsigned lumi;
};

class MyMuon {
 public : 
  TLorentzVector p4;
  float trkRelChi2,  staRelChi2;
  float Itrk;
  unsigned isGlobal;
  unsigned isTracker;
  unsigned isSA;

};

class MyPhoton {
 public : 
  TLorentzVector p4;
  float Iecal, Ihcal, Itrk;
  float hOe;
  float sigieie;
  unsigned isEB, isEE;
  unsigned hasPixSeed;
};

class MyPU {
 public : 
  int bunchCrossing;
  int nInteractions;
};
