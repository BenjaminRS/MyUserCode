#ifndef MYOBJS_H
#define MYOBJS_H

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TObject.h"

namespace brs{
	class MyGenParticle : public TObject{
	 public:
		MyGenParticle(){}
		~MyGenParticle(){}
	  int pdg;
	  TLorentzVector p4;
	  unsigned label;
	  ClassDef(MyGenParticle,1)
	};

	class MyRecoElectron : public TObject{
	 public:
		MyRecoElectron(){}
		~MyRecoElectron(){}
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
	  double E;
	  double p_in;
	  double trkDxy;
	  double trkDz;
	  double convR;
	  double mva_epi; // for PF
	  unsigned isEB, isEE;
	  ClassDef(MyRecoElectron,1)
	};

	class MyHLTObj : public TObject{
	 public:
		MyHLTObj(){}
		~MyHLTObj(){}
	  TVector3 p3;
	  unsigned label;
	  std::string triggername;
	  std::string filtername;
	  ClassDef(MyHLTObj,1)
	};

	class MyHLT : public TObject{
	 public:
		MyHLT(){}
		~MyHLT(){}
	  unsigned long long bits;
		  // 0x1   HLT_Jet30U
		  // 0x2   HLT_Jet50U
		  // 0x4   HLT_Jet70U
		  // 0x8   HLT_Jet100U
	  ClassDef(MyHLT,1)
	};

	class MyMET : public TObject{
	 public :
		MyMET(){}
		~MyMET(){}
	  TVector2 pf;
	  TVector2 calo;
	  ClassDef(MyMET,1)
	};


	class MyGsfTrack : public TObject{
	 public:
		MyGsfTrack(){}
		~MyGsfTrack(){}
	  unsigned lostHits;
	  unsigned npixHits;
	  unsigned pixHitMask;
	  ClassDef(MyGsfTrack,1)
	};

	class MyPFElectron : public TObject{
	 public:
		MyPFElectron(){}
		~MyPFElectron(){}
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
	  ClassDef(MyPFElectron,1)
	};


	class MyPFJet : public TObject{
	 public:
		MyPFJet(){}
		~MyPFJet(){}
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
	  double met;
	  TLorentzVector metp4;
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
	  ClassDef(MyPFJet,1)
	};


	class MyCaloJet : public TObject{
	 public:
		MyCaloJet(){}
		~MyCaloJet(){}
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
	  unsigned nele;
	  unsigned neleFixed;
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
	  ClassDef(MyCaloJet,1)
	};

	class MyEvent : public TObject{
	 public :
		MyEvent(){}
		~MyEvent(){}
	  unsigned run;
	  unsigned event;
	  unsigned lumi;
	  bool goodVertex;
	  ClassDef(MyEvent,1)
	};

	class MyMuon : public TObject{
	 public :
		MyMuon(){}
		virtual ~MyMuon(){}
	  TLorentzVector p4;
	  float trkRelChi2,  staRelChi2;
	  float Itrk;
	  unsigned isGlobal;
	  unsigned isTracker;
	  unsigned isSA;
	  double globalTrkNormChi2;
	  int globalTrkMuHits;
	  int 	numMatchedStations;
	  double trkIPxy;
	  double trkIPz;
	  int numPixelHits;
	  int numTrkHits;
	  double innerTrkNormChi2;
	  int numPixelLayers;
	  double combIso;
	  double trkIso;
	  bool globalTrkIsNonNull;
	  bool innerTrkIsNonNull;
	  ClassDef(MyMuon,1)
	};

	class MyPhoton : public TObject{
	 public :
		MyPhoton(){}
		~MyPhoton(){}
	  TLorentzVector p4;
	  float Iecal, Ihcal, Itrk;
	  float hOe;
	  float sigieie;
	  unsigned isEB, isEE;
	  unsigned hasPixSeed;
	  ClassDef(MyPhoton,1)
	};

	class MySC : public TObject{
	 public :
		MySC(){}
		~MySC(){}
	  float energy;
	  TVector3 point;
	  float Epreshower;
	  float phiWidth;
	  float etaWidth;
	  ClassDef(MySC,1)
	};

	class MyPU : public TObject{
	 public :
		MyPU(){}
		~MyPU(){}
	  int bunchCrossing;
	  int nInteractions;
	  ClassDef(MyPU,1)
	};
}

ClassImp(brs::MyGenParticle)
ClassImp(brs::MyRecoElectron)
ClassImp(brs::MyHLTObj)
ClassImp(brs::MyHLT)
ClassImp(brs::MyMET)
ClassImp(brs::MyGsfTrack)
ClassImp(brs::MyPFElectron)
ClassImp(brs::MyPFJet)
ClassImp(brs::MyCaloJet)
ClassImp(brs::MyEvent)
ClassImp(brs::MyMuon)
ClassImp(brs::MyPhoton)
ClassImp(brs::MySC)
ClassImp(brs::MyPU)
#endif
