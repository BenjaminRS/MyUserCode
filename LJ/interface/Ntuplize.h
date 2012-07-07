// system include files
#include <memory>
#include <vector>

//#include "TVector2.h"
#include "TVector3.h"
#include "TSystem.h"
#include "TMath.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

// #include "SimDataFormats/Track/interface/SimTrack.h"
// #include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateElectronExtra.h"
#include "RecoParticleFlow/PFProducer/interface/PFElectronExtraEqual.h"

#include "DataFormats/MuonReco/interface/Muon.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/TrackJet.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "CommonTools/TriggerUtils/interface/GenericTriggerEventFlag.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "TH1.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "MyObjs.h"

class Ntuplize: public edm::EDAnalyzer {
public:
	explicit Ntuplize(const edm::ParameterSet&);
	~Ntuplize();
	TFile * outf;
	TTree * tree;
	TH1F * hevt;
	TNtuple * nt;
	std::string tupfile;
	std::string run_scheme;
	std::string trigger_scheme;
	bool do_data, do_signalMC, do_zeeMC, do_wenuMC, do_wwMC, do_wjetMC;
	bool do_jet_trigger, do_electron_trigger, do_muon_trigger, do_photon_trigger;
	std::vector<std::pair<reco::GenParticle, unsigned> > elevec;
	std::vector<std::pair<reco::GenParticle, unsigned> > genpvec;
	std::vector<MyHLTObj*> * hltobjVector;
	std::vector<MyGenParticle*> * genpVector;
	std::vector<MyGenParticle*> * genelectronVector;
	std::vector<MyGenParticle*> * genmuonVector;
	std::vector<MyRecoElectron*> * recoelectronVector;
//	std::vector<MyCaloJet*> * calogenjetVector;
	std::vector<MyPFJet*> * genjetVector;
	std::vector<MyPFJet*> * pfjetVector;
	std::vector<MyPFElectron*> * pfeleVector;
	std::vector<MyCaloJet*> * ak5calojetVector;
	std::vector<MyCaloJet*> * ak7calojetVector;
	std::vector<MyMuon*> * muonVector;
	std::vector<MyPhoton*> * photonVector;
	std::vector<MyPU*> * puVector;
	MyHLT hlt;
	MyMET met;
	MyEvent evtdata;
	bool _hlt_initialized;
	HLTConfigProvider hltConfig_;
	std::string processname_;
	std::vector<unsigned> triggerIndices_;
	std::vector<std::string> triggerNames_;
	std::vector<std::string> lastFilterLabels_;
private:
	void findSignalElectrons(const reco::GenParticle &p, std::vector<std::pair<reco::GenParticle, unsigned> > &vec, unsigned& label);
	void findZElectronsMC(const reco::GenParticle &p, std::vector<std::pair<reco::GenParticle, unsigned> > &vec, unsigned& label);
	void findWElectronsMC(const reco::GenParticle &p, std::vector<std::pair<reco::GenParticle, unsigned> > &vec, unsigned& label);
	void findWDaughters(const reco::GenParticle &p, std::vector<std::pair<reco::GenParticle, unsigned> > &vec, unsigned& label);
	virtual void beginJob();
	virtual void beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob();
	unsigned long long checkTrigger(const edm::TriggerResults & HLTR, std::vector<unsigned> & triggerIndices);	// BRS: not sure how 32bit compilers behave here
	bool verbose;
};
