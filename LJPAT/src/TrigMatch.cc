// Original Author:  Benjamin Radburn-Smith
//         Created:  Wed Nov 16 18:15:13 GMT 2011
#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Provenance/interface/ParameterSetID.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <DataFormats/MuonReco/interface/Muon.h>
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "TrigToolsFuncs.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "TH1F.h"

//#include "TVector3.h"

//namespace trigger{
//  class TriggerEvent;
//}
//namespace edm{
//  class Event;
//  class EventSetup;
//  class ParameterSet;
//}

class TrigMatch : public edm::EDAnalyzer {
	public:
		explicit TrigMatch(const edm::ParameterSet&);
		~TrigMatch();
	private:
		virtual void beginJob() ;
		virtual void beginRun(const edm::Run&,const edm::EventSetup&);
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		TH1F* AllMu, *PassedMu;
		TH1F* AllE, *PassedE;
		TH1F* dist;
		TH1F* numObjsTriggered;
		TH1F* numEJTriggered;
		edm::InputTag trigEventTag_;
		std::string filterName_;
		std::string pathName_;
		std::vector<std::string> filtersOfPath_; //the filters of the path
		HLTConfigProvider hltConfig_; //to translate path names to filter names
};

TrigMatch::TrigMatch(const edm::ParameterSet& iPara){
	edm::Service<TFileService> fs;
	AllMu=fs->make<TH1F>("AllMu","MC Muons Pasing Triggers",1000,0.0,1000.0);
	PassedMu=fs->make<TH1F>("PassedMu","MC Muons Pasing Triggers",1000,0.0,1000.0);
	AllE=fs->make<TH1F>("AllE","MC Electrons Pasing Triggers",1000,0.0,1000.0);
	PassedE=fs->make<TH1F>("PassedE","MC Electrons Pasing Triggers",1000,0.0,1000.0);
	dist=fs->make<TH1F>("dist","MC electron dist from trigged obj",1000,0.0,1.0);
	numObjsTriggered=fs->make<TH1F>("numObjsTriggered","Number of objects in event that passes the trigger",100,0.0,100.0);
	numEJTriggered=fs->make<TH1F>("numEJTriggered","Number of EJs that passes the trigger per event",10,0.0,10.0);
	trigEventTag_ = iPara.getParameter<edm::InputTag>("trigEventTag");
	filterName_ = iPara.getParameter<std::string>("filterName");
	pathName_ = iPara.getParameter<std::string>("pathName");
}

TrigMatch::~TrigMatch(){}

void TrigMatch::beginRun(const edm::Run& run,const edm::EventSetup& setup){
//	std::cerr << "Attempting to match p4's against the filter: "<<filterName_<<std::endl;
//	std::cout << "BeginRun\n";
	bool changed=false;
//	bool foundFilter=false;
	hltConfig_.init(run,setup,trigEventTag_.process(),changed);
//	std::cerr << "Launching change check:\n";
//	if(changed) {
//		if(hltConfig_.triggerIndex(pathName_)<hltConfig_.size()){ //checks trigger exists
//			std::cout <<"Trigger Path: " <<pathName_ << " exists! #filters before=" << filtersOfPath_.size() << " ";
//			filtersOfPath_ =hltConfig_.saveTagsModules(pathName_); //hlt config has changed, load the filters used by the path into our vector
//			std::cout <<"and #filters after="<<filtersOfPath_.size()<<" moduleLabels.size()="<< hltConfig_.moduleLabels(pathName_).size()<<std::endl;
//			std::cerr << "\thas changed\n";
//		for (unsigned int i(0);i<hltConfig_.moduleLabels(pathName_).size();++i){
//			if(hltConfig_.moduleLabels(pathName_).at(i) == filterName_) foundFilter=true;
//		}
//		std::vector<std::string> pathNames;
//		pathNames = hltConfig_.triggerNames();
//		for (std::vector<std::string>::iterator name=pathNames.begin(); name!=pathNames.end();++name){
//			std::cerr << "Path: (PS="<<  <<") " << *name << std::endl;
//		}
//			std::cout << "hltConfig_.moduleLabels("<<pathName_<<").at("<<i<<"): " <<  hltConfig_.moduleLabels(pathName_).at(i) <<std::endl;
//		}
//		else std::cout <<"Trigger Path: " <<pathName_ << " does not exist!\n";
//	}
//	std::cout << "foundFilter?: " << foundFilter << std::endl;
	std::cout << "hltConfig_.tableName: " << hltConfig_.tableName() << std::endl;
}

void TrigMatch::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
//	std::cout << "Event[" << iEvent.id().event() << "] of run: " << iEvent.id().run() <<"\n";
	edm::Handle<trigger::TriggerEvent> trigEvent;	iEvent.getByLabel(trigEventTag_,trigEvent);
	edm::Handle<std::vector<reco::GenParticle> > genParticleCollection; iEvent.getByLabel("genParticles", genParticleCollection);
	edm::Handle<std::vector<reco::Muon> > muonCollection;	iEvent.getByLabel("muons", muonCollection);

	if (!trigEvent.isValid()){
		edm::LogError("TrigMatch") << "Failed to get trigger event, " << trigEventTag_;
		return;
	}

	const reco::Candidate* pMuFromW;
	const reco::Candidate* pEFromW;
	bool IsWmunuEvent = false;
	bool IsWenuEvent = false;
	for(std::vector<reco::GenParticle>::const_iterator MCparticle=genParticleCollection->begin(); MCparticle!=genParticleCollection->end(); ++MCparticle){
		if (MCparticle->pdgId() == 24 || MCparticle->pdgId() == -24){
			if ((MCparticle->daughter(0))->pdgId() != 24 && (MCparticle->daughter(0))->pdgId() != -24){
				std::cerr << "Found a W - which decays to: " ;//<< MCparticle->numberOfDaughters();
				for (unsigned int daug(0);daug < MCparticle->numberOfDaughters(); ++daug){
					if ((MCparticle->daughter(daug))->pdgId() == 13 || (MCparticle->daughter(daug))->pdgId() == -13){
						pMuFromW = MCparticle->daughter(daug);	//Assumes only one W per event
						IsWmunuEvent = true;
					}
					if ((MCparticle->daughter(daug))->pdgId() == 11 || (MCparticle->daughter(daug))->pdgId() == -11){
						pEFromW = MCparticle->daughter(daug);	//Assumes only one W per event
						IsWenuEvent = true;
					}
					std::cerr << (MCparticle->daughter(daug))->pdgId() << " with et: " << (MCparticle->daughter(daug))->et() << "\t";
				}
				std::cerr << std::endl;
			}
		}
	}

	if (IsWmunuEvent){
//		using namespace muon;
		std::cerr << "This is a Wmunu Event\n";
		for(std::vector<reco::Muon>::const_iterator Mu=muonCollection->begin(); Mu!=muonCollection->end(); ++Mu){
//			if (Mu->isGlobalMuon() && )
			bool muPassedId =false;
//			enum muon::SelectionType type=muon::All;
			muPassedId = muon::isGoodMuon(*Mu, muon::AllGlobalMuons);
			if (muPassedId && (deltaR(*Mu,*pMuFromW)<0.3)){	// ToDo: Possibilty for duplicates here -> need to study
				std::cerr << "Muon passed Id cuts!\n";
				std::vector<math::XYZTLorentzVector> trigObjP4s1, trigObjP4s2, trigObjP4s3, trigObjP4s4, trigObjP4s5, trigObjP4s6;
		//		trigtools::getP4sOfObsPassingFilter(trigObjP4s,*trigEvent,"hltSingleMuIsoL3IsoFiltered17",trigEventTag_.process()); //HLT_IsoMu17_v5
		//		trigtools::getP4sOfObsPassingFilter(trigObjP4s2,*trigEvent,"",trigEventTag_.process()); //
		//		trigtools::getP4sOfObsPassingFilter(trigObjP4s3,*trigEvent,"hltSingleMu30L3Filtered30",trigEventTag_.process()); //HLT_Mu30_v1
				trigtools::getP4sOfObsPassingFilter(trigObjP4s1,*trigEvent,"hltSingleMuL2QualIsoL3IsoFiltered20",trigEventTag_.process()); //HLT_IsoMu20_v9 - 3e33
				trigtools::getP4sOfObsPassingFilter(trigObjP4s2,*trigEvent,"hltSingleMuL2QualIsoL3IsoFiltered24",trigEventTag_.process()); //HLT_IsoMu24_v9 - 3e33
				trigtools::getP4sOfObsPassingFilter(trigObjP4s3,*trigEvent,"hltL3IsoL1sMu14Eta2p1L1f0L2f14QL2IsoL3f24L3IsoFiltered",trigEventTag_.process()); //HLT_IsoMu24_eta2p1_v3 - 3e33
				trigtools::getP4sOfObsPassingFilter(trigObjP4s4,*trigEvent,"hltL3IsoL1sMu14Eta2p1L1f0L2f14QL2IsoL3f30L3IsoFiltered",trigEventTag_.process()); //HLT_IsoMu30_eta2p1_v3 - 3e33
				trigtools::getP4sOfObsPassingFilter(trigObjP4s5,*trigEvent,"hltSingleMu40L2QualL3Filtered40",trigEventTag_.process()); //HLT_Mu40_v6 - 3e33
				trigtools::getP4sOfObsPassingFilter(trigObjP4s6,*trigEvent,"hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40",trigEventTag_.process()); //HLT_Mu40_eta2p1_v1 - 3e33
		//		AllMu->Fill(pMuFromW->pt());
				AllMu->Fill(Mu->pt());
				bool passedPath=false;
				for(size_t objNr=0;objNr<trigObjP4s1.size();objNr++){
		//			if (deltaR(trigObjP4s1.at(objNr),*pMuFromW)<0.3){
					if (deltaR(trigObjP4s1.at(objNr),*Mu)<0.3){
						passedPath=true;
		//				std::cerr << "\tFilter1 Matched: " << pMuFromW->pdgId()
		//						<< " to the obj passing the trig at dist= " << deltaR(trigObjP4s.at(objNr),*pMuFromW)
		//						<< "\twith et= " << pMuFromW->et()
		//						<< "\tand parent= " << pMuFromW->mother()->pdgId()
		//						<< std::endl;
					}
				}
				for(size_t objNr=0;objNr<trigObjP4s2.size();objNr++){
		//			if (deltaR(trigObjP4s2.at(objNr),*pMuFromW)<0.3){
					if (deltaR(trigObjP4s2.at(objNr),*Mu)<0.3)	passedPath=true;
				}
				for(size_t objNr=0;objNr<trigObjP4s3.size();objNr++){
		//			if (deltaR(trigObjP4s3.at(objNr),*pMuFromW)<0.3) passedPath=true;
					if (deltaR(trigObjP4s3.at(objNr),*Mu)<0.3) passedPath=true;
				}
				for(size_t objNr=0;objNr<trigObjP4s4.size();objNr++){
		//			if (deltaR(trigObjP4s4.at(objNr),*pMuFromW)<0.3) passedPath=true;
					if (deltaR(trigObjP4s4.at(objNr),*Mu)<0.3) passedPath=true;
				}
				for(size_t objNr=0;objNr<trigObjP4s5.size();objNr++){
		//			if (deltaR(trigObjP4s5.at(objNr),*pMuFromW)<0.3) passedPath=true;
					if (deltaR(trigObjP4s5.at(objNr),*Mu)<0.3) passedPath=true;
				}
				for(size_t objNr=0;objNr<trigObjP4s6.size();objNr++){
		//			if (deltaR(trigObjP4s6.at(objNr),*pMuFromW)<0.3) passedPath=true;
					if (deltaR(trigObjP4s6.at(objNr),*Mu)<0.3) passedPath=true;
				}
		//		if (passedPath) PassedMu->Fill(pMuFromW->pt());
				if (passedPath) PassedMu->Fill(Mu->pt());
			}
		}
	}

	if (IsWenuEvent){
		std::cerr << "This is a Wenu Event\n";
		std::vector<math::XYZTLorentzVector> trigObjP4sE, trigObjP4sE2, trigObjP4sE3;
//		trigtools::getP4sOfObsPassingFilter(trigObjP4sE,*trigEvent,"hltEle27CaloIdTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",trigEventTag_.process()); //HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2
		trigtools::getP4sOfObsPassingFilter(trigObjP4sE,*trigEvent,"hltEle32CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",trigEventTag_.process()); //HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v7 - 3e33 PS60
		trigtools::getP4sOfObsPassingFilter(trigObjP4sE2,*trigEvent,"",trigEventTag_.process()); //HLT_Ele30_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_DiCentralJet30_PFMHT25_v1 - 3e33 PS1
//		trigtools::getP4sOfObsPassingFilter(trigObjP4sE3,*trigEvent,"",trigEventTag_.process()); //
		AllE->Fill(pEFromW->pt());
		bool passedPathE=false;
		for(size_t objNr=0;objNr<trigObjP4sE.size();objNr++){
			if (deltaR(trigObjP4sE.at(objNr),*pEFromW)<0.3){
				passedPathE=true;
//				std::cerr << "\tFilter1 Matched: " << pMuFromW->pdgId()
//						<< " to the obj passing the trig at dist= " << deltaR(trigObjP4s.at(objNr),*pMuFromW)
//						<< "\twith et= " << pMuFromW->et()
//						<< "\tand parent= " << pMuFromW->mother()->pdgId()
//						<< std::endl;
			}
		}
//		for(size_t objNr=0;objNr<trigObjP4s3.size();objNr++){
//			if (deltaR(trigObjP4sE3.at(objNr),*pEFromW)<0.3){
//				passedPath=true;
////						std::cerr << "\tFilter3 Matched: " << pMuFromW->pdgId()
////								<< " to the obj passing the trig at dist= " << deltaR(trigObjP4s3.at(objNr),*pMuFromW)
////								<< "\twith et= " << pMuFromW->et()
////								<< "\tand parent= " << pMuFromW->mother()->pdgId()
////								<< std::endl;
//			}
//		}
		if (passedPathE) PassedE->Fill(pEFromW->pt());
	}






//
//	std::string filter2="HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2";
//	std::string filter3="HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_v2";
//
//	std::vector<math::XYZTLorentzVector> trigObjP4s, trigObjP4s2, trigObjP4s3;
//	trigtools::getP4sOfObsPassingFilter(trigObjP4s,*trigEvent,filterName_,trigEventTag_.process());
////	std::cerr << "trigObjP4s.size(): " << trigObjP4s.size() << std::endl;
////	std::cerr << "trigEventTag_.process(): " << trigEventTag_.process() << std::endl;
////	we get the process name from the trigger event process as they have to be the same...
//
////  do something with them...
//	std::cerr <<"for filter "<<filterName_<<" nr objects passing = "<<trigObjP4s.size()<<std::endl;
//	numObjsTriggered->Fill(trigObjP4s.size());
//	unsigned int NumEJsPassing(0);
//	for(size_t objNr=0;objNr<trigObjP4s.size();objNr++){
//		bool EJtriggered=false;
//		std::cerr <<"     et obj "<<objNr<<" : "<<trigObjP4s[objNr].Et()<<std::endl;
//		for(std::vector<reco::GenParticle>::const_iterator MCparticle=genParticleCollection->begin(); MCparticle!=genParticleCollection->end(); ++MCparticle){
//			if (deltaR(trigObjP4s.at(objNr),*MCparticle)<0.3){
//				if (MCparticle->pdgId() == 11 || MCparticle->pdgId() == -11){
//					std::cerr << "\tMatched: " << MCparticle->pdgId()
//							<< " to the obj passing the trig at dist= " << deltaR(trigObjP4s.at(objNr),*MCparticle)
//							<< "\twith et= " << MCparticle->et()
//							<< "\tand parent= " << MCparticle->mother()->pdgId()
//							<< std::endl;
//					dist->Fill(deltaR(trigObjP4s.at(objNr),*MCparticle));
//					if (MCparticle->mother()->pdgId() == 3000001) EJtriggered=true;
//				}
//			}
//		}
//		if (EJtriggered) ++NumEJsPassing;
//	}
//	numEJTriggered->Fill(NumEJsPassing);
//	for(std::vector<reco::GenParticle>::const_iterator MCparticle=genParticleCollection->begin(); MCparticle!=genParticleCollection->end(); ++MCparticle){
//		if (MCparticle->pdgId() == 3000006) std::cerr << "hd2 et = " << MCparticle->et() << std::endl;
//	}
//
////	//now we are going to do the same again for the paths filters
////	for(size_t filterNr=0;filterNr<filtersOfPath_.size();filterNr++){
////		trigtools::getP4sOfObsPassingFilter(trigObjP4s,*trigEvent,filtersOfPath_[filterNr],trigEventTag_.process());
////		std::cout <<"\tfilter "<<filtersOfPath_[filterNr]<<" nr objects passing "<<trigObjP4s.size()<<std::endl;
////		for(size_t objNr=0;objNr<trigObjP4s.size();objNr++){
////			std::cout <<"\tet obj "<<objNr<<" : "<<trigObjP4s[objNr].Et()
////					<< " " << trigObjP4s[objNr].X()
////					<< " " << trigObjP4s[objNr].Y()
////					<< " " << trigObjP4s[objNr].Z()
////					<< "\t" << trigObjP4s[objNr].Px()
////					<< " " << trigObjP4s[objNr].Py()
////					<< " " << trigObjP4s[objNr].Pz()
////					<<std::endl;
////		}
////	}

}

void TrigMatch::beginJob(){
}
//define this as a plug-in
DEFINE_FWK_MODULE(TrigMatch);
