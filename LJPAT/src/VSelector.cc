#include <map>
#include <string>

#include "TH1.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TLorentzVector.h"
//#include "PhysicsTools/SelectorUtils/interface/SimpleCutBasedElectronIDSelectionFunctor.h"
//#include <DataFormats/FWLite/interface/ESHandle.h>


class VSelector : public edm::EDAnalyzer {
	public:
		explicit VSelector(const edm::ParameterSet&);
		~VSelector();
	private:
		virtual void beginJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;
		// simple map to contain all histograms; histograms are booked in the beginJob() method
		std::map<std::string,TH1F*> histContainer_;
		// plot number of towers per jet
//		TH1F* jetTowers_;
//		edm::InputTag photonSrc_;
		edm::InputTag elecSrc_;
//		edm::InputTag muonSrc_;
//		edm::InputTag tauSrc_;
//		edm::InputTag jetSrc_;
//		edm::InputTag metSrc_;
//		TH1F* SigmaIetaIetaB;
//		TH1F* SigmaEtaEtaB;
//		TH1F* SigmaIetaIetaE;
//		TH1F* SigmaEtaEtaE;
//		TH1F* SCpass;
//		TH1F* SCall;
//		TH1F* numElectrons;
//		TH1F* numMCElectrons;
//		TH1F* numGSFElectrons;
//		TH1F* MCElePt;
//		TH1F* PatElePt;
//		TH1F* ElePt;
//		TH1F* NumPFCands;
//		TH1F* PFElePt;
//		TH1F* PFCandObjs;
//		TH1F* MCParticles;
		TH1F* AllPatElePt;
		TH1F* SelectedPatElePt;
		TH1F* TMP;
		unsigned int EventNum;
		unsigned int numMCele;
		unsigned int numPFMatched;
		unsigned int tracksNextToEle;
};

VSelector::VSelector(const edm::ParameterSet& iConfig):
	histContainer_(),
//	photonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("photonSrc")),
	elecSrc_(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc"))
//	muonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc")),
//	tauSrc_(iConfig.getUntrackedParameter<edm::InputTag>("tauSrc" )),
//	jetSrc_(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc" ))
//	metSrc_(iConfig.getUntrackedParameter<edm::InputTag>("metSrc" ))
{

}

VSelector::~VSelector(){}

void VSelector::beginJob(){
	EventNum=0;
	edm::Service<TFileService> fs;
  // book histograms: uncomment the following line to book the jetTowers_ histogram
  //jetTowers_= fs->make<TH1F>("jetTowers", "towers per jet",   90, 0,  90);
//  histContainer_["photons"]=fs->make<TH1F>("photons", "photon multiplicity",   10, 0,  10);
//	histContainer_["elecs"  ]=fs->make<TH1F>("elecs",   "electron multiplicity", 10, 0,  10);
//	histContainer_["muons"  ]=fs->make<TH1F>("muons",   "muon multiplicity",     10, 0,  10);
//	histContainer_["taus"   ]=fs->make<TH1F>("taus",    "tau multiplicity",      10, 0,  10);
//	histContainer_["jets"   ]=fs->make<TH1F>("jets",    "jet multiplicity",      10, 0,  10);
//	histContainer_["met"    ]=fs->make<TH1F>("met",     "missing E_{T}",         20, 0, 100);
//	SigmaIetaIetaB=fs->make<TH1F>("SigmaIetaIetaB","SigmaIetaIeta Barrel",	200,0.0,0.05);
//	SigmaEtaEtaB=fs->make<TH1F>("SigmaEtaEtaB","SigmaEtaEta Barrel",		200,0.0,0.05);
//	SigmaIetaIetaE=fs->make<TH1F>("SigmaIetaIetaE","SigmaIetaIeta Endcap",	200,0.0,0.05);
//	SigmaEtaEtaE=fs->make<TH1F>("SigmaEtaEtaE","SigmaEtaEta Endcap",		200,0.0,0.05);
//	SCpass=fs->make<TH1F>("SCPass","LJets which pass Scheme C",	100,0.0,500.0);
//	SCall=fs->make<TH1F>("SCAll","All LJets",					100,0.0,500.0);
//	numElectrons=fs->make<TH1F>("numElectrons","Number of PAT::Electrons in LJ", 20,0.0,20.0);
//	numMCElectrons=fs->make<TH1F>("numMCElectrons","Number of MC Electrons in LJ", 20,0.0,20.0);
//	numGSFElectrons=fs->make<TH1F>("numGSFElectrons","Number of Gsf Electrons in LJ", 20,0.0,20.0);
//	MCElePt=fs->make<TH1F>("MCElectronsPt","Pt of MC Electrons in LJ", 100,0.0,100.0);
//	PatElePt=fs->make<TH1F>("PatElectronsPt","Pt of Pat Electrons in LJ", 100,0.0,100.0);
	AllPatElePt=fs->make<TH1F>("AllPatElectronsPt","Pt of All Pat Electrons in Event", 200,0.0,200.0);
	SelectedPatElePt=fs->make<TH1F>("SelectedPatElectronsPt","Pt of Selected Pat Electrons in Event", 200,0.0,200.0);
//	ElePt=fs->make<TH1F>("GsfElectronsPt","Pt of GsfElectrons in LJ", 100,0.0,100.0);
//	NumPFCands=fs->make<TH1F>("numPFCands","Number of PFCandidates in LJ", 50,0.0,50.0);
//	PFElePt=fs->make<TH1F>("PFElectronsPt","Pt of PF Electrons (RECO) in LJ", 100,0.0,100.0);
//	PFCandObjs=fs->make<TH1F>("PFCandsInLJ","ID of PFCandidates in LJ", 10,0.0,10.0);
//	MCParticles=fs->make<TH1F>("MCParticlesInLJ","ID of MCParticles in LJ", 4000000,0.0,4000000.0);
//	TMP=fs->make<TH1F>("tmp","tmp", 200,-5.0,5.0);
	TMP=fs->make<TH1F>("tmp","tmp", 10,0.0,10.0);
//	numMCele=0;
//	numPFMatched=0;
//	tracksNextToEle=0;
}

void VSelector::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
	edm::Handle<edm::View<pat::Electron> > electrons;	iEvent.getByLabel(elecSrc_,electrons);
//	edm::Handle<edm::View<pat::Muon> > muons;	iEvent.getByLabel(muonSrc_,muons);
//	edm::Handle<edm::View<pat::Tau> > taus;	iEvent.getByLabel(tauSrc_,taus);
//	edm::Handle<edm::View<pat::Jet> > jets;	iEvent.getByLabel(jetSrc_,jets);
//	edm::Handle<edm::View<pat::MET> > mets;	iEvent.getByLabel(metSrc_,mets);
//	edm::Handle<edm::View<pat::Photon> > photons; iEvent.getByLabel(photonSrc_,photons);
//	edm::Handle<std::vector<reco::GenParticle> > genParticleCollection; iEvent.getByLabel("genParticles", genParticleCollection);
	edm::Handle<std::vector<reco::PFCandidate> > pfCandCollection; iEvent.getByLabel("pfAllElectronsPFlow", pfCandCollection);
//	edm::Handle<std::vector<reco::PFCandidate> > pfCandCollection; iEvent.getByLabel("particleFlow","electrons",pfCandCollection );
//	edm::Handle<std::vector<reco::PFCandidate> > pfCandCollection; iEvent.getByLabel("particleFlow","",pfCandCollection );
//	edm::Handle<std::vector<reco::Track> > trackCollection; iEvent.getByLabel("generalTracks","",trackCollection );
//	edm::Handle<std::vector<reco::GsfElectron> > electronCollection; iEvent.getByLabel("gsfElectrons","",electronCollection );

	std::cerr << "Event["<< EventNum << "]\n";

//    for (unsigned int MCpos(0); MCpos < genParticleCollection->size(); ++MCpos){
//		const reco::GenParticle* MCparticle = &(genParticleCollection->at(MCpos));
//		if ( MCparticle->pdgId() == 11 || MCparticle->pdgId() == -11 ){ //electron
//			if (MCparticle->mother()->pdgId() == 23){  //From a Z
//				for(edm::View<pat::Electron>::const_iterator patele=electrons->begin(); patele!=electrons->end(); ++patele){
//					if (deltaR(*patele,*MCparticle)<0.3){
//						std::cerr << "PAT electron matched to a MC gen particle electron\n";
//						if (patele->isEE()){
//							SigmaIetaIetaE->Fill(patele->sigmaIetaIeta());
//							SigmaEtaEtaE->Fill(patele->sigmaEtaEta());
//						}
//						else{
//							SigmaIetaIetaB->Fill(patele->sigmaIetaIeta());
//							SigmaEtaEtaB->Fill(patele->sigmaEtaEta());
//						}
//					}
//				}
//			}
//		}
//    }

//	################ Z selection ######################

//	 // get the track collection
//	 edm::Handle<reco::TrackCollection> ctfTracks;
//	 iEvent.getByLabel("generalTracks",ctfTracks);
//	 // the magnetic field: Warning! this is only for MC, for data check this twiki: https://twiki.cern.ch/twiki/bin/view/CMS/ConversionBackgroundRejection
//	 edm::ESHandle<MagneticField> magneticField;
//	 esetup.get<IdealMagneticFieldRecord>().get(magneticField);
//	 double bfield = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
////	  define your function for  WP95 using relative isolations
//	 SimpleCutBasedElectronIDSelectionFunctor patSele95(SimpleCutBasedElectronIDSelectionFunctor::relIso95, bfield, ctfTracks);


//	for(std::vector<reco::PFCandidate>::const_iterator PFele=pfCandCollection->begin(); PFele!=pfCandCollection->end(); ++PFele){
//		if ((PFele->eta()>1.566 && PFele->eta()<2.5) || (PFele->eta()<-1.566 && PFele->eta()>-2.5)){ //Endcap
////			TMP->Fill(PFele->eta());
//			if (PFele->)
//		}
//		if (PFele->eta()<1.4442 && PFele->eta()>-1.4442){ //Barrel
////			TMP->Fill(PFele->eta());
//		}
//	}
//	std::vector<const pat::Electron*> eleCands;
	std::vector<edm::View<pat::Electron>::const_iterator> eleCands;
	eleCands.clear();
	unsigned int numEle(0);
	for(edm::View<pat::Electron>::const_iterator patele=electrons->begin(); patele!=electrons->end(); ++patele){
//		bool pass=false;
//		pass = patSele95(patele);
//		std::cout << "Pass?: " << pass << std::endl;
		AllPatElePt->Fill(patele->pt());
		if (patele->isEE()){ //Endcap - does this avoid gap?
//			TMP->Fill(patele->eta());
			// Now for WP95:
			if ((patele->sigmaIetaIeta() < 0.03) &&
				(patele->deltaPhiSuperClusterTrackAtVtx() < 0.7) &&
				(patele->deltaEtaSuperClusterTrackAtVtx() < 0.01) &&
				(patele->hadronicOverEm() < 0.15)){
				SelectedPatElePt->Fill(patele->pt());
				if (patele->pt() > 20){
//					eleCands.push_back(&*patele);
					eleCands.push_back(patele);
					++numEle;
				}
			}

		}
		if (patele->isEB()){ //Barrel - does this avoid gap?
			// Now for WP95:
			if ((patele->sigmaIetaIeta() < 0.01) &&
				(patele->deltaPhiSuperClusterTrackAtVtx() < 0.8) &&
				(patele->deltaEtaSuperClusterTrackAtVtx() < 0.007) &&
				(patele->hadronicOverEm() < 0.15)){
				SelectedPatElePt->Fill(patele->pt());
				if (patele->pt() > 20){
//					eleCands.push_back(&*patele);
					eleCands.push_back(patele);
					++numEle;
				}
			}

		}
	}

//	for (size_t i(0); i < eleCands.size(); ++i){
//		std::cout << "eleCands.at("<<i<<"): " << *eleCands.at(i) << std::endl;
//		std::cout << "\tID size: " << eleCands.at(i)->electronIDs().size() << std::endl;
//		for (size_t j(0);j<eleCands.at(i)->electronIDs().size(); ++j){
//			std::cout << "\tIDs: " << eleCands.at(i)->electronIDs().at(j).first << " " << eleCands.at(i)->electronIDs().at(j).second << std::endl;
//		}
////		std::cout << "\tID: " << eleCands.at(i)->electronID("simpleEleId95relIso") << std::endl;
//	}

	for (size_t firstPos(0); firstPos < eleCands.size(); ++firstPos){
		for (size_t secondPos((firstPos+1)); secondPos < eleCands.size(); ++secondPos){
			TLorentzVector e1(eleCands.at(firstPos)->px(), eleCands.at(firstPos)->py(), eleCands.at(firstPos)->pz(), eleCands.at(firstPos)->energy());
			TLorentzVector e2(eleCands.at(secondPos)->px(),eleCands.at(secondPos)->py(),eleCands.at(secondPos)->pz(),eleCands.at(secondPos)->energy());
			TLorentzVector Z = e1 + e2;
			if (Z.M()>75.0 && Z.M()<105.0){
				std::cerr << "Found a Z Cand!: " << Z.M() << std::endl;
				std::cerr << "\tele1: " << *eleCands.at(firstPos) << std::endl;
				std::cerr << "\tele2: " << *eleCands.at(secondPos) << std::endl;
			}
		}
	}


//	std::cout << "numEle: " << numEle << std::endl;
	TMP->Fill(numEle);
//	################ Z selection END ##################
//
//
//
//	const pat::Jet* jet1=NULL;
//	const pat::Jet* jet2=NULL;
////	std::cout << "trackCollection->size(): " << trackCollection->size() << std::endl;
//
//	for(std::vector<reco::GenParticle>::const_iterator MCparticle=genParticleCollection->begin(); MCparticle!=genParticleCollection->end(); ++MCparticle){
////		if ( MCparticle->pdgId() == 11 || MCparticle->pdgId() == -11){	//e
//		if ( MCparticle->pdgId() == 3000006){	//hd2
////			if (MCparticle->mother()->pdgId() == 3000001){	//gamma_d
//			if (MCparticle->mother()->pdgId() == 25){  // Protect from loopback
////				++numMCele;
//			std::cerr << "MCparticle->pdgId(): " << MCparticle->pdgId() << " ("<< MCparticle->eta() << "," << MCparticle->phi() << ")" << std::endl;
////			if (MCparticle->mother()) std::cout << " has mother: " << MCparticle->mother()->pdgId() << std::endl;
////			else std::cout << std::endl;
////				for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
//				for( size_t k=0; k<jets->size(); ++k ) {
//					const pat::Jet* jet = &(jets->at(k));
////				std::cout << " genjet pt: " << jet->pt() << " at eta: " << jet->eta() << " phi: " << jet->phi() << std::endl;
//					if (deltaR(*jet,*MCparticle)<0.3){	//Match PatJet to the GenParticle = hd2
//						std::cerr << "\tpatjet matched to hd2; jet1= " << jet1 << " jet2= " << jet2 << std::endl;
//						if (!jet1) jet1=jet;
//						else jet2=jet;
//						std::cerr << "\tafter: jet1= " << jet1 << " jet2= " << jet2 << std::endl;
////						++numMatched;
////				std::cout << " genjet matched to hd3 eta: "
////					<< MCparticle->eta() << " phi: " << MCparticle->phi()  << " with dR: " << deltaR(*jet,*MCparticle) << "\n";
//					}
//				}
////			tracksNextToEle=0;
////			for(std::vector<reco::Track>::const_iterator tracks=trackCollection->begin(); tracks!=trackCollection->end(); ++tracks){
////				if (deltaR(*MCparticle,*tracks) < 0.5){
////					++tracksNextToEle;
////				}
////			}
////			std::cout << "tracksNextToEle " << tracksNextToEle << std::endl;
//
//
////				for(std::vector<reco::PFCandidate>::const_iterator PFparticle=pfCandCollection->begin(); PFparticle!=pfCandCollection->end(); ++PFparticle){
////					if (deltaR(*MCparticle,*PFparticle) < 0.002){
////						if (PFparticle->pt()>1.0)
////						std::cerr << "\tFound PFCand next to MC electron: " << PFparticle->particleId()
////								<< " at dist: " << deltaR(*MCparticle,*PFparticle)
////								<< " ("<< PFparticle->eta() << "," << PFparticle->phi() << ")"
////								<< std::endl;
////						PFCandObjs->Fill(PFparticle->particleId());
////						++numPFMatched;
////					}
////				}
//
//			}
////			std::cout << "Found hd2 == lepton jet particle\n";
//		}
//	}
//
//	std::vector<const pat::Jet*> PATJets;
//	if (jet1) PATJets.push_back(jet1);
//	if (jet2) PATJets.push_back(jet2);
//
//	for(std::vector<const pat::Jet*>::const_iterator jet=PATJets.begin(); jet!=PATJets.end(); ++jet){
//		double pt = (*jet)->pt();
//		SCall->Fill(pt);
//		std::cerr<< "Jet pt: " << pt << std::endl;
//		double ECALenergy(0), HCALenergy(0);
//		for (std::vector<reco::PFCandidatePtr>::const_iterator pfcand=(*jet)->getPFConstituents().begin(); pfcand!=(*jet)->getPFConstituents().end(); ++pfcand){
//			ECALenergy += (*pfcand)->ecalEnergy();
//			HCALenergy += (*pfcand)->hcalEnergy();
//		}
//////		for (unsigned int cons(0); cons < jet->getPFConstituents().size(); ++cons){
//////			reco::PFCandidatePtr pfcand = pfjet->getPFConstituent(cons);
//////			ECALenergy += pfcand->ecalEnergy();
//////			HCALenergy += pfcand->hcalEnergy();
//////		}
//		double EMF = (ECALenergy/(ECALenergy+HCALenergy));
//		std::cerr<< "EMF: " << EMF << std::endl;
//
////		unsigned int numEle(0);
////		for(edm::View<pat::Electron>::const_iterator patele=electrons->begin(); patele!=electrons->end(); ++patele){
////			if (deltaR(*(*jet),*patele) < 0.5){
////				PatElePt->Fill(patele->pt());
////				++numEle;
////			}
////		}
////		std::cerr << "Number of Pat electrons in PatJet: " << numEle << std::endl;
////		numElectrons->Fill(numEle);
////
//		unsigned int gsfNextToEle=0;
//		for(std::vector<reco::GsfElectron>::const_iterator ele=electronCollection->begin(); ele!=electronCollection->end(); ++ele){
//			if (deltaR(*(*jet),*ele) < 0.5){
//				++gsfNextToEle;
//				ElePt->Fill(ele->pt());
//			}
//		}
//		numGSFElectrons->Fill(gsfNextToEle);
//		std::cerr << "gsfNextToEle " << gsfNextToEle << std::endl;
//
//
////		unsigned int numMCEle(0);
////		std::cout << "In Jet\n";
////		for(std::vector<reco::GenParticle>::const_iterator MCparticle=genParticleCollection->begin(); MCparticle!=genParticleCollection->end(); ++MCparticle){
//////			if ( MCparticle->pdgId() == 11 || MCparticle->pdgId() == -11){	//e
////			if (deltaR(*(*jet),*MCparticle) < 0.5){
//////					MCElePt->Fill(MCparticle->pt());
//////					++numMCEle;
////				std::cout << "\tMC Particle: " << MCparticle->pdgId()
////						<< "\tat: " << MCparticle->p4().X() << " "
////									<< MCparticle->p4().Y() << " "
////									<< MCparticle->p4().Z()
////						<< "\tpt: " << MCparticle->pt()
////						<< std::endl;
//////				MCParticles->Fill(MCparticle->pdgId());
////			}
//////			}
////		}
////		std::cerr << "Number of MC electrons in PatJet = " << numMCEle << std::endl;
////		numMCElectrons->Fill(numMCEle);
//
//
////		unsigned int numPFEle(0);
////		for(std::vector<reco::PFCandidate>::const_iterator PFparticle=pfCandCollection->begin(); PFparticle!=pfCandCollection->end(); ++PFparticle){
////			if (deltaR(*(*jet),*PFparticle) < 0.5){
//////				PFElePt->Fill(PFparticle->pt());
//////				++numPFEle;
//////				if (PFparticle->particleId() == 0)	PFCandObjs->Fill("X");
////				std::cout << "\tPF Particle: " << PFparticle->particleId()
////								<< "\tat: " << PFparticle->p4().X() << " "
////											<< PFparticle->p4().Y() << " "
////											<< PFparticle->p4().Z()
////								<< "\tpt: " << PFparticle->pt()
////						<< std::endl;
////			}
////		}
////		std::cerr << "Number of PF electrons in PatJet = " << numPFEle << std::endl;
////		NumPFCands->Fill(numPFEle);
//
//		double ntrks = ((*jet)->associatedTracks().size());
//////					fout  << "1\t" << EMF  << "\t" << ntrks << "\t" << numEle << "\t" << pt << "\n"; //EMF too high needs to be down at 90
//////				hist->Fill(EMF);
////		std::cerr<< "ntrks: " << ntrks << std::endl;
//		if (EMF>0.9 && ntrks>2 && gsfNextToEle>0){
////					SCpass->Fill(pt);
//			if (1==gsfNextToEle && ntrks/pt<0.45){
////						SCall->Fill(pt);
//				if ( (pt<70 && (*jet)->phiphiMoment()>0.005) || (*jet)->phiphiMoment()>0.005-(0.0001*(pt-70)) )
////							fout  << "1\t" << EMF  << "\t" << ntrks << "\t" << numEle << "\t" << pt << "\n";
//					SCpass->Fill(pt);  //SC
//			}
//			else if(gsfNextToEle>1) SCpass->Fill(pt);
//		}
//	}
//
//
////	size_t nJets=0;
////	for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
////		if(jet->pt()>50){
////			++nJets;
////		}
////		// uncomment the following line to fill the jetTowers_ histogram
////		// jetTowers_->Fill(jet->getCaloConstituents().size());
////	}
//
////	histContainer_["jets"]->Fill(nJets);
////	// do something similar for the other candidates
//////	histContainer_["photons"]->Fill(photons->size() );
////	histContainer_["elecs" ]->Fill(electrons->size());
////	histContainer_["muons"]->Fill(muons->size() );
////	histContainer_["taus" ]->Fill(taus->size()  );
////	histContainer_["met"  ]->Fill(mets->empty() ? 0 : (*mets)[0].et());
	++EventNum;
}



void VSelector::endJob(){
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(VSelector);
