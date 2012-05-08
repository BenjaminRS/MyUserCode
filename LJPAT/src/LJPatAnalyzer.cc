#include <map>
#include <string>
#include <fstream>

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

class PatBasicAnalyzer : public edm::EDAnalyzer {
	public:
		explicit PatBasicAnalyzer(const edm::ParameterSet&);
		~PatBasicAnalyzer();
	private:
		virtual void beginJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;
		// simple map to contain all histograms; histograms are booked in the beginJob() method
		std::map<std::string,TH1F*> histContainer_;
		// plot number of towers per jet
		TH1F* jetTowers_;
//		edm::InputTag photonSrc_;
		edm::InputTag elecSrc_;
//		edm::InputTag muonSrc_;
//		edm::InputTag tauSrc_;
		edm::InputTag jetSrc_;
//		edm::InputTag metSrc_;
		TH1F* SigmaIetaIetaB;
		TH1F* SigmaEtaEtaB;
		TH1F* SigmaIetaIetaE;
		TH1F* SigmaEtaEtaE;
		TH1F* SCpass;
		TH1F* SCall;
		TH1F* numElectrons;
		TH1F* numMCElectrons;
		TH1F* numGSFElectrons;
		TH1F* MCElePt;
		TH1F* PatElePt;
		TH1F* ElePt;
		TH1F* NumPFCands;
		TH1F* PFElePt;
		TH1F* PFCandObjs;
		TH1F* MCParticles;
		TH1F* ETA;
		TH1F* PT;
		TH1F* JetPT;
		unsigned int EventNum;
		unsigned int numMCele;
		unsigned int numPFMatched;
		unsigned int tracksNextToEle;
		unsigned int numMCparts;
		unsigned int nummatchedJets;
		unsigned int numMultiGenMatch;
		unsigned int numMultiRecoMatch;
		unsigned int numMultiPatMatch;
		unsigned int numRecoJetMatched;
		unsigned int numPATJetMatched;
		unsigned int numGenJetMatched;
		std::string filename_;
		unsigned int step_;
		unsigned int mass_;
};

PatBasicAnalyzer::PatBasicAnalyzer(const edm::ParameterSet& iConfig):
	histContainer_(),
	filename_(iConfig.getUntrackedParameter<std::string>("filename")),
	mass_(iConfig.getUntrackedParameter<unsigned int>("mass")),
	step_(iConfig.getUntrackedParameter<unsigned int>("step"))
//	photonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("photonSrc")),
//	elecSrc_(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc")),
//	muonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc")),
//	tauSrc_(iConfig.getUntrackedParameter<edm::InputTag>("tauSrc" )),
//	jetSrc_(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc" ))
//	metSrc_(iConfig.getUntrackedParameter<edm::InputTag>("metSrc" ))
{

}

PatBasicAnalyzer::~PatBasicAnalyzer(){}

void PatBasicAnalyzer::beginJob(){
	edm::Service<TFileService> fs;
  // book histograms: uncomment the following line to book the jetTowers_ histogram
  //jetTowers_= fs->make<TH1F>("jetTowers", "towers per jet",   90, 0,  90);
//  histContainer_["photons"]=fs->make<TH1F>("photons", "photon multiplicity",   10, 0,  10);
	histContainer_["elecs"  ]=fs->make<TH1F>("elecs",   "electron multiplicity", 10, 0,  10);
	histContainer_["muons"  ]=fs->make<TH1F>("muons",   "muon multiplicity",     10, 0,  10);
	histContainer_["taus"   ]=fs->make<TH1F>("taus",    "tau multiplicity",      10, 0,  10);
	histContainer_["jets"   ]=fs->make<TH1F>("jets",    "jet multiplicity",      10, 0,  10);
	histContainer_["met"    ]=fs->make<TH1F>("met",     "missing E_{T}",         20, 0, 100);
	SigmaIetaIetaB=fs->make<TH1F>("SigmaIetaIetaB","SigmaIetaIeta Barrel",	200,0.0,0.05);
	SigmaEtaEtaB=fs->make<TH1F>("SigmaEtaEtaB","SigmaEtaEta Barrel",		200,0.0,0.05);
	SigmaIetaIetaE=fs->make<TH1F>("SigmaIetaIetaE","SigmaIetaIeta Endcap",	200,0.0,0.05);
	SigmaEtaEtaE=fs->make<TH1F>("SigmaEtaEtaE","SigmaEtaEta Endcap",		200,0.0,0.05);
	SCpass=fs->make<TH1F>("SCPass","LJets which pass Scheme C",	100,0.0,500.0);
	SCall=fs->make<TH1F>("SCAll","All LJets",					100,0.0,500.0);
	numElectrons=fs->make<TH1F>("numElectrons","Number of PAT::Electrons in LJ", 20,0.0,20.0);
	numMCElectrons=fs->make<TH1F>("numMCElectrons","Number of MC Electrons in LJ", 20,0.0,20.0);
	numGSFElectrons=fs->make<TH1F>("numGSFElectrons","Number of Gsf Electrons in LJ", 20,0.0,20.0);
	MCElePt=fs->make<TH1F>("MCElectronsPt","Pt of MC Electrons in LJ", 100,0.0,100.0);
	PatElePt=fs->make<TH1F>("PatElectronsPt","Pt of Pat Electrons in LJ", 100,0.0,100.0);
	ElePt=fs->make<TH1F>("GsfElectronsPt","Pt of GsfElectrons in LJ", 100,0.0,100.0);
	NumPFCands=fs->make<TH1F>("numPFCands","Number of PFCandidates in LJ", 50,0.0,50.0);
	PFElePt=fs->make<TH1F>("PFElectronsPt","Pt of PF Electrons (RECO) in LJ", 100,0.0,100.0);
	PFCandObjs=fs->make<TH1F>("PFCandsInLJ","ID of PFCandidates in LJ", 10,0.0,10.0);
	MCParticles=fs->make<TH1F>("MCParticlesInLJ","ID of MCParticles in LJ", 4000000,0.0,4000000.0);
	ETA=fs->make<TH1F>("dR","dR", 500,0.0,5.0);
	PT=fs->make<TH1F>("pt","MC particle Pt", 200,0.0,400.0);
	JetPT=fs->make<TH1F>("Jetpt","Jet Pt", 200,0.0,400.0);
	numMCele=0;
	numPFMatched=0;
	tracksNextToEle=0;
	numMCparts=0;
	nummatchedJets=0;
	numRecoJetMatched=0;
	numPATJetMatched=0;
	numMultiGenMatch=numMultiRecoMatch=numMultiPatMatch=0;
	numGenJetMatched=0;
}

void PatBasicAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
//	edm::Handle<edm::View<pat::Electron> > electrons;	iEvent.getByLabel(elecSrc_,electrons);
//	edm::Handle<edm::View<pat::Muon> > muons;	iEvent.getByLabel(muonSrc_,muons);
//	edm::Handle<edm::View<pat::Tau> > taus;	iEvent.getByLabel(tauSrc_,taus);
//	edm::Handle<edm::View<pat::Jet> > jets;	iEvent.getByLabel(jetSrc_,jets);
//	edm::Handle<edm::View<pat::MET> > mets;	iEvent.getByLabel(metSrc_,mets);
//	edm::Handle<edm::View<pat::Photon> > photons; iEvent.getByLabel(photonSrc_,photons);
	edm::Handle<std::vector<reco::GenParticle> > genParticleCollection; iEvent.getByLabel("genParticles", genParticleCollection);
//	edm::Handle<std::vector<reco::PFCandidate> > pfCandCollection; iEvent.getByLabel("pfAllElectronsPFlow", pfCandCollection);
//	edm::Handle<std::vector<reco::PFCandidate> > pfCandCollection; iEvent.getByLabel("particleFlow","electrons",pfCandCollection );
	edm::Handle<std::vector<reco::PFCandidate> > pfCandCollection; iEvent.getByLabel("particleFlow","",pfCandCollection );
//	edm::Handle<std::vector<reco::Track> > trackCollection; iEvent.getByLabel("generalTracks","",trackCollection );
//	edm::Handle<std::vector<reco::GsfElectron> > electronCollection; iEvent.getByLabel("gsfElectrons","",electronCollection );

	edm::Handle<std::vector<reco::GsfElectron> > electronCollection;	iEvent.getByLabel("gsfElectrons", electronCollection);
	edm::Handle<std::vector<reco::GenJet> > jetCollection;	iEvent.getByLabel("ak5GenJets", jetCollection);
	edm::Handle<std::vector<reco::PFJet> > pfJetCollection;	iEvent.getByLabel("ak5PFJets", pfJetCollection);
	edm::Handle<edm::View<pat::Jet> > PATjets;	iEvent.getByLabel("selectedPatJetsPFlow",PATjets);

	std::ofstream fout(filename_.c_str(), std::ios_base::app);
//	std::ofstream fout("test.csv", std::ios_base::app);

	std::cerr << "Event["<< iEvent.id() << "]\n";
//	std::cout << "Number of electrons in collection:\t" << electrons->size() << std::endl;
//	std::cout << "Number of PFCandidates in collection:\t" << pfCandCollection->size() << std::endl;
//	for (std::vector<reco::PFCandidate>::const_iterator PFCand=pfCandCollection->begin(); PFCand!=pfCandCollection->end(); ++PFCand){
//		enum reco::PFCandidate::ParticleType type=PFCand->particleId();
////		std::cout << "\t" << type <<std::endl;
//		NumPFCands->Fill(type);
//	}
//	NumPFCands->Fill(pfCandCollection->size());

//	##############################################
//	############### Signal START #################
//	##############################################
//	const pat::Jet* jet1=NULL;
//	const pat::Jet* jet2=NULL;

//	for(unsigned int i(0); i<genParticleCollection->size(); ++i){
//		const reco::GenParticle* MCparticle = &genParticleCollection->at(i);
//		if ( MCparticle->pdgId() == 3000006){	//hd2
//			if (MCparticle->mother()->pdgId() == 25){  // Protect from loopback
//
////				std::cout << "eta: " << MCparticle->eta()
////						<< " phi: " << MCparticle->phi()
////						<< " pt: " << MCparticle->pt()
////						<< " px: " << MCparticle->px()
////						<< " py: " << MCparticle->py()
////						<< " pz: " << MCparticle->pz()
////						<< " decays to:" << std::endl;
//				for(unsigned int j(0); j<MCparticle->numberOfDaughters(); ++j){
////					std::cout << "\t" << MCparticle->daughterRef(j)->pdgId()
////								<< " eta: " << MCparticle->daughterRef(j)->eta()
////								<< " phi: " << MCparticle->daughterRef(j)->phi()
////								<< " pt: " << MCparticle->daughterRef(j)->pt()
////								<< " px: " << MCparticle->daughterRef(j)->px()
////								<< " py: " << MCparticle->daughterRef(j)->py()
////								<< " pz: " << MCparticle->daughterRef(j)->pz()
////								<< std::endl;
//					for(unsigned int k(0); k<MCparticle->daughterRef(j)->numberOfDaughters(); ++k){
////						std::cout << "\t\t" << MCparticle->daughterRef(j)->daughterRef(k)->pdgId()
////									<< " eta: " << MCparticle->daughterRef(j)->daughterRef(k)->eta()
////									<< " phi: " << MCparticle->daughterRef(j)->daughterRef(k)->phi()
////									<< " pt: " << MCparticle->daughterRef(j)->daughterRef(k)->pt()
////									<< " px: " << MCparticle->daughterRef(j)->daughterRef(k)->px()
////									<< " py: " << MCparticle->daughterRef(j)->daughterRef(k)->py()
////									<< " pz: " << MCparticle->daughterRef(j)->daughterRef(k)->pz()
////									<< std::endl;
//						if (MCparticle->daughterRef(j)->daughterRef(k)->pdgId() == 3000001){
//							ETA->Fill(deltaR(*MCparticle,*MCparticle->daughterRef(j)->daughterRef(k)));
//						}
//					}
//				}
//			}
//		}
//	}

	for(std::vector<reco::GenParticle>::const_iterator MCparticle=genParticleCollection->begin(); MCparticle!=genParticleCollection->end(); ++MCparticle){
		if ( MCparticle->pdgId() == 3000006){	//hd2
			if (MCparticle->mother()->pdgId() == 25){  // Protect from loopback
				std::cerr << "MCparticle->pdgId(): " << MCparticle->pdgId() << " ("<< MCparticle->eta() << "," << MCparticle->phi() << ") with pt: " << MCparticle->pt() << std::endl;
				fout << iEvent.id().event() << "\t" << mass_ << "\t" << step_ << "\t0\t"<< MCparticle->eta() << "\t" << MCparticle->phi() << "\t" << MCparticle->pt() << std::endl;
//				for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
				++numMCparts;
				ETA->Fill(MCparticle->eta());
				PT->Fill(MCparticle->pt());
//				std::cerr << "MCparticle->pt(): " << MCparticle->pt() << std::endl;
				int x=numGenJetMatched;
				for(std::vector<reco::GenJet>::const_iterator genjet=jetCollection->begin(); genjet!=jetCollection->end(); ++genjet){
					if (deltaR(*genjet,*MCparticle)<0.3){
						std::cerr << "\tMatched GenJet: (" << genjet->eta() << "," << genjet->phi() << ") with pt: " << genjet->pt() << std::endl;
						++numGenJetMatched;
						fout << iEvent.id().event() << "\t" << mass_ << "\t" << step_ << "\t1\t"<< genjet->eta() << "\t" << genjet->phi() << "\t" << genjet->pt() << std::endl;
					}
				}
				int y=numGenJetMatched;
				if (y-x > 1){
					std::cerr << "\tMatched more than one GenJet!\n";
					++numMultiGenMatch;
				}
				x=numRecoJetMatched;
				for(std::vector<reco::PFJet>::const_iterator jet=pfJetCollection->begin(); jet!=pfJetCollection->end(); ++jet){
					if (deltaR(*jet,*MCparticle)<0.3){
						++numRecoJetMatched;
						std::cerr << "\tMatched RecoJet: (" << jet->eta() << "," << jet->phi() << ") with pt: " << jet->pt() << std::endl;
						fout << iEvent.id().event() << "\t" << mass_ << "\t" << step_ << "\t2\t"<< jet->eta() << "\t" << jet->phi() << "\t" << jet->pt() << std::endl;
					}
				}
				y=numRecoJetMatched;
				if (y-x > 1){
					std::cerr << "\tMatched more than one RecoJet!\n";
					++numMultiRecoMatch;
				}
				x=numPATJetMatched;
				for(edm::View<pat::Jet>::const_iterator pjet=PATjets->begin(); pjet!=PATjets->end(); ++pjet){
					if (deltaR(*pjet,*MCparticle)<0.3){
						++numPATJetMatched;
						std::cerr << "\tMatched PATJet: (" << pjet->eta() << "," << pjet->phi() << ") with pt: " << pjet->pt() << std::endl;
						fout << iEvent.id().event() << "\t" << mass_ << "\t" << step_ << "\t3\t"<< pjet->eta() << "\t" << pjet->phi() << "\t" << pjet->pt() << std::endl;
					}
				}
				y=numPATJetMatched;
				if (y-x > 1){
					std::cerr << "\tMatched more than one PatJet!\n";
					++numMultiPatMatch;
				}

//				for( size_t k=0; k<jets->size(); ++k ) {
//					const pat::Jet* jet = &(jets->at(k));
//					if (deltaR(*jet,*MCparticle)<0.3){	//Match PatJet to the GenParticle = hd2
//						++nummatchedJets;
//						std::cerr << "\tMatched jet: (" << jet->eta() << "," << jet->phi() << ") with pt: " << jet->pt() << std::endl;
//						JetPT->Fill(jet->pt());
////						std::cerr << "\t\tpatjet matched to hd2; jet1= " << jet1 << " jet2= " << jet2 << std::endl;
//						if (!jet1) jet1=jet;
//						else jet2=jet;
////						std::cerr << "\tafter: jet1= " << jet1 << " jet2= " << jet2 << std::endl;
//					}
//				}

			}
		}
	}


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
////		for (unsigned int cons(0); cons < jet->getPFConstituents().size(); ++cons){
////			reco::PFCandidatePtr pfcand = pfjet->getPFConstituent(cons);
////			ECALenergy += pfcand->ecalEnergy();
////			HCALenergy += pfcand->hcalEnergy();
////		}
//		double EMF = (ECALenergy/(ECALenergy+HCALenergy));
//		std::cerr<< "EMF: " << EMF << std::endl;
//
//		edm::View<pat::Electron>::const_iterator ThePatEle;
//		unsigned int numEle(0);
//		for(edm::View<pat::Electron>::const_iterator patele=electrons->begin(); patele!=electrons->end(); ++patele){
//			if (deltaR(*(*jet),*patele) < 0.5){
//				PatElePt->Fill(patele->pt());
//				++numEle;
//				ThePatEle=patele;
//			}
//		}
//		std::cerr << "Number of Pat electrons in PatJet: " << numEle << std::endl;
//		numElectrons->Fill(numEle);
////
////		unsigned int gsfNextToEle=0;
////		for(std::vector<reco::GsfElectron>::const_iterator ele=electronCollection->begin(); ele!=electronCollection->end(); ++ele){
////			if (deltaR(*(*jet),*ele) < 0.5){
////				++gsfNextToEle;
////				ElePt->Fill(ele->pt());
////			}
////		}
////		numGSFElectrons->Fill(gsfNextToEle);
////		std::cerr << "gsfNextToEle " << gsfNextToEle << std::endl;
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
////		if (EMF>0.9 && ntrks>2 && gsfNextToEle>0){
//		if (EMF>0.9 && ntrks>2 && numEle>0){
////					SCpass->Fill(pt);
////			if (1==gsfNextToEle && ntrks/pt<0.45){
//			if (1==numEle && ntrks/pt<0.45){
////						SCall->Fill(pt);
////				if ( (pt<70 && (*jet)->phiphiMoment()>0.005) || (*jet)->phiphiMoment()>0.005-(0.0001*(pt-70)) )		SCpass->Fill(pt);  //SC
//				if ( ThePatEle->dr04EcalRecHitSumEt()/ThePatEle->p4().Pt()>0.04 )		SCpass->Fill(pt);  //SC
//			}
////			else if(gsfNextToEle>1) SCpass->Fill(pt);
//			else if(numEle>1) SCpass->Fill(pt);
//		}
//	}
//	##############################################
//	############### Signal END ###################
//	##############################################

//	##############################################
//	########### Jet BACKGROUND START #############
//	##############################################
//	for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
//		double pt = jet->pt();
//		SCall->Fill(pt);
//		std::cerr<< "Jet pt: " << pt << std::endl;
//		double ECALenergy(0), HCALenergy(0);
//		for (std::vector<reco::PFCandidatePtr>::const_iterator pfcand=jet->getPFConstituents().begin(); pfcand!=jet->getPFConstituents().end(); ++pfcand){
//			ECALenergy += (*pfcand)->ecalEnergy();
//			HCALenergy += (*pfcand)->hcalEnergy();
//		}
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
//		for(edm::View<pat::Electron>::const_iterator ele=electrons->begin(); ele!=electrons->end(); ++ele){
//			if (deltaR(*jet,*ele) < 0.7){
//				++gsfNextToEle;
//				ElePt->Fill(ele->pt());
//			}
//		}
//		numGSFElectrons->Fill(gsfNextToEle);
//		std::cerr << "gsfNextToEle " << gsfNextToEle << std::endl;
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
//		double ntrks = (jet->associatedTracks().size());
//////					fout  << "1\t" << EMF  << "\t" << ntrks << "\t" << numEle << "\t" << pt << "\n"; //EMF too high needs to be down at 90
//////				hist->Fill(EMF);
////		std::cerr<< "ntrks: " << ntrks << std::endl;
//		if (EMF>0.9 && ntrks>2 && gsfNextToEle>0){
////					SCpass->Fill(pt);
//			if (1==gsfNextToEle && ntrks/pt<0.45){
////						SCall->Fill(pt);
//				if ( (pt<70 && jet->phiphiMoment()>0.005) || jet->phiphiMoment()>0.005-(0.0001*(pt-70)) )
////							fout  << "1\t" << EMF  << "\t" << ntrks << "\t" << numEle << "\t" << pt << "\n";
//					SCpass->Fill(pt);  //SC
//			}
//			else if(gsfNextToEle>1) SCpass->Fill(pt);
//		}
//	}
//	##############################################
//	########### Jet BACKGROUND END ###############
//	##############################################

//	##############################################
//	########### ELE BACKGROUND START #############
//	##############################################
//	for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
//		for(std::vector<reco::GenParticle>::const_iterator MCparticle=genParticleCollection->begin(); MCparticle!=genParticleCollection->end(); ++MCparticle){
//			if ( MCparticle->pdgId() == 11 || MCparticle->pdgId() == -11){	//e
//				if (MCparticle->mother()->pdgId() == 23){  // From Z
//					if (deltaR(*jet,*MCparticle) < 0.3){
//						std::cerr << "\tFound matched jet: " << " ("<< jet->eta() << "," << jet->phi() << ")"
//								 << " to MC electron: " << " ("<< MCparticle->eta() << "," << MCparticle->phi() << ")"
//								<< std::endl;
//						double pt = jet->pt();
//						SCall->Fill(pt);
//						std::cerr<< "Jet pt: " << pt << std::endl;
//						double ECALenergy(0), HCALenergy(0);
//						for (std::vector<reco::PFCandidatePtr>::const_iterator pfcand=jet->getPFConstituents().begin(); pfcand!=jet->getPFConstituents().end(); ++pfcand){
//							ECALenergy += (*pfcand)->ecalEnergy();
//							HCALenergy += (*pfcand)->hcalEnergy();
//						}
//						double EMF = (ECALenergy/(ECALenergy+HCALenergy));
//						unsigned int gsfNextToEle=0;
//						for(std::vector<reco::GsfElectron>::const_iterator electron=electronCollection->begin(); electron!=electronCollection->end(); ++electron){
//							if (deltaR(*jet,*electron) < 0.5){
//								++gsfNextToEle;
//								ElePt->Fill(electron->pt());
//							}
//						}
//						double ntrks = (jet->associatedTracks().size());
//						if (EMF>0.9 && ntrks>2 && gsfNextToEle>0){
//							if (1==gsfNextToEle && ntrks/pt<0.45){
//								if ( (pt<70 && jet->phiphiMoment()>0.005) || jet->phiphiMoment()>0.005-(0.0001*(pt-70)) )
//									SCpass->Fill(pt);  //SC
//							}
//							else if(gsfNextToEle>1){
//								SCpass->Fill(pt);
//								std::cout << "WARNING: I heard you liked electrons, so I put and electron in your electron\n";
//							}
//						}
//
//					}
//				}
//			}
//		}
//	}
//
//	for(std::vector<reco::GenParticle>::const_iterator MCparticle=genParticleCollection->begin(); MCparticle!=genParticleCollection->end(); ++MCparticle){
//		if ( MCparticle->pdgId() == 11 || MCparticle->pdgId() == -11){	//e
//			if (MCparticle->mother()->pdgId() == 23){  // From Z
////				++numMCele;
//				std::cerr << "MCparticle->pdgId(): " << MCparticle->pdgId() << " ("<< MCparticle->eta() << "," << MCparticle->phi() << ")" << std::endl;
//				for(std::vector<reco::GsfElectron>::const_iterator ele=electronCollection->begin(); ele!=electronCollection->end(); ++ele){
//					if (deltaR(*MCparticle,*ele) < 0.5){
//						std::cerr << "\tFound matched electron: " << " ("<< ele->eta() << "," << ele->phi() << ")" << std::endl;
//						for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
//							if (deltaR(*jet,*ele) < 0.5){
//								std::cerr << "\tFound matched jet: " << " ("<< jet->eta() << "," << jet->phi() << ")" << std::endl;
//								double pt = jet->pt();
//								SCall->Fill(pt);
//								std::cerr<< "Jet pt: " << pt << std::endl;
//								double ECALenergy(0), HCALenergy(0);
//								for (std::vector<reco::PFCandidatePtr>::const_iterator pfcand=jet->getPFConstituents().begin(); pfcand!=jet->getPFConstituents().end(); ++pfcand){
//									ECALenergy += (*pfcand)->ecalEnergy();
//									HCALenergy += (*pfcand)->hcalEnergy();
//								}
//								double EMF = (ECALenergy/(ECALenergy+HCALenergy));
//								unsigned int gsfNextToEle=0;
//								for(std::vector<reco::GsfElectron>::const_iterator electron=electronCollection->begin(); electron!=electronCollection->end(); ++electron){
//									if (deltaR(*jet,*electron) < 0.5){
//										++gsfNextToEle;
//										ElePt->Fill(electron->pt());
//									}
//								}
//								double ntrks = (jet->associatedTracks().size());
//								if (EMF>0.9 && ntrks>2 && gsfNextToEle>0){
//									if (1==gsfNextToEle && ntrks/pt<0.45){
//										if ( (pt<70 && jet->phiphiMoment()>0.005) || jet->phiphiMoment()>0.005-(0.0001*(pt-70)) )
//				//							fout  << "1\t" << EMF  << "\t" << ntrks << "\t" << numEle << "\t" << pt << "\n";
//											SCpass->Fill(pt);  //SC
//									}
//									else if(gsfNextToEle>1) SCpass->Fill(pt);
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//	##############################################
//	########### ELE BACKGROUND END ###############
//	##############################################


}



void PatBasicAnalyzer::endJob(){
//	std::cout << "numMCele: " << numMCele << std::endl;
//	std::cout << "numPFMatched: " << numPFMatched << std::endl;
	std::cerr << "numMCparts: " << numMCparts << std::endl;
	std::cerr << "numGenJetMatched: " << numGenJetMatched << std::endl;
	std::cerr << "numRecoJetMatched: " << numRecoJetMatched << std::endl;
	std::cerr << "numPATJetMatched: " << numPATJetMatched << std::endl;
//	std::cout << "nummatchedJets: " << nummatchedJets << std::endl;
	std::cerr << "numMultiGenMatch: " << numMultiGenMatch << std::endl;
	std::cerr << "numMultiRecoMatch: " << numMultiRecoMatch << std::endl;
	std::cerr << "numMultiPatMatch: " << numMultiPatMatch << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PatBasicAnalyzer);
