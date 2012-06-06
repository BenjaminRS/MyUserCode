// Original Author:  benjamin radburn-smith
//         Created:  Wed Jun 22 09:48:57 CDT 2011
#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include <fstream>
//#include "PhysicsTools/SelectorUtils/interface/SimpleCutBasedElectronIDSelectionFunctor.h"

class LJ : public edm::EDAnalyzer {
	public:
		explicit LJ(const edm::ParameterSet&);
		~LJ();
	private:
		virtual void beginJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;
		virtual void RunEJSelection(std::vector<const reco::PFJet*>&);
		virtual void RunEJSelection(std::vector<const reco::CaloJet*>&);
		virtual void FindEJs();
		TH1F* SCpass;
		TH1F* SCall;
		TH1F* PassEMF;
		TH1F* PassNtrks;
		TH1F* PassNEle;
		TH1F* MC;
		TH2F* JetPts;
		int numRecoJetMatched;
		int numGenJetMatched;
		int numCaloRecoJetMatched;
		int numMCparts;
		TH1F* DeltaR;
		TH1F* HiggsMass, *HiggsMass2;
		TH1F* PFJetPtComp, *GenJetPtComp;
		std::vector<double> PFJetPtCompStore;
		std::vector<double> GenJetPtCompStore;
		std::vector<double> CaloJetPtCompStore;
		TH1F* NumGoodEJ;
		unsigned int numOfGoodElectronJets;
		const reco::PFJet* PFJet1;
		const reco::PFJet* PFJet2;
		const reco::GenJet* GenJet1;
		const reco::GenJet* GenJet2;
		const reco::CaloJet* CaloJet1;
		const reco::CaloJet* CaloJet2;
		edm::Handle<std::vector<reco::GenParticle> > genParticleCollection;
		edm::Handle<std::vector<reco::GsfElectron> > electronCollection;
		edm::Handle<std::vector<reco::GenJet> > genjetCollection;
		edm::Handle<std::vector<reco::PFJet> > pfJetCollection;
//		edm::Handle<std::vector<reco::Muon> > muonCollection;
		edm::Handle<std::vector<reco::CaloJet> > caloJetCollection;
		edm::Handle < reco::JetTracksAssociationCollection > assocCaloJetTracks;
};

LJ::LJ(const edm::ParameterSet& iConfig){
	edm::Service<TFileService> fs;
	SCpass=fs->make<TH1F>("SCPass","Scheme C Pass - ak7",500,0.0,500.0);
	SCall=fs->make<TH1F>("SCAll","Scheme C All - ak7",500,0.0,500.0);
	PassEMF=fs->make<TH1F>("PassEMF","PassEMF - ak7",500,0.0,500.0);
	PassNtrks=fs->make<TH1F>("PassNtrks","PassNtrks - ak7",500,0.0,500.0);
	PassNEle=fs->make<TH1F>("PassNEle","PassNEle - ak7",500,0.0,500.0);
	MC=fs->make<TH1F>("MC","MCJets",100,0.0,500.0);
	PFJetPtComp=fs->make<TH1F>("PFJetPtComp","ak7 PFJet Pt Compared to Hd2",800,-2.0,2.0);
	GenJetPtComp=fs->make<TH1F>("GenJetPtComp","ak7 GenJet Pt Compared to Hd2",800,-2.0,2.0);
	HiggsMass=fs->make<TH1F>("HiggsMass","Higgs Mass from ak7 GenJet",1500,0.0,1500.0);
	HiggsMass2=fs->make<TH1F>("HiggsMass2","Higgs Mass from ak7 PFJet",1500,0.0,1500.0);
	JetPts=fs->make<TH2F>("JetPts","ak7 PFJet Vs ak7 GenJet Pts",1000,0.0,1000.0,1000,0.0,1000.0);
	NumGoodEJ=fs->make<TH1F>("NumGoodEJ","Number of EJ per Event which Pass Scheme C",4,0.0,4.0);
	numRecoJetMatched=0;
	numGenJetMatched=0;
	numCaloRecoJetMatched=0;
	numMCparts=0;
	DeltaR=fs->make<TH1F>("dR","dR", 500,0.0,5.0);
	numOfGoodElectronJets=0;
}

LJ::~LJ(){}

void LJ::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
//	std::ofstream fout("1110ElePropertiesZeeBkg.csv", std::ios_base::app);
	using namespace edm;
//	std::ofstream fout("040811Values.csv", std::ios_base::app);
	std::cerr << "Event["<< iEvent.id() << "]\n";
	iEvent.getByLabel("gsfElectrons", electronCollection);
	iEvent.getByLabel("ak7GenJets", genjetCollection);
//	iEvent.getByLabel("ak7GenJetsWithHd0", genjetCollection);
	iEvent.getByLabel("genParticles", genParticleCollection);
	iEvent.getByLabel("ak7PFJets", pfJetCollection);
//	iEvent.getByLabel("muons", muonCollection);
	iEvent.getByLabel("ak7CaloJets", caloJetCollection);
	iEvent.getByLabel("ak7JetTracksAssociatorAtVertex", assocCaloJetTracks);


//	SimpleCutBasedElectronIDSelectionFunctor patSele95(SimpleCutBasedElectronIDSelectionFunctor::relIso95);

// ############################################## Zee BACKGROUND START
//	std::cout << "# of electrons: " << electronCollection->size() <<std::endl;
//	if (electronCollection->size()>2){
//		const reco::GsfElectron* reconele = &(electronCollection->at(0));
//		const reco::GsfElectron* reconele2 = &(electronCollection->at(1));
//		TLorentzVector e1(reconele->px(), reconele->py(), reconele->pz(), reconele->energy()), e2(reconele2->px(),reconele2->py(),reconele2->pz(),reconele2->energy());
//		TLorentzVector Zmass = e1 + e2;
//
//		if (deltaR(*reconele,*reconele2) > 2.8 && Zmass.M() < 100 && Zmass.M() > 80){   // likely a Z event
//			for (unsigned int pos(0); pos<electronCollection->size(); ++pos){
//				const reco::GsfElectron* reconele = &(electronCollection->at(pos));
//				All->Fill(reconele->pt());
//				if (reconele->isEE()){
//					if (reconele->sigmaIetaIeta() < 0.01 &&
//							reconele->deltaEtaSuperClusterTrackAtVtx() < 0.004 &&
//							reconele->deltaPhiSuperClusterTrackAtVtx() < 0.06 &&
//							reconele->hadronicOverEm() < 0.04) Pass->Fill(reconele->pt());
//				}
//				if (reconele->isEB()){
//					if (reconele->sigmaIetaIeta() < 0.03 &&
//							reconele->deltaEtaSuperClusterTrackAtVtx() < 0.007 &&
//							reconele->deltaPhiSuperClusterTrackAtVtx() < 0.03 &&
//							reconele->hadronicOverEm() < 0.15) Pass->Fill(reconele->pt());
//				}
//
//
//				All->Fill(reconele2->pt());
//				if (reconele2->isEE()){
//					if (reconele2->sigmaIetaIeta() < 0.01 &&
//							reconele2->deltaEtaSuperClusterTrackAtVtx() < 0.004 &&
//							reconele2->deltaPhiSuperClusterTrackAtVtx() < 0.06 &&
//							reconele2->hadronicOverEm() < 0.04) Pass->Fill(reconele2->pt());
//				}
//				if (reconele2->isEB()){
//					if (reconele2->sigmaIetaIeta() < 0.03 &&
//							reconele2->deltaEtaSuperClusterTrackAtVtx() < 0.007 &&
//							reconele2->deltaPhiSuperClusterTrackAtVtx() < 0.03 &&
//							reconele2->hadronicOverEm() < 0.15) Pass->Fill(reconele2->pt());
//				}
//			}
//		}
//	}
//		for (unsigned int x(0); x<electronCollection->size(); ++x){
//			if (electronCollection->at(x).pt()>10){
//				for (unsigned int y(x+1); y<electronCollection->size(); ++y){
//					if (electronCollection->at(y).pt()>10 && x!=y){
//                				const reco::GsfElectron* reconele = &(electronCollection->at(x));
//				                const reco::GsfElectron* reconele2 = &(electronCollection->at(y));
//				                TLorentzVector e1(reconele->px(), reconele->py(), reconele->pz(), reconele->energy()), e2(reconele2->px(),reconele2->py(),reconele2->pz(),reconele2->energy());
//				                TLorentzVector Zmass = e1 + e2;
//				                if (deltaR(*reconele,*reconele2) > 2.8 && Zmass.M() < 100 && Zmass.M() > 80){
////				                	bool pass = patSele95(*reconele);
////				                	std::cout << "pass?: " << pass << std::endl;
//
////							pt->Fill(reconele->pt());
////							pt->Fill(reconele2->pt());
///*
//							OneByFive->Fill(reconele->e1x5());
//                                			TwoByFive->Fill(reconele->e2x5Max());
//                                			FiveByFive->Fill(reconele->e5x5());
//                                			SigmaIetaIeta->Fill(reconele->sigmaIetaIeta());
//                                			SigmaEtaEta->Fill(reconele->sigmaEtaEta());
//                                                      OneByFive->Fill(reconele2->e1x5());
//                                                      TwoByFive->Fill(reconele2->e2x5Max());
//                                                      FiveByFive->Fill(reconele2->e5x5());
//                                                      SigmaIetaIeta->Fill(reconele2->sigmaIetaIeta());
//                                                      SigmaEtaEta->Fill(reconele2->sigmaEtaEta());
//
//*/
//
////                                        if (reconele->eta()<1.5 && reconele->eta()>-1.5){
//////                                                SigmaIetaIetaB->Fill(reconele->sigmaIetaIeta());
////  //                                              SigmaEtaEtaB->Fill(reconele->sigmaEtaEta());
////						fout << "0\t" << reconele->e1x5() << "\t"
////                                                                << reconele->e2x5Max() << "\t"
////                                                                << reconele->e5x5() << "\t"
////                                                                << reconele->sigmaIetaIeta() << "\t"
////                                                                << reconele->sigmaEtaEta() << "\t"
////								<< reconele->p4().Pt()
////								<< std::endl;
////                                        }
//////                                        else{
//////                                                SigmaIetaIetaE->Fill(reconele->sigmaIetaIeta());
//////                                                SigmaEtaEtaE->Fill(reconele->sigmaEtaEta());
//////                                        }
////					if (reconele2->eta()<1.5 && reconele2->eta()>-1.5){
//////                                                SigmaIetaIetaB->Fill(reconele2->sigmaIetaIeta());
//////                                                SigmaEtaEtaB->Fill(reconele2->sigmaEtaEta());
////						fout << "0\t" << reconele2->e1x5() << "\t"
////                                                                << reconele2->e2x5Max() << "\t"
////                                                                << reconele2->e5x5() << "\t"
////                                                                << reconele2->sigmaIetaIeta() << "\t"
////                                                                << reconele2->sigmaEtaEta() << "\t"
////								<< reconele2->p4().Pt()
////								<< std::endl;
////                                        }
////                                        else{
////                                                SigmaIetaIetaE->Fill(reconele2->sigmaIetaIeta());
////                                                SigmaEtaEtaE->Fill(reconele2->sigmaEtaEta());
////                                        }
//
///*
//                                        if (reconele->eta()<1.5 && reconele->eta()>-1.5){
//                                                for( unsigned  k=0; k<pfJetCollection->size(); k++ ) { // Find corresponding jet
//                                                        const reco::PFJet* pfjet = &(pfJetCollection->at(k));
//                                                        if (deltaR(*reconele,*pfjet)<0.3){
//                                                                double pt = pfjet->pt();
//                                                                if (reconele->sigmaIetaIeta()<0.01){
//                                                                        SCpass->Fill(pt);
//                                                                }
//                                                                SCall->Fill(pt);
//                                                        }
//                                                }
//                                        }
//                                        if (reconele2->eta()<1.5 && reconele2->eta()>-1.5){
//                                                for( unsigned  k=0; k<pfJetCollection->size(); k++ ) { // Find corresponding jet
//                                                        const reco::PFJet* pfjet = &(pfJetCollection->at(k));
//                                                        if (deltaR(*reconele2,*pfjet)<0.3){
//                                                                double pt = pfjet->pt();
//                                                                if (reconele2->sigmaIetaIeta()<0.01){
//                                                                        SCpass->Fill(pt);
//                                                                }
//                                                                SCall->Fill(pt);
//                                                        }
//                                                }
//                                        }
//*/
///*	fout << "0\t" << reconele->e1x5() << "\t"
//                                                                << reconele->e2x5Max() << "\t"
//                                                                << reconele->e5x5() << "\t"
//                                                                << reconele->sigmaIetaIeta() << "\t"
//                                                                << reconele->sigmaEtaEta() << std::endl;
//                                                        fout << "0\t" << reconele2->e1x5() << "\t"
//                                                                << reconele2->e2x5Max() << "\t"
//                                                                << reconele2->e5x5() << "\t"
//                                                                << reconele2->sigmaIetaIeta() << "\t"
//                                                                << reconele2->sigmaEtaEta() << std::endl;
//*/
//						}
//					}
//				}
//			}
//		}
//	}
//                        for (unsigned int re(0); re < 2; ++re){
//                                const reco::GsfElectron* reconele = &(electronCollection->at(re));
//                              std::cout << "Iecal o ePt: " << reconele->Iecal/reconele->p4.Vect().Pt() << std::endl;
//                              IECAL->Fill(reconele->dr03EcalRecHitSumEt());
//				IECAL->Fill(reconele->dr04EcalRecHitSumEt()/reconele->p4().Pt());
//				std::cout << "Iecal: " << reconele->dr03EcalRecHitSumEt() << std::endl;
//                                IECALOeT->Fill(reconele->Iecal/reconele->p4.Vect().Et());
//                              fout << "0\t" << (reconele->Iecal)/(reconele->p4.Vect().Pt()) << "\t"
//                                              << (reconele->Iecal)<< "\t" << (reconele->p4.Vect().Pt()) << std::endl;
//				for( unsigned  k=0; k<pfJetCollection->size(); k++ ) { // Find corresponding jet
//					const reco::PFJet* pfjet = &(pfJetCollection->at(k));
//					if (deltaR(*reconele,*pfjet)<0.3){
//						double pt = pfjet->pt();
//						PtDist->Fill(pt);
//						SCall->Fill(pt);
///*						double ECALenergy(0), HCALenergy(0);
//						for (unsigned int cons(0); cons < pfjet->getPFConstituents().size(); ++cons){
//							reco::PFCandidatePtr pfcand = pfjet->getPFConstituent(cons);
//							ECALenergy += pfcand->ecalEnergy();
//							HCALenergy += pfcand->hcalEnergy();
//						}
//						double EMF=(ECALenergy/(ECALenergy+HCALenergy));
//						double ntrks = (pfjet->getTrackRefs().size());
////			 fout << "0\t" << pt << "\t" << EMF << "\t" << ntrks  << "\t" << ntrks/pt  << "\t" << reconele->p4().Pt() << "\t" << reconele->dr03EcalRecHitSumEt()/reconele->p4().Pt() << std::endl;
//
//						if(EMF>0.90 && ntrks>2){
////							if (ntrks/pt<0.45 && ((pt<70 && pfjet->phiphiMoment()>0.005) || pfjet->phiphiMoment()>0.005-(0.0001*(pt-70))) ) SCpass->Fill(pt);  //SC
//							if (ntrks/pt<0.45 && reconele->dr04EcalRecHitSumEt()/reconele->p4().Pt()>0.04 ) SCpass->Fill(pt); //SD
//						}
//		//                                                IECALOpT->Fill(reconele->Iecal/pt);
//
//						fout << "0\t" << reconele->e1x5() << "\t"
//                                                                << reconele->e2x5Max() << "\t"
//                                                                << reconele->e5x5() << "\t"
//                                                                << reconele->sigmaIetaIeta() << "\t"
//                                                                << reconele->sigmaEtaEta() << std::endl;
////					}
////				}
//
//
//        		}
//		}
//	}
//*/
// ############################################## Zee BACKGROUND END

/*
// ############################################## Jet BACKGROUND START
	for( unsigned  p=0; p<pfJetCollection->size(); p++ ) {
		const reco::PFJet* pfjet = &(pfJetCollection->at(p));
		double pt = pfjet->pt();
		SCall->Fill(pt);
		double ECALenergy(0), HCALenergy(0);
		for (unsigned int cons(0); cons < pfjet->getPFConstituents().size(); ++cons){
			reco::PFCandidatePtr pfcand = pfjet->getPFConstituent(cons);
			ECALenergy += pfcand->ecalEnergy();
			HCALenergy += pfcand->hcalEnergy();
		}
		double EMF = (ECALenergy/(ECALenergy+HCALenergy));
		int numEle(0);
		const reco::GsfElectron* reconele = NULL;
		for ( unsigned int e(0);e<electronCollection->size();++e){
			const reco::GsfElectron* electron = &(electronCollection->at(e));
			if (deltaR(*pfjet,*electron) < 0.7){
				++numEle;   // electron in Jet
				reconele=electron;
			}
		}
		double ntrks = (pfjet->getTrackRefs().size());
		if (EMF>0.9 && ntrks>2 && numEle>0){
			if (1==numEle && ntrks/pt<0.45){
//				if (reconele->dr03EcalRecHitSumEt()/reconele->p4().Pt()>0.04) SCpass->Fill(pt);  //SD
				if ( (pt<70 && pfjet->phiphiMoment()>0.005) || pfjet->phiphiMoment()>0.005-(0.0001*(pt-70)) ) SCpass->Fill(pt);
			}
			else if(numEle>1) SCpass->Fill(pt);
		}
	}
// ############################################## Jet BACKGROUND END
*/

// ############################################## SIGNAL START
//	for(std::vector<reco::GenParticle>::const_iterator MCparticle=genParticleCollection->begin(); MCparticle!=genParticleCollection->end(); ++MCparticle){
//		if ( MCparticle->pdgId() == 3000006){	//hd2
//			if (MCparticle->mother()->pdgId() == 25){  // Protect from loopback
//				for(unsigned int j(0); j<MCparticle->numberOfDaughters(); ++j){
//					for(unsigned int k(0); k<MCparticle->daughterRef(j)->numberOfDaughters(); ++k){
//						if (MCparticle->daughterRef(j)->daughterRef(k)->pdgId() == 3000001){
//							DeltaR->Fill(deltaR(*MCparticle,*MCparticle->daughterRef(j)->daughterRef(k)));
//						}
//					}
//				}
//			}
//		}
//	}

////	##### SELECTION #####
//	std::cerr << "Beginning Selection\n";
//	unsigned int numOfGoodJets(0);
//	for(std::vector<reco::PFJet>::const_iterator jet=pfJetCollection->begin(); jet!=pfJetCollection->end(); ++jet){
//		if (jet->pt() > 20){
//			std::cerr << "\tjet (eta,phi,pt): " << jet->eta() << "," << jet->phi() <<","<<jet->pt() << std::endl;
//			++numOfGoodJets;
//		}
//	}
//	std::cerr << "\tnumOfGoodJets: " << numOfGoodJets << std::endl;

//	##### FIND EJ's in MC #####
	PFJet1 = NULL;
	PFJet2 = NULL;
	GenJet1 = NULL;
	GenJet2 = NULL;
	CaloJet1 = NULL;
	CaloJet2 = NULL;

	PFJetPtCompStore.clear();
	GenJetPtCompStore.clear();
	CaloJetPtCompStore.clear();
	this->FindEJs();

	std::vector<const reco::GenJet*> GenJets;
	GenJets.push_back(GenJet1);
	GenJets.push_back(GenJet2);

	std::vector<const reco::PFJet*> PFJets;
	PFJets.push_back(PFJet1);
	PFJets.push_back(PFJet2);

	std::vector<const reco::CaloJet*> CaloJets;
	CaloJets.push_back(CaloJet1);
	CaloJets.push_back(CaloJet2);

//	##### RUN SELECTION ON EJ's #####
	numOfGoodElectronJets=0;
	this->RunEJSelection(PFJets);
	std::cerr << "numOfGoodElectronJets: " << numOfGoodElectronJets << std::endl;
	NumGoodEJ->Fill(numOfGoodElectronJets);
	if (numOfGoodElectronJets == 2){
		std::cerr << "\t2 EJ found\n";
		if (GenJet1 && GenJet2){	// Does it make sense to have this here as this if depends on PF not Gen
			TLorentzVector gjet1(GenJet1->px(), GenJet1->py(), GenJet1->pz(), GenJet1->energy());
			TLorentzVector gjet2(GenJet2->px(), GenJet2->py(), GenJet2->pz(), GenJet2->energy());
			TLorentzVector Higgs = gjet1 + gjet2;
			HiggsMass->Fill(Higgs.M());
		}
		if (PFJet1 && PFJet2){
			TLorentzVector pfjet1(PFJet1->px(), PFJet1->py(), PFJet1->pz(), PFJet1->energy());
			TLorentzVector pfjet2(PFJet2->px(), PFJet2->py(), PFJet2->pz(), PFJet2->energy());
			TLorentzVector Higgs2 = pfjet1 + pfjet2;
			HiggsMass2->Fill(Higgs2.M());
		}
		for (size_t i(0); i<PFJetPtCompStore.size();++i){	// ie print only the jets which passed the event selection
			PFJetPtComp->Fill(PFJetPtCompStore.at(i));
			std::cerr << "PFJetPtComp: " << PFJetPtCompStore.at(i) << std::endl;
		}
		for (size_t j(0); j<GenJetPtCompStore.size();++j){
			GenJetPtComp->Fill(GenJetPtCompStore.at(j));
			std::cerr << "GenJetPtComp: " << GenJetPtCompStore.at(j) << std::endl;
		}
	}

// ############################################## SIGNAL END
}

void LJ::FindEJs(){
	for (unsigned int MCpos(0); MCpos < genParticleCollection->size(); ++MCpos){	// Loop over genParticles to find hd2
		const reco::GenParticle* MCparticle = &(genParticleCollection->at(MCpos));
		if ( MCparticle->pdgId() == 3000006){	// Hd2
			if (MCparticle->mother()->pdgId() == 25){  // Protect from loopback
//			std::cout << "MCparticle->pdgId(): " << MCparticle->pdgId();
//			if (MCparticle->mother()) std::cout << " has mother: " << MCparticle->mother()->pdgId() << std::endl;
//			else std::cout << std::endl;
				float foundHd0=false;
				for(unsigned int j(0); j<MCparticle->numberOfDaughters(); ++j){
					for(unsigned int k(0); k<MCparticle->daughterRef(j)->numberOfDaughters(); ++k){
						if (MCparticle->daughterRef(j)->daughterRef(k)->pdgId() == 3000004) foundHd0=true;
					}
				}
				if (!foundHd0){	// = No MET in EJ
					float PFPt(0), GenJetPt(0), CaloPt(0);
					bool PFmatched(false), GenJetMatched(false), Calomatched(false);
					std::cerr << "MCparticle->pdgId(): " << MCparticle->pdgId() << " ("<< MCparticle->eta() << "," << MCparticle->phi() << ") with pt: " << MCparticle->pt() << std::endl;
					++numMCparts;
	//				for( unsigned  k=0; k<genjetCollection->size(); k++ ) {
	//					const reco::PFJet* jet = &(pfJetCollection->at(k));
	//				for(std::vector<reco::PFJet>::const_iterator jet=pfJetCollection->begin(); jet!=pfJetCollection->end(); ++jet){
					for( unsigned int p=0; p<pfJetCollection->size(); p++ ) {
						const reco::PFJet* pfjet = &(pfJetCollection->at(p));
	//				std::cout << " genjet pt: " << jet->pt() << " at eta: " << jet->eta() << " phi: " << jet->phi() << std::endl;
						if (deltaR(*pfjet,*MCparticle)<0.3){	// ## Jet is 'matched' if within dR<0.3 of MCParticle->EJ
	//						&& fabs(jet->pt()-MCparticle->pt())/MCparticle->pt() <0.25){	//Match GenJet to the GenParticle = hd2
	//						std::cout << "genjet matched to hd3; jet1= " << jet1 << " jet2= " << jet2 << std::endl;
							++numRecoJetMatched;
							PFmatched=true;
							std::cerr << "\tMatched RecoJet: (" << pfjet->eta() << "," << pfjet->phi() << ") with pt: " << pfjet->pt()
											<< "\trel diff: " << fabs(pfjet->pt()-MCparticle->pt())/MCparticle->pt()
											<< std::endl;
							PFPt=pfjet->pt();
							PFJetPtCompStore.push_back((pfjet->pt()-MCparticle->pt())/MCparticle->pt());
							if (!PFJet1) PFJet1=pfjet;
							else PFJet2=pfjet;
							std::cerr << "\tafter: PFJet1= " << PFJet1 << " PFJet2= " << PFJet2 << std::endl;
	//std::cout << " genjet matched to hd3 eta: "
	//					<< MCparticle->eta() << " phi: " << MCparticle->phi()  << " with dR: " << deltaR(*jet,*MCparticle) << "\n";
						}
					}
					for( unsigned int c=0; c<caloJetCollection->size(); c++ ) {
						const reco::CaloJet* calojet = &(caloJetCollection->at(c));
						if (deltaR(*calojet,*MCparticle)<0.3){	// ## Jet is 'matched' if within dR<0.3 of MCParticle->
							++numCaloRecoJetMatched;
							Calomatched=true;
							std::cerr << "\tMatched calojet: (" << calojet->eta() << "," << calojet->phi() << ") with pt: " << calojet->pt()
											<< "\trel diff: " << fabs(calojet->pt()-MCparticle->pt())/MCparticle->pt()
											<< std::endl;
							CaloPt=calojet->pt();
							CaloJetPtCompStore.push_back((calojet->pt()-MCparticle->pt())/MCparticle->pt());
							if (!CaloJet1) CaloJet1=calojet;
							else CaloJet2=calojet;
							std::cerr << "\tafter: CaloJet1= " << CaloJet1 << " CaloJet2= " << CaloJet2 << std::endl;
						}
					}
	//				for(std::vector<reco::GenJet>::const_iterator genjet=genjetCollection->begin(); genjet!=genjetCollection->end(); ++genjet){
					for( unsigned int g=0; g<genjetCollection->size(); g++ ) {
						const reco::GenJet* genjet = &(genjetCollection->at(g));
						if (deltaR(*genjet,*MCparticle)<0.3){	// ## Jet is 'matched' if within dR<0.3 of MCParticle->EJ
	//						&& fabs(genjet->pt()-MCparticle->pt())/MCparticle->pt() <0.25 ){
							std::cerr << "\tMatched GenJet: (" << genjet->eta() << "," << genjet->phi() << ") with pt: " << genjet->pt()
									<< "\trel diff: " << fabs(genjet->pt()-MCparticle->pt())/MCparticle->pt()
											<< std::endl;
							++numGenJetMatched;
							GenJetPtCompStore.push_back((genjet->pt()-MCparticle->pt())/MCparticle->pt());
							GenJetMatched=true;
							GenJetPt=genjet->pt();
							if (!GenJet1) GenJet1=genjet;
							else GenJet2=genjet;
						}
					}
					if (PFmatched && GenJetMatched) JetPts->Fill(GenJetPt,PFPt);
				} // no MET in EJ
			}	// if EJ from H
		}	// if EJ
	}	// loop genParticles
}

void LJ::RunEJSelection(std::vector<const reco::PFJet*>& PFJets){
	for( unsigned  p=0; p<PFJets.size(); p++ ) {
		const reco::PFJet* pfjet;	// Can put these outside loop to stop mutliple constructions
		pfjet = PFJets.at(p);
		if (pfjet){
			double pt = pfjet->pt();
			SCall->Fill(pt);
			double ECALenergy(0), HCALenergy(0);
			for (unsigned int cons(0); cons < pfjet->getPFConstituents().size(); ++cons){
				reco::PFCandidatePtr pfcand = pfjet->getPFConstituent(cons);
				ECALenergy += pfcand->ecalEnergy();
				HCALenergy += pfcand->hcalEnergy();
			}
			double EMF = (ECALenergy/(ECALenergy+HCALenergy));
			int numEle(0);
			for ( unsigned int e(0);e<electronCollection->size();++e){
				const reco::GsfElectron* electron = &(electronCollection->at(e));
				if (deltaR(*pfjet,*electron) < 0.7) ++numEle;   // electron in Jet ==> ak7
			}
			double ntrks = (pfjet->getTrackRefs().size());
			if (EMF>0.95){
				PassEMF->Fill(pt);
				if (ntrks>2){
					PassNtrks->Fill(pt);
					if(numEle>0){
						PassNEle->Fill(pt);
						if (1==numEle && ntrks/pt<0.45){
							if ( (pt<70 && pfjet->phiphiMoment()>0.005) || pfjet->phiphiMoment()>0.005-(0.0001*(pt-70)) ){
								SCpass->Fill(pt);
								++numOfGoodElectronJets;
							}
						}
						else if(numEle>1){
							SCpass->Fill(pt);
							++numOfGoodElectronJets;
						}
					}
				}
			}
		} //pfjet exists
	} //pfJets
}

void LJ::RunEJSelection(std::vector<const reco::CaloJet*>& CaloJets){
	for( unsigned  p=0; p<CaloJets.size(); p++ ) {
		const reco::CaloJet* calojet;	// Can put these outside loop to stop mutliple constructions
		calojet = CaloJets.at(p);
		if (calojet){
			double pt = calojet->pt();
			SCall->Fill(pt);
			double EMF = calojet->emEnergyFraction();
			int numEle(0);
			for ( unsigned int e(0);e<electronCollection->size();++e){
				const reco::GsfElectron* electron = &(electronCollection->at(e));
				if (deltaR(*calojet,*electron) < 0.7) ++numEle;   // electron in Jet ==> ak7
			}
			double ntrks = (reco::JetTracksAssociation::tracksNumber(*assocCaloJetTracks, calojet));
			if (EMF>0.95){
				PassEMF->Fill(pt);
				if (ntrks>2){
					PassNtrks->Fill(pt);
					if(numEle>0){
						PassNEle->Fill(pt);
						if (1==numEle && ntrks/pt<0.45){
							if ( (pt<70 && calojet->phiphiMoment()>0.005) || calojet->phiphiMoment()>0.005-(0.0001*(pt-70)) ){
								SCpass->Fill(pt);
								++numOfGoodElectronJets;
							}
						}
						else if(numEle>1){
							SCpass->Fill(pt);
							++numOfGoodElectronJets;
						}
					}
				}
			}
		} //calojet exists
	} //caloJets
}


void LJ::beginJob(){}
void LJ::endJob() {
	std::cout << "numMCparts: " << numMCparts << std::endl;
	std::cout << "numRecoJetMatched: " << numRecoJetMatched << std::endl;
	std::cout << "numGenJetMatched: " << numGenJetMatched << std::endl;
}

DEFINE_FWK_MODULE(LJ);
