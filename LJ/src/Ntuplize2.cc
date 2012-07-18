// Package:    Ntuplize2
// Original Author:  Kristian Hahn	Created:  Tue Aug  4 08:56:23 CEST 2009
// Modified: Benjamin R-S, Chris M

#include "../interface/Ntuplize2.h"
#include "../interface/MCfind2.h"

Ntuplize2::Ntuplize2(const edm::ParameterSet& iConfig){
	verbose= iConfig.getUntrackedParameter<bool>("verbosity",false);
	if (verbose) std::cerr << "Ntuplize2 created ... " << std::endl;
	tupfile = iConfig.getUntrackedParameter <std::string>("TupFileName", "tup.root");
	do_data = do_signalMC = do_zeeMC = do_zmumuMC = do_wenuMC = do_wmunuMC = do_wwMC = do_zzMC = do_wjetMC = do_zjetMC = false;
	run_scheme = iConfig.getUntrackedParameter <std::string>("RunningScheme", "DATA");
	if (run_scheme == std::string("DATA")) do_data = true;
	else if (run_scheme == std::string("WENU_MC")) do_wenuMC = true;
	else if (run_scheme == std::string("WMUNU_MC")) do_wmunuMC = true;
	else if (run_scheme == std::string("WW_MC")) do_wwMC = true;
	else if (run_scheme == std::string("ZZ_MC")) do_zzMC = true;
	else if (run_scheme == std::string("WJET_MC")) do_wjetMC = true;
	else if (run_scheme == std::string("ZJET_MC")) do_zjetMC = true;
	else if (run_scheme == std::string("ZEE_MC")) do_zeeMC = true;
	else if (run_scheme == std::string("ZMUMU_MC")) do_zmumuMC = true;
	else if (run_scheme == std::string("SIGNAL_MC")) do_signalMC = true;
	do_electron_trigger = do_muon_trigger = do_jet_trigger = do_photon_trigger = false;
	trigger_scheme = iConfig.getUntrackedParameter <std::string>("TriggerScheme", "");
	if (trigger_scheme == std::string("jet")) do_jet_trigger = true;
	else if (trigger_scheme == std::string("electron")) do_electron_trigger = true;
	else if (trigger_scheme == std::string("muon")) do_muon_trigger = true;
	else if (trigger_scheme == std::string("photon")) do_photon_trigger = true;
	eventNr=0;
}

Ntuplize2::~Ntuplize2() {}

unsigned long long Ntuplize2::checkTrigger(const edm::TriggerResults& HLTR,std::vector<unsigned>& triggerIndices) {
	unsigned long bits = 0x0;
	if (verbose) std::cout << "sizeof(bits) : " << sizeof(bits) << std::endl;
	if (verbose) std::cout << "looping over trigger indices ... HLTR.size() " << HLTR.size() << std::endl;
	for (unsigned i = 0; i < triggerIndices.size(); i++) {
		unsigned triggerIndex = triggerIndices[i];
		if (verbose) std::cout << "triggerIndex=" << triggerIndices[i] << std::endl;
		// triggerIndex must be less than the size of HLTR or you get a CMSException: _M_range_check
		if (triggerIndex < HLTR.size()) {
			if (HLTR.accept(triggerIndex)) bits |= (1ULL << i);
		}
	}
	return bits;
}

void Ntuplize2::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
//	if (!do_data) return;
	triggerNames_.clear();
	triggerIndices_.clear();
	lastFilterLabels_.clear();
	//--- htlConfig_
	processname_ = "HLT";
	bool changed(true);
	_hlt_initialized = hltConfig_.init(iRun, iSetup, processname_, changed);
	if (!_hlt_initialized) {
		//if (!hltConfig_.init(iRun,iSetup,processname_,changed)) {
		processname_ = "FU";
		_hlt_initialized = hltConfig_.init(iRun, iSetup, processname_, changed);
		if (!_hlt_initialized) {
			//if (!hltConfig_.init(iRun,iSetup,processname_,changed)){
			LogDebug("FakeRate") << "HLTConfigProvider failed to initialize.";
		}
	}
	std::vector < std::string > tmptrignames;
	for (unsigned i = 0; i < hltConfig_.triggerNames().size(); ++i) {
		if (verbose) std::cout << "trigger: " << i << "\tname: "	<< hltConfig_.triggerNames()[i] << std::endl;
		std::string str(hltConfig_.triggerNames()[i]);
		if (do_data){
			if (do_electron_trigger) {
				if (strstr(str.c_str(), "HLT_Ele") != NULL && strstr(str.c_str(), "Jet") == NULL && strstr(str.c_str(), "Tau") == NULL) {
					if (verbose) std::cout << "\t\t found an ele trigger!!! " << std::endl;
					tmptrignames.push_back(str);
				}
				else if (strstr(str.c_str(), "HLT_DoubleEle") != NULL) {
					if (verbose) std::cout << "\t\t found an double ele trigger!!! "	<< std::endl;
					tmptrignames.push_back(str);
				}
			}
			if (do_muon_trigger) {
				if (strstr(str.c_str(), "HLT_Mu") != NULL && strstr(str.c_str(), "Jet") == NULL && strstr(str.c_str(), "Jpsi") == NULL && strstr(str.c_str(), "CaloId") == NULL) {
					if (verbose) std::cout << "\t\t found an mu trigger!!! " << std::endl;
					tmptrignames.push_back(str);
				}
				else if (strstr(str.c_str(), "HLT_DoubleMu") != NULL) {
					if (verbose) std::cout << "\t\t found an double mu trigger!!! " << std::endl;
					tmptrignames.push_back(str);
				}
			}
			if (do_photon_trigger) {
				if (strstr(str.c_str(), "HLT_Photon") != NULL && strstr(str.c_str(), "HT") == NULL && strstr(str.c_str(), "HLT_Photon2") == NULL) {
					if (verbose) std::cout << "\t\t found an photon trigger!!! " << std::endl;
					tmptrignames.push_back(str);
				}
			}
			if (do_jet_trigger) {
				if (strstr(str.c_str(), "HLT_Jet") != NULL && strstr(str.c_str(), "NoiseFiltered") == NULL) {
					if (verbose) std::cout << "\t\t found an jet trigger!!! " << std::endl;
					tmptrignames.push_back(str);
				}
			}
		}
		if (!do_data){
			if (strstr(str.c_str(), "HLT_Ele") != NULL && strstr(str.c_str(), "Jet") == NULL && strstr(str.c_str(), "Tau") == NULL) {
				if (verbose) std::cout << "\t\t found an ele trigger!!! " << std::endl;
				tmptrignames.push_back(str);
			}
			if (strstr(str.c_str(), "HLT_DoubleEle") != NULL) {
				if (verbose) std::cout << "\t\t found an double ele trigger!!! "	<< std::endl;
				tmptrignames.push_back(str);
			}
			if (strstr(str.c_str(), "HLT_Mu") != NULL && strstr(str.c_str(), "Jet") == NULL && strstr(str.c_str(), "Jpsi") == NULL && strstr(str.c_str(), "CaloId") == NULL) {
				if (verbose) std::cout << "\t\t found an mu trigger!!! " << std::endl;
				tmptrignames.push_back(str);
			}
			if (strstr(str.c_str(), "HLT_DoubleMu") != NULL) {
				if (verbose) std::cout << "\t\t found an double mu trigger!!! " << std::endl;
				tmptrignames.push_back(str);
			}
			if (strstr(str.c_str(), "HLT_Photon") != NULL && strstr(str.c_str(), "HT") == NULL && strstr(str.c_str(), "HLT_Photon2") == NULL) {
				if (verbose) std::cout << "\t\t found an photon trigger!!! " << std::endl;
				tmptrignames.push_back(str);
			}
			if (strstr(str.c_str(), "HLT_Jet") != NULL && strstr(str.c_str(), "NoiseFiltered") == NULL) {
				if (verbose) std::cout << "\t\t found an jet trigger!!! " << std::endl;
				tmptrignames.push_back(str);
			}
		}
	}

	int index;
	int limit = hltConfig_.triggerNames().size() - 1;
	for (unsigned i = 0; i < tmptrignames.size(); i++) {
		index = hltConfig_.triggerIndex(tmptrignames[i].c_str());
		if (index < limit) {
			triggerNames_.push_back(tmptrignames[i]);
			triggerIndices_.push_back((unsigned) index);
		}
	}
	for (unsigned i = 0; i < triggerIndices_.size(); i++) {
		if (verbose) std::cerr << "trigger index: " << triggerIndices_[i] << " ";
		unsigned modindex = hltConfig_.moduleLabels(triggerIndices_[i]).size() - 2;
		lastFilterLabels_.push_back(hltConfig_.moduleLabels(triggerIndices_.at(i)).at(modindex));
		if (verbose){
			std::cerr << "last filter : " << hltConfig_.moduleLabels(triggerIndices_[i])[modindex]	<< std::endl;
			std::cerr << "hltConfig_.moduleLabels(triggerIndices_.at(i)) of size "<< hltConfig_.moduleLabels(triggerIndices_.at(i)).size() << " ="  << std::endl;
			for (unsigned int modLabs=0; modLabs < hltConfig_.moduleLabels(triggerIndices_.at(i)).size(); ++modLabs ){
				std::cerr << "\t" << hltConfig_.moduleLabels(triggerIndices_.at(i)).at(modLabs) << std::endl;
			}
		}
		allFilterLabels_.push_back(hltConfig_.moduleLabels(triggerIndices_.at(i)));
	}
}

// ------------ method called once each job just before starting event loop  ------------
void Ntuplize2::beginJob() {
	if (verbose) std::cerr << "output file is " << tupfile << std::endl;
	gSystem->Load("./MyObjs_C.so");

	genpVector = new std::vector<MyGenParticle*>();
	genmuonVector = new std::vector<MyGenParticle*>();
	genelectronVector = new std::vector<MyGenParticle*>();
	recoelectronVector = new std::vector<MyRecoElectron*>();
	pfjetVector = new std::vector<MyPFJet*>();
	pfeleVector = new std::vector<MyPFElectron*>();
	ak5calojetVector = new std::vector<MyCaloJet*>();
	ak7calojetVector = new std::vector<MyCaloJet*>();
	genjetVector = new std::vector<MyPFJet*>();
	muonVector = new std::vector<MyMuon*>();
	SCVector = new std::vector<MySC*>();
	photonVector = new std::vector<MyPhoton*>();
	hltobjVector = new std::vector<MyHLTObj*>();
	puVector = new std::vector<MyPU*>();
	genSigMetVector = new std::vector<MyGenParticle*>();

	outf = new TFile(tupfile.c_str(), "RECREATE");
	hevt = new TH1F("hevt", "hev", 3, 0, 2);
	tree = new TTree("tree", "thetree");
	nt = new TNtuple("nt", "nt", "mva_e_pi:mva_gamma:pix");

	tree->Branch("gen_particles", "vector<MyGenParticle*>", &genpVector);
	tree->Branch("reco_electrons", "vector<MyRecoElectron*>", &recoelectronVector);
	tree->Branch("pf_jets", "vector<MyPFJet*>", &pfjetVector);
	tree->Branch("pf_electrons", "vector<MyPFElectron*>", &pfeleVector);
	tree->Branch("ak5calo_jets", "vector<MyCaloJet*>", &ak5calojetVector);
	tree->Branch("ak7calo_jets", "vector<MyCaloJet*>", &ak7calojetVector);
	tree->Branch("mygen_jets", "vector<MyPFJet*>", &genjetVector);
	tree->Branch("reco_muons", "vector<MyMuon*>", &muonVector);
	tree->Branch("hlt_objs", "vector<MyHLTObj*>", &hltobjVector);
	tree->Branch("gen_electrons", "vector<MyGenParticle*>", &genelectronVector);
	tree->Branch("gen_muons", "vector<MyGenParticle*>", &genmuonVector);
	tree->Branch("superclusters", "vector<MySC*>", &SCVector);
	tree->Branch("hlt", &hlt);
	tree->Branch("met", &met);
	tree->Branch("event_data", &evtdata);
	tree->Branch("photons", &photonVector);
	if (!do_data) tree->Branch("pileup", &puVector);
	if (verbose) tree->Print();
}


void Ntuplize2::analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup) {
	++eventNr;
	if (verbose) std::cout << "EVENT count = " << eventNr << std::endl;
	if (verbose) std::cerr << "Inside analyze ..." << std::endl;
	using namespace edm;
	edm::Handle <reco::TrackCollection> hGeneralTracks;				iEvent.getByLabel("generalTracks", hGeneralTracks);
	edm::Handle <reco::GenParticleCollection> genp;					iEvent.getByLabel("genParticles", genp);
	edm::Handle <reco::GsfElectronCollection> gsfElectrons;				iEvent.getByLabel("gsfElectrons", gsfElectrons);
	edm::Handle <reco::PFJetCollection> pfJets;					iEvent.getByLabel("ak5PFJets", pfJets);
	edm::Handle <reco::PFCandidateCollection> pfCandidates;				iEvent.getByLabel("particleFlow", pfCandidates);
	edm::Handle <reco::CaloJetCollection> ak5caloJets;				iEvent.getByLabel("ak5CaloJets", ak5caloJets);
	edm::Handle <reco::CaloJetCollection> ak7caloJets;				iEvent.getByLabel("ak7CaloJets", ak7caloJets);
	edm::Handle <reco::JetTracksAssociationCollection> ak5assocCaloJetTracks;	iEvent.getByLabel("ak5JetTracksAssociatorAtVertex", ak5assocCaloJetTracks);
	edm::Handle <reco::JetTracksAssociationCollection> ak7assocCaloJetTracks;	iEvent.getByLabel("ak7JetTracksAssociatorAtVertex", ak7assocCaloJetTracks);
	edm::Handle <reco::JetTagCollection> bTag_secvtx_Handle;			iEvent.getByLabel("combinedSecondaryVertexBJetTags", bTag_secvtx_Handle);
	edm::Handle <reco::JetTagCollection> bTag_softele_Handle;			iEvent.getByLabel("softElectronByIP3dBJetTags", bTag_softele_Handle);
	edm::Handle <reco::JetTagCollection> bTag_jetprob_Handle;			iEvent.getByLabel("jetBProbabilityBJetTags", bTag_jetprob_Handle);
	edm::Handle <reco::PFMETCollection> pfMetCollection;				iEvent.getByLabel("pfMet", pfMetCollection);
	edm::Handle <reco::CaloMETCollection> caloMetCollection;			iEvent.getByLabel("corMetGlobalMuons", caloMetCollection);
	edm::Handle <reco::SuperClusterCollection> SCCollectionEE;			iEvent.getByLabel("correctedMulti5x5SuperClustersWithPreshower", SCCollectionEE);
	edm::Handle <reco::SuperClusterCollection> SCCollectionEB;			iEvent.getByLabel("correctedHybridSuperClusters", SCCollectionEB);
	edm::Handle <reco::MuonCollection> recoMuons;					iEvent.getByLabel("muons", recoMuons);
	edm::Handle <reco::PhotonCollection> Photons;					iEvent.getByLabel("photons", Photons);
	edm::Handle <edm::TriggerResults> HLTR;						iEvent.getByLabel(edm::InputTag("TriggerResults", "", "HLT"), HLTR);
	edm::Handle <trigger::TriggerEvent> trgEvent;					iEvent.getByLabel(edm::InputTag("hltTriggerSummaryAOD", "", "HLT"), trgEvent);
//	edm::Handle <reco::BeamSpot> beamSpot;						iEvent.getByLabel ("offlineBeamSpot",beamSpot);
	edm::Handle <reco::VertexCollection> primaryVertex;					iEvent.getByLabel("offlinePrimaryVertices",primaryVertex);
	std::vector <trigger::TriggerObject> trigObjs;		// is this used?

	if (!do_data) {
		// PU
		if (verbose) std::cerr << "gonna do PU .. " << std::endl;
		puVector->clear();
		edm::Handle < std::vector<PileupSummaryInfo> > PupInfo;			iEvent.getByLabel("addPileupInfo", PupInfo);
		std::vector<PileupSummaryInfo>::const_iterator PVI;
		for (PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
			MyPU* pu = new MyPU();
			pu->bunchCrossing = PVI->getBunchCrossing();
			pu->nInteractions = PVI->getPU_NumInteractions();
			puVector->push_back(pu);
		}
		if (verbose) std::cerr << "done w/ PU .. " << std::endl;
	}
	if (verbose) std::cerr << "================================================" << std::endl;
	hevt->Fill(1);

//########################### HLT INFO ###########################
// BRS: store all HLT info even if MC
	hltobjVector->clear();
	hlt.bits = 0;
	if (verbose) std::cerr << "getting TriggerResults ... " << std::endl;
	if (verbose) std::cout << "checking trigger  ... " << std::endl;
	unsigned long bits = checkTrigger((*HLTR), triggerIndices_);
	// Something of interest fired, store HLT stuff
	if (verbose) std::cout << "storing HLT info ... " << std::endl;
	hlt.bits = bits;
	const trigger::TriggerObjectCollection& TOC(trgEvent->getObjects());
	if (verbose) std::cerr << "looping over " << triggerIndices_.size() << " indices" << std::endl;

	for (unsigned i = 0; i < triggerIndices_.size(); i++) {
		unsigned long long mask = 0x0;
		mask |= (1ULL << i);
		if (!(bits & mask)) continue;
		if (verbose) std::cerr << "fired bit " << i << "\triggername: " << triggerNames_.at(i) << std::endl;

//##### Find the p4's of objects only matching the last filter in a path:
//		std::cerr << "EVENT["<<iEvent.id().event()<<"] HLT MODIFICATIONS:.. For HLTPATH: " << triggerNames_.at(i) << std::endl;
//		edm::InputTag myLastFilter = edm::InputTag(lastFilterLabels_[i].c_str(), "", "HLT");
////		if (verbose)
//			std::cerr << "checking filterindex vs size of filters : " << trgEvent->filterIndex(myLastFilter) << " vs " << trgEvent->sizeFilters() << std::endl;
//		if (trgEvent->filterIndex(myLastFilter) < trgEvent->sizeFilters()){
////			if (verbose)
//				std::cerr << "going to fetch keys for " << myLastFilter << " ... " << std::endl;
//			const trigger::Keys& keys( trgEvent->filterKeys( trgEvent->filterIndex(myLastFilter) ) );
////			if (verbose)
//				std::cerr << "now looping over " << keys.size() << " keys ... " << std::endl;
//			for (unsigned hlto = 0; hlto < keys.size(); hlto++) {
//				unsigned hltf = keys[hlto];
//				const trigger::TriggerObject& L3obj(TOC[hltf]);
//				TVector3 tvec(0, 0, 0);
//				tvec.SetPtEtaPhi(L3obj.pt(), L3obj.eta(), L3obj.phi());
//				MyHLTObj* o = new MyHLTObj();
//				o->p3 = TVector3(tvec);
//				o->label = i;
//				o->triggername = triggerNames_[i];
//				o->filtername = lastFilterLabels_[i];
//				hltobjVector->push_back(o);
////				if (verbose)
//					std::cerr << "obj pt: " << o->p3.Pt() << std::endl;
//			}
//		}

//##### Find the p4's of objects matching all filters in a path --> as some paths will be double legs:
		if (verbose) std::cerr << "EVENT["<<iEvent.id().event()<<"] HLT MODIFICATIONS:.. For HLTPATH: " << triggerNames_.at(i)
				<< "\tallFilterLabels_.at(i).size(): " << allFilterLabels_.at(i).size() << std::endl;
		for (unsigned int filterLab=0; filterLab < allFilterLabels_.at(i).size(); ++filterLab ){
			if (verbose) std::cerr << "\tallFilterLabels_.at(i).at(filterLab).c_str(): " << allFilterLabels_.at(i).at(filterLab).c_str() << std::endl;
			edm::InputTag myFilter = edm::InputTag(allFilterLabels_.at(i).at(filterLab).c_str(), "", "HLT");
			if (verbose) std::cerr << "\ttrgEvent->filterIndex(myFilter): " << trgEvent->filterIndex(myFilter) << "\ttrgEvent->sizeFilters(): " << trgEvent->sizeFilters() << std::endl;
			if (trgEvent->filterIndex(myFilter) < trgEvent->sizeFilters()){	// not sure if this comp will now work
				if (verbose) std::cerr << "going to fetch keys for " << myFilter << " ... " << std::endl;
				const trigger::Keys& keys( trgEvent->filterKeys( trgEvent->filterIndex(myFilter) ) );
				if (verbose) std::cerr << "now looping over " << keys.size() << " keys ... " << std::endl;
				for (unsigned hlto = 0; hlto < keys.size(); hlto++) {
					unsigned hltf = keys[hlto];
					const trigger::TriggerObject& L3obj(TOC[hltf]);
					TVector3 tvec(0, 0, 0);
					tvec.SetPtEtaPhi(L3obj.pt(), L3obj.eta(), L3obj.phi());
					MyHLTObj* o = new MyHLTObj();
					o->p3 = TVector3(tvec);
					o->label = i;
					o->triggername = triggerNames_.at(i);
					o->filtername = allFilterLabels_.at(i).at(filterLab);
					hltobjVector->push_back(o);
					if (verbose) std::cerr << "obj pt: " << o->p3.Pt() << std::endl;
				}
			}
		}
	}
	if (verbose) std::cout << "HLT info stored ... " << std::endl;
	// end HLT


//########################### Primary Vertex ###########################
	bool goodPV = false;
	std::vector<const reco::Vertex*> goodVertices;
	if (verbose) std::cerr << "Number of Primary Vertices = " << primaryVertex->size() << std::endl;
	for (size_t vIndex=0; vIndex<primaryVertex->size(); ++vIndex){
		const reco::Vertex* vtx=&(primaryVertex->at(vIndex));
		if (fabs(vtx->z()<24.0) &&	//z<24cm
			(fabs(vtx->position().rho()) < 2.0) &&	//radially<2cm
			(vtx->ndof()>4.0)){	//ndof>4
			goodPV=true;
			goodVertices.push_back(vtx);
		}
	}
	if (verbose) std::cerr << "Num Good Primary Vertices = " << goodVertices.size() << std::endl;



//#################################################################
//########################## Start of MC ##########################

//########################### WW MC ###########################
	if (do_wwMC) {
		if (verbose) std::cerr << "clearing gen vecs ... " << std::endl;
		genpvec.clear();
		genpVector->clear();
		if (verbose) std::cerr << "gonna loop ... " << std::endl;
		unsigned nlabels = 0, old_label = 999;
		unsigned Wcount = 0;
		for (unsigned i = 0; i < genp->size(); i++) {
			if (verbose) std::cerr << "i: " << i << "\tpdg: " << ((*genp)[i]).pdgId()
					<< "\tstatus: " << ((*genp)[i]).status() << "\tpt: "
					<< ((*genp)[i]).pt() << "\teta: " << ((*genp)[i]).eta()
					<< std::endl;
			if (abs(((*genp)[i]).pdgId()) == 24 && ((*genp)[i]).status() == 3) {
				findWDaughters(((*genp)[i]), genpvec, Wcount);
				Wcount++;
			}
		}
		if (verbose) std::cerr << "how many gen particles : " << genpvec.size() << std::endl;
		for (unsigned j = 0; j < genpvec.size(); j++) {
			if (verbose) std::cerr << "j: " << j
					<< "\tpdg: " << genpvec[j].first.pdgId()
					<< "\tpt: " << elevec[j].first.pt()
					<< "\teta: " << elevec[j].first.eta()
					<< "\tphi: " << elevec[j].first.phi()
					<< "\tlabel: " << elevec[j].second
					<< std::endl;
			MyGenParticle* g = new MyGenParticle();
			TLorentzVector tvec(0, 0, 0, 0);
			tvec.SetPtEtaPhiE(genpvec[j].first.pt(), genpvec[j].first.eta(),genpvec[j].first.phi(), genpvec[j].first.energy());
			g->pdg = genpvec[j].first.pdgId();
			g->p4 = tvec;
			g->label = genpvec[j].second;
			genpVector->push_back(g);
			if (g->label != old_label) {
				old_label = g->label;
				nlabels++;
			}
		}
		if (verbose) std::cerr << "checking again ... " << std::endl;
		for (unsigned y = 0; y < genpVector->size(); y++) {
			if (verbose) std::cerr << "y: " << y << "\tpdg: " << (*genpVector)[y]->pdg
					<< "\tpt: " << (*genpVector)[y]->p4.Pt() << std::endl;
		}
	} // do_WW

//########################### ZZ MC ###########################
	if (do_zzMC) {
		if (verbose) std::cerr << "clearing gen vecs ... " << std::endl;
		genpvec.clear();
		genpVector->clear();
		if (verbose) std::cerr << "gonna loop ... " << std::endl;
		unsigned nlabels = 0, old_label = 999;
		unsigned Zcount = 0;
		for (unsigned i = 0; i < genp->size(); i++) {
			if (verbose) std::cerr << "i: " << i << "\tpdg: " << ((*genp)[i]).pdgId()
					<< "\tstatus: " << ((*genp)[i]).status() << "\tpt: "
					<< ((*genp)[i]).pt() << "\teta: " << ((*genp)[i]).eta()
					<< std::endl;
			if (abs(((*genp)[i]).pdgId()) == 23 && ((*genp)[i]).status() == 3) {
				findZDaughters(((*genp)[i]), genpvec, Zcount);
				Zcount++;
			}
		}
		if (verbose) std::cerr << "how many gen particles : " << genpvec.size() << std::endl;
		for (unsigned j = 0; j < genpvec.size(); j++) {
			if (verbose) std::cerr << "j: " << j
					<< "\tpdg: " << genpvec[j].first.pdgId()
					<< "\tpt: " << elevec[j].first.pt()
					<< "\teta: " << elevec[j].first.eta()
					<< "\tphi: " << elevec[j].first.phi()
					<< "\tlabel: " << elevec[j].second
					<< std::endl;
			MyGenParticle* g = new MyGenParticle();
			TLorentzVector tvec(0, 0, 0, 0);
			tvec.SetPtEtaPhiE(genpvec[j].first.pt(), genpvec[j].first.eta(),genpvec[j].first.phi(), genpvec[j].first.energy());
			g->pdg = genpvec[j].first.pdgId();
			g->p4 = tvec;
			g->label = genpvec[j].second;
			genpVector->push_back(g);
			if (g->label != old_label) {
				old_label = g->label;
				nlabels++;
			}
		}
		if (verbose) std::cerr << "checking again ... " << std::endl;
		for (unsigned y = 0; y < genpVector->size(); y++) {
			if (verbose) std::cerr << "y: " << y << "\tpdg: " << (*genpVector)[y]->pdg
					<< "\tpt: " << (*genpVector)[y]->p4.Pt() << std::endl;
		}
	} // do_ZZ

//########################### WJet MC ###########################
	if (do_wjetMC) {
		if (verbose) if (verbose) std::cerr << "clearing gen vecs ... " << std::endl;
		genpvec.clear();
		genpVector->clear();
		if (verbose) std::cerr << "gonna loop ... " << std::endl;
		unsigned nlabels = 0, old_label = 999;
		unsigned Wcount = 0;
		for (unsigned i = 0; i < genp->size(); i++) {
			if (verbose) std::cerr << "i: " << i << "\tpdg: " << ((*genp)[i]).pdgId()
							<< "\tstatus: " << ((*genp)[i]).status() << "\tpt: "
							<< ((*genp)[i]).pt() << "\teta: " << ((*genp)[i]).eta()
							<< std::endl;
			if (abs(((*genp)[i]).pdgId()) == 24 && ((*genp)[i]).status() == 3) {
				findWDaughters(((*genp)[i]), genpvec, Wcount);
				Wcount++;
			}
		}
		if (verbose) std::cerr << "how many gen particles : " << genpvec.size() << std::endl;
		for (unsigned j = 0; j < genpvec.size(); j++) {
			if (verbose) std::cerr << "j: " << j
					<< "\tpdg: " << genpvec[j].first.pdgId()
					<< "\tpt: " << elevec[j].first.pt()
					<< "\teta: " << elevec[j].first.eta()
					<< "\tphi: " << elevec[j].first.phi()
					<< "\tlabel: " << elevec[j].second
					<< std::endl;
			MyGenParticle* g = new MyGenParticle();
			TLorentzVector tvec(0, 0, 0, 0);
			tvec.SetPtEtaPhiE(genpvec[j].first.pt(), genpvec[j].first.eta(),genpvec[j].first.phi(), genpvec[j].first.energy());
			g->pdg = genpvec[j].first.pdgId();
			g->p4 = tvec;
			g->label = genpvec[j].second;
			genpVector->push_back(g);
			if (g->label != old_label) {
				old_label = g->label;
				nlabels++;
			}
		}
		if (verbose) std::cerr << "checking again ... " << std::endl;
		for (unsigned y = 0; y < genpVector->size(); y++) {
			if (verbose) std::cerr << "y: " << y << "\tpdg: " << (*genpVector)[y]->pdg << "\tpt: " << (*genpVector)[y]->p4.Pt() << std::endl;
		}
	} // do_Wjet

//########################### ZJet MC ###########################
	if (do_zjetMC) {
		if (verbose) if (verbose) std::cerr << "clearing gen vecs ... " << std::endl;
		genpvec.clear();
		genpVector->clear();
		if (verbose) std::cerr << "gonna loop ... " << std::endl;
		unsigned nlabels = 0, old_label = 999;
		unsigned Zcount = 0;
		for (unsigned i = 0; i < genp->size(); i++) {
			if (verbose) std::cerr << "i: " << i << "\tpdg: " << ((*genp)[i]).pdgId()
							<< "\tstatus: " << ((*genp)[i]).status() << "\tpt: "
							<< ((*genp)[i]).pt() << "\teta: " << ((*genp)[i]).eta()
							<< std::endl;
			if (abs(((*genp)[i]).pdgId()) == 23 && ((*genp)[i]).status() == 3) {
				findZDaughters(((*genp)[i]), genpvec, Zcount);
				Zcount++;
			}
		}
		if (verbose) std::cerr << "how many gen particles : " << genpvec.size() << std::endl;
		for (unsigned j = 0; j < genpvec.size(); j++) {
			if (verbose) std::cerr << "j: " << j
					<< "\tpdg: " << genpvec[j].first.pdgId()
					<< "\tpt: " << elevec[j].first.pt()
					<< "\teta: " << elevec[j].first.eta()
					<< "\tphi: " << elevec[j].first.phi()
					<< "\tlabel: " << elevec[j].second
					<< std::endl;
			MyGenParticle* g = new MyGenParticle();
			TLorentzVector tvec(0, 0, 0, 0);
			tvec.SetPtEtaPhiE(genpvec[j].first.pt(), genpvec[j].first.eta(),genpvec[j].first.phi(), genpvec[j].first.energy());
			g->pdg = genpvec[j].first.pdgId();
			g->p4 = tvec;
			g->label = genpvec[j].second;
			genpVector->push_back(g);
			if (g->label != old_label) {
				old_label = g->label;
				nlabels++;
			}
		}
		if (verbose) std::cerr << "checking again ... " << std::endl;
		for (unsigned y = 0; y < genpVector->size(); y++) {
			if (verbose) std::cerr << "y: " << y << "\tpdg: " << (*genpVector)[y]->pdg << "\tpt: " << (*genpVector)[y]->p4.Pt() << std::endl;
		}
	} // do_Zjet

//########################### Wenu MC ###########################
	if (do_wenuMC) {
		elevec.clear();
		genelectronVector->clear();
		unsigned nlabels = 0, old_label = 999;
		for (unsigned i = 0; i < genp->size(); i++) {
			if (abs(((*genp)[i]).pdgId()) == 24 && ((*genp)[i]).status() == 3) {
				unsigned label;
				findWElectronsMC(((*genp)[i]), elevec, label);
				if (verbose) std::cerr << "how many electrons : " << elevec.size() << std::endl;
				for (unsigned j = 0; j < elevec.size(); j++) {
					if (verbose) std::cerr << "j: " << j
							<< "\tpt: " << elevec[j].first.pt()
							<< "\teta: " << elevec[j].first.eta()
							<< "\tphi: " << elevec[j].first.phi()
							<< "\tlabel: " << elevec[j].second
							<< std::endl;
					MyGenParticle* gele = new MyGenParticle();
					TLorentzVector tvec(0, 0, 0, 0);
					tvec.SetPtEtaPhiE(elevec[j].first.pt(),elevec[j].first.eta(), elevec[j].first.phi(),elevec[j].first.energy());
					gele->p4 = tvec;
					gele->label = elevec[j].second;
					genelectronVector->push_back(gele);
					if (gele->label != old_label) {
						old_label = gele->label;
						nlabels++;
					}
				}
			}
		}
	} // do_wenu

//########################### Wmunu MC ###########################
	if (do_wenuMC) {
		muvec.clear();
		genmuonVector->clear();
		unsigned nlabels = 0, old_label = 999;
		for (unsigned i = 0; i < genp->size(); i++) {
			if (abs(((*genp)[i]).pdgId()) == 24 && ((*genp)[i]).status() == 3) {
				unsigned label;
				findWMuonsMC(((*genp)[i]), muvec, label);
				if (verbose) std::cerr << "how many muons : " << muvec.size() << std::endl;
				for (unsigned j = 0; j < muvec.size(); j++) {
					if (verbose) std::cerr << "j: " << j
							<< "\tpt: " << muvec[j].first.pt()
							<< "\teta: " << muvec[j].first.eta()
							<< "\tphi: " << muvec[j].first.phi()
							<< "\tlabel: " << muvec[j].second
							<< std::endl;
					MyGenParticle* gmu = new MyGenParticle();
					TLorentzVector tvec(0, 0, 0, 0);
					tvec.SetPtEtaPhiE(muvec[j].first.pt(),muvec[j].first.eta(), muvec[j].first.phi(),muvec[j].first.energy());
					gmu->p4 = tvec;
					gmu->label = muvec[j].second;
					genelectronVector->push_back(gmu);
					if (gmu->label != old_label) {
						old_label = gmu->label;
						nlabels++;
					}
				}
			}
		}
	} // do_wmunu

//########################### Zee MC ###########################
	if (do_zeeMC) {
		elevec.clear();
		genelectronVector->clear();
		unsigned nlabels = 0, old_label = 999;
		for (unsigned i = 0; i < genp->size(); i++) {
			if (abs(((*genp)[i]).pdgId()) == 23 && ((*genp)[i]).status() == 3) {
				unsigned label;
				findZElectronsMC(((*genp)[i]), elevec, label);
				if (verbose) std::cerr << "how many electrons : " << elevec.size() << std::endl;
				for (unsigned j = 0; j < elevec.size(); j++) {
					if (verbose) std::cerr << "j: " << j
							<< "\tpt: " << elevec[j].first.pt()
							<< "\teta: " << elevec[j].first.eta()
							<< "\tphi: " << elevec[j].first.phi()
							<< "\tlabel: " << elevec[j].second
							<< std::endl;
					MyGenParticle* gele = new MyGenParticle();
					TLorentzVector tvec(0, 0, 0, 0);
					tvec.SetPtEtaPhiE(elevec[j].first.pt(),elevec[j].first.eta(), elevec[j].first.phi(),elevec[j].first.energy());
					gele->p4 = tvec;
					gele->label = elevec[j].second;
					genelectronVector->push_back(gele);
					if (gele->label != old_label) {
						old_label = gele->label;
						nlabels++;
					}
				}
			}
		}
	} // do_zeeMC

//########################### Zmumu MC ###########################
	if (do_zmumuMC) {
		muvec.clear();
		genmuonVector->clear();
		unsigned nlabels = 0, old_label = 999;
		for (unsigned i = 0; i < genp->size(); i++) {
			if (abs(((*genp)[i]).pdgId()) == 23 && ((*genp)[i]).status() == 3) {
				unsigned label;
				findZMuonsMC(((*genp)[i]), muvec, label);
				if (verbose) std::cerr << "how many muons : " << muvec.size() << std::endl;
				for (unsigned j = 0; j < muvec.size(); j++) {
					if (verbose) std::cerr << "j: " << j
							<< "\tpt: " << muvec[j].first.pt()
							<< "\teta: " << muvec[j].first.eta()
							<< "\tphi: " << muvec[j].first.phi()
							<< "\tlabel: " << muvec[j].second
							<< std::endl;
					MyGenParticle* gmu = new MyGenParticle();
					TLorentzVector tvec(0, 0, 0, 0);
					tvec.SetPtEtaPhiE(muvec[j].first.pt(),muvec[j].first.eta(), muvec[j].first.phi(),muvec[j].first.energy());
					gmu->p4 = tvec;
					gmu->label = muvec[j].second;
					genmuonVector->push_back(gmu);
					if (gmu->label != old_label) {
						old_label = gmu->label;
						nlabels++;
					}
				}
			}
		}
	} // do_zmumuMC

//########################### H2EJ MC ###########################
	if (do_signalMC) {
		if (verbose) std::cerr << "IN SIGNAL SECTION\n";
		if (verbose) std::cerr << "clearing gen vecs ... " << std::endl;
		elevec.clear();
		hd0vec.clear();
		genpvec.clear();
		genSigMetVector->clear();
		genpVector->clear();
		genelectronVector->clear();
		unsigned nlabels = 0, nWlabels = 0, old_label = 999;
//		unsigned Wcount(0);
		unsigned int Zcount(0);
		if (verbose) std::cerr << "Finding Particles" << std::endl;
		for (unsigned i = 0; i < genp->size(); i++) {
			// change for pythia8
//			if (abs(((*genp)[i]).pdgId()) == 24	&& ((*genp)[i]).status() == 62) {	// Loop over W's  -- not sure what status 62 is?  ===> WH not looked at for now
//				findWDaughters(((*genp)[i]), genpvec, Wcount);
//				Wcount++;
//			}
//			else
			if (((*genp)[i]).pdgId() == 23	&& ((*genp)[i]).status() == 62) {	// Loop over Z's
				findZDaughters(((*genp)[i]), genpvec, Zcount);
				Zcount++;
			}
			// change for pythia8
			if (abs(((*genp)[i]).pdgId()) == 25	&& ((*genp)[i]).status() == 62) { // If Higgs
//				std::cerr << "Finding Higgs properties" << std::endl;
				unsigned label;
				findSignalElectrons(((*genp)[i]), elevec, label); // Going down through the decay chain to find e's
				findSignalMET(((*genp)[i]), hd0vec, label);
				if (verbose) std::cerr << "how many hd0 : " << hd0vec.size() << std::endl;
				if (verbose) std::cerr << "how many electrons : " << elevec.size() << std::endl;
				for (unsigned j = 0; j < elevec.size(); j++) {	//Signal Electrons in EJ
					if (verbose){
						std::cerr << "j: " << j
								<< "\tpt: " << elevec[j].first.pt()
								<< "\teta: " << elevec[j].first.eta()
								<< "\tphi: " << elevec[j].first.phi()
								<< "\tlabel: " << elevec[j].second
								<< std::endl;
					}
					MyGenParticle* gele = new MyGenParticle();
					TLorentzVector tvec(0, 0, 0, 0);
					tvec.SetPtEtaPhiE(elevec[j].first.pt(),elevec[j].first.eta(), elevec[j].first.phi(),elevec[j].first.energy());
					gele->p4 = tvec;
					gele->label = elevec[j].second;
					genelectronVector->push_back(gele);
					if (gele->label != old_label) {	//we have switched to a different EJ
						old_label = gele->label;
						nlabels++;
					}
				}
				for (unsigned h = 0; h < hd0vec.size(); h++) {	//Signal MET in EJ
					MyGenParticle* ghd0 = new MyGenParticle();
					TLorentzVector tvec(0, 0, 0, 0);
					tvec.SetPtEtaPhiE(hd0vec[h].first.pt(),hd0vec[h].first.eta(), hd0vec[h].first.phi(),hd0vec[h].first.energy());
					ghd0->p4 = tvec;
					ghd0->label = hd0vec[h].second;
					genSigMetVector->push_back(ghd0);
				}
			}
		}

		if (verbose) std::cerr << "how many gen particles : " << genpvec.size() << std::endl;
		for (unsigned j = 0; j < genpvec.size(); j++) {	// Filled from findVDaughters
			if (verbose){
				std::cerr << "j: " << j
						<< "\tpdg: " << genpvec[j].first.pdgId()
						<< "\tpt: " << genpvec[j].first.pt()
						<< "\teta: " << genpvec[j].first.eta()
						<< "\tphi: " << genpvec[j].first.phi()
						<< "\tlabel: " << genpvec[j].second
						<< std::endl;
			}
			MyGenParticle* g = new MyGenParticle();
			TLorentzVector tvec(0, 0, 0, 0);
			tvec.SetPtEtaPhiE(genpvec[j].first.pt(), genpvec[j].first.eta(),genpvec[j].first.phi(), genpvec[j].first.energy());
			g->pdg = genpvec[j].first.pdgId();
			g->p4 = tvec;
			g->label = genpvec[j].second;
			genpVector->push_back(g);
			if (g->label != old_label) {
				old_label = g->label;
				nWlabels++;
			}
		}

		if (verbose) std::cerr << "checking again ... " << std::endl;
		for (unsigned y = 0; y < genpVector->size(); y++) {
			if (verbose) std::cerr << "y: " << y << "\tpdg: " << (*genpVector)[y]->pdg << "\tpt: " << (*genpVector)[y]->p4.Pt() << std::endl;
		}
		if (nlabels != 2) {
			if (verbose) std::cerr << "!!!!!!!!!! nlabels != 2 !!!!!!!!!!!\nnlabels : " << nlabels << std::endl;
			return;
		}

		// find smallest cone that will fit the electrons ...
		genjetVector->clear();
//		calogenjetVector->clear();
		if (elevec.size()) {	// The (genparticle,int) vector filled by findSignalElectrons
			std::vector<double> energyFromEles;
			unsigned jetnum = -1;
			std::vector < math::XYZTLorentzVector > centroids;
			std::vector<unsigned> neles;
			if (verbose) std::cout << "looping over gele's ... " << std::endl;
			for (unsigned j = 0; j < elevec.size(); j++) {
				if (verbose) std::cout << "\tele: " << j << "\tlabel: " << elevec[j].second << std::endl;
				if (elevec[j].second != jetnum) {
					centroids.push_back(math::XYZTLorentzVector(0, 0, 0, 0));
					neles.push_back(0);
					energyFromEles.push_back(0);
					jetnum = elevec[j].second;	//to only find the two EJs
				}
				centroids[jetnum] += math::XYZTLorentzVector(elevec[j].first.px(), elevec[j].first.py(),elevec[j].first.pz(), elevec[j].first.p()); // massless hypothesis ;
				neles[jetnum] += 1;
				energyFromEles[jetnum]+=elevec[j].first.energy();
			} // loop over gen eles ...

			if (verbose) std::cout << "N centroids : " << centroids.size() << std::endl;
			for (unsigned j = 0; j < centroids.size(); j++) {
				if (verbose) std::cerr << "centroid " << j << "\tpt: " << centroids[j].Pt() << "\teta: " << centroids[j].Eta() << "\tphi: " << centroids[j].Phi() << std::endl;
				if (verbose) std::cerr << centroids[j] << std::endl;
				TVector3 axis;
				axis.SetPtEtaPhi(centroids[j].Pt(), centroids[j].Eta(), centroids[j].Phi());
				double maxdR = -1;
				unsigned maxindex(0);
				double maxdRminpT = -1;
				unsigned maxindexminpT;
				for (unsigned k = 0; k < elevec.size(); k++) {
					if (verbose) std::cerr << "\tlooping over eles ... " << std::endl;
					if (elevec[k].second == j) {
						TVector3 e(elevec[k].first.px(), elevec[k].first.py(), elevec[k].first.pz());
						double dR = axis.DrEtaPhi(e);
						if (dR > maxdR) {
							maxdR = dR;
							maxindex = k;
						}
						if (elevec[k].first.pt() >= 2.) {
							if (dR > maxdRminpT) {
								maxdRminpT = dR;
								maxindexminpT = k;
							}
						}
					}
				} // loop over gen eles ...
				if (verbose) std::cout << "\t--> max dist: " << maxdR << "\t--> max dist min pT: " << maxdRminpT << "\t index: " << maxindex << std::endl;
				// save ...
//				if (verbose) std::cout << "energyFromEles: " << energyFromEles[j] << std::endl;
				double sigMET(0);
				TLorentzVector sigMETp4(0,0,0,0);
				if (verbose) std::cout << "genSigMetVector.size(): " << genSigMetVector->size() << std::endl;
				for (unsigned metIndex = 0; metIndex < genSigMetVector->size(); metIndex++) {
					MyGenParticle* m = genSigMetVector->at(metIndex);
					if (m->label == j){
						sigMET+=m->p4.Pt();
						sigMETp4 += m->p4;
//						sigMETp4 += math::XYZTLorentzVector(m->p4.Px(),m->p4.Py(),m->p4.Pz(),m->p4.Energy());
					}
				}
				if (verbose) std::cout << "sigMET: " << sigMET << std::endl;
				MyPFJet* genjet = new MyPFJet();
				genjet->p4 = TLorentzVector(axis, energyFromEles[j]);
				genjet->dR = maxdR;
				genjet->dRminpT = maxdRminpT;
				genjet->nele = neles[j];
				genjet->met = sigMET;
				genjet->metp4 = sigMETp4;
//				genjet->metp4 = math::XYZTLorentzVector(sigMETp4.p4.Px(),sigMETp4.p4.Py(),sigMETp4.p4.Pz(),sigMETp4.p4.Energy());
				genjet->index = j;
				genjetVector->push_back(genjet);
				if (verbose) std::cout << "genjet->p4.Energy(): " << genjet->p4.Energy() << std::endl;
				if (verbose) std::cout << "genjetVector->size(): " << genjetVector->size() << std::endl;
//				MyCaloJet* calogenjet = new MyCaloJet();
//				calogenjet->p4 = TLorentzVector(axis, 0);
//				calogenjet->dR = maxdR;
//				calogenjet->dRminpT = maxdRminpT;
//				calogenjet->nele = neles[j];
//				calogenjet->index = j;
//				calogenjetVector->push_back(calogenjet);
			} // loop over centroids
		} // if we found some eles
		if (verbose) std::cout << "--- gjetvec size: " << genjetVector->size() << std::endl;
	} // do_signal
//########################### End of MC ###########################
//#################################################################


	// save evt info
	evtdata.run = iEvent.id().run();
	evtdata.event = iEvent.id().event();
	evtdata.lumi = iEvent.id().luminosityBlock();

//########################### B-tagging ###########################
	const reco::JetTagCollection & bTags_secvtx = *(bTag_secvtx_Handle.product());
	const reco::JetTagCollection & bTags_softele = *(bTag_softele_Handle.product());
	const reco::JetTagCollection & bTags_jetprob = *(bTag_jetprob_Handle.product());

//########################### Photons ###########################
	photonVector->clear();
	for (unsigned i = 0; i < Photons->size(); i++) {
		MyPhoton* photon = new MyPhoton;
		photon->p4 = TLorentzVector(0, 0, 0, 0);
		photon->p4.SetPtEtaPhiE((*Photons)[i].pt(), (*Photons)[i].eta(),(*Photons)[i].phi(), (*Photons)[i].energy());
		photon->Iecal = (*Photons)[i].ecalRecHitSumEtConeDR04();
		photon->Ihcal = (*Photons)[i].hcalTowerSumEtConeDR04();
		photon->Itrk = (*Photons)[i].trkSumPtHollowConeDR04();
		photon->hOe = (*Photons)[i].hadronicOverEm();
		photon->sigieie = (*Photons)[i].sigmaIetaIeta();
		photon->hasPixSeed = (*Photons)[i].hasPixelSeed();
		photon->isEB = (*Photons)[i].isEB();
		photon->isEE = (*Photons)[i].isEE();
		photonVector->push_back(photon);
	}

//########################### MET ###########################
	const reco::PFMET pfMet = (*pfMetCollection)[0];
	const reco::CaloMET caloMet = (*caloMetCollection)[0];
	met.pf = TVector2(pfMet.px(), pfMet.py());
	met.calo = TVector2(caloMet.px(), caloMet.py());

//########################### ak5CaloJets ###########################
	if (verbose) std::cout << "--- ak5caloJets: " << ak5caloJets->size() << std::endl;
	ak5calojetVector->clear();
	for (unsigned i = 0; i < ak5caloJets->size(); i++) {
		int ntrks = reco::JetTracksAssociation::tracksNumber(*ak5assocCaloJetTracks, (*ak5caloJets)[i]);
		float sumPt = reco::JetTracksAssociation::tracksP4(*ak5assocCaloJetTracks,(*ak5caloJets)[i]).pt();
		reco::TrackRefVector calotracks = reco::JetTracksAssociation::getValue(*ak5assocCaloJetTracks, (*ak5caloJets)[i]);
		reco::CaloJet correctedJet = (*ak5caloJets)[i];
//		double scale = corrector->correction((*ak5caloJets)[i].p4());  //calculate the correction - BRS: comment out, add next line
		double scale = 1.0;
		correctedJet.scaleEnergy(scale); // apply correction

		if (verbose) std::cout << "\tpt: " << (*ak5caloJets)[i].pt()
								<< "\tpt (cor): " << correctedJet.pt()
								<< "\teta: " << (*ak5caloJets)[i].eta()
								<< "\tphi: " << (*ak5caloJets)[i].phi()
								<< "\temf: " << (*ak5caloJets)[i].emEnergyFraction()
								<< "\tcr: " << sumPt / (*ak5caloJets)[i].pt()
								<< std::endl;

		unsigned int numEle(0), numFixedEle(0);
		if (verbose) std::cout << "\t###gsfElectrons->size(): " << gsfElectrons->size() << std::endl;
		for (size_t ele = 0; ele < gsfElectrons->size(); ++ele) {
			double deltR = deltaR((*ak5caloJets)[i].p4(),(*gsfElectrons)[ele].p4());
			if (verbose) std::cout << "\t\t###deltR: " << deltR << std::endl;
			if (deltR < (*ak5caloJets)[i].maxDistance()) ++numEle;   // electron in Jet ==> within jet's max dist
			if (deltR < 0.5) ++numFixedEle;   // electron in Jet ==> within fixed cone size
		}
		if (verbose) std::cout << "\t###finished filling numEle: " << numEle << " and numFixedEle: " << numFixedEle << std::endl;

		MyCaloJet* calojet = new MyCaloJet;
		calojet->p4 = TLorentzVector(0, 0, 0, 0);
		calojet->p4.SetPtEtaPhiE((*ak5caloJets)[i].pt(), (*ak5caloJets)[i].eta(), (*ak5caloJets)[i].phi(), (*ak5caloJets)[i].energy());
		calojet->p4cor.SetPtEtaPhiE(correctedJet.pt(), correctedJet.eta(), correctedJet.phi(), correctedJet.energy());
		calojet->energy = (*ak5caloJets)[i].energy();
		calojet->et = (*ak5caloJets)[i].et();
		calojet->area = (*ak5caloJets)[i].jetArea();
		calojet->etaeta = (*ak5caloJets)[i].etaetaMoment();
		calojet->phiphi = (*ak5caloJets)[i].phiphiMoment();
		calojet->etaphi = (*ak5caloJets)[i].etaphiMoment();
		calojet->nele = numEle;
		calojet->neleFixed = numFixedEle;
		calojet->emf = (*ak5caloJets)[i].emEnergyFraction();
		calojet->ntrks = ntrks;
		calojet->dR = (*ak5caloJets)[i].maxDistance();
		calojet->crtrans = sumPt / ((*ak5caloJets)[i].et() * (*ak5caloJets)[i].emEnergyFraction());
		calojet->cr = sumPt / ((*ak5caloJets)[i].energy() * (*ak5caloJets)[i].emEnergyFraction());
		calojet->Ical = (*ak5caloJets)[i].etInAnnulus(0.2, 0.4);
		calojet->Et_0_2 = (*ak5caloJets)[i].etInAnnulus(0.0, 0.2);
		double Itrk = 0.0;
		unsigned ntrkAbove1 = 0, ntrkAbove2_5 = 0, ntrkAbove5 = 0;
		for (unsigned t = 0; t < calotracks.size(); ++t) {
			const reco::Track& track = *(calotracks[t]);
			if (track.pt() < 0.5) continue;
			if (track.pt() >= 1.0) ntrkAbove1++;
			if (track.pt() >= 2.5) ntrkAbove2_5++;
			if (track.pt() >= 5.0) ntrkAbove5++;
			TVector3 tvec(0, 0, 0);
			tvec.SetPtEtaPhi(track.pt(), track.eta(), track.phi());
			double dR = calojet->p4.Vect().DrEtaPhi(tvec);
			if (dR < 0.4 && dR > 0.2) Itrk += track.pt();	// BRS: where do these numbers come from? Ical Itrk
		}
		calojet->Itrk = Itrk;
		calojet->ntrkAbove1 = ntrkAbove1;
		calojet->ntrkAbove2_5 = ntrkAbove2_5;
		calojet->ntrkAbove5 = ntrkAbove5;

		// Loop over btags
		bool found_btag_info = false;
		for (unsigned int b = 0; b != bTags_secvtx.size(); ++b) {
			if ((*ak5caloJets)[i].pt() == bTags_secvtx[b].first->pt()) {
				//	 if (verbose) std::cout << "-->found btag info ... discrim: " << bTags_secvtx[b].second << std::endl;
				calojet->btag_secvtx = bTags_secvtx[b].second;
				found_btag_info = true;
				break;
			}
		}
		// Loop over btags
		found_btag_info = false;
		for (unsigned int b = 0; b != bTags_softele.size(); ++b) {
			if ((*ak5caloJets)[i].pt() == bTags_softele[b].first->pt()) {
				//	 if (verbose) std::cout << "-->found btag info ... discrim: " << bTags_softele[b].second << std::endl;
				calojet->btag_softele = bTags_softele[b].second;
				found_btag_info = true;
				break;
			}
		}
		// Loop over btags
		found_btag_info = false;
		for (unsigned int b = 0; b != bTags_jetprob.size(); ++b) {
			if ((*ak5caloJets)[i].pt() == bTags_jetprob[b].first->pt()) {
				//	 if (verbose) std::cout << "-->found btag info ... discrim: " << bTags_jetprob[b].second << std::endl;
				calojet->btag_jetprob = bTags_jetprob[b].second;
				found_btag_info = true;
				break;
			}
		}
		if (!(found_btag_info)) {
			if (verbose) std::cout << "no btag info found for this jet ... " << std::endl;
		}
		ak5calojetVector->push_back(calojet);
		if (verbose) std::cout << "calo EM: " << (*ak5caloJets)[i].energy() * (*ak5caloJets)[i].emEnergyFraction() << std::endl;
		if (verbose) std::cout << "calo HAD: " << (*ak5caloJets)[i].energy() * (1 - (*ak5caloJets)[i].emEnergyFraction()) << std::endl;
	}

//########################### ak7CaloJets ###########################
	if (verbose) std::cout << "--- ak7caloJets: " << ak7caloJets->size() << std::endl;
	ak7calojetVector->clear();
	for (unsigned i = 0; i < ak7caloJets->size(); i++) {
		int ntrks = reco::JetTracksAssociation::tracksNumber(*ak7assocCaloJetTracks, (*ak7caloJets)[i]);
		float sumPt = reco::JetTracksAssociation::tracksP4(*ak7assocCaloJetTracks,(*ak7caloJets)[i]).pt();
		reco::TrackRefVector calotracks = reco::JetTracksAssociation::getValue(*ak7assocCaloJetTracks, (*ak7caloJets)[i]);
		reco::CaloJet correctedJet = (*ak7caloJets)[i];
//		double scale = corrector->correction((*ak7caloJets)[i].p4());  //calculate the correction - BRS: comment out, add next line
		double scale = 1.0;
		correctedJet.scaleEnergy(scale); // apply correction

		if (verbose) std::cout << "\tpt: " << (*ak7caloJets)[i].pt()
								<< "\tpt (cor): " << correctedJet.pt()
								<< "\teta: " << (*ak7caloJets)[i].eta()
								<< "\tphi: " << (*ak7caloJets)[i].phi()
								<< "\temf: " << (*ak7caloJets)[i].emEnergyFraction()
								<< "\tcr: " << sumPt / (*ak7caloJets)[i].pt()
								<< std::endl;

		unsigned int numEle(0), numFixedEle(0);
		if (verbose) std::cout << "\t###gsfElectrons->size(): " << gsfElectrons->size() << std::endl;
		for (size_t ele = 0; ele < gsfElectrons->size(); ++ele) {
			double deltR = deltaR((*ak7caloJets)[i].p4(),(*gsfElectrons)[ele].p4());
			if (deltR < (*ak7caloJets)[i].maxDistance()) ++numEle;   // electron in Jet ==> within jet's max dist
			if (deltR < 0.7) ++numFixedEle;   // electron in Jet ==> within fixed cone size
		}
		if (verbose) std::cout << "\t###finished filling numEle: " << numEle << " and numFixedEle: " << numFixedEle << std::endl;

		MyCaloJet* calojet = new MyCaloJet;
		calojet->p4 = TLorentzVector(0, 0, 0, 0);
		calojet->p4.SetPtEtaPhiE((*ak7caloJets)[i].pt(), (*ak7caloJets)[i].eta(), (*ak7caloJets)[i].phi(), (*ak7caloJets)[i].energy());
		calojet->p4cor.SetPtEtaPhiE(correctedJet.pt(), correctedJet.eta(), correctedJet.phi(), correctedJet.energy());
		calojet->energy = (*ak7caloJets)[i].energy();
		calojet->et = (*ak7caloJets)[i].et();
		calojet->area = (*ak7caloJets)[i].jetArea();
		calojet->etaeta = (*ak7caloJets)[i].etaetaMoment();
		calojet->phiphi = (*ak7caloJets)[i].phiphiMoment();
		calojet->etaphi = (*ak7caloJets)[i].etaphiMoment();
		calojet->nele = numEle;
		calojet->neleFixed = numFixedEle;
		calojet->emf = (*ak7caloJets)[i].emEnergyFraction();
		calojet->ntrks = ntrks;
		calojet->dR = (*ak7caloJets)[i].maxDistance();
		calojet->crtrans = sumPt / ((*ak7caloJets)[i].et() * (*ak7caloJets)[i].emEnergyFraction());
		calojet->cr = sumPt / ((*ak7caloJets)[i].energy() * (*ak7caloJets)[i].emEnergyFraction());
		calojet->Ical = (*ak7caloJets)[i].etInAnnulus(0.2, 0.4);
		calojet->Et_0_2 = (*ak7caloJets)[i].etInAnnulus(0.0, 0.2);
		double Itrk = 0.0;
		unsigned ntrkAbove1 = 0, ntrkAbove2_5 = 0, ntrkAbove5 = 0;
		for (unsigned t = 0; t < calotracks.size(); ++t) {
			const reco::Track& track = *(calotracks[t]);
			if (track.pt() < 0.5) continue;
			if (track.pt() >= 1.0) ntrkAbove1++;
			if (track.pt() >= 2.5) ntrkAbove2_5++;
			if (track.pt() >= 5.0) ntrkAbove5++;
			TVector3 tvec(0, 0, 0);
			tvec.SetPtEtaPhi(track.pt(), track.eta(), track.phi());
			double dR = calojet->p4.Vect().DrEtaPhi(tvec);
			if (dR < 0.4 && dR > 0.2) Itrk += track.pt();	// BRS: where do these numbers come from? Ical Itrk
		}
		calojet->Itrk = Itrk;
		calojet->ntrkAbove1 = ntrkAbove1;
		calojet->ntrkAbove2_5 = ntrkAbove2_5;
		calojet->ntrkAbove5 = ntrkAbove5;

		// Loop over btags
		bool found_btag_info = false;
		for (unsigned int b = 0; b != bTags_secvtx.size(); ++b) {
			if ((*ak7caloJets)[i].pt() == bTags_secvtx[b].first->pt()) {
				//	 if (verbose) std::cout << "-->found btag info ... discrim: " << bTags_secvtx[b].second << std::endl;
				calojet->btag_secvtx = bTags_secvtx[b].second;
				found_btag_info = true;
				break;
			}
		}
		// Loop over btags
		found_btag_info = false;
		for (unsigned int b = 0; b != bTags_softele.size(); ++b) {
			if ((*ak7caloJets)[i].pt() == bTags_softele[b].first->pt()) {
				//	 if (verbose) std::cout << "-->found btag info ... discrim: " << bTags_softele[b].second << std::endl;
				calojet->btag_softele = bTags_softele[b].second;
				found_btag_info = true;
				break;
			}
		}
		// Loop over btags
		found_btag_info = false;
		for (unsigned int b = 0; b != bTags_jetprob.size(); ++b) {
			if ((*ak7caloJets)[i].pt() == bTags_jetprob[b].first->pt()) {
				//	 if (verbose) std::cout << "-->found btag info ... discrim: " << bTags_jetprob[b].second << std::endl;
				calojet->btag_jetprob = bTags_jetprob[b].second;
				found_btag_info = true;
				break;
			}
		}
		if (!(found_btag_info)) {
			if (verbose) std::cout << "no btag info found for this jet ... " << std::endl;
		}
		ak7calojetVector->push_back(calojet);
		if (verbose) std::cout << "calo EM: " << (*ak7caloJets)[i].energy() * (*ak7caloJets)[i].emEnergyFraction() << std::endl;
		if (verbose) std::cout << "calo HAD: " << (*ak7caloJets)[i].energy() * (1 - (*ak7caloJets)[i].emEnergyFraction()) << std::endl;
	}

//########################### PFJets ###########################
	if (verbose) std::cout << "--- pfjets: " << pfJets->size() << std::endl;
	pfjetVector->clear();
	pfeleVector->clear();
	for (unsigned i = 0; i < pfJets->size(); i++) {
		reco::TrackRefVector pftracks = (*pfJets)[i].getTrackRefs();
		math::XYZTLorentzVector result(0, 0, 0, 0);
		for (unsigned t = 0; t < pftracks.size(); ++t) {
			// vector sum for axis
			const reco::Track& track = *(pftracks[t]);
			result += math::XYZTLorentzVector(track.px(), track.py(),
					track.pz(), track.p()); // massless hypothesis
		}
		float sumPt = result.pt();
		if (verbose) std::cout << "\tPF:: nele: " << (*pfJets)[i].electronMultiplicity()
				<< "\tnphoton: " << (*pfJets)[i].photonMultiplicity()
				<< "\teleFrac: " << (*pfJets)[i].electronEnergyFraction()
				<< std::endl;

		MyPFJet* pfjet = new MyPFJet;
		pfjet->p4 = TLorentzVector(0, 0, 0, 0);
		pfjet->p4.SetPtEtaPhiE((*pfJets)[i].pt(), (*pfJets)[i].eta(),(*pfJets)[i].phi(), (*pfJets)[i].energy());
		pfjet->area = (*pfJets)[i].jetArea();
		pfjet->etaeta = (*pfJets)[i].etaetaMoment();
		pfjet->phiphi = (*pfJets)[i].phiphiMoment();
		pfjet->etaphi = (*pfJets)[i].etaphiMoment();
		pfjet->et = (*pfJets)[i].et();
		pfjet->ntrks = pftracks.size(); //(*pfJets)[i].electronMultiplicity();
		pfjet->nele = (*pfJets)[i].electronMultiplicity();
		pfjet->nphoton = (*pfJets)[i].photonMultiplicity();
		pfjet->electron_energy = (*pfJets)[i].electronEnergy();
		pfjet->photon_energy = (*pfJets)[i].photonEnergy();
		pfjet->electron_energyFrac = (*pfJets)[i].electronEnergyFraction();
		pfjet->photon_energyFrac = (*pfJets)[i].photonEnergyFraction();
		pfjet->neutralHadFrac = (*pfJets)[i].neutralHadronEnergyFraction();
		pfjet->chargedHadFrac = (*pfJets)[i].chargedHadronEnergyFraction();
		pfjet->neutralEMFrac = (*pfJets)[i].neutralEmEnergyFraction();
		pfjet->chargedEMFrac = (*pfJets)[i].chargedEmEnergyFraction();
		pfjet->dR = (*pfJets)[i].maxDistance();
		pfjet->emf = (*pfJets)[i].chargedEmEnergyFraction()	+ (*pfJets)[i].neutralEmEnergyFraction();
//		pfjet->cr = sumPt/((*pfJets)[i].chargedEmEnergy() + (*pfJets)[i].neutralEmEnergy());
		pfjet->cr = sumPt / ((*pfJets)[i].chargedEmEnergy()	+ (*pfJets)[i].neutralEmEnergy());
		pfjet->crtrans = sumPt / ((*pfJets)[i].et() * (*pfJets)[i].chargedEmEnergyFraction()+(*pfJets)[i].neutralEmEnergyFraction());
		pfjet->Ical = (*pfJets)[i].etInAnnulus(0.2, 0.4);
		pfjet->Et_0_2 = (*pfJets)[i].etInAnnulus(0.0, 0.2);

		pfjet->gsfeleindices = 0x0;
		unsigned ntrks = 0, mynele = 0;
		unsigned save_pfele_vecsize = pfeleVector->size();
		float ecalConstitEnergy = 0, hcalConstitEnergy = 0;
		for (unsigned c = 0; c < (*pfJets)[i].getPFConstituents().size(); c++) {
			reco::PFCandidatePtr ptr = (*pfJets)[i].getPFConstituents()[c];
			ecalConstitEnergy += ptr->ecalEnergy();
			hcalConstitEnergy += ptr->hcalEnergy();
			if (abs(ptr->pdgId()) == 11) {
				mynele++;
				if (verbose) std::cerr << "-------------------------" << std::endl;
				if (verbose) std::cerr << "\tconstit is an ele ... " << std::endl;

				reco::GsfTrackRef gtrk = ptr->gsfTrackRef();
				if (gtrk.isNull()) {
					if (verbose) std::cerr << "\tpf gsftrk ref not valid, skipping ... "	<< std::endl;
					continue;
				}
				ntrks++;

				reco::GsfElectronRef geleref = ptr->gsfElectronRef();
				if (geleref.isNull()){	if (verbose) std::cerr << "\tpf gele ref not valid, skipping ... " << std::endl;}
					//	   continue;
				else{	if (verbose) std::cerr << "\tpf gele IS  valid, skipping ... " << std::endl;}

				bool matchedEle = false;
				unsigned gsfIndex(0);
				reco::GsfElectron gele;
				for (unsigned e = 0; e < gsfElectrons->size(); e++) {
					if (gtrk == (*gsfElectrons)[e].gsfTrack()) {
						if (verbose) std::cerr << "\t\tmatched to gsfEle " << e << std::endl;
						pfjet->gsfeleindices |= (1 << e );
						gsfIndex = e;
						matchedEle = true;
						gele = (*gsfElectrons)[e];
						if (verbose) std::cerr << "\t\tPFstatus: " << (*gsfElectrons)[e].mvaOutput().status << std::endl;
						if (verbose) std::cerr << "\t\tECAL driven: " << (*gsfElectrons)[e].ecalDriven() << std::endl;
						if (verbose) std::cerr << "\t\tECAL driven seed: " << (*gsfElectrons)[e].ecalDrivenSeed() << std::endl;
						if (verbose) std::cerr << "\t\tTK driven seed: " << (*gsfElectrons)[e].trackerDrivenSeed() << std::endl;
						if (verbose) std::cerr << "\t\tECAL presel: " << (*gsfElectrons)[e].passingCutBasedPreselection() << std::endl;
						if (verbose) std::cerr << "\t\tPF presel: " << (*gsfElectrons)[e].passingMvaPreselection() << std::endl;
						//LorentzVector elep4 = (*gsfElectrons)[e].p4(reco::GsfElectron::P4_COMBINATION);
						if (verbose) std::cerr << "\t\tmva: " << (*gsfElectrons)[e].mvaOutput().mva << std::endl;
						if (verbose) std::cerr << "\t\tpt: " << (*gsfElectrons)[e].p4().Pt() << std::endl;
						if (verbose) std::cerr << "\t\tconvR: " << (*gsfElectrons)[e].convRadius() << std::endl;
					}
				}
				if (!matchedEle) {
					if (verbose) std::cerr << "\t no matched GSF! " << std::endl;
					continue;
				}
				else if (verbose) std::cerr << "\tmatched to GSF ele! " << std::endl;

				MyPFElectron* pfele = new MyPFElectron;
				pfele->gsfIndex = gsfIndex;
				pfele->ecalDriven = gele.ecalDriven();
				pfele->ecalDrivenSeed = gele.ecalDrivenSeed();
				pfele->tkDrivenSeed = gele.trackerDrivenSeed();
				pfele->passMvaPreselection = gele.passingMvaPreselection();
				pfele->passCutBasedPreselection = gele.passingCutBasedPreselection();
				pfele->pt = gele.p4().Pt();
				pfele->mva_e_pi = gele.mvaOutput().mva;
				pfele->convR = gele.convRadius();

				unsigned hitmask = 0x0;
				unsigned nhits = 0;
				const reco::HitPattern& hp = gtrk->hitPattern();
				pfele->npixHits = hp.pixelLayersWithMeasurement();
				pfele->lostHits = gtrk->lost();
				for (int h = 0; h < hp.numberOfHits(); h++) {
					uint32_t hit = hp.getHitPattern(h);
					if (hp.validHitFilter(hit) && hp.pixelHitFilter(hit)) {
						nhits++;
						hitmask |= (1 << (hp.getLayer(hit) - 1));
						if (verbose) std::cerr << "\tvalid pix layer: " << hp.getLayer(hit) << std::endl;
					}
				}
				if (verbose) std::cerr << "\ttotal valid pix :  " << nhits << std::endl;
				if (verbose) std::cerr << "-------------------------" << std::endl;
				pfele->pixHitMask = hitmask;
				if (verbose) std::cout << "pushing back track ... " << std::endl;
				pfeleVector->push_back(pfele);
			}
		}

		if (pfjet->nele != mynele) {
			if (verbose) std::cerr << "NELE DIFFERS: " << std::endl
					<< "pf: " << pfjet->nele << std::endl
					<< "my: " << mynele << std::endl
					<< "quitting ... " << std::endl;
			exit(1);
		}
		if (pfjet->nele > 0) {
			pfjet->pfele_min_index = save_pfele_vecsize;
			pfjet->pfele_max_index = pfeleVector->size() - 1;
		}
		else {
			pfjet->pfele_min_index = -1;
			pfjet->pfele_max_index = -1;
		}
		if (verbose) std::cerr << "nele for jet " << i << " : " << pfjet->nele << std::endl
				<< "min index: " << pfjet->pfele_min_index << std::endl
				<< "max index: " << pfjet->pfele_max_index << std::endl;

		double Itrk = 0.;
		for (unsigned t = 0; t < pftracks.size(); ++t) {
			const reco::Track& track = *(pftracks[t]);
			if (track.pt() < 0.5) continue;
			TVector3 tvec(0, 0, 0);
			tvec.SetPtEtaPhi(track.pt(), track.eta(), track.phi());
			double dR = pfjet->p4.Vect().DrEtaPhi(tvec);
			if (dR < 0.4 && dR > 0.2) Itrk += track.pt(); // scalar sum for isolation
		}
		pfjet->Itrk = Itrk;
		pfjet->ecalConstitEnergy = ecalConstitEnergy;
		pfjet->hcalConstitEnergy = hcalConstitEnergy;
		pfjetVector->push_back(pfjet);

		if (verbose) std::cout << "pt: " << pfjet->p4.Vect().Pt()
				<< "\teta: " << pfjet->p4.Vect().Eta()
				<< "\tphi: " << pfjet->p4.Vect().Phi()
				<< "\temf: " << pfjet->emf
				<< "\tcr: " << pfjet->cr
				<< "\tECAL: " << pfjet->ecalConstitEnergy
				<< std::endl;
		if (verbose) std::cout << "\tHCAL: " << pfjet->hcalConstitEnergy << std::endl;
		if (verbose) std::cout << std::endl;
	}
	if (verbose) std::cout << "total pfele size: " << pfeleVector->size() << std::endl;

//########################### SuperClusters ###########################
	if (verbose) std::cout << "--- superclusters: " << SCCollectionEB->size()+SCCollectionEE->size() << std::endl;

	SCVector->clear();
	//Barrel Superclusters (Corrected)
	for( unsigned i=0; i<SCCollectionEB->size(); i++ ) {
		MySC* SC = new MySC;
		SC->energy = (*SCCollectionEB)[i].rawEnergy();
		SC->Epreshower = (*SCCollectionEB)[i].preshowerEnergy();
		SC->point.SetXYZ((*SCCollectionEB)[i].position().x(),(*SCCollectionEB)[i].position().y(),(*SCCollectionEB)[i].position().z());
		SC->phiWidth = (*SCCollectionEB)[i].phiWidth();
		SC->etaWidth = (*SCCollectionEB)[i].etaWidth();

		SCVector->push_back(SC);
	}

	//Endcap Superclusters (Corrected 5x5)
	for( unsigned i=0; i<SCCollectionEE->size(); i++ ) {
		MySC* SC = new MySC;
		SC->energy = (*SCCollectionEE)[i].rawEnergy();
		SC->Epreshower = (*SCCollectionEE)[i].preshowerEnergy();
		SC->point.SetXYZ((*SCCollectionEE)[i].position().x(),(*SCCollectionEE)[i].position().y(),(*SCCollectionEE)[i].position().z());
		SC->phiWidth = (*SCCollectionEE)[i].phiWidth();
		SC->etaWidth = (*SCCollectionEE)[i].etaWidth();

		SCVector->push_back(SC);
	}  

//########################### Muons ###########################
	if (verbose) std::cout << "--- muons: " << recoMuons->size() << std::endl;
	muonVector->clear();
	for (unsigned i = 0; i < recoMuons->size(); i++) {
		MyMuon * muon = new MyMuon;
		muon->p4 = TLorentzVector((*recoMuons)[i].px(), (*recoMuons)[i].py(),(*recoMuons)[i].pz(), (*recoMuons)[i].energy());
		muon->isGlobal = (*recoMuons)[i].isGlobalMuon();
		muon->isTracker = (*recoMuons)[i].isTrackerMuon();
		muon->isSA = (*recoMuons)[i].isStandAloneMuon();
		muon->Itrk = (*recoMuons)[i].isolationR03().sumPt;	//sum of pt tracks
		muon->trkRelChi2 = (*recoMuons)[i].combinedQuality().trkRelChi2;
		muon->staRelChi2 = (*recoMuons)[i].combinedQuality().staRelChi2;

		if (verbose) std::cout << "\tlooking for muon ID variables ..." << std::endl;	// can probably use automated funcs here if you trust them
		if ((*recoMuons)[i].globalTrack().isNonnull()){
			muon->globalTrkNormChi2 = (*recoMuons)[i].globalTrack()->normalizedChi2();
			muon->globalTrkMuHits = (*recoMuons)[i].globalTrack()->hitPattern().numberOfValidMuonHits();
		}
		else{
			muon->globalTrkNormChi2 = -1;
			muon->globalTrkMuHits = -1;
		}
		muon->numMatchedStations = (*recoMuons)[i].numberOfMatchedStations();
		if ((*recoMuons)[i].innerTrack().isNonnull()){
			if (goodVertices.at(0)){
				muon->trkIPxy = fabs((*recoMuons)[i].innerTrack()->dxy(goodVertices.at(0)->position()));
				muon->trkIPz = fabs((*recoMuons)[i].innerTrack()->dz(goodVertices.at(0)->position()));
			}
			else{
				muon->trkIPxy = fabs((*recoMuons)[i].innerTrack()->dxy());
				muon->trkIPz = fabs((*recoMuons)[i].innerTrack()->dz());
			}
			muon->numPixelHits = (*recoMuons)[i].innerTrack()->hitPattern().numberOfValidPixelHits();
			muon->numTrkHits = (*recoMuons)[i].innerTrack()->hitPattern().numberOfValidTrackerHits(); //Older selection = same as VHbb; still acceptable
			muon->innerTrkNormChi2 = (*recoMuons)[i].innerTrack()->normalizedChi2();
			muon->numPixelLayers = (*recoMuons)[i].innerTrack()->hitPattern().pixelLayersWithMeasurement();
		}
		else{
			muon->trkIPxy = -1.0;
			muon->trkIPz = -1.0;
			muon->numPixelHits = -1;
			muon->numTrkHits = -1;
			muon->innerTrkNormChi2 = -1;
			muon->numPixelLayers = -1;
		}

		if (verbose) std::cout << "\tlooking for muon isolation variables ..." << std::endl;
		muon->combIso = ((*recoMuons)[i].isolationR03().sumPt + (*recoMuons)[i].isolationR03().emEt + (*recoMuons)[i].isolationR03().hadEt)/(*recoMuons)[i].pt();
		muon->trkIso = (*recoMuons)[i].isolationR03().sumPt/(*recoMuons)[i].pt();

		if (verbose) std::cout << "\tsaving muon" << std::endl;
		muonVector->push_back(muon);
	}

//########################### Electrons ###########################
	if (verbose) std::cout << "--- electrons: " << gsfElectrons->size() << std::endl;
	recoelectronVector->clear();
	for (unsigned i = 0; i < gsfElectrons->size(); i++) {
		if (verbose) std::cout << "pt: " << (*gsfElectrons)[i].pt() << "\teta: " << (*gsfElectrons)[i].eta();
//		if (verbose) std::cout << "\tfflags: " << std::hex << (*gsfElectrons)[i].fiducialFlags() << std::dec;
		if (verbose) std::cout << std::endl;
		MyRecoElectron* ele = new MyRecoElectron;
		ele->p4 = TLorentzVector(0, 0, 0, 0);
		ele->p4.SetPtEtaPhiE((*gsfElectrons)[i].pt(), (*gsfElectrons)[i].eta(),	(*gsfElectrons)[i].phi(), (*gsfElectrons)[i].energy());

		if (verbose) std::cout << "\tfinding electron ID variables\n";	// can probably use automated funcs here if you trust them
		ele->sigieie = (*gsfElectrons)[i].sigmaIetaIeta();
		ele->deta = (*gsfElectrons)[i].deltaEtaSuperClusterTrackAtVtx();
		ele->dphi = (*gsfElectrons)[i].deltaPhiSuperClusterTrackAtVtx();
		ele->HoverE = (*gsfElectrons)[i].hadronicOverEm();
		ele->eSCoverP = (*gsfElectrons)[i].eSuperClusterOverP();
		ele->Itrk = (*gsfElectrons)[i].dr03TkSumPt();
		ele->Iecal = (*gsfElectrons)[i].dr03EcalRecHitSumEt();	// switched from dr04 to dr03
		ele->Ihcal = (*gsfElectrons)[i].dr03HcalTowerSumEt();
//		ele->convDist = (*gsfElectrons)[i].convDist();
//		ele->convDcot = (*gsfElectrons)[i].convDcot();
		ele->misHits = (*gsfElectrons)[i].gsfTrack()->trackerExpectedHitsInner().numberOfHits();
		ele->fBrem = (*gsfElectrons)[i].fbrem();
		ele->E = (*gsfElectrons)[i].ecalEnergy();
		ele->p_in = (*gsfElectrons)[i].ecalEnergy()/(*gsfElectrons)[i].eSuperClusterOverP();
//		ele->p_in = (*gsfElectrons)[i].trackMomentumAtVtx().p();
		if (goodVertices.at(0)){
			ele->trkDxy = (*gsfElectrons)[i].gsfTrack()->dxy(goodVertices.at(0)->position());
			ele->trkDz = (*gsfElectrons)[i].gsfTrack()->dz(goodVertices.at(0)->position());
		}
		else{
			ele->trkDxy = (*gsfElectrons)[i].gsfTrack()->dxy();
			ele->trkDz = (*gsfElectrons)[i].gsfTrack()->dz();
//			ele->trkDxy = -1.0;
//			ele->trkDz = -1.0;
		}
		ele->isEB = (*gsfElectrons)[i].isEB();
		ele->isEE = (*gsfElectrons)[i].isEE();
		((*gsfElectrons)[i].ecalDrivenSeed() != 0 ) ? ele->ecalDriven = 1 : ele->ecalDriven = 0;
		((*gsfElectrons)[i].trackerDrivenSeed() != 0 ) ? ele->trkDriven = 1 : ele->trkDriven = 0;
		ele->SCenergy = (*gsfElectrons)[i].superCluster()->energy();
		TVector3 scvec((*gsfElectrons)[i].superCluster()->position().x(),(*gsfElectrons)[i].superCluster()->position().y(),(*gsfElectrons)[i].superCluster()->position().z());
		ele->SCeta = scvec.Eta();
		ele->SCphi = scvec.Phi();

		//fill conversion partner track info
		ConversionFinder convFinder;
		//const double bfield = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
		const double bfield = 3.8;
		ConversionInfo convInfo = convFinder.getConversionInfo((*gsfElectrons)[i], hGeneralTracks, bfield);
		ele->convDcot = convInfo.dcot();
		ele->convDist = convInfo.dist();
		ele->convR = convInfo.radiusOfConversion();

		/*
		 reco::TrackRef convTrackRef = convInfo.conversionPartnerTk();
		 if (trackerTrackMap_ && convTrackRef.isNonnull()) {
		 outElectron->SetConvPartnerTrk(trackerTrackMap_->GetMit(refToPtr(convTrackRef)));
		 }
		 //fill additional conversion flag
		 outElectron->SetMatchesVertexConversion(convMatcher.matchesGoodConversion(*iM,hConversions));

		 */
		if (verbose) std::cout << "\tSaving electron\n";
		recoelectronVector->push_back(ele);
	}

	if (verbose) std::cout << "Filling the tree ...\n";
	tree->Fill();
	if (verbose) std::cout << "================================================"	<< std::endl;
	if (verbose) std::cout << "Freeing memory ...\n";

	if (verbose) std::cout << "Freeing memory puVector...\n";
	for (unsigned x = 0; x < puVector->size(); x++) {
		delete puVector->at(x);
		puVector->at(x) = NULL;
	}
	if (verbose) std::cout << "Freeing memory genmuonVector...\n";
	for (unsigned x = 0; x < genmuonVector->size(); x++) {
		delete genmuonVector->at(x);
		genmuonVector->at(x) = NULL;
	}
	if (verbose) std::cout << "Freeing memory genelectronVector...\n";
	for (unsigned x = 0; x < genelectronVector->size(); x++) {
		delete genelectronVector->at(x);
		genelectronVector->at(x) = NULL;
	}
	if (verbose) std::cout << "Freeing memory genSigMetVector..." << std::endl;
	for (unsigned x = 0; x < genSigMetVector->size(); x++) {
		delete genSigMetVector->at(x);
		genSigMetVector->at(x) = NULL;
	}
	if (verbose) std::cout << "Freeing memory genpVector...\n";
	for (unsigned x = 0; x < genpVector->size(); x++) {
		delete genpVector->at(x);
		genpVector->at(x) = NULL;
	}
	if (verbose) std::cout << "Freeing memory genjetVector...\n";
	for (unsigned x = 0; x < genjetVector->size(); x++) {
		delete genjetVector->at(x);
		genjetVector->at(x) = NULL;
	}
	if (verbose) std::cout << "Freeing memory photonVector...\n";
	for (unsigned x = 0; x < photonVector->size(); x++) {
		delete photonVector->at(x);
		photonVector->at(x) = NULL;
	}
	if (verbose) std::cout << "Freeing memory ak5calojetVector...\n";
	for (unsigned x = 0; x < ak5calojetVector->size(); x++) {
		delete ak5calojetVector->at(x);
		ak5calojetVector->at(x) = NULL;
	}
	if (verbose) std::cout << "Freeing memory ak7calojetVector...\n";
	for (unsigned x = 0; x < ak7calojetVector->size(); x++) {
		delete ak7calojetVector->at(x);
		ak7calojetVector->at(x) = NULL;
	}
	if (verbose) std::cout << "Freeing memory pfjetVector...\n";
	for (unsigned x = 0; x < pfjetVector->size(); x++) {
		delete pfjetVector->at(x);
		pfjetVector->at(x) = NULL;
	}
	if (verbose) std::cout << "Freeing memory SCVector...\n";
	for (unsigned x = 0; x < SCVector->size(); x++) {
		delete SCVector->at(x);
		SCVector->at(x) = NULL;
	}
	if (verbose) std::cout << "Freeing memory muonVector...\n";
	for (unsigned x = 0; x < muonVector->size(); x++) {
		delete muonVector->at(x);
		muonVector->at(x) = NULL;
	}
	if (verbose) std::cout << "Freeing memory recoelectronVector...\n";
	for (unsigned x = 0; x < recoelectronVector->size(); x++) {
		delete recoelectronVector->at(x);
		recoelectronVector->at(x) = NULL;
	}

	if (verbose) std::cout << "================================================"	<< std::endl;
	if (verbose) std::cout << std::endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void Ntuplize2::endJob() {
	//  tree->Write();
	outf->cd();
	nt->Write();
	hevt->Write();
	tree->Print();
	outf->Write();
	outf->Close();

//	delete genpVector; // crashes if try to delete -> so just leave them
//	delete genmuonVector;
//	delete genelectronVector;
//	delete recoelectronVector;
//	delete pfjetVector;
//	delete pfeleVector;
//	delete ak5calojetVector;
//	delete ak7calojetVector;
//	delete genjetVector;
//	delete muonVector;
//	delete SCVector;
//	delete photonVector;
//	delete hltobjVector;
//	delete puVector;
//	delete genSigMetVector;
//	delete outf;
//	delete hevt;
//	delete tree;
//	delete nt;

}

DEFINE_FWK_MODULE(Ntuplize2);
