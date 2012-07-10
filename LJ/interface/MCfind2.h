void Ntuplize2::findSignalElectrons( const reco::GenParticle &p, std::vector<std::pair<reco::GenParticle, unsigned> > &vec, unsigned& label ) {
  //std::cerr << "starting daughter loop for a " << p.pdgId() << std::endl;
  for( reco::GenParticle::const_iterator it = p.begin(); it<p.end(); it++ ) {
    if( p.pdgId() == 25 ) { 
    	if (verbose) std::cerr << std::endl << "found a Higgs!, mass:" << p.mass() << std::endl << std::endl;
      if ( it->pdgId() == 25 )  continue;
      label = it - p.begin(); 
    }
    if (verbose){
    	std::cerr << "\t daughter " << it->pdgId()
	      << "\tpx: " << it->px() 
	      << "\teta: " << it->eta() 
	      << "\tphi: " << it->phi() 
	      << "\tstatus: " << it->status()
	      << std::endl;
    }
    reco::GenParticle thegen = *((reco::GenParticle*)&(*it));
    if( abs(it->pdgId()) == 11 && it->status() == 1 ) { 
      vec.push_back(std::pair<reco::GenParticle, unsigned> (thegen, label));
      if (verbose) std::cerr << "\t\t pushing back ele" << std::endl;
    }
    else {
    	if (verbose) std::cerr << "\t\t looking some more ..." << std::endl;
      findSignalElectrons( thegen, vec, label );
    }
  }
  //  std::cerr << "returned" << std::endl;
}


void Ntuplize2::findWDaughters( const reco::GenParticle &p, std::vector<std::pair<reco::GenParticle, unsigned> > &vec, unsigned &label ) {
  for( reco::GenParticle::const_iterator it = p.begin(); it<p.end(); it++ ) {
    if( abs(p.pdgId()) == 24 ) {	//W
    	if (verbose){
    		std::cerr << std::endl;
    		std::cerr << "found a W!, mass:" << p.mass() << std::endl;
    		std::cerr << std::endl;
    	}
      if ( abs(it->pdgId()) == 24 )  continue;
    }
    if (verbose) std::cerr << "\t daughter " << it->pdgId() << "\tpx: " << it->px() << "\teta: " << it->eta() << "\tphi: " << it->phi()<< "\tstatus: " << it->status() << std::endl;
    reco::GenParticle thegen = *((reco::GenParticle*)&(*it));
    if( ((abs(it->pdgId()) == 11 || abs(it->pdgId()) == 13) && it->status() == 1 ) || (abs(it->pdgId()) == 15 && it->status() == 2 ) || (abs(it->pdgId() <=6 ) && it->status() == 3) ) {
      vec.push_back(std::pair<reco::GenParticle, unsigned> (thegen, label));
      if (verbose) std::cerr << "\t\t pushing back gen particle" << std::endl;
    } else {
    	if (verbose) std::cerr << "\t\t looking some more ..." << std::endl;
      findWDaughters( thegen, vec, label );
    }
  }
}

void Ntuplize2::findZDaughters( const reco::GenParticle &p, std::vector<std::pair<reco::GenParticle, unsigned> > &vec, unsigned &label ) {
  for( reco::GenParticle::const_iterator it = p.begin(); it<p.end(); it++ ) {
    if( p.pdgId() == 23 ) {	//Z
    	if (verbose){
    		std::cerr << std::endl;
    		std::cerr << "found a W!, mass:" << p.mass() << std::endl;
    		std::cerr << std::endl;
    	}
      if ( it->pdgId() == 23 )  continue;
    }
    if (verbose) std::cerr << "\t daughter " << it->pdgId() << "\tpx: " << it->px() << "\teta: " << it->eta() << "\tphi: " << it->phi()<< "\tstatus: " << it->status() << std::endl;
    reco::GenParticle thegen = *((reco::GenParticle*)&(*it));
    if( ((abs(it->pdgId()) == 11 || abs(it->pdgId()) == 13) && it->status() == 1 ) || (abs(it->pdgId()) == 15 && it->status() == 2 ) || (abs(it->pdgId() <=6 ) && it->status() == 3) ) {
      vec.push_back(std::pair<reco::GenParticle, unsigned> (thegen, label));
      if (verbose) std::cerr << "\t\t pushing back gen particle" << std::endl;
    } else { 
    	if (verbose) std::cerr << "\t\t looking some more ..." << std::endl;
      findZDaughters( thegen, vec, label );
    }
  }
}



void Ntuplize2::findZElectronsMC( const reco::GenParticle &p, std::vector<std::pair<reco::GenParticle, unsigned> > &vec, unsigned& label ) {
  //std::cerr << "starting daughter loop for a " << p.pdgId() << std::endl;
  for( reco::GenParticle::const_iterator it = p.begin(); it<p.end(); it++ ) {
    if( p.pdgId() == 23 ) { 
    	if (verbose) std::cerr << std::endl << "found a Z!, mass:" << p.mass() << std::endl << std::endl;
      if ( it->pdgId() == 23 )  continue;
      label = it - p.begin(); 
    }
    if (verbose){
		std::cerr << "\t daughter " << it->pdgId()
			  << "\tpx: " << it->px()
			  << "\teta: " << it->eta()
			  << "\tphi: " << it->phi()
			  << "\tstatus: " << it->status()
			  << std::endl;
    }
    reco::GenParticle thegen = *((reco::GenParticle*)&(*it));
    if( abs(it->pdgId()) == 11 && it->status() == 1 ) { 
      vec.push_back(std::pair<reco::GenParticle, unsigned> (thegen, label));
      if (verbose)  std::cerr << "\t\t pushing back ele" << std::endl;
    } else { 
    	if (verbose) std::cerr << "\t\t looking some more ..." << std::endl;
      findZElectronsMC( thegen, vec, label ); }
  }
  //  std::cerr << "returned" << std::endl;
}


void Ntuplize2::findWElectronsMC( const reco::GenParticle &p, std::vector<std::pair<reco::GenParticle, unsigned> > &vec, unsigned& label ) {
  //std::cerr << "starting daughter loop for a " << p.pdgId() << std::endl;
  for( reco::GenParticle::const_iterator it = p.begin(); it<p.end(); it++ ) {
    if( p.pdgId() == 24 ) { 
    	if (verbose){
		  std::cerr << std::endl;
		  std::cerr << "found a W!, mass:" << p.mass() << std::endl;
		  std::cerr << std::endl;
    	}
      if ( it->pdgId() == 24 )  continue;
      label = it - p.begin(); 
    }
    if (verbose){
		std::cerr << "\t daughter " << it->pdgId()
			  << "\tpx: " << it->px()
			  << "\teta: " << it->eta()
			  << "\tphi: " << it->phi()
			  << "\tstatus: " << it->status()
			  << std::endl;
    }
    reco::GenParticle thegen = *((reco::GenParticle*)&(*it));
    if( abs(it->pdgId()) == 11 && it->status() == 1 ) { 
      vec.push_back(std::pair<reco::GenParticle, unsigned> (thegen, label));
      if (verbose) std::cerr << "\t\t pushing back ele" << std::endl;
    } else { 
    	if (verbose) std::cerr << "\t\t looking some more ..." << std::endl;
      findWElectronsMC( thegen, vec, label ); }
  }
  //  std::cerr << "returned" << std::endl;
}
