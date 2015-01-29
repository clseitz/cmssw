// -*- C++ -*-
//
// Package:    RecoJets/JetProducers
// Class:      addSubjetProducer
// 
/**\class addSubjetProducer addSubjetProducer.cc RecoJets/JetProducers/plugins/addSubjetProducer.cc

 Description: Simple subjet producer. Stores subjets with at least a 3% of the total pt of the jet.

*/
//
// Original Author:  Alejandro Gomez Espinosa
//         Created:  Sun, 07 Dec 2014 23:45:39 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "RecoJets/JetProducers/plugins/CompoundJetProducer.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
//
// class declaration
//

class addSubjetProducer : public edm::EDProducer {
   public:
      explicit addSubjetProducer(const edm::ParameterSet&);
      ~addSubjetProducer();

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      
      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::PFJetCollection> jets_;
      edm::EDGetTokenT<pat::JetCollection> patjets_;
      double subJetRadius_, ptFracSubjet_;
};

addSubjetProducer::addSubjetProducer(const edm::ParameterSet& iConfig):
	jets_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
	patjets_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("patjets"))),
	subJetRadius_(iConfig.getUntrackedParameter<double>("subJetRadius", 0.2)),
	ptFracSubjet_(iConfig.getUntrackedParameter<double>("ptFracSubjet", 0.03))

{
	produces< reco::BasicJetCollection >("Subjets"); 
	produces< pat::JetCollection >(); 
}


addSubjetProducer::~addSubjetProducer()
{

}


void
addSubjetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   	using namespace edm;

	edm::Handle<reco::PFJetCollection> jets;
	iEvent.getByToken(jets_, jets);

	edm::Handle<pat::JetCollection> patjets;
	iEvent.getByToken(patjets_, patjets);

	std::vector<fastjet::PseudoJet> subJets;      
	std::vector<fastjet::PseudoJet> constituents;      
	fastjet::JetDefinition aktp2(fastjet::antikt_algorithm, subJetRadius_ );

	std::auto_ptr< reco::BasicJetCollection > subjetCollection ( new reco::BasicJetCollection () );
	std::auto_ptr< pat::JetCollection > outputs( new pat::JetCollection () );
  	
	// This will store the handle for the subjets after we write them
	edm::OrphanHandle< std::vector<reco::BasicJet> > subjetHandleAfterPut;
	// this is the mapping of subjet to hard jet
	std::vector< std::vector<int> > indices;
	indices.resize( jets->size() );
	// this is the jet areas
	std::vector<double> areaJets;
	//reco::PFJetCollection subjetsPerJet; 

	int jetIndex = -1;
	for( const reco::PFJet &iJet : *jets ){

		jetIndex++;
		areaJets.push_back( iJet.jetArea() );
		
		// get jet constituents
		constituents.clear();      
		std::vector<const reco::Candidate * > jetConst = iJet.getJetConstituentsQuick();

		// loop over jet constituents for iJet
		// for each jet constituent, add to vector< pseudojet > for reclustering
		for ( unsigned int  iJetConst = 0 ; iJetConst < jetConst.size() ; iJetConst++ ){

			// get pseudojets to cluster
			constituents.push_back( fastjet::PseudoJet( jetConst[ iJetConst ]->px(), jetConst[ iJetConst ]->py(), jetConst[ iJetConst ]->pz(), jetConst[ iJetConst ]->energy() ) );
			//LogWarning("Constituents") << "jet const : " << jetConst[ iJetConst ]->px() << " " << jetConst[ iJetConst ]->py() << " " << jetConst[ iJetConst ]->pz() << " " << jetConst[ iJetConst ]->energy();

		}

		// recluster 
		fastjet::ClusterSequence cs_aktp2(constituents, aktp2);
		subJets = sorted_by_pt(cs_aktp2.inclusive_jets());

		std::vector<CompoundPseudoSubJet>  subjetsOutput;
		for( fastjet::PseudoJet & iSubJet : subJets ){

			if( iSubJet.pt() > ptFracSubjet_ * iJet.pt() ) {

				reco::Candidate::LorentzVector p4SubJet( iSubJet.px(), iSubJet.py(), iSubJet.pz(), iSubJet.e() );
				reco::Particle::Point point(0,0,0);

				std::vector<reco::CandidatePtr> subjetConstituents;
				std::vector<fastjet::PseudoJet > subjetConst = iSubJet.constituents();
				std::vector<int> constituents;
				std::vector<fastjet::PseudoJet>::const_iterator subConstIt = subjetConst.begin();
				for ( ; subConstIt != subjetConst.end(); ++subConstIt ) {
					if (subConstIt->user_index() >= 0) constituents.push_back( subConstIt->user_index() );
				}

				double subJetArea = 0.0;
				if ( iSubJet.has_area() ) subJetArea = iSubJet.area();

				// Make a CompoundPseudoSubJet object to hold this subjet and the indices of its constituents
				subjetsOutput.push_back( CompoundPseudoSubJet( iSubJet, subJetArea, constituents ) );

				// This holds the subjet-to-hardjet mapping
				indices[ jetIndex ].push_back( subjetCollection->size() );      

				// Add the concrete subjet type to the subjet list to write to event record
				reco::BasicJet jet;
				reco::writeSpecific( jet, p4SubJet, point, subjetConstituents, iSetup);
				jet.setJetArea( subJetArea );
				subjetCollection->push_back( jet );
			}
		}
	}
	// put subjets into event record
	subjetHandleAfterPut = iEvent.put( subjetCollection, "Subjets" );

	int dummy = -1;
	for( const pat::Jet & ip4 : *patjets ) {

		outputs->push_back( ip4 );
		dummy++;
		std::vector<int> & ind = indices[ dummy ];
		std::vector<reco::CandidatePtr> iJetConstituents;
	
		for( std::vector<int>::const_iterator isub = ind.begin(); isub != ind.end(); ++isub ) {

			reco::CandidatePtr candPtr( subjetHandleAfterPut, *isub, false );
			iJetConstituents.push_back( candPtr );
		}   
		outputs->back().clearDaughters();
		for( const reco::CandidatePtr & subj : iJetConstituents ){
			outputs->back().addDaughter( subj );
		}
	}

	// put jets into event record
	iEvent.put( outputs );
 
}

//define this as a plug-in
DEFINE_FWK_MODULE(addSubjetProducer);
