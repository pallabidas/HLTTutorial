// -*- C++ -*-
//
// Package:    HLTAnalysis/TriggerAnalyzerRAWMiniAOD
// Class:      TriggerAnalyzerRAWMiniAOD
// 
/**\class TriggerAnalyzerRAWMiniAOD TriggerAnalyzerRAWMiniAOD.cc HLTAnalysis/TriggerAnalyzerRAWMiniAOD/plugins/TriggerAnalyzerRAWMiniAOD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Laurent Thomas
//         Created:  Fri, 24 Mar 2017 04:09:55 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/AssociationMap.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TLorentzVector.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class TriggerAnalyzerRAWMiniAOD : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TriggerAnalyzerRAWMiniAOD(const edm::ParameterSet&);
      ~TriggerAnalyzerRAWMiniAOD();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------


  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigobjectsMINIAODToken_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsORIGToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> trigobjectsRAWToken_;
  edm::EDGetTokenT<edm::TriggerResults>  trgresultsHLT2Token_;

  edm::EDGetTokenT<std::vector<pat::Jet> > jet_token;
  edm::EDGetTokenT<std::vector<reco::Vertex> > PV_token;

  edm::Service<TFileService> fs;
  
  TH1F* h_PFJet60_den;
  TH1F* h_PFJet60_num;

  TH1F* h_PFJet110_den;
  TH1F* h_PFJet110_num;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TriggerAnalyzerRAWMiniAOD::TriggerAnalyzerRAWMiniAOD(const edm::ParameterSet& iConfig)

{
  trigobjectsMINIAODToken_ = consumes<pat::TriggerObjectStandAloneCollection>( edm::InputTag("slimmedPatTrigger"));
  trigobjectsRAWToken_=consumes<trigger::TriggerEvent>(edm::InputTag("hltTriggerSummaryAOD::HLT2"));  

  trgresultsORIGToken_= consumes<edm::TriggerResults>( edm::InputTag("TriggerResults::HLT") );
  trgresultsHLT2Token_= consumes<edm::TriggerResults>( edm::InputTag("TriggerResults::HLT2") );
 
  jet_token = consumes< std::vector<pat::Jet> >(edm::InputTag("slimmedJets") );
  PV_token = consumes<std::vector<reco::Vertex> > (edm::InputTag("offlineSlimmedPrimaryVertices"));
  

  //now do what ever initialization is needed
  //   usesResource("TFileService");

  h_PFJet60_den= fs->make<TH1F>("h_PFJet60_den","",50,0,200);
  h_PFJet60_num= fs->make<TH1F>("h_PFJet60_num","",50,0,200);

  h_PFJet110_den= fs->make<TH1F>("h_PFJet110_den","",50,0,200);
  h_PFJet110_num= fs->make<TH1F>("h_PFJet110_num","",50,0,200);

}


TriggerAnalyzerRAWMiniAOD::~TriggerAnalyzerRAWMiniAOD()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TriggerAnalyzerRAWMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;
   using namespace reco;
   using namespace std;



   // ****************Part 1. Accessing some trigger information ************* 
   bool passHLT_PFJet60(false);
   bool passHLT_PFJet110(false);

   //Accessing trigger bits:
   //This works in both RAW, AOD or MINIAOD 
   //Here we access the decision provided by the HLT (i.e. original trigger step). 
   edm::Handle<edm::TriggerResults> trigResults;
   iEvent.getByToken(trgresultsORIGToken_, trigResults);
   if( !trigResults.failedToGet() ) {
     int N_Triggers = trigResults->size();
     const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);

     for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
       if (trigResults.product()->accept(i_Trig)) {
	 TString TrigPath =trigName.triggerName(i_Trig);
	 	 //cout << "Passed path: " << TrigPath<<endl;
	 if(TrigPath.Index("HLT_PFJet60_v") >=0) passHLT_PFJet60=true; 
	 //Notice the special syntax: since the path version can change during data taking one only looks for the string "HLT_IsoMu24_v"
       }
     }
   }

   edm::Handle<edm::TriggerResults> trigResults2;
   iEvent.getByToken(trgresultsHLT2Token_, trigResults2);
   if( !trigResults2.failedToGet() ) {
     int N_Triggers = trigResults2->size();
     const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults2);

     for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
       if (trigResults2.product()->accept(i_Trig)) {
         TString TrigPath =trigName.triggerName(i_Trig);
         if(TrigPath.Index("HLT_PFJet110_v") >=0) passHLT_PFJet110=true;
       }
     }
   }

   //Accessing the trigger objects in MINIAOD
   //This recipe works for MINIAOD only
   //edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   //iEvent.getByToken(trigobjectsMINIAODToken_, triggerObjects);

   //const edm::TriggerNames &names = iEvent.triggerNames(*trigResults);
   //for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
   //  obj.unpackFilterLabels(iEvent,*trigResults);
   //  obj.unpackPathNames(names);
   //  for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
   //    string myfillabl=obj.filterLabels()[h];
   //    cout << "Trigger object name, pt, eta, phi: "
   //    	    << myfillabl<<", " << obj.pt()<<", "<<obj.eta()<<", "<<obj.phi() << endl;
   //  }
   //}

   
   // **************** Part 2. Accessing some offline information ************** 
   
   //What you really want to do is to assess the trigger performances on top of an offline selection. 
    
   
   //Offline jets
   edm::Handle< std::vector<pat::Jet> > jets;
   iEvent.getByToken(jet_token,jets );

   double leadingjetpt(-100);
   //, leadingjeteta(-100), leadingjetphi(-100); 

   if((*jets).size() > 0){
     leadingjetpt = (*jets)[0].pt();
     //leadingjeteta = (*jets)[0].eta();
     //leadingjetphi = (*jets)[0].phi();
   }
   //for( std::vector<pat::Jet>::const_iterator jet = (*jets).begin(); jet != (*jets).end(); jet++ ) {
   //  leadingjetpt = jet->pt();
   //   = jet->eta();
   //  double phijet = jet->phi();     
   //}

   
   //We can now fill some histograms (numerator and denominator to study the efficiency of our favourite path).
   //Here we factorize the muon and jet legs and measure their efficiencies separately
   h_PFJet60_den->Fill(leadingjetpt); 
   if(passHLT_PFJet60) h_PFJet60_num->Fill(leadingjetpt);

   h_PFJet110_den->Fill(leadingjetpt);
   if(passHLT_PFJet110) h_PFJet110_num->Fill(leadingjetpt);

}


// ------------ method called once each job just before starting event loop  ------------
void 
TriggerAnalyzerRAWMiniAOD::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerAnalyzerRAWMiniAOD::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerAnalyzerRAWMiniAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerAnalyzerRAWMiniAOD);
