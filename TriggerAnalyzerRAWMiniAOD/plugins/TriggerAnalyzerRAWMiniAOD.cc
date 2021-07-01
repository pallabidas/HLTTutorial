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
#include "TTree.h"

#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

bool compareByPt (l1extra::L1JetParticle i, l1extra::L1JetParticle j) { return(i.pt()>j.pt()); };
using namespace l1extra;
using namespace std;

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
      bool RecoHLTMatching(const edm::Event&,double recoeta, double recophi, std::string filtername, double dRmatching = 0.3);

      // ----------member data ---------------------------

  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigobjectsMINIAODToken_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsORIGToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> trigobjectsRAWToken_;
  edm::EDGetTokenT<edm::TriggerResults>  trgresultsHLT2Token_;

  edm::EDGetTokenT<std::vector<pat::Jet> > jet_token;
  edm::EDGetTokenT<std::vector<pat::Jet> > jetak8_token;
  edm::EDGetTokenT<BXVector<l1t::Jet>> stage2Jet_token;
  edm::EDGetTokenT<std::vector<l1extra::L1JetParticle> > l1Boosted_token;
  edm::EDGetTokenT<std::vector<l1extra::L1JetParticle> > central_token;
  edm::EDGetTokenT<std::vector<l1extra::L1JetParticle> > forward_token;
  edm::EDGetTokenT<reco::GenParticleCollection> gen_token;

  edm::Service<TFileService> fs;
 
  TH1F* h_leadingJetPt;
  TH1F* h_L1SingleJet_den;

  TH1F* h_L1SingleJet180_num1;
  TH1F* h_L1SingleJet90_num1;
  TH1F* h_L1boosted100_num1;
  TH1F* h_L1boosted120_num1;

  TH1F* h_L1SingleJet180_num2;
  TH1F* h_L1SingleJet90_num2;
  TH1F* h_L1boosted100_num2;
  TH1F* h_L1boosted120_num2;
  TH1F* h_L1boosted100_250_num2;
  TH1F* h_L1boosted120_250_num2;

  TH1F* h_L1SingleJet180_num3;
  TH1F* h_L1SingleJet90_num3;
  TH1F* h_L1boosted100_num3;
  TH1F* h_L1boosted120_num3;
  TH1F* h_L1boosted100_250_num3;
  TH1F* h_L1boosted120_250_num3;

  TH1F* h_L1SingleJet180_num4;
  TH1F* h_L1SingleJet90_num4;
  TH1F* h_L1boosted100_num4;
  TH1F* h_L1boosted120_num4;
  TH1F* h_L1boosted100_250_num4;
  TH1F* h_L1boosted120_250_num4;

  TH1F* h_goodJetPt;
  TH1F* h_offlinePt;
  TH1F* h_l1Pt35;
  TH1F* h_l1Pt90;
  TH1F* h_l1Pt180;
  TH1F* h_l1boosted100;
  TH1F* h_l1boosted120;
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
  jetak8_token = consumes< std::vector<pat::Jet> >(edm::InputTag("slimmedJetsAK8") );
  stage2Jet_token = consumes<BXVector<l1t::Jet>>( edm::InputTag("caloStage2Digis","Jet","RECO") );
  l1Boosted_token = consumes<vector<l1extra::L1JetParticle>>( edm::InputTag("uct2016EmulatorDigis","Boosted","") );
  central_token = consumes< std::vector<l1extra::L1JetParticle> >(edm::InputTag("l1extraParticles","Central","RECO") );
  forward_token = consumes< std::vector<l1extra::L1JetParticle> >(edm::InputTag("l1extraParticles","Forward","RECO") );
  gen_token = consumes<reco::GenParticleCollection> (edm::InputTag("genParticles", "", "HLT"));

  h_leadingJetPt = fs->make<TH1F>("h_leadingJetPt","",40,0,500);
  h_L1SingleJet_den = fs->make<TH1F>("h_L1SingleJet_den","",40,0,500);

  h_L1SingleJet180_num1 = fs->make<TH1F>("h_L1SingleJet180_num1","",40,0,500);
  h_L1SingleJet90_num1 = fs->make<TH1F>("h_L1SingleJet90_num1","",40,0,500);
  h_L1boosted100_num1 = fs->make<TH1F>("h_L1boosted100_num1","",40,0,500);
  h_L1boosted120_num1 = fs->make<TH1F>("h_L1boosted120_num1","",40,0,500);

  h_L1SingleJet180_num2 = fs->make<TH1F>("h_L1SingleJet180_num2","",40,0,500);
  h_L1SingleJet90_num2 = fs->make<TH1F>("h_L1SingleJet90_num2","",40,0,500);
  h_L1boosted100_num2 = fs->make<TH1F>("h_L1boosted100_num2","",40,0,500);
  h_L1boosted120_num2 = fs->make<TH1F>("h_L1boosted120_num2","",40,0,500);
  h_L1boosted100_250_num2 = fs->make<TH1F>("h_L1boosted100_250_num2","",40,0,500);
  h_L1boosted120_250_num2 = fs->make<TH1F>("h_L1boosted120_250_num2","",40,0,500);

  h_L1SingleJet180_num3 = fs->make<TH1F>("h_L1SingleJet180_num3","",40,0,500);
  h_L1SingleJet90_num3 = fs->make<TH1F>("h_L1SingleJet90_num3","",40,0,500);
  h_L1boosted100_num3 = fs->make<TH1F>("h_L1boosted100_num3","",40,0,500);
  h_L1boosted120_num3 = fs->make<TH1F>("h_L1boosted120_num3","",40,0,500);
  h_L1boosted100_250_num3 = fs->make<TH1F>("h_L1boosted100_250_num3","",40,0,500);
  h_L1boosted120_250_num3 = fs->make<TH1F>("h_L1boosted120_250_num3","",40,0,500);

  h_L1SingleJet180_num4 = fs->make<TH1F>("h_L1SingleJet180_num4","",40,0,500);
  h_L1SingleJet90_num4 = fs->make<TH1F>("h_L1SingleJet90_num4","",40,0,500);
  h_L1boosted100_num4 = fs->make<TH1F>("h_L1boosted100_num4","",40,0,500);
  h_L1boosted120_num4 = fs->make<TH1F>("h_L1boosted120_num4","",40,0,500);
  h_L1boosted100_250_num4 = fs->make<TH1F>("h_L1boosted100_250_num4","",40,0,500);
  h_L1boosted120_250_num4 = fs->make<TH1F>("h_L1boosted120_250_num4","",40,0,500);

  h_goodJetPt = fs->make<TH1F>("h_goodJetPt","",40,0,500);
  h_offlinePt = fs->make<TH1F>("h_offlinePt","",40,0,500);
  h_l1Pt35 = fs->make<TH1F>("h_l1Pt35","",40,0,500);
  h_l1Pt90 = fs->make<TH1F>("h_l1Pt90","",40,0,500);
  h_l1Pt180 = fs->make<TH1F>("h_l1Pt180","",40,0,500);
  h_l1boosted100 = fs->make<TH1F>("h_l1boosted100","",40,0,500);
  h_l1boosted120 = fs->make<TH1F>("h_l1boosted120","",40,0,500);
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

   bool passHLT_AK8PFJet400_TrimMass30(false);
   bool passHLT_AK8PFJet400_TrimMass30_L190(false);
   bool passHLT_AK8PFJet400_TrimMass30_L1boosted(false);
   bool passHLT_AK8PFJet360_TrimMass30(false);
   bool passHLT_AK8PFJet360_TrimMass30_L190(false);
   bool passHLT_AK8PFJet360_TrimMass30_L1boosted(false);
   bool passHLT_AK8PFJet250_TrimMass30_L1boosted(false);
   bool passHLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB(false);
   bool passHLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_L190(false);
   bool passHLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_L1boosted(false);
   bool passHLT_AK8PFJet250_TrimMass30_PFAK8BoostedDoubleB_L1boosted(false);
   bool passHLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV(false);
   bool passHLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_L190(false);
   bool passHLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_L1boosted(false);
   bool passHLT_AK8PFJet250_TrimMass30_PFAK8BTagDeepCSV_L1boosted(false);
   bool passL1SingleJet35(false);
   bool passL1SingleJet90(false);
   bool passL1SingleJet180(false);
   bool passL1boosted100(false);
   bool passL1boosted120(false);

   edm::Handle<reco::GenParticleCollection> genParticles;
   if(!iEvent.getByToken(gen_token, genParticles)) cout<<"ERROR GETTING GEN PARTICLES"<<std::endl;
   iEvent.getByToken(gen_token, genParticles);
   reco::GenParticle genHiggs;
   for( reco::GenParticleCollection::const_iterator genparticle = genParticles->begin(); genparticle != genParticles->end(); genparticle++){
      if(genparticle->status() > 21 && genparticle->status() < 41 && genparticle->pdgId() == 25){
         genHiggs = *genparticle;
         //cout<<"gen particle pt: "<<genHiggs.pt()<<"\t"<<"eta: "<<genHiggs.eta()<<"\t"<<"phi: "<<genHiggs.phi()<<"\t"<<"mass: "<<genHiggs.mass()<<endl;
      }
   }

   edm::Handle<BXVector<l1t::Jet>> stage2Jets;
   if(!iEvent.getByToken(stage2Jet_token, stage2Jets)) cout<<"ERROR GETTING THE STAGE 2 JETS"<<std::endl;
   iEvent.getByToken(stage2Jet_token, stage2Jets);
   vector<l1t::Jet> seeds;
   seeds.clear();
   const BXVector<l1t::Jet> &s2j = *stage2Jets;
   for(auto obj : s2j) {
      seeds.push_back(obj);
   }

   edm::Handle<vector<l1extra::L1JetParticle>> l1Boosted;
   if(!iEvent.getByToken(l1Boosted_token, l1Boosted)) cout<<"ERROR GETTING THE L1BOOSTED JETS"<<std::endl;
   iEvent.getByToken(l1Boosted_token, l1Boosted);
   vector<l1extra::L1JetParticle> l1JetsSorted;
   l1JetsSorted.clear();
   for( vector<l1extra::L1JetParticle>::const_iterator l1Jet = l1Boosted->begin(); l1Jet != l1Boosted->end(); l1Jet++ ){
      l1JetsSorted.push_back(*l1Jet);
   }
   if(l1JetsSorted.size() > 1){  std::sort(l1JetsSorted.begin(),l1JetsSorted.end(),compareByPt);}

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
	 if(TrigPath.Index("HLT_AK8PFJet400_TrimMass30_v") >= 0) passHLT_AK8PFJet400_TrimMass30 = true;
	 if(TrigPath.Index("HLT_AK8PFJet360_TrimMass30_v") >= 0) passHLT_AK8PFJet360_TrimMass30 = true;
         if(TrigPath.Index("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_v") >= 0) passHLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB = true;
         if(TrigPath.Index("HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_v") >= 0) passHLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV = true;
	 //Notice the special syntax: since the path version can change during data taking one only looks for the string "HLT_IsoMu24_v"
       }
     }
   }

   edm::Handle<edm::TriggerResults> trigResultsHLT2;
   iEvent.getByToken(trgresultsHLT2Token_, trigResultsHLT2);
   if( !trigResultsHLT2.failedToGet() ) {
     int N_Triggers = trigResultsHLT2->size();
     const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResultsHLT2);

     for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
       if (trigResultsHLT2.product()->accept(i_Trig)) {
	 TString TrigPath =trigName.triggerName(i_Trig);
	 if(TrigPath.Index("HLT_AK8PFJet400_TrimMass30_L190_v") >= 0) passHLT_AK8PFJet400_TrimMass30_L190 = true;
         if(TrigPath.Index("HLT_AK8PFJet400_TrimMass30_L1boosted_v") >= 0) passHLT_AK8PFJet400_TrimMass30_L1boosted = true;
         if(TrigPath.Index("HLT_AK8PFJet360_TrimMass30_L190_v") >= 0) passHLT_AK8PFJet360_TrimMass30_L190 = true;
         if(TrigPath.Index("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_L190_v") >= 0) passHLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_L190 = true;
         if(TrigPath.Index("HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_L190_v") >= 0) passHLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_L190 = true;
         if(TrigPath.Index("HLT_AK8PFJet360_TrimMass30_L1boosted_v") >= 0) passHLT_AK8PFJet360_TrimMass30_L1boosted = true;
         if(TrigPath.Index("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_L1boosted_v") >= 0) passHLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_L1boosted = true;
         if(TrigPath.Index("HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_L1boosted_v") >= 0) passHLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_L1boosted = true;
         if(TrigPath.Index("HLT_AK8PFJet250_TrimMass30_L1boosted_v") >= 0) passHLT_AK8PFJet250_TrimMass30_L1boosted = true;
         if(TrigPath.Index("HLT_AK8PFJet250_TrimMass30_PFAK8BoostedDoubleB_np2_L1boosted_v") >= 0) passHLT_AK8PFJet250_TrimMass30_PFAK8BoostedDoubleB_L1boosted = true;
         if(TrigPath.Index("HLT_AK8PFJet250_TrimMass30_PFAK8BTagDeepCSV_p17_L1boosted_v") >= 0) passHLT_AK8PFJet250_TrimMass30_PFAK8BTagDeepCSV_L1boosted = true;
       }
     }
   }

   //Accessing the trigger objects in MINIAOD
   //This recipe works for MINIAOD only
   edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   iEvent.getByToken(trigobjectsMINIAODToken_, triggerObjects);

   const edm::TriggerNames &names = iEvent.triggerNames(*trigResults);
   for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
     obj.unpackFilterLabels(iEvent,*trigResults);
     obj.unpackPathNames(names);
     for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
       string myfillabl=obj.filterLabels()[h];
       //cout << "Trigger object name, pt, eta, phi: "
       //	    << myfillabl<<", " << obj.pt()<<", "<<obj.eta()<<", "<<obj.phi() << endl;
       if(!passL1SingleJet180 && myfillabl=="hltL1sSingleJet180") { h_l1Pt180->Fill(obj.pt()); passL1SingleJet180 = true; }
       if(!passL1SingleJet90 && myfillabl=="hltL1sSingleJet90") { h_l1Pt90->Fill(obj.pt()); passL1SingleJet90 = true; }
       if(!passL1SingleJet35 && myfillabl=="hltL1sSingleJet35") { h_l1Pt35->Fill(obj.pt()); passL1SingleJet35 = true; }
       if(myfillabl=="hltL1sSingleJet35"){
         for(auto jet : l1JetsSorted){
           if(abs(jet.eta()) < 2.5 && deltaR(jet.eta(), jet.phi(), obj.eta(), obj.phi()) < 0.4){
             if(jet.pt()*1.25 > 100.) { passL1boosted100 = true; h_l1boosted100->Fill(jet.pt()*1.25); } // a scale factor for L1boosted jet pT
             if(jet.pt()*1.25 > 120.) { passL1boosted120 = true; h_l1boosted120->Fill(jet.pt()*1.25); }
           }
         }
       }
     }
   }
 
   edm::Handle< std::vector<pat::Jet> > jets;
   iEvent.getByToken(jet_token,jets);

   edm::Handle< std::vector<pat::Jet> > jetsak8;
   iEvent.getByToken(jetak8_token,jetsak8);

   std::vector<pat::Jet> goodJetsAK8;
   goodJetsAK8.clear();

   double leadingjetpt = -100.; // leading AK4 jet
   double leadingjetak8pt = -100.; // leading AK8 jet
   double goodjetpt = -100.; // leading AK8 jet with exactly two subjets and at least one b hadron
   bool hasSoftDropMass40(false);

   if(iEvent.getByToken(jetak8_token, jetsak8)){
     for (const pat::Jet &jet : *jetsak8) {
       double softDropMass = jet.userFloat("ak8PFJetsPuppiSoftDropMass");
       if(!hasSoftDropMass40 && softDropMass > 40. && abs(jet.eta()) < 2.5 && reco::deltaR(jet, genHiggs) < 0.4) { hasSoftDropMass40 = true;  leadingjetak8pt = jet.pt(); } // leading AK8 jet within central region matched to gen Higgs
       if(jet.subjets("SoftDropPuppi").size() ==  2 && jet.jetFlavourInfo().getbHadrons().size() > 1){
         goodJetsAK8.push_back(jet);
       }
     }
   }

   if((*jets).size() > 0) { leadingjetpt = (*jets).at(0).pt(); } // leading AK4 jet
   if(goodJetsAK8.size() > 0) { goodjetpt = goodJetsAK8.at(0).pt(); } // leading AK8 jet having 2 subjets and at leat 2 b hadrons


//	   //Accessing the trigger objects in RAW/AOD
//	   edm::Handle<trigger::TriggerEvent> triggerObjectsSummary;
//	   iEvent.getByToken(trigobjectsRAWToken_ ,triggerObjectsSummary);
//	
//	   trigger::TriggerObjectCollection selectedObjects;
//	   if (triggerObjectsSummary.isValid()) {
//	     size_t filterIndex = (*triggerObjectsSummary).filterIndex( edm::InputTag("hltL1sSingleJet180") );
//	     trigger::TriggerObjectCollection allTriggerObjects = triggerObjectsSummary->getObjects();
//	     if (filterIndex < (*triggerObjectsSummary).sizeFilters()) { 
//	       const trigger::Keys &keys = (*triggerObjectsSummary).filterKeys(filterIndex);
//	       for (size_t j = 0; j < keys.size(); j++) {
//	         trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
//	         h_l1Pt180->Fill(foundObject.pt());
//	       }
//	     }
//	     size_t filterIndex1 = (*triggerObjectsSummary).filterIndex( edm::InputTag("hltL1sSingleJet90") );
//	     if (filterIndex1 < (*triggerObjectsSummary).sizeFilters()) {
//	       const trigger::Keys &keys = (*triggerObjectsSummary).filterKeys(filterIndex1);
//	       for (size_t j = 0; j < keys.size(); j++) {
//	         trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
//	         h_l1Pt90->Fill(foundObject.pt());
//	         for(auto jet : l1JetsSorted){
//	           if(reco::deltaR(jet, foundObject)<0.4){
//             if(jet.pt()*1.25 > 100.) passL1boosted100 = true;
//             if(jet.pt()*1.25 > 120.) passL1boosted120 = true;
//           }
//         }
//       }
//     }
//   }

   if(leadingjetpt > 0) h_leadingJetPt->Fill(leadingjetpt);
   if(goodjetpt > 0) { h_goodJetPt->Fill(goodjetpt); }

//// For finding efficiency of HLT paths by changing L1 seed

   if(leadingjetak8pt > 0) {

     h_L1SingleJet_den->Fill(leadingjetak8pt); // this is the common denominator

//// HLT_AK8PFJet400_TrimMass30

     if(passHLT_AK8PFJet400_TrimMass30) { h_L1SingleJet180_num1->Fill(leadingjetak8pt); }

     if(passHLT_AK8PFJet400_TrimMass30_L190) { h_L1SingleJet90_num1->Fill(leadingjetak8pt); }

     if(passHLT_AK8PFJet400_TrimMass30_L1boosted && passL1boosted100) { h_L1boosted100_num1->Fill(leadingjetak8pt); }

     if(passHLT_AK8PFJet400_TrimMass30_L1boosted && passL1boosted120) { h_L1boosted120_num1->Fill(leadingjetak8pt); }

//// HLT_AK8PFJet360_TrimMass30

     if(passHLT_AK8PFJet360_TrimMass30) { h_L1SingleJet180_num2->Fill(leadingjetak8pt); }

     if(passHLT_AK8PFJet360_TrimMass30_L190) { h_L1SingleJet90_num2->Fill(leadingjetak8pt); }

     if(passHLT_AK8PFJet250_TrimMass30_L1boosted && passL1boosted100) { h_L1boosted100_250_num2->Fill(leadingjetak8pt); }

     if(passHLT_AK8PFJet360_TrimMass30_L1boosted && passL1boosted100) { h_L1boosted100_num2->Fill(leadingjetak8pt); }

     if(passHLT_AK8PFJet250_TrimMass30_L1boosted && passL1boosted120) { h_L1boosted120_250_num2->Fill(leadingjetak8pt); }

     if(passHLT_AK8PFJet360_TrimMass30_L1boosted && passL1boosted120) { h_L1boosted120_num2->Fill(leadingjetak8pt); }


//// HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2

     if(passHLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB) { h_L1SingleJet180_num3->Fill(leadingjetak8pt); }

     if(passHLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_L190) { h_L1SingleJet90_num3->Fill(leadingjetak8pt); }

     if(passHLT_AK8PFJet250_TrimMass30_PFAK8BoostedDoubleB_L1boosted && passL1boosted100) { h_L1boosted100_250_num3->Fill(leadingjetak8pt); }

     if(passHLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_L1boosted && passL1boosted100) { h_L1boosted100_num3->Fill(leadingjetak8pt); }

     if(passHLT_AK8PFJet250_TrimMass30_PFAK8BoostedDoubleB_L1boosted && passL1boosted120) { h_L1boosted120_250_num3->Fill(leadingjetak8pt); }

     if(passHLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_L1boosted && passL1boosted120) { h_L1boosted120_num3->Fill(leadingjetak8pt); }


//// HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17

     if(passHLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV) { h_L1SingleJet180_num4->Fill(leadingjetak8pt); }

     if(passHLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_L190) { h_L1SingleJet90_num4->Fill(leadingjetak8pt); }

     if(passHLT_AK8PFJet250_TrimMass30_PFAK8BTagDeepCSV_L1boosted && passL1boosted100) { h_L1boosted100_250_num4->Fill(leadingjetak8pt); }

     if(passHLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_L1boosted && passL1boosted100) { h_L1boosted100_num4->Fill(leadingjetak8pt); }

     if(passHLT_AK8PFJet250_TrimMass30_PFAK8BTagDeepCSV_L1boosted && passL1boosted120) { h_L1boosted120_250_num4->Fill(leadingjetak8pt); }

     if(passHLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_L1boosted && passL1boosted120) { h_L1boosted120_num4->Fill(leadingjetak8pt); }

   }
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

bool TriggerAnalyzerRAWMiniAOD::RecoHLTMatching(const edm::Event& iEvent, double recoeta, double recophi, std::string filtername, double dRmatching){
  //In the next few lines one loops over all the trigger objects (corresponding to a given filter) and check whether one of them matches the reco object under study                                       
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsORIGToken_, trigResults);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(trigobjectsMINIAODToken_, triggerObjects);

  const edm::TriggerNames &names = iEvent.triggerNames(*trigResults);
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    obj.unpackFilterLabels(iEvent,*trigResults);
    obj.unpackPathNames(names);
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
      std::string myfillabl=obj.filterLabels()[h];
      if( myfillabl.find(filtername)!=std::string::npos   && deltaR(recoeta,recophi, obj.eta(),obj.phi())<dRmatching ) return true;
    }
  }

  return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerAnalyzerRAWMiniAOD);
