/*******************************************************************************
 *
 * Simple algorithm that reading and processing xAOD input.
 * Kinematic and other variables of the jets and objects are written to an
 * output ROOT ntuple, with the following structure:
 *   - ROOT file
 *    |- TFile object
 *      |- TTrees for jets, objects, etc.
 *        |- leaves containing the variables that are saved
 *
 * yes it's weird that the top level directory in the file is a TFile object
 * but hey it works so whatever
 *
 *******************************************************************************/

// c++ include(s):
#include <iostream>

// EL include(s):
#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>

// EDM include(s):
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODBase/IParticleHelpers.h"
#include "xAODBase/IParticleContainer.h"
#include "xAODBase/IParticle.h"
#include "AthContainers/ConstDataVector.h"
#include "AthContainers/DataVector.h"
#include "xAODCore/ShallowCopy.h"
#include "AssociationUtils/OverlapRemovalInit.h"
#include "AssociationUtils/ToolBox.h"
#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"

// package include(s):
#include "xAODAnaHelpers/HelperFunctions.h"
#include "xAODAnaHelpers/Jet.h"
#include "xAODAnaHelpers/mySimpleAlg.h"

// ROOT includes:
#include "TFile.h"
#include "TSystem.h"

// this is needed to distribute the algorithm to the workers
ClassImp(mySimpleAlg)

mySimpleAlg :: mySimpleAlg () :
    Algorithm("mySimpleAlg")
{
}

// free up memory in destructor
mySimpleAlg :: ~mySimpleAlg () {}

EL::StatusCode mySimpleAlg :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.

  ANA_MSG_INFO( "Calling setupJob" );

  job.useXAOD ();
  xAOD::Init( "mySimpleAlg" ).ignore(); // call before opening first file

  // create an output stream to save ntuples to
  EL::OutputStream myOutputNtuples("ntuples");
  job.outputAdd( myOutputNtuples );

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode mySimpleAlg :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  ANA_CHECK( xAH::Algorithm::algInitialize() );

  // create output file
  TFile *m_file = wk()->getOutputFile( "ntuples" );
  m_file->mkdir( m_outFileName.c_str() );
  m_file->cd( m_outFileName.c_str() );

  // Here we specify all the trees that will be saved to the output file:
  // event info
  m_infoTree = new TTree(m_infoTreeName.c_str(),"infoTree");
  m_infoTree->Branch("runNumber", &m_runNumber, "runNumber/I");
  m_infoTree->Branch("eventNumber", &m_evtNumber, "eventNumber/I");
  m_infoTree->Branch("nVtx", &m_nv, "nVtx/I");
  m_infoTree->Branch("nPriVtx", &m_npv, "nPriVtx/I");
  m_infoTree->Branch("nPileUpVtx", &m_npuv, "nPileUpVtx/I");
  m_infoTree->Branch("averageInteractionsPerCrossing", &m_mu, "averageInteractionsPerCrossing/I");
  m_infoTree->Branch("mcChannelNumber", &m_mcChannelNumber, "mcChannelNumber/I");
  m_infoTree->Branch("mcEventWeights", &m_mcEventWeights, "mcEventWeights/D");
  //m_infoTree->Branch("initialSumW", &m_MD_initialSumW, "initialSumW/D"); //weighted number of events

  // PFlow jets
  if( !m_pflowJetContainerName.empty() ) {
    m_pflowJetTree = new TTree(m_pflowJetTreeName.c_str(),"pflowJetTree");
    m_pflowJetTree->Branch("pt", &m_pt);
    m_pflowJetTree->Branch("eta", &m_eta);
    m_pflowJetTree->Branch("phi", &m_phi);
    m_pflowJetTree->Branch("e", &m_e);
    // m_pflowJetTree->Branch("isTag", &m_isTag);
    // m_pflowJetTree->Branch("nTrk", &m_nTrk);
    // m_pflowJetTree->Branch("passJVT", &m_passJVT);
    // m_pflowJetTree->Branch("pt_constit", &m_pt_constit);
    // m_pflowJetTree->Branch("eta_constit", &m_eta_constit);
    // m_pflowJetTree->Branch("phi_constit", &m_phi_constit);
    // m_pflowJetTree->Branch("e_constit", &m_e_constit);
    m_pflowJetTree->Branch("isIsoJetDR0p6", &m_isIsoJetDR0p6);
    // m_pflowJetTree->Branch("detEta", &m_detEta);
    m_pflowJetTree->Branch("neutralPFOindex", &m_neutralPFOindex);
    m_pflowJetTree->Branch("chargedPFOindex", &m_chargedPFOindex);
    m_pflowJetTree->Branch("neutralPFOenergy", &m_neutralPFOenergy);
    m_pflowJetTree->Branch("chargedPFOenergy", &m_chargedPFOenergy);
    m_pflowJetTree->Branch("calpfo_NLeadingTruthParticlePdgId", &m_calpfo_NLeadingTruthParticlePdgId);
    m_pflowJetTree->Branch("calpfo_NLeadingTruthParticleBarcode", &m_calpfo_NLeadingTruthParticleBarcode);
    m_pflowJetTree->Branch("calpfo_NLeadingTruthParticleEnergy", &m_calpfo_NLeadingTruthParticleEnergy);
  }

  // EMTopo jets
  if( !m_topoJetContainerName.empty() ) {
    m_topoJetTree = new TTree(m_topoJetTreeName.c_str(),"topoJetTree");
    m_topoJetTree->Branch("pt", &m_pt);
    m_topoJetTree->Branch("eta", &m_eta);
    m_topoJetTree->Branch("phi", &m_phi);
    m_topoJetTree->Branch("e", &m_e);
    // m_topoJetTree->Branch("isTag", &m_isTag);
    // m_topoJetTree->Branch("nTrk", &m_nTrk);
    // m_topoJetTree->Branch("passJVT", &m_passJVT);
    // m_topoJetTree->Branch("detEta", &m_detEta);
  }

  // electrons
  if( !m_electronContainerName.empty() ) {
    m_electronTree = new TTree(m_electronTreeName.c_str(),"electronTree");
    m_electronTree->Branch("pt", &m_pt);
    m_electronTree->Branch("eta", &m_eta);
    m_electronTree->Branch("phi", &m_phi);
    m_electronTree->Branch("e", &m_e);
    // m_electronTree->Branch("passSel", &m_passSel);
    m_electronTree->Branch("neutralPFOindex", &m_neutralPFOindex);
    m_electronTree->Branch("chargedPFOindex", &m_chargedPFOindex);
    m_electronTree->Branch("neutralPFOenergy", &m_neutralPFOenergy);
    m_electronTree->Branch("chargedPFOenergy", &m_chargedPFOenergy);
    m_electronTree->Branch("calpfo_NLeadingTruthParticlePdgId", &m_calpfo_NLeadingTruthParticlePdgId);
    m_electronTree->Branch("calpfo_NLeadingTruthParticleBarcode", &m_calpfo_NLeadingTruthParticleBarcode);
    m_electronTree->Branch("calpfo_NLeadingTruthParticleEnergy", &m_calpfo_NLeadingTruthParticleEnergy);
    if( m_calclusIsAvail ) {
      m_electronTree->Branch("calclus_NLeadingTruthParticlePdgId", &m_calclus_NLeadingTruthParticlePdgId);
      m_electronTree->Branch("calclus_NLeadingTruthParticleBarcode", &m_calclus_NLeadingTruthParticleBarcode);
      m_electronTree->Branch("calclus_NLeadingTruthParticleEnergy", &m_calclus_NLeadingTruthParticleEnergy);
    }
  }

  // photons
  if( !m_photonContainerName.empty() ) {
    m_photonTree = new TTree(m_photonTreeName.c_str(),"photonTree");
    m_photonTree->Branch("pt", &m_pt);
    m_photonTree->Branch("eta", &m_eta);
    m_photonTree->Branch("phi", &m_phi);
    m_photonTree->Branch("e", &m_e);
    // m_photonTree->Branch("passSel", &m_passSel);
    m_photonTree->Branch("neutralPFOindex", &m_neutralPFOindex);
    m_photonTree->Branch("chargedPFOindex", &m_chargedPFOindex);
    m_photonTree->Branch("neutralPFOenergy", &m_neutralPFOenergy);
    m_photonTree->Branch("chargedPFOenergy", &m_chargedPFOenergy);
  }

  // muons
  if( !m_muonContainerName.empty() ) {
    m_muonTree = new TTree(m_muonTreeName.c_str(),"muonTree");
    m_muonTree->Branch("pt", &m_pt);
    m_muonTree->Branch("eta", &m_eta);
    m_muonTree->Branch("phi", &m_phi);
    m_muonTree->Branch("e", &m_e);
    // m_muonTree->Branch("passSel", &m_passSel);
    m_muonTree->Branch("neutralPFOindex", &m_neutralPFOindex);
    m_muonTree->Branch("chargedPFOindex", &m_chargedPFOindex);
    m_muonTree->Branch("neutralPFOenergy", &m_neutralPFOenergy);
    m_muonTree->Branch("chargedPFOenergy", &m_chargedPFOenergy);
  }

  // truth objects (electrons/muons/etc. depending on topology)
  if( !m_truthContainerName.empty() ) {
    m_truthTree = new TTree(m_truthTreeName.c_str(),"truthTree");
    m_truthTree->Branch("barcode", &m_truthBarcode);
    m_truthTree->Branch("pdgId", &m_truthID);
    m_truthTree->Branch("pt", &m_pt);
    m_truthTree->Branch("eta", &m_eta);
    m_truthTree->Branch("phi", &m_phi);
    m_truthTree->Branch("e", &m_e);
  }

  // truth jets (intended to be AntiKt4TruthJets)
  if( !m_truthJetContainerName.empty() ) {
    m_truthJetTree = new TTree(m_truthJetTreeName.c_str(),"truthJetTree");
    m_truthJetTree->Branch("pt", &m_pt);
    m_truthJetTree->Branch("eta", &m_eta);
    m_truthJetTree->Branch("phi", &m_phi);
    m_truthJetTree->Branch("e", &m_e);
    m_truthJetTree->Branch("isIsoJetDR1p0", &m_isIsoJetDR1p0);
  }

  // truth WZ jets (AntiKt4Truth(Dressed)WZJets)
  if( !m_truthWZJetContainerName.empty() ) {
    m_truthWZJetTree = new TTree(m_truthWZJetTreeName.c_str(),"truthWZJetTree");
    m_truthWZJetTree->Branch("pt", &m_pt);
    m_truthWZJetTree->Branch("eta", &m_eta);
    m_truthWZJetTree->Branch("phi", &m_phi);
    m_truthWZJetTree->Branch("e", &m_e);
    m_truthWZJetTree->Branch("isIsoJetDR1p0", &m_isIsoJetDR1p0);
  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode mySimpleAlg :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode mySimpleAlg :: changeInput (bool /*firstFile*/)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode mySimpleAlg :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  ANA_MSG_INFO( "Initializing mySimpleAlg Interface... ");
  ANA_MSG_INFO( "mySimpleAlg Interface successfully initialized!" );

  // set conversion factor if we want to save pT and energy in GeV
  if( m_convertMeVToGeV ) m_conversionFactor = 0.001;
  else m_conversionFactor = 1;

  // initialise decoration keys for links between objects and FlowElements
  ANA_CHECK(m_electronNeutralFEReadDecorKey.assign("Electrons.neutralGlobalFELinks"));
  ANA_CHECK(m_electronChargedFEReadDecorKey.assign("Electrons.chargedGlobalFELinks"));
  ANA_CHECK(m_photonNeutralFEReadDecorKey.assign("Photons.neutralGlobalFELinks"));
  ANA_CHECK(m_muonNeutralFEReadDecorKey.assign("Muons.neutralGlobalFELinks"));
  ANA_CHECK(m_muonChargedFEReadDecorKey.assign("Muons.chargedGlobalFELinks"));
  ANA_CHECK(m_electronNeutralFEReadDecorKey.initialize());
  ANA_CHECK(m_electronChargedFEReadDecorKey.initialize());
  ANA_CHECK(m_photonNeutralFEReadDecorKey.initialize());
  ANA_CHECK(m_muonNeutralFEReadDecorKey.initialize());
  ANA_CHECK(m_muonChargedFEReadDecorKey.initialize());

  // initialise containers
  if( !m_pflowJetContainerName.empty() ) {
    ANA_CHECK(m_jetContKey.assign(m_pflowJetContainerName));
    ANA_CHECK(m_jetContKey.initialize());
  }

  if( !m_topoJetContainerName.empty() ) {
    ANA_CHECK(m_jetTopoContKey.assign(m_topoJetContainerName));
    ANA_CHECK(m_jetTopoContKey.initialize());
  }

  if( !m_truthJetContainerName.empty() ) {
    ANA_CHECK(m_jetTruthContKey.assign(m_truthJetContainerName));
    ANA_CHECK(m_jetTruthContKey.initialize());
  }

  if( !m_truthWZJetContainerName.empty() ) {
    ANA_CHECK(m_jetTruthWZContKey.assign(m_truthWZJetContainerName));
    ANA_CHECK(m_jetTruthWZContKey.initialize());
  }

  if( !m_truthContainerName.empty() ) {
    ANA_CHECK(m_truthContKey.assign(m_truthContainerName));
    ANA_CHECK(m_truthContKey.initialize());
  }

  if( !m_electronContainerName.empty() ) {
    ANA_CHECK(m_eleContKey.assign(m_electronContainerName));
    ANA_CHECK(m_eleContKey.initialize());
  }

  if( !m_photonContainerName.empty() ) {
    ANA_CHECK(m_phoContKey.assign(m_photonContainerName));
    ANA_CHECK(m_phoContKey.initialize());
  }

  if( !m_muonContainerName.empty() ) {
    ANA_CHECK(m_muContKey.assign(m_muonContainerName));
    ANA_CHECK(m_muContKey.initialize());
  }

  ANA_CHECK(m_truthParticlesContKey.assign("TruthParticles"));
  ANA_CHECK(m_truthParticlesContKey.initialize());

  // initialise jet constituent container(s)
  std::string inputContainerBaseGlobal = "Global";
  std::string inputContainerBaseCHSG = "CHSG";
  // std::string inputContainerBaseCSSKG = "CSSKG";
  std::string chargedContainerName = "ChargedParticleFlowObjects";
  std::string neutralContainerName = "NeutralParticleFlowObjects";
  m_inGlobalChargedFEKey = inputContainerBaseGlobal + chargedContainerName;
  m_inGlobalNeutralFEKey = inputContainerBaseGlobal + neutralContainerName;
  m_inCHSGChargedFEKey = inputContainerBaseCHSG + chargedContainerName;
  m_inCHSGNeutralFEKey = inputContainerBaseCHSG + neutralContainerName;
  // m_inCSSKGChargedFEKey = inputContainerBaseCSSKG + chargedContainerName;
  // m_inCSSKGNeutralFEKey = inputContainerBaseCSSKG + neutralContainerName;
  ANA_CHECK(m_inGlobalChargedFEKey.initialize());
  ANA_CHECK(m_inGlobalNeutralFEKey.initialize());
  ANA_CHECK(m_inCHSGChargedFEKey.initialize());
  ANA_CHECK(m_inCHSGNeutralFEKey.initialize());
  // ANA_CHECK(m_inCSSKGChargedFEKey.initialize());
  // ANA_CHECK(m_inCSSKGNeutralFEKey.initialize());

  // // more decorations
  // ANA_CHECK(m_neutralFECellsRemovedReadDecorKey.assign("GlobalNeutralParticleFlowObjects.cellsRemovedFromNeutralFE_barcodeEnergyPair"));
  // ANA_CHECK(m_neutralFECellsRemovedReadDecorKey.initialize());

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode mySimpleAlg :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  // get output file (required for the method for running on grid)
  TFile *m_file = wk()->getOutputFile( "ntuples" );

  /*********************
   * NTUPLE FILLER
   *********************/
  // // fetch primary vertex information for this event for later reference (needed to find number of jet tracks)
  // size_t vtxIdx = 0;
  // m_nv = 0;
  // m_npv = 0;
  // m_npuv = 0;
  // const xAOD::Vertex *priVtx = nullptr;
  // const xAOD::VertexContainer* vertices = nullptr;
  // if(evtStore()->retrieve(vertices, "PrimaryVertices").isSuccess()) {
  //   for(const xAOD::Vertex* vtx : *vertices) {
  //     m_nv++;

  //     if(vtx->vertexType() == xAOD::VxType::PriVtx) {
  // 	if( !priVtx ) priVtx = vtx;
  // 	m_npv++;
  //     }

  //     if(vtx->vertexType() == xAOD::VxType::PileUp) {
  // 	m_npuv++;
  //     }
  //   }
  // } else {
  //   ANA_MSG_WARNING( "Failed to retrieve primary vertex container" );
  // }
  // vtxIdx = priVtx->index();

  // EVENT INFO
  const xAOD::EventInfo *ei = nullptr;
  ANA_CHECK( evtStore()->retrieve( ei , "EventInfo" ) );
  m_runNumber = ei->runNumber();
  m_evtNumber = ei->eventNumber();
  m_mu = ei->averageInteractionsPerCrossing();
  m_mcChannelNumber = ei->mcChannelNumber();
  m_mcEventWeights = ei->mcEventWeights().at(0);
  m_infoTree->Fill();

  // PFLOW JETS
  if( !m_pflowJetContainerName.empty() ) {
    // tell the tree to go to the file
    // doesn't actually seem to be needed (at least when running locally - maybe necessary on grid? leaving the line here as a comment for now.)
    //m_pflowJetTree->SetDirectory( m_file->GetDirectory( m_outFileName.c_str() ) );

    SG::ReadHandle<xAOD::FlowElementContainer> inGlobalNeutralFEHandle = makeHandle(m_inGlobalNeutralFEKey);
    SG::ReadHandle<xAOD::FlowElementContainer> inGlobalChargedFEHandle = makeHandle(m_inGlobalChargedFEKey);
    SG::ReadHandle<xAOD::FlowElementContainer> inCHSGNeutralFEHandle = makeHandle(m_inCHSGNeutralFEKey);
    SG::ReadHandle<xAOD::FlowElementContainer> inCHSGChargedFEHandle = makeHandle(m_inCHSGChargedFEKey);
    SG::ReadHandle<xAOD::JetContainer> inContainer(m_jetContKey);
    SG::ReadHandle<xAOD::TruthParticleContainer> truthParticles(m_truthParticlesContKey);

    // placeholder if no scale factor information
    static const std::vector<float> junk(1,-999);

    // loop over jets
    for( const xAOD::Jet *jet : *inContainer ) {
      // vectors to store index and energy of constituents
      std::vector<int> chargedPFOindex;
      std::vector<int> neutralPFOindex;
      std::vector<double> chargedPFOenergy;
      std::vector<double> neutralPFOenergy;

      // vectors to store truth info from calpfo decoration
      std::vector<std::vector<int>> calpfo_NLeadingTruthParticlePdgId;
      std::vector<std::vector<int>> calpfo_NLeadingTruthParticleBarcode;
      std::vector<std::vector<double>> calpfo_NLeadingTruthParticleEnergy;

      // vectors to store truth info from calclus decoration
      std::vector<std::vector<int>> calclus_NLeadingTruthParticlePdgId;
      std::vector<std::vector<int>> calclus_NLeadingTruthParticleBarcode;
      std::vector<std::vector<double>> calclus_NLeadingTruthParticleEnergy;

      // loop over jet constituents
      ANA_MSG_VERBOSE( "Number of jet constituents: " << jet->numConstituents() );
      for( size_t consti = 0; consti < jet->numConstituents(); consti++) {
	// constituents point to the CHSG PFlow containers!!! Note that only the Global containers has calpfo decoration
      	const xAOD::FlowElement *fe = static_cast<const xAOD::FlowElement*>(jet->rawConstituent(consti));

	if( inGlobalChargedFEHandle->size() != inCHSGChargedFEHandle->size() || inGlobalNeutralFEHandle->size() != inCHSGNeutralFEHandle->size() )
	  ANA_MSG_INFO( "Ruh roh... Global and CHSG containers have different sizes..." );

	//if( fe->isCharged() ) {
	if( inCHSGChargedFEHandle->size() > fe->index() && inCHSGChargedFEHandle->at(fe->index()) == fe ) { //use this line if isCharged() doesn't work...
	  const xAOD::FlowElement *fe_global = inGlobalChargedFEHandle->at(fe->index());
	  chargedPFOindex.push_back(fe_global->index());
	  chargedPFOenergy.push_back(fe_global->e()*m_conversionFactor);
      	  ANA_MSG_VERBOSE( "Charged Flow Element index: " << fe_global->index() );
      	}
	//else if( !fe->isCharged() ) {
      	else if( inCHSGNeutralFEHandle->size() > fe->index() && inCHSGNeutralFEHandle->at(fe->index()) == fe ){ //use this line if isCharged() doesn't work...
	  const xAOD::FlowElement *fe_global = inGlobalNeutralFEHandle->at(fe->index());
	  neutralPFOindex.push_back(fe_global->index());
	  neutralPFOenergy.push_back(fe_global->e()*m_conversionFactor);
      	  ANA_MSG_VERBOSE( "Neutral Flow Element index: " << fe_global->index() );

	  // for neutral PFOs, also save calibration hit information
	  // by default, (up to) the three largest contributions to each neutral FE is saved to the calpfo vector
	  SG::AuxElement::ConstAccessor< std::vector<std::pair<unsigned int,double>> > calpfoVec("calpfo_100LeadingTruthParticleBarcodeEnergyPairs");
	  std::vector<std::pair<unsigned int,double>> barcodeEnergyPair = calpfoVec(*fe_global);
	  std::vector<Int_t> truthIDs;
	  std::vector<Int_t> truthBarcodes;
	  std::vector<Double_t> truthEnergies; //note: this is the amount of energy deposited in the nPFO by the truth particle

      	  ANA_MSG_VERBOSE( "Got the calpfo vector for this neutral constituent! Here are its details:" );

	  for( Size_t truthContrib = 0; truthContrib < barcodeEnergyPair.size(); truthContrib++ ) {
	    // find truth particle with matching barcode
	    bool foundMatchingBarcode = false;
	    for( const xAOD::TruthParticle *truthParticle : *truthParticles ) {
	      if( barcodeEnergyPair.at(truthContrib).first != truthParticle->barcode() ) continue;
	      foundMatchingBarcode = true;
	      truthIDs.push_back( truthParticle->pdgId() );
	      truthBarcodes.push_back( truthParticle->barcode() );
	      truthEnergies.push_back( barcodeEnergyPair.at(truthContrib).second * m_conversionFactor );

	      ANA_MSG_VERBOSE( "[neutral jet constituent]   PDGID: " << truthParticle->pdgId() << "   Barcode: " << truthParticle->barcode() << "   energy deposited: " << barcodeEnergyPair.at(truthContrib).second * m_conversionFactor << " GeV" );
	      break;
	    }
	    if( !foundMatchingBarcode ) ANA_MSG_VERBOSE( "Huh..? Didn't find any matching barcodes... This is unexpected." );
	  } //end loop over calpfo vector
	  calpfo_NLeadingTruthParticlePdgId.push_back(truthIDs);
	  calpfo_NLeadingTruthParticleBarcode.push_back(truthBarcodes);
	  calpfo_NLeadingTruthParticleEnergy.push_back(truthEnergies);

	  // save calibration hit information for the topo-cluster corresponding to this neutral FE (if decoration is available)
	  if( m_calclusIsAvail ) {
	    xAOD::CaloCluster *linkedCluster = (xAOD::CaloCluster*) fe_global->otherObjects().at(0);
	    // const xAOD::IParticle* FE_Iparticle=thisFE->otherObjects().at(0);
	    // const xAOD::CaloCluster* thisCaloCluster = dynamic_cast<const xAOD::CaloCluster*>(FE_Iparticle);

	    ANA_MSG_VERBOSE( "Fetching calclus vector for the topo-cluster linked to the neutral FE that is linked to this electron..." );
	    SG::AuxElement::ConstAccessor< std::vector<std::pair<unsigned int,double>> > calpfoVec("calclus_100LeadingTruthParticleBarcodeEnergyPairs");
	    barcodeEnergyPair = calpfoVec(*linkedCluster);
	    ANA_MSG_VERBOSE( "Got the calclus vector! Here are its details:" );

	    // find truth particle with matching barcode and check its PDG ID
	    truthIDs.clear();
	    truthBarcodes.clear();
	    truthEnergies.clear();
	    for( Size_t truthContrib = 0; truthContrib < barcodeEnergyPair.size(); truthContrib++ ) {
	      bool foundMatchingBarcode = false;
	      for( const xAOD::TruthParticle *truthParticle : *truthParticles ) {
		if( barcodeEnergyPair.at(truthContrib).first != truthParticle->barcode() ) continue;
		foundMatchingBarcode = true;
		truthIDs.push_back( truthParticle->pdgId() );
		truthBarcodes.push_back( truthParticle->barcode() );
		truthEnergies.push_back( barcodeEnergyPair.at(truthContrib).second * m_conversionFactor );

		ANA_MSG_VERBOSE( "[topo-cluster linked to the electron-linked neutral FE]   PDGID: " << truthParticle->pdgId() << "   Barcode: " << truthParticle->barcode() << "   energy deposited: " << barcodeEnergyPair.at(truthContrib).second * m_conversionFactor << " GeV" );
		break;
	      }
	      if( !foundMatchingBarcode ) ANA_MSG_VERBOSE( "Huh..? Didn't find any matching barcodes... This is unexpected." );
	    } //end loop over calpfo vector
	    calclus_NLeadingTruthParticlePdgId.push_back(truthIDs);
	    calclus_NLeadingTruthParticleBarcode.push_back(truthBarcodes);
	    calclus_NLeadingTruthParticleEnergy.push_back(truthEnergies);
	  }
      	}
      } //end loop over jet constituents

      m_pt.push_back( jet->pt() * m_conversionFactor );
      m_eta.push_back( jet->eta() );
      m_phi.push_back( jet->phi() );
      m_e.push_back( jet->e() * m_conversionFactor );
      m_neutralPFOindex.push_back( neutralPFOindex );
      m_chargedPFOindex.push_back( chargedPFOindex );
      m_neutralPFOenergy.push_back( neutralPFOenergy );
      m_chargedPFOenergy.push_back( chargedPFOenergy );
      m_calpfo_NLeadingTruthParticlePdgId.push_back( calpfo_NLeadingTruthParticlePdgId );
      m_calpfo_NLeadingTruthParticleBarcode.push_back( calpfo_NLeadingTruthParticleBarcode );
      m_calpfo_NLeadingTruthParticleEnergy.push_back( calpfo_NLeadingTruthParticleEnergy );
      m_calclus_NLeadingTruthParticlePdgId.push_back( calclus_NLeadingTruthParticlePdgId );
      m_calclus_NLeadingTruthParticleBarcode.push_back( calclus_NLeadingTruthParticleBarcode );
      m_calclus_NLeadingTruthParticleEnergy.push_back( calclus_NLeadingTruthParticleEnergy );

      // // b-tag and scale factor
      // SG::AuxElement::ConstAccessor< char > isTag("BTag_DL1dv01_FixedCutBEff_77"); // actually an int: 0 = not b-tagged, 1 = tagged
      // SG::AuxElement::ConstAccessor< std::vector<float> > sf("BTag_SF_DL1dv01_FixedCutBEff_77");
      // m_isTag.push_back( isTag.isAvailable(*jet) ? isTag(*jet) : -1 );
      // m_sf.push_back( sf.isAvailable(*jet) ? sf(*jet) : junk);
      // if( isTag.isAvailable(*jet) && msgLvl(MSG::VERBOSE) ) std::cout << "Jet b-tag flag: " << (int)isTag(*jet) << "\n";

      // // number of tracks associated with jet
      // static const SG::AuxElement::ConstAccessor< std::vector<int> > acc("NumTrkPt500");
      // int numTracks = acc(*jet).at(vtxIdx);
      // m_nTrk.push_back( numTracks );
      // ANA_MSG_VERBOSE( "Number of jet tracks: " << numTracks );

      // // JVT
      // static const SG::AuxElement::ConstAccessor< char > passJVT("JetJVT_Passed_Tight");
      // m_passJVT.push_back( passJVT.isAvailable(*jet) ? passJVT(*jet) : -1 );
      // ANA_MSG_VERBOSE( "Pass JVT cut: " << (int)passJVT(*jet) );

      // // constituent-scale quantities
      // xAOD::JetFourMom_t constit4vec = jet->getAttribute<xAOD::JetFourMom_t>("JetConstitScaleMomentum");
      // m_pt_constit.push_back( constit4vec.Pt() * m_conversionFactor );
      // m_eta_constit.push_back( constit4vec.Eta() );
      // m_phi_constit.push_back( constit4vec.Phi() );
      // m_e_constit.push_back( constit4vec.E() * m_conversionFactor );

      // isolation check: see if there are any other jets with pt > 7 GeV within deltaR < 0.6
      bool foundClosebyJet = false;
      for( const xAOD::Jet *jet2 : *inContainer ) {
	if( jet2 == jet || jet2->pt() < 7000 ) {
	  continue;
	}

	double jetPhiDiff = TVector2::Phi_mpi_pi(jet->phi() - jet2->phi());
	double jetEtaDiff = jet->eta() - jet2->eta();
	double jetDeltaR = sqrt( jetPhiDiff*jetPhiDiff + jetEtaDiff*jetEtaDiff );

	if( jetDeltaR < 0.6 ) {
	  foundClosebyJet = true;
	  break;
	}
      } //end inner loop over jets

      if( !foundClosebyJet ) m_isIsoJetDR0p6.push_back( 1 );
      else m_isIsoJetDR0p6.push_back( 0 );

      // // detector eta
      // m_detEta.push_back( jet->getAttribute<float>("DetectorEta") );

    } //end loop over jets

    // fill tree
    m_pflowJetTree->Fill();

    // clear vectors before going to next event
    m_pt.clear();
    m_eta.clear();
    m_phi.clear();
    m_e.clear();
    // m_isTag.clear();
    // m_sf.clear();
    // m_nTrk.clear();
    // m_passJVT.clear();
    // m_pt_constit.clear();
    // m_eta_constit.clear();
    // m_phi_constit.clear();
    // m_e_constit.clear();
    m_isIsoJetDR0p6.clear();
    // m_detEta.clear();
    m_neutralPFOindex.clear();
    m_chargedPFOindex.clear();
    m_neutralPFOenergy.clear();
    m_chargedPFOenergy.clear();
    m_calpfo_NLeadingTruthParticlePdgId.clear();
    m_calpfo_NLeadingTruthParticleBarcode.clear();
    m_calpfo_NLeadingTruthParticleEnergy.clear();
    m_calclus_NLeadingTruthParticlePdgId.clear();
    m_calclus_NLeadingTruthParticleBarcode.clear();
    m_calclus_NLeadingTruthParticleEnergy.clear();
  }

  // EMTOPO JETS
  if( !m_topoJetContainerName.empty() ) {
    SG::ReadHandle<xAOD::JetContainer> inContainer(m_jetTopoContKey);

    // placeholder if no scale factor information
    static const std::vector<float> junk(1,-999);

    // loop over jets
    for( const xAOD::Jet *jet : *inContainer ) {
      m_pt.push_back( jet->pt() * m_conversionFactor );
      m_eta.push_back( jet->eta() );
      m_phi.push_back( jet->phi() );
      m_e.push_back( jet->e() * m_conversionFactor );

    //   // b-tag and scale factor -- actually nothing for EMTopo jets available...
    //   SG::AuxElement::ConstAccessor< char > isTag("BTag_DL1dv01_FixedCutBEff_77"); // actually an int: 0 = not b-tagged, 1 = tagged
    //   SG::AuxElement::ConstAccessor< std::vector<float> > sf("BTag_SF_DL1dv01_FixedCutBEff_77");
    //   m_isTag.push_back( isTag.isAvailable(*jet) ? isTag(*jet) : -1 );
    //   m_sf.push_back( sf.isAvailable(*jet) ? sf(*jet) : junk);
    //   if( isTag.isAvailable(*jet) && msgLvl(MSG::VERBOSE) ) std::cout << "Jet b-tag flag: " << (int)isTag(*jet) << "\n";

    //   // number of tracks associated with jet
    //   static const SG::AuxElement::ConstAccessor< std::vector<int> > acc("NumTrkPt500");
    //   int numTracks = acc(*jet).at(vtxIdx);
    //   m_nTrk.push_back( numTracks );
    //   ANA_MSG_VERBOSE( "Number of jet tracks: " << numTracks );

    //   // JVT
    //   static const SG::AuxElement::ConstAccessor< char > passJVT("JetJVT_Passed_Tight");
    //   m_passJVT.push_back( passJVT.isAvailable(*jet) ? passJVT(*jet) : -1 );
    //   ANA_MSG_VERBOSE( "Pass JVT cut: " << (int)passJVT(*jet) );

    //   // detector eta
    //   m_detEta.push_back( jet->getAttribute<float>("DetectorEta") );
    }

    // fill tree
    m_topoJetTree->Fill();

    // clear vectors before going to next event
    m_pt.clear();
    m_eta.clear();
    m_phi.clear();
    m_e.clear();
    // m_isTag.clear();
    // m_sf.clear();
    // m_nTrk.clear();
    // m_passJVT.clear();
    // m_detEta.clear();
  }

  /// ELECTRONS
  if( !m_electronContainerName.empty() ) {
    SG::ReadHandle<xAOD::ElectronContainer> inContainer(m_eleContKey);
    SG::ReadHandle<xAOD::TruthParticleContainer> truthParticles(m_truthParticlesContKey);

    Int_t countElectron = 0;
    for( const xAOD::Electron *electron : *inContainer ) {
      // // identification selection check
      // SG::AuxElement::ConstAccessor< char > passSel("passSel");

      // vectors to store FE indices and energies for current electron
      std::vector<int> neutralPFOindex;
      std::vector<int> chargedPFOindex;
      std::vector<double> neutralPFOenergy;
      std::vector<double> chargedPFOenergy;

      // FE links to Global PFlow containers
      SG::ReadDecorHandle<xAOD::ElectronContainer, std::vector<ElementLink<xAOD::FlowElementContainer>>> neutralFEReadDecorHandle(m_electronNeutralFEReadDecorKey);
      SG::ReadDecorHandle<xAOD::ElectronContainer, std::vector<ElementLink<xAOD::FlowElementContainer>>> chargedFEReadDecorHandle(m_electronChargedFEReadDecorKey);
      std::vector<ElementLink<xAOD::FlowElementContainer>> electronNFELinks = neutralFEReadDecorHandle(*electron);
      std::vector<ElementLink<xAOD::FlowElementContainer>> electronCFELinks = chargedFEReadDecorHandle(*electron);

      // vectors to store calibration hit information for each neutral FE
      // by default, (up to) the three largest contributions is saved to the calpfo vector
      std::vector<std::vector<int>> calpfo_NLeadingTruthParticlePdgId;
      std::vector<std::vector<int>> calpfo_NLeadingTruthParticleBarcode;
      std::vector<std::vector<double>> calpfo_NLeadingTruthParticleEnergy;

      // also prepare vectors for CaloCalTopoCluster calibration hits
      std::vector<std::vector<int>> calclus_NLeadingTruthParticlePdgId;
      std::vector<std::vector<int>> calclus_NLeadingTruthParticleBarcode;
      std::vector<std::vector<double>> calclus_NLeadingTruthParticleEnergy;

      // loop over each neutral FE linked to the electron
      for( ElementLink<xAOD::FlowElementContainer> feLink : electronNFELinks ) {
	if( feLink.isValid() ) {
	  const xAOD::FlowElement *electronNeutralGlobalFlowElement = *feLink;
	  neutralPFOindex.push_back(electronNeutralGlobalFlowElement->index());
	  neutralPFOenergy.push_back(electronNeutralGlobalFlowElement->e() * m_conversionFactor);

	  // for debug
	  ANA_MSG_VERBOSE( "Electron number " << countElectron << " is linked to a neutral FE with index: " << electronNeutralGlobalFlowElement->index() << " with energy: " << electronNeutralGlobalFlowElement->e() * m_conversionFactor << " GeV." );

	  // save calibration hit information
	  ANA_MSG_VERBOSE( "Fetching calpfo vector for the neutral FE linked to this electron..." );
	  SG::AuxElement::ConstAccessor< std::vector<std::pair<unsigned int,double>> > calpfoVec("calpfo_100LeadingTruthParticleBarcodeEnergyPairs");
	  std::vector<std::pair<unsigned int,double>> barcodeEnergyPair = calpfoVec(*electronNeutralGlobalFlowElement);
	  ANA_MSG_VERBOSE( "Got the calpfo vector! Here are its details:" );

	  // find truth particle with matching barcode and check its PDG ID
	  std::vector<Int_t> truthIDs;
	  std::vector<Int_t> truthBarcodes;
	  std::vector<Double_t> truthEnergies;
	  for( Size_t truthContrib = 0; truthContrib < barcodeEnergyPair.size(); truthContrib++ ) {
	    bool foundMatchingBarcode = false;
	    for( const xAOD::TruthParticle *truthParticle : *truthParticles ) {
	      if( barcodeEnergyPair.at(truthContrib).first != truthParticle->barcode() ) continue;
	      foundMatchingBarcode = true;
	      truthIDs.push_back( truthParticle->pdgId() );
	      truthBarcodes.push_back( truthParticle->barcode() );
	      truthEnergies.push_back( barcodeEnergyPair.at(truthContrib).second * m_conversionFactor );

	      ANA_MSG_VERBOSE( "[electron-linked neutral FE]   PDGID: " << truthParticle->pdgId() << "   Barcode: " << truthParticle->barcode() << "   energy deposited: " << barcodeEnergyPair.at(truthContrib).second * m_conversionFactor << " GeV" );
	      break;
	    }
	    if( !foundMatchingBarcode ) ANA_MSG_VERBOSE( "Huh..? Didn't find any matching barcodes... This is unexpected." );
	  } //end loop over calpfo vector
	  calpfo_NLeadingTruthParticlePdgId.push_back(truthIDs);
	  calpfo_NLeadingTruthParticleBarcode.push_back(truthBarcodes);
	  calpfo_NLeadingTruthParticleEnergy.push_back(truthEnergies);

	  // save calibration hit information for the topo-cluster corresponding to this neutral FE (if decoration is available)
	  if( m_calclusIsAvail ) {
	    xAOD::CaloCluster *linkedCluster = (xAOD::CaloCluster*) electronNeutralGlobalFlowElement->otherObjects().at(0);

	    ANA_MSG_VERBOSE( "Fetching calclus vector for the topo-cluster linked to the neutral FE that is linked to this electron..." );
	    SG::AuxElement::ConstAccessor< std::vector<std::pair<unsigned int,double>> > calpfoVec("calclus_100LeadingTruthParticleBarcodeEnergyPairs");
	    barcodeEnergyPair = calpfoVec(*linkedCluster);
	    ANA_MSG_VERBOSE( "Got the calclus vector! Here are its details:" );

	    // find truth particle with matching barcode and check its PDG ID
	    truthIDs.clear();
	    truthBarcodes.clear();
	    truthEnergies.clear();
	    for( Size_t truthContrib = 0; truthContrib < barcodeEnergyPair.size(); truthContrib++ ) {
	      bool foundMatchingBarcode = false;
	      for( const xAOD::TruthParticle *truthParticle : *truthParticles ) {
		if( barcodeEnergyPair.at(truthContrib).first != truthParticle->barcode() ) continue;
		foundMatchingBarcode = true;
		truthIDs.push_back( truthParticle->pdgId() );
		truthBarcodes.push_back( truthParticle->barcode() );
		truthEnergies.push_back( barcodeEnergyPair.at(truthContrib).second * m_conversionFactor );

		ANA_MSG_VERBOSE( "[topo-cluster linked to the electron-linked neutral FE]   PDGID: " << truthParticle->pdgId() << "   Barcode: " << truthParticle->barcode() << "   energy deposited: " << barcodeEnergyPair.at(truthContrib).second * m_conversionFactor << " GeV" );
		break;
	      }
	      if( !foundMatchingBarcode ) ANA_MSG_VERBOSE( "Huh..? Didn't find any matching barcodes... This is unexpected." );
	    } //end loop over calpfo vector
	    calclus_NLeadingTruthParticlePdgId.push_back(truthIDs);
	    calclus_NLeadingTruthParticleBarcode.push_back(truthBarcodes);
	    calclus_NLeadingTruthParticleEnergy.push_back(truthEnergies);
	  }
	}
      } //end loop over linked neutral FEs

      // loop over each charged FE linked to the electron
      for( ElementLink<xAOD::FlowElementContainer> feLink : electronCFELinks ) {
	if( feLink.isValid() ) {
	  const xAOD::FlowElement *electronChargedGlobalFlowElement = *feLink;
	  chargedPFOindex.push_back(electronChargedGlobalFlowElement->index());
	  chargedPFOenergy.push_back(electronChargedGlobalFlowElement->e() * m_conversionFactor);

	  // for debug
	  ANA_MSG_VERBOSE( "Electron number " << countElectron << ", charged Global FE link has index: " << electronChargedGlobalFlowElement->index() << " and energy: " << electronChargedGlobalFlowElement->e() * m_conversionFactor << " GeV." );
	}
      } //end loop over linked charged FEs

      m_pt.push_back( electron->pt() * m_conversionFactor );
      m_eta.push_back( electron->eta() );
      m_phi.push_back( electron->phi() );
      m_e.push_back( electron->e() * m_conversionFactor );
      // m_passSel.push_back( passSel.isAvailable(*electron) ? passSel(*electron) : -1 );
      m_neutralPFOindex.push_back( neutralPFOindex );
      m_chargedPFOindex.push_back( chargedPFOindex );
      m_neutralPFOenergy.push_back( neutralPFOenergy );
      m_chargedPFOenergy.push_back( chargedPFOenergy );
      m_calpfo_NLeadingTruthParticlePdgId.push_back( calpfo_NLeadingTruthParticlePdgId );
      m_calpfo_NLeadingTruthParticleBarcode.push_back( calpfo_NLeadingTruthParticleBarcode );
      m_calpfo_NLeadingTruthParticleEnergy.push_back( calpfo_NLeadingTruthParticleEnergy );
      m_calclus_NLeadingTruthParticlePdgId.push_back( calclus_NLeadingTruthParticlePdgId );
      m_calclus_NLeadingTruthParticleBarcode.push_back( calclus_NLeadingTruthParticleBarcode );
      m_calclus_NLeadingTruthParticleEnergy.push_back( calclus_NLeadingTruthParticleEnergy );

      countElectron++;
    } //end loop over electrons

    // fill tree
    m_electronTree->Fill();

    // clear vectors before going to next event
    m_pt.clear();
    m_eta.clear();
    m_phi.clear();
    m_e.clear();
    // m_passSel.clear();
    m_neutralPFOindex.clear();
    m_chargedPFOindex.clear();
    m_neutralPFOenergy.clear();
    m_chargedPFOenergy.clear();
    m_calpfo_NLeadingTruthParticlePdgId.clear();
    m_calpfo_NLeadingTruthParticleBarcode.clear();
    m_calpfo_NLeadingTruthParticleEnergy.clear();
    m_calclus_NLeadingTruthParticlePdgId.clear();
    m_calclus_NLeadingTruthParticleBarcode.clear();
    m_calclus_NLeadingTruthParticleEnergy.clear();
  }

  /// PHOTONS
  if( !m_photonContainerName.empty() ) {
    SG::ReadHandle<xAOD::PhotonContainer> inContainer(m_phoContKey);

    Int_t countPhoton = 0;
    for( const xAOD::Photon *photon : *inContainer ) {
      // // identification selection check
      // SG::AuxElement::ConstAccessor< char > passSel("passSel");

      // vectors to store FE indices and energies for current photon
      std::vector<int> neutralPFOindex;
      std::vector<double> neutralPFOenergy;

      // FE links to Global PFlow containers
      SG::ReadDecorHandle<xAOD::PhotonContainer, std::vector<ElementLink<xAOD::FlowElementContainer>>> neutralFEReadDecorHandle(m_photonNeutralFEReadDecorKey);
      std::vector<ElementLink<xAOD::FlowElementContainer>> photonNFELinks = neutralFEReadDecorHandle(*photon);

      // loop over each neutral FE linked to the photon
      for( ElementLink<xAOD::FlowElementContainer> feLink : photonNFELinks ) {
	if( feLink.isValid() ) {
	  const xAOD::FlowElement *photonNFE = *feLink;
	  neutralPFOindex.push_back(photonNFE->index());
	  neutralPFOenergy.push_back(photonNFE->e() * m_conversionFactor);
	  ANA_MSG_VERBOSE( "Photon number " << countPhoton << ", link has index: " << photonNFE->index() );
	}
      }

      m_pt.push_back( photon->pt() * m_conversionFactor );
      m_eta.push_back( photon->eta() );
      m_phi.push_back( photon->phi() );
      m_e.push_back( photon->e() * m_conversionFactor );
      // m_passSel.push_back( passSel.isAvailable(*photon) ? passSel(*photon) : -1 );
      m_neutralPFOindex.push_back( neutralPFOindex );
      m_neutralPFOenergy.push_back( neutralPFOenergy );

      countPhoton++;
    } //end loop over photons

    // fill tree
    m_photonTree->Fill();

    // clear vectors before going to next event
    m_pt.clear();
    m_eta.clear();
    m_phi.clear();
    m_e.clear();
    // m_passSel.clear();
    m_neutralPFOindex.clear();
    m_neutralPFOenergy.clear();
  }

  /// MUONS
  if( !m_muonContainerName.empty() ) {
    SG::ReadHandle<xAOD::MuonContainer> inContainer(m_muContKey);

    Int_t countMuon = 0;
    for( const xAOD::Muon *muon : *inContainer ) {
      // // identification selection check
      // SG::AuxElement::ConstAccessor< char > passSel("passSel");

      // vectors to store FE indices and energies for current muon
      std::vector<int> neutralPFOindex;
      std::vector<int> chargedPFOindex;
      std::vector<double> neutralPFOenergy;
      std::vector<double> chargedPFOenergy;

      // FE links to Global PFlow containers
      SG::ReadDecorHandle<xAOD::MuonContainer, std::vector<ElementLink<xAOD::FlowElementContainer>>> neutralFEReadDecorHandle(m_muonNeutralFEReadDecorKey);
      SG::ReadDecorHandle<xAOD::MuonContainer, std::vector<ElementLink<xAOD::FlowElementContainer>>> chargedFEReadDecorHandle(m_muonChargedFEReadDecorKey);
      std::vector<ElementLink<xAOD::FlowElementContainer>> muonNFELinks = neutralFEReadDecorHandle(*muon);
      std::vector<ElementLink<xAOD::FlowElementContainer>> muonCFELinks = chargedFEReadDecorHandle(*muon);

      // loop over each neutral FE linked to the muon
      for( ElementLink<xAOD::FlowElementContainer> feLink : muonNFELinks ) {
	if( feLink.isValid() ) {
	  const xAOD::FlowElement *muonNFE = *feLink;
	  neutralPFOindex.push_back(muonNFE->index());
	  neutralPFOenergy.push_back(muonNFE->e() * m_conversionFactor);
	  ANA_MSG_VERBOSE( "Muon number " << countMuon << ", neutral FE link has index: " << muonNFE->index() );
	}
      }

      // loop over each charged FE linked to the muon
      for( ElementLink<xAOD::FlowElementContainer> feLink : muonCFELinks ) {
	if( feLink.isValid() ) {
	  const xAOD::FlowElement *muonCFE = *feLink;
	  chargedPFOindex.push_back(muonCFE->index());
	  chargedPFOenergy.push_back(muonCFE->e() * m_conversionFactor);
	  ANA_MSG_VERBOSE( "Muon number " << countMuon << ", charged FE link has index: " << muonCFE->index() );
	}
      }

      m_pt.push_back( muon->pt() * m_conversionFactor );
      m_eta.push_back( muon->eta() );
      m_phi.push_back( muon->phi() );
      m_e.push_back( muon->e() * m_conversionFactor );
      // m_passSel.push_back( passSel.isAvailable(*muon) ? passSel(*muon) : -1 );
      m_neutralPFOindex.push_back( neutralPFOindex );
      m_chargedPFOindex.push_back( chargedPFOindex );
      m_neutralPFOenergy.push_back( neutralPFOenergy );
      m_chargedPFOenergy.push_back( chargedPFOenergy );

      countMuon++;
    }

    // fill tree
    m_muonTree->Fill();

    // clear vectors before going to next event
    m_pt.clear();
    m_eta.clear();
    m_phi.clear();
    m_e.clear();
    // m_passSel.clear();
    m_neutralPFOindex.clear();
    m_chargedPFOindex.clear();
    m_neutralPFOenergy.clear();
    m_chargedPFOenergy.clear();
  }

  // // TRUTH PARTICLES (for Z+jets samples; expect 2 truth electrons/muons)
  // // side note: for ttbar topology, we already saved truth particles during the preselection step
  // if( !m_truthContainerName.empty() && (m_isZee || m_isZmumu) ) {
  //   if( m_isZee ) m_truthPDGID = 11;
  //   else m_truthPDGID = 13;

  //   SG::ReadHandle<xAOD::TruthParticleContainer> inContainer(m_truthContKey);

  //   // save pT-ordered truth particles
  //   int truthCount = 0;
  //   double leadingTruthPt = 0;
  //   double leadingTruthEta = 0;
  //   double leadingTruthPhi = 0;
  //   double leadingTruthE = 0;
  //   int leadingTruthID = 0;
  //   int leadingTruthBarcode = 0;
  //   double subleadingTruthPt = 0;
  //   double subleadingTruthEta = 0;
  //   double subleadingTruthPhi = 0;
  //   double subleadingTruthE = 0;
  //   int subleadingTruthID = 0;
  //   int subleadingTruthBarcode = 0;
  //   for( const xAOD::TruthParticle* truth : *inContainer ) {
  //     if( std::abs( truth->pdgId() ) != m_truthPDGID ) continue;

  // 	truthCount++;
  // 	if( truthCount == 1 ) {
  // 	  leadingTruthPt = truth->pt() * m_conversionFactor;
  // 	  leadingTruthEta = truth->eta();
  // 	  leadingTruthPhi = truth->phi();
  // 	  leadingTruthE = truth->e() * m_conversionFactor;
  // 	  leadingTruthID = truth->pdgId();
  // 	  leadingTruthBarcode = truth->barcode();
  // 	} else if( truthCount > 1 ) {
  // 	  if( truth->pt() * m_conversionFactor > leadingTruthPt ) {
  // 	    // set previous truth as subleading, and current truth as leading
  // 	    subleadingTruthPt = leadingTruthPt;
  // 	    subleadingTruthEta = leadingTruthEta;
  // 	    subleadingTruthPhi = leadingTruthPhi;
  // 	    subleadingTruthE = leadingTruthE;
  // 	    subleadingTruthID = leadingTruthID;
  // 	    subleadingTruthBarcode = leadingTruthBarcode;
  // 	    leadingTruthPt = truth->pt() * m_conversionFactor;
  // 	    leadingTruthEta = truth->eta();
  // 	    leadingTruthPhi = truth->phi();
  // 	    leadingTruthE = truth->e() * m_conversionFactor;
  // 	    leadingTruthID = truth->pdgId();
  // 	    leadingTruthBarcode = truth->barcode();
  // 	  } else {
  // 	    // set current truth as subleading
  // 	    subleadingTruthPt = truth->pt() * m_conversionFactor;
  // 	    subleadingTruthEta = truth->eta();
  // 	    subleadingTruthPhi = truth->phi();
  // 	    subleadingTruthE = truth->e() * m_conversionFactor;
  // 	    subleadingTruthID = truth->pdgId();
  // 	    subleadingTruthBarcode = truth->barcode();
  // 	  }

  // 	  // stop iterating once we have two truths
  // 	  break;
  // 	} // end if truthCount > 1
  //   } // end loop over truth container

  //   m_pt.push_back( leadingTruthPt );
  //   m_pt.push_back( subleadingTruthPt );
  //   m_eta.push_back( leadingTruthEta );
  //   m_eta.push_back( subleadingTruthEta );
  //   m_phi.push_back( leadingTruthPhi );
  //   m_phi.push_back( subleadingTruthPhi );
  //   m_e.push_back( leadingTruthE );
  //   m_e.push_back( subleadingTruthE );
  //   m_truthID.push_back( leadingTruthID );
  //   m_truthID.push_back( subleadingTruthID );
  //   m_truthBarcode.push_back( leadingTruthBarcode );
  //   m_truthBarcode.push_back( subleadingTruthBarcode );

  //   // fill tree
  //   m_truthTree->Fill();

  //   // clear vectors before going to next event
  //   m_pt.clear();
  //   m_eta.clear();
  //   m_phi.clear();
  //   m_e.clear();
  //   m_truthID.clear();
  //   m_truthBarcode.clear();
  // }

  // // TRUTH PARTICLES (for SinglePhoton samples; expect 1 truth photon)
  // if( !m_truthContainerName.empty() && m_isSinglePhoton ) {
  //   SG::ReadHandle<xAOD::TruthParticleContainer> inContainer(m_truthContKey);

  //   m_truthPDGID = 22;
  //   double leadingTruthPt = 0;
  //   double leadingTruthEta = 0;
  //   double leadingTruthPhi = 0;
  //   double leadingTruthE = 0;
  //   int leadingTruthID = 0;
  //   int leadingTruthBarcode = 0;
  //   for( const xAOD::TruthParticle* truth : *inContainer ) {
  //     if( std::abs( truth->pdgId() ) != m_truthPDGID ) continue;

  //     leadingTruthPt = truth->pt() * m_conversionFactor;
  //     leadingTruthEta = truth->eta();
  //     leadingTruthPhi = truth->phi();
  //     leadingTruthE = truth->e() * m_conversionFactor;
  //     leadingTruthID = truth->pdgId();
  //     leadingTruthBarcode = truth->barcode();

  //     // stop iterating once we have a photon
  //     break;
  //   } // end loop over truth container

  //   m_pt.push_back( leadingTruthPt );
  //   m_eta.push_back( leadingTruthEta );
  //   m_phi.push_back( leadingTruthPhi );
  //   m_e.push_back( leadingTruthE );
  //   m_truthID.push_back( leadingTruthID );
  //   m_truthBarcode.push_back( leadingTruthBarcode );

  //   // fill tree
  //   m_truthTree->Fill();

  //   // clear vectors before going to next event
  //   m_pt.clear();
  //   m_eta.clear();
  //   m_phi.clear();
  //   m_e.clear();
  //   m_truthID.clear();
  //   m_truthBarcode.clear();
  // }

  // TRUTH JETS
  if( !m_truthJetContainerName.empty() ) {
    SG::ReadHandle<xAOD::JetContainer> inContainer(m_jetTruthContKey);

    for( const xAOD::Jet *jet : *inContainer ) {
      m_pt.push_back( jet->pt() * m_conversionFactor );
      m_eta.push_back( jet->eta() );
      m_phi.push_back( jet->phi() );
      m_e.push_back( jet->e() * m_conversionFactor );

      // isolation check: see if there are any other truth jets with pT > 7 GeV within DeltaR < 1
      bool foundClosebyJet = false;
      for( const xAOD::Jet *jet2 : *inContainer ) {
	if( jet2 == jet || jet2->pt() < 7000 ) {
	  continue;
	}

	double jetPhiDiff = TVector2::Phi_mpi_pi(jet->phi() - jet2->phi());
	double jetEtaDiff = jet->eta() - jet2->eta();
	double jetDeltaR = sqrt( jetPhiDiff*jetPhiDiff + jetEtaDiff*jetEtaDiff );

	if( jetDeltaR < 1.0 ) {
	  foundClosebyJet = true;
	  break;
	}
      } //end inner loop over truth jets

      if( !foundClosebyJet ) m_isIsoJetDR1p0.push_back( 1 );
      else m_isIsoJetDR1p0.push_back( 0 );

    } //end loop over truth jets

    // fill tree
    m_truthJetTree->Fill();

    // clear vectors before going to next event
    m_pt.clear();
    m_eta.clear();
    m_phi.clear();
    m_e.clear();
    m_isIsoJetDR1p0.clear();
  }

  // TRUTH WZ JETS
  if( !m_truthWZJetContainerName.empty() ) {
    SG::ReadHandle<xAOD::JetContainer> inContainer(m_jetTruthWZContKey);

    for( const xAOD::Jet *jet : *inContainer ) {
      m_pt.push_back( jet->pt() * m_conversionFactor );
      m_eta.push_back( jet->eta() );
      m_phi.push_back( jet->phi() );
      m_e.push_back( jet->e() * m_conversionFactor );

      // isolation check: see if there are any other truth jets with pT > 7 GeV within DeltaR < 1
      bool foundClosebyJet = false;
      for( const xAOD::Jet *jet2 : *inContainer ) {
	if( jet2 == jet || jet2->pt() < 7000 ) {
	  continue;
	}

	double jetPhiDiff = TVector2::Phi_mpi_pi(jet->phi() - jet2->phi());
	double jetEtaDiff = jet->eta() - jet2->eta();
	double jetDeltaR = sqrt( jetPhiDiff*jetPhiDiff + jetEtaDiff*jetEtaDiff );

	if( jetDeltaR < 1.0 ) {
	  foundClosebyJet = true;
	  break;
	}
      } //end inner loop over truth jets

      if( !foundClosebyJet ) m_isIsoJetDR1p0.push_back( 1 );
      else m_isIsoJetDR1p0.push_back( 0 );

    } //end loop over truth jets

    // fill tree
    m_truthWZJetTree->Fill();

    // clear vectors before going to next event
    m_pt.clear();
    m_eta.clear();
    m_phi.clear();
    m_e.clear();
    m_isIsoJetDR1p0.clear();
  }
  
  /*********************
   * FOR DEBUG STUDIES
   *********************/
  // SG::ReadHandle<xAOD::FlowElementContainer> inGlobalNeutralFEHandle = makeHandle(m_inGlobalNeutralFEKey);
  // SG::ReadHandle<xAOD::FlowElementContainer> inGlobalChargedFEHandle = makeHandle(m_inGlobalChargedFEKey);
  // SG::ReadHandle<xAOD::FlowElementContainer> inCHSGNeutralFEHandle = makeHandle(m_inCHSGNeutralFEKey);
  // SG::ReadHandle<xAOD::FlowElementContainer> inCHSGChargedFEHandle = makeHandle(m_inCHSGChargedFEKey);
  // SG::ReadHandle<xAOD::FlowElementContainer> inCSSKGNeutralFEHandle = makeHandle(m_inCSSKGNeutralFEKey);
  // SG::ReadHandle<xAOD::FlowElementContainer> inCSSKGChargedFEHandle = makeHandle(m_inCSSKGChargedFEKey);
  // SG::ReadHandle<xAOD::ElectronContainer> inContainer(m_eleContKey);
  // SG::ReadHandle<xAOD::TruthParticleContainer> truthParticles(m_truthParticleContKey);
  // SG::ReadHandle<xAOD::JetContainer> ufoJets("AntiKt10UFOCSSKJets");
  // SG::ReadHandle<xAOD::JetContainer> pflowJets("AntiKt4EMPFlowJets");

  // // this block tests fetching UFO jet constituents and the otherObjects() and chargedObjects() of those constituents
  // if( false ) {
  //   for( const xAOD::Jet *jet : *ufoJets ) {
  //     for( size_t consti = 0; consti < jet->numConstituents(); consti++) {
  // 	const xAOD::FlowElement *constit = static_cast<const xAOD::FlowElement*>(jet->rawConstituent(consti)); //assumption: this returns constituents in UFOCSSK
  // 	//ANA_MSG_INFO( "Trying to access UFO constituent energy... " << constit->e() << " Did it work? If so, rawConstituent is valid." );

  // 	if( constit->isCharged() ) {
  // 	  std::vector<const xAOD::IParticle*> chargedObjects = constit->chargedObjects(); //assumption: this returns charged objects in CSSKGChargedParticleFlowObjects
  // 	  Int_t chargedObjSize = chargedObjects.size();
  // 	  //ANA_MSG_INFO( "Size of UFO jet constituent chargedObjects() vector: " << chargedObjSize );
	
  // 	  if( chargedObjSize > 0 && chargedObjects.at(0) ) ANA_MSG_INFO( "Valid charged object..." );
  // 	  else ANA_MSG_INFO( "Charged object is a nullptr." );
  // 	}

  // 	if( !constit->isCharged() ) {
  // 	  std::vector<const xAOD::IParticle*> otherObjects = constit->otherObjects(); //assumption: this returns other objects in CSSKGNeutralParticleFlowObjects
  // 	  Int_t otherObjSize = otherObjects.size();
  // 	  //ANA_MSG_INFO( "Size of UFO jet constituent otherObjects() vector: " << otherObjSize );
	
  // 	  if( otherObjSize > 0 && otherObjects.at(0) ) ANA_MSG_INFO( "Valid other object..." );
  // 	  else ANA_MSG_INFO( "Other object is a nullptr." );
  // 	}
  //     }
  //   }
  // }

  // // this block tests the isCharged() accessor for FEs
  // if( false ) {
  //   for( const xAOD::FlowElement *fe : *inCHSGNeutralFEHandle ) {
  //     const xAOD::FlowElement *fe_orig = static_cast<const xAOD::FlowElement*>(getOriginalObject(*fe));
  //     if( fe_orig->isCharged() ) ANA_MSG_INFO( "This FE is charged!" );
  //     else ANA_MSG_INFO( "This FE is neutral!" );
  //   }
  // }

  // // this block checks the size of the calpfo vector
  // if( true ) {
  //   for( const xAOD::FlowElement *fe : *inGlobalNeutralFEHandle ) {
  //     SG::AuxElement::ConstAccessor< std::vector<std::pair<unsigned int,double>> > calpfoVec("calpfo_100LeadingTruthParticleBarcodeEnergyPairs");
  //     std::vector<std::pair<unsigned int,double>> barcodeEnergyPair = calpfoVec(*fe);
  //     ANA_MSG_INFO( "Got the calpfo vector! It has size: " << barcodeEnergyPair.size() );
  //   }
  // }

  // // this block tests the calibration hits of neutral FEs linked to electrons
  // if( false ) {
  //   Int_t countElectron = 0;
  //   for( const xAOD::Electron *electron : *inContainer ) {

  //     // decorations
  //     SG::ReadDecorHandle<xAOD::ElectronContainer, std::vector<ElementLink<xAOD::FlowElementContainer>>> neutralFEReadDecorHandle(m_electronNeutralFEReadDecorKey);
  //     SG::ReadDecorHandle<xAOD::ElectronContainer, std::vector<ElementLink<xAOD::FlowElementContainer>>> chargedFEReadDecorHandle(m_electronChargedFEReadDecorKey);
  //     std::vector<ElementLink<xAOD::FlowElementContainer>> electronNFELinks = neutralFEReadDecorHandle(*electron);
  //     std::vector<ElementLink<xAOD::FlowElementContainer>> electronCFELinks = chargedFEReadDecorHandle(*electron);

  //     // vectors to store calibration hit information for each neutral FE
  //     // by default, (up to) the three largest contributions is saved to the calpfo vector
  //     std::vector<std::vector<int>> calpfo_NLeadingTruthParticlePdgId;
  //     std::vector<std::vector<int>> calpfo_NLeadingTruthParticleBarcode;
  //     std::vector<std::vector<double>> calpfo_NLeadingTruthParticleEnergy;

  //     // // also prepare vectors for CaloCalTopoCluster calibration hits
  //     // std::vector<std::vector<int>> calclus_NLeadingTruthParticlePdgId;
  //     // std::vector<std::vector<int>> calclus_NLeadingTruthParticleBarcode;
  //     // std::vector<std::vector<double>> calclus_NLeadingTruthParticleEnergy;

  //     // loop over each neutral FE linked to the electron
  //     for( ElementLink<xAOD::FlowElementContainer> feLink : electronNFELinks ) {
  // 	if( feLink.isValid() ) {
  // 	  const xAOD::FlowElement *electronNeutralGlobalFlowElement = *feLink;

  // 	  ANA_MSG_INFO( "Electron number " << countElectron << " is linked to a neutral FE with index: " << electronNeutralGlobalFlowElement->index() << " with energy: " << electronNeutralGlobalFlowElement->e() * 0.001 << " GeV." );

  // 	  // cells removed from the neutral PFO
  // 	  SG::ReadDecorHandle<xAOD::FlowElementContainer, std::vector<std::pair<unsigned int, double>>> neutralFECellsRemovedReadDecorHandle(m_neutralFECellsRemovedReadDecorKey);
  // 	  std::vector<std::pair<unsigned int, double>> caloCellsRemoved = neutralFECellsRemovedReadDecorHandle(*electronNeutralGlobalFlowElement);
  // 	  ANA_MSG_INFO( "Size of the cellsRemovedFromNeutralFE decoration: " << caloCellsRemoved.size() );

  // 	  // calibration hit information
  // 	  ANA_MSG_INFO( "Fetching calpfo vector for the neutral FE linked to this electron..." );
  // 	  SG::AuxElement::ConstAccessor< std::vector<std::pair<unsigned int,double>> > calpfoVec("calpfo_100LeadingTruthParticleBarcodeEnergyPairs");
  // 	  std::vector<std::pair<unsigned int,double>> barcodeEnergyPair = calpfoVec(*electronNeutralGlobalFlowElement);
  // 	  ANA_MSG_INFO( "Got the calpfo vector! Here are its details:" );

  // 	  // find truth particle with matching barcode and check its PDG ID
  // 	  std::vector<Int_t> truthIDs;
  // 	  std::vector<Int_t> truthBarcodes;
  // 	  std::vector<Double_t> truthEnergies;
  // 	  for( Size_t truthContrib = 0; truthContrib < barcodeEnergyPair.size(); truthContrib++ ) {
  // 	    bool foundMatchingBarcode = false;
  // 	    for( const xAOD::TruthParticle *truthParticle : *truthParticles ) {
  // 	      if( barcodeEnergyPair.at(truthContrib).first != truthParticle->barcode() ) continue;
  // 	      foundMatchingBarcode = true;
  // 	      truthIDs.push_back( truthParticle->pdgId() );
  // 	      truthBarcodes.push_back( truthParticle->barcode() );
  // 	      truthEnergies.push_back( barcodeEnergyPair.at(truthContrib).second * 0.001 );

  // 	      ANA_MSG_INFO( "[electron-linked neutral FE]   PDGID: " << truthParticle->pdgId() << "   Barcode: " << truthParticle->barcode() << "   energy deposited: " << barcodeEnergyPair.at(truthContrib).second * 0.001 << " GeV" );
  // 	      break;
  // 	    }
  // 	    if( !foundMatchingBarcode ) ANA_MSG_INFO( "Huh..? Didn't find any matching barcodes... This is unexpected." );
  // 	  } //end loop over calpfo vector
  // 	  calpfo_NLeadingTruthParticlePdgId.push_back(truthIDs);
  // 	  calpfo_NLeadingTruthParticleBarcode.push_back(truthBarcodes);
  // 	  calpfo_NLeadingTruthParticleEnergy.push_back(truthEnergies);

  // 	  // // save calibration hit information for the topo-cluster corresponding to this neutral FE (if decoration is available)
  // 	  // if( m_calclusIsAvail ) {
  // 	  //   xAOD::CaloCluster *linkedCluster = (xAOD::CaloCluster*) electronNeutralGlobalFlowElement->otherObjects().at(0);

  // 	  //   ANA_MSG_INFO( "Fetching calclus vector for the topo-cluster linked to the neutral FE that is linked to this electron..." );
  // 	  //   SG::AuxElement::ConstAccessor< std::vector<std::pair<unsigned int,double>> > calpfoVec("calclus_NLeadingTruthParticleBarcodeEnergyPairs");
  // 	  //   barcodeEnergyPair = calpfoVec(*linkedCluster);
  // 	  //   ANA_MSG_INFO( "Got the calclus vector! Here are its details:" );

  // 	  //   // find truth particle with matching barcode and check its PDG ID
  // 	  //   truthIDs.clear();
  // 	  //   truthBarcodes.clear();
  // 	  //   truthEnergies.clear();
  // 	  //   for( Size_t truthContrib = 0; truthContrib < barcodeEnergyPair.size(); truthContrib++ ) {
  // 	  //     bool foundMatchingBarcode = false;
  // 	  //     for( const xAOD::TruthParticle *truthParticle : *truthParticles ) {
  // 	  //       if( barcodeEnergyPair.at(truthContrib).first != truthParticle->barcode() ) continue;
  // 	  //       foundMatchingBarcode = true;
  // 	  //       truthIDs.push_back( truthParticle->pdgId() );
  // 	  //       truthBarcodes.push_back( truthParticle->barcode() );
  // 	  //       truthEnergies.push_back( barcodeEnergyPair.at(truthContrib).second * 0.001 );

  // 	  //       ANA_MSG_INFO( "[topo-cluster linked to the electron-linked neutral FE]   PDGID: " << truthParticle->pdgId() << "   Barcode: " << truthParticle->barcode() << "   energy deposited: " << barcodeEnergyPair.at(truthContrib).second * 0.001 << " GeV" );
  // 	  //       break;
  // 	  //     }
  // 	  //     if( !foundMatchingBarcode ) ANA_MSG_INFO( "Huh..? Didn't find any matching barcodes... This is unexpected." );
  // 	  //   } //end loop over calpfo vector
  // 	  //   calclus_NLeadingTruthParticlePdgId.push_back(truthIDs);
  // 	  //   calclus_NLeadingTruthParticleBarcode.push_back(truthBarcodes);
  // 	  //   calclus_NLeadingTruthParticleEnergy.push_back(truthEnergies);
  // 	  // }
  // 	}
  //     } //end loop over linked neutral FEs
  //   } //end loop over electrons
  // }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode mySimpleAlg :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.

  ANA_MSG_DEBUG("Calling postExecute");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode mySimpleAlg :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode mySimpleAlg :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.

  ANA_MSG_INFO( "Calling histFinalize" );
  ANA_CHECK( xAH::Algorithm::algFinalize() );

  return EL::StatusCode::SUCCESS;
}
