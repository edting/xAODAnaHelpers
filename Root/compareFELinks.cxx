/******************************************
 *
 * compareFELinks:
 *
 * produces ntuples containing electron, muon, photon and jet kinematic variables (pT, eta, etc.)
 * as well as a vector of values corresponding to the index of the Flow Elements they each link to
 * n.b. for jets, these values are extracted from constituents
 *
 * Edmund Ting
 * Sep 2021
 *
 ******************************************/

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
#include "xAODAnaHelpers/compareFELinks.h"

// ROOT includes:
#include "TFile.h"
#include "TSystem.h"

// this is needed to distribute the algorithm to the workers
ClassImp(compareFELinks)

compareFELinks :: compareFELinks () :
    Algorithm("compareFELinks")
{
}

// free up memory in destructor
compareFELinks :: ~compareFELinks () {}

EL::StatusCode compareFELinks :: setupJob (EL::Job& job)
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
  xAOD::Init( "compareFELinks" ).ignore(); // call before opening first file

  // create an output stream to save ntuples to
  EL::OutputStream myOutputNtuples("ntuples");
  job.outputAdd( myOutputNtuples );

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode compareFELinks :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  ANA_CHECK( xAH::Algorithm::algInitialize() );

  // create output file
  //  m_file = new TFile(m_outFileName.c_str(), "RECREATE"); //this method will not work on the grid; use the three lines below instead
  TFile *m_file = wk()->getOutputFile( "ntuples" );
  m_file->mkdir( m_outFileName.c_str() );
  m_file->cd( m_outFileName.c_str() );

  // store run and event number
  m_infoTree = new TTree(m_infoTreeName.c_str(),"infoTree");
  m_infoTree->Branch("runNumber", &m_runNumber, "runNumber/I");
  m_infoTree->Branch("eventNumber", &m_evtNumber, "eventNumber/I");
  m_infoTree->Branch("nVtx", &m_nv, "nVtx/I");
  m_infoTree->Branch("nPriVtx", &m_npv, "nPriVtx/I");
  m_infoTree->Branch("nPileUpVtx", &m_npuv, "nPileUpVtx/I");
  m_infoTree->Branch("averageInteractionsPerCrossing", &m_mu, "averageInteractionsPerCrossing/I");
  m_infoTree->Branch("mcChannelNumber", &m_mcChannelNumber, "mcChannelNumber/I");
  m_infoTree->Branch("mcEventWeights", &m_mcEventWeights, "mcEventWeights/D");
  m_infoTree->Branch("initialSumW", &m_MD_initialSumW, "initialSumW/D"); //weighted number of events

  // create tree(s) and branches, for containers whose strings aren't empty
  if( !m_pflowJetContainerName.empty() ) {
    m_pflowJetTree = new TTree(m_pflowJetTreeName.c_str(),"pflowJetTree");
    m_pflowJetTree->Branch("pt", &m_pt);
    m_pflowJetTree->Branch("eta", &m_eta);
    m_pflowJetTree->Branch("phi", &m_phi);
    m_pflowJetTree->Branch("e", &m_e);
    m_pflowJetTree->Branch("isTag", &m_isTag);
    //    m_pflowJetTree->Branch("sf", &m_sf);
    m_pflowJetTree->Branch("nTrk", &m_nTrk);
    m_pflowJetTree->Branch("passJVT", &m_passJVT);
    m_pflowJetTree->Branch("pt_constit", &m_pt_constit);
    m_pflowJetTree->Branch("eta_constit", &m_eta_constit);
    m_pflowJetTree->Branch("phi_constit", &m_phi_constit);
    m_pflowJetTree->Branch("e_constit", &m_e_constit);
    m_pflowJetTree->Branch("isIsoJetDR0p6", &m_isIsoJetDR0p6);
    m_pflowJetTree->Branch("detEta", &m_detEta);
    m_pflowJetTree->Branch("neutralPFOindex", &m_neutralPFOindex);
    m_pflowJetTree->Branch("chargedPFOindex", &m_chargedPFOindex);
    m_pflowJetTree->Branch("neutralPFOenergy", &m_neutralPFOenergy);
    m_pflowJetTree->Branch("chargedPFOenergy", &m_chargedPFOenergy);
    m_pflowJetTree->Branch("calpfo_NLeadingTruthParticlePdgId", &m_calpfo_NLeadingTruthParticlePdgId);
    m_pflowJetTree->Branch("calpfo_NLeadingTruthParticleBarcode", &m_calpfo_NLeadingTruthParticleBarcode);
    m_pflowJetTree->Branch("calpfo_NLeadingTruthParticleEnergy", &m_calpfo_NLeadingTruthParticleEnergy);
  }

  if( !m_topoJetContainerName.empty() ) {
    m_topoJetTree = new TTree(m_topoJetTreeName.c_str(),"topoJetTree");
    m_topoJetTree->Branch("pt", &m_pt);
    m_topoJetTree->Branch("eta", &m_eta);
    m_topoJetTree->Branch("phi", &m_phi);
    m_topoJetTree->Branch("e", &m_e);
    m_topoJetTree->Branch("isTag", &m_isTag);
    //    m_topoJetTree->Branch("sf", &m_sf);
    m_topoJetTree->Branch("nTrk", &m_nTrk);
    m_topoJetTree->Branch("passJVT", &m_passJVT);
    m_topoJetTree->Branch("detEta", &m_detEta);
  }

  if( !m_electronContainerName.empty() ) {
    m_electronTree = new TTree(m_electronTreeName.c_str(),"electronTree");
    m_electronTree->Branch("pt", &m_pt);
    m_electronTree->Branch("eta", &m_eta);
    m_electronTree->Branch("phi", &m_phi);
    m_electronTree->Branch("e", &m_e);
    m_electronTree->Branch("passSel", &m_passSel);
    m_electronTree->Branch("neutralPFOindex", &m_neutralPFOindex);
    m_electronTree->Branch("chargedPFOindex", &m_chargedPFOindex);
    m_electronTree->Branch("neutralPFOenergy", &m_neutralPFOenergy);
    m_electronTree->Branch("chargedPFOenergy", &m_chargedPFOenergy);
    m_electronTree->Branch("calpfo_NLeadingTruthParticlePdgId", &m_calpfo_NLeadingTruthParticlePdgId);
    m_electronTree->Branch("calpfo_NLeadingTruthParticleBarcode", &m_calpfo_NLeadingTruthParticleBarcode);
    m_electronTree->Branch("calpfo_NLeadingTruthParticleEnergy", &m_calpfo_NLeadingTruthParticleEnergy);
  }

  if( !m_photonContainerName.empty() ) {
    m_photonTree = new TTree(m_photonTreeName.c_str(),"photonTree");
    m_photonTree->Branch("pt", &m_pt);
    m_photonTree->Branch("eta", &m_eta);
    m_photonTree->Branch("phi", &m_phi);
    m_photonTree->Branch("e", &m_e);
    m_photonTree->Branch("passSel", &m_passSel);
    m_photonTree->Branch("neutralPFOindex", &m_neutralPFOindex);
    m_photonTree->Branch("chargedPFOindex", &m_chargedPFOindex);
    m_photonTree->Branch("neutralPFOenergy", &m_neutralPFOenergy);
    m_photonTree->Branch("chargedPFOenergy", &m_chargedPFOenergy);
  }

  if( !m_muonContainerName.empty() ) {
    m_muonTree = new TTree(m_muonTreeName.c_str(),"muonTree");
    m_muonTree->Branch("pt", &m_pt);
    m_muonTree->Branch("eta", &m_eta);
    m_muonTree->Branch("phi", &m_phi);
    m_muonTree->Branch("e", &m_e);
    m_muonTree->Branch("passSel", &m_passSel);
    m_muonTree->Branch("neutralPFOindex", &m_neutralPFOindex);
    m_muonTree->Branch("chargedPFOindex", &m_chargedPFOindex);
    m_muonTree->Branch("neutralPFOenergy", &m_neutralPFOenergy);
    m_muonTree->Branch("chargedPFOenergy", &m_chargedPFOenergy);
  }

  if( !m_truthContainerName.empty() ) {
    m_truthTree = new TTree(m_truthTreeName.c_str(),"truthTree");
    m_truthTree->Branch("barcode", &m_truthBarcode);
    m_truthTree->Branch("pdgId", &m_truthID);
    m_truthTree->Branch("pt", &m_pt);
    m_truthTree->Branch("eta", &m_eta);
    m_truthTree->Branch("phi", &m_phi);
    m_truthTree->Branch("e", &m_e);
  }

  if( !m_truthJetContainerName.empty() ) {
    m_truthJetTree = new TTree(m_truthJetTreeName.c_str(),"truthJetTree");
    m_truthJetTree->Branch("pt", &m_pt);
    m_truthJetTree->Branch("eta", &m_eta);
    m_truthJetTree->Branch("phi", &m_phi);
    m_truthJetTree->Branch("e", &m_e);
    m_truthJetTree->Branch("isIsoJetDR1p0", &m_isIsoJetDR1p0);
  }

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



EL::StatusCode compareFELinks :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed

  // Below was copied from BasicEventSelection.cxx
  // get TEvent and TStore - must be done here b/c we need to retrieve CutBookkeepers container from TEvent!
  //
  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();

  // get the MetaData tree once a new file is opened, with
  //
  TTree* MetaData = dynamic_cast<TTree*>( wk()->inputFile()->Get("MetaData") );
  if ( !MetaData ) {
    ANA_MSG_ERROR( "MetaData tree not found! Exiting.");
    return EL::StatusCode::FAILURE;
  }
  MetaData->LoadTree(0);

  //---------------------------
  // Meta data - CutBookkepers
  //---------------------------
  //
  // Metadata for intial N (weighted) events are used to correctly normalise MC
  // if running on a MC DAOD which had some skimming applied at the derivation stage

  //check if file is from a DxAOD
  bool m_isDerivation = !MetaData->GetBranch("StreamAOD");

  if( m_useMetaData ) {

    // Check for potential file corruption
    //
    // If there are some Incomplete CBK, throw a WARNING,
    // unless ALL of them have inputStream == "unknownStream"
    //
    const xAOD::CutBookkeeperContainer* incompleteCBC(nullptr);
    if ( !m_event->retrieveMetaInput(incompleteCBC, "IncompleteCutBookkeepers").isSuccess() ) {
      ANA_MSG_ERROR("Failed to retrieve IncompleteCutBookkeepers from MetaData! Exiting.");
      return EL::StatusCode::FAILURE;
    }
    bool allFromUnknownStream(true);
    if ( incompleteCBC->size() != 0 ) {

      std::string stream("");
      for ( auto cbk : *incompleteCBC ) {
	ANA_MSG_INFO("Incomplete cbk name: " << cbk->name() << " - stream: " << cbk->inputStream());
	if ( cbk->inputStream() != "unknownStream" ) {
	  allFromUnknownStream = false;
	  stream = cbk->inputStream();
	  break;
	}
      }
      if ( !allFromUnknownStream ) { ANA_MSG_WARNING("Found incomplete CBK from stream: " << stream << ". This is not necessarily a sign of file corruption (incomplete CBK appear when 'maxevents' is set in the AOD jo, for instance), but you may still want to check input file for potential corruption..." ); }

    }

    // Now, let's find the actual information
    //
    const xAOD::CutBookkeeperContainer* completeCBC(nullptr);
    if ( !m_event->retrieveMetaInput(completeCBC, "CutBookkeepers").isSuccess() ) {
      ANA_MSG_ERROR("Failed to retrieve CutBookkeepers from MetaData! Exiting.");
      return EL::StatusCode::FAILURE;
    }

    // Find the smallest cycle number, the original first processing step/cycle
    int minCycle(10000);
    for ( auto cbk : *completeCBC ) {
      if ( !( cbk->name().empty() )  && ( minCycle > cbk->cycle() ) ){ minCycle = cbk->cycle(); }
    }

    // Now, let's actually find the right one that contains all the needed info...
    const xAOD::CutBookkeeper* allEventsCBK(nullptr);
    const xAOD::CutBookkeeper* DxAODEventsCBK(nullptr);

    if ( m_isDerivation ) {
      if(m_derivationName != ""){
	ANA_MSG_INFO("Override auto config to look at DAOD made by Derivation Algorithm: " << m_derivationName);
      }else{
	ANA_MSG_INFO("Will autoconfig to look at DAOD made by Derivation Algorithm.");
      }
    }

    int maxCycle(-1);
    for ( const xAOD::CutBookkeeper* cbk: *completeCBC )
      {
	// Find initial book keeper
	ANA_MSG_INFO("Complete cbk name: " << cbk->name() << " - stream: " << cbk->inputStream() );
	if( cbk->cycle() > maxCycle && cbk->name() == "AllExecutedEvents" && cbk->inputStream() == "StreamAOD" )
	  {
	    allEventsCBK = cbk;
	    maxCycle = cbk->cycle();
	  }

	// Find derivation book keeper
	if ( m_isDerivation )
	  {

	    if(m_derivationName != "")
	      {
		if ( cbk->name() == m_derivationName )
		  DxAODEventsCBK = cbk;
	      }
	    else if( cbk->name().find("Kernel") != std::string::npos )
	      {
		ANA_MSG_INFO("Auto config found DAOD made by Derivation Algorithm: " << cbk->name());
		DxAODEventsCBK = cbk;
	      }
	  } // is derivation

	// // Find and record AOD-level sumW information for all alternate weights
	// //  The nominal AllExecutedEvents will be filled later, due to posibility of multiple CBK entries
	// if(cbk->name().substr(0,37) == "AllExecutedEvents_NonNominalMCWeight_" && cbk->name().length()!=17 && cbk->inputStream() == "StreamAOD")
	//   {
	//     // Extract weight index from name
	//     int32_t idx=std::stoi(cbk->name().substr(37));

	//     // Fill histogram, making sure that there is space
	//     // Note will fill bin index = idx+1
	//     while(idx>=m_histSumW->GetNbinsX())
	//       m_histSumW->LabelsInflate("X");
	//     m_histSumW->Fill(idx, cbk->sumOfEventWeights());
	//   }
      }

    if(allEventsCBK == nullptr) {
      ANA_MSG_WARNING("No allEventsCBK found (this is expected for DataScouting, otherwise not). Event numbers set to 0.");
      //m_MD_initialNevents     = 0;
      m_MD_initialSumW        = 0;
      //m_MD_initialSumWSquared = 0;
    }
    else {
      //m_MD_initialNevents     = allEventsCBK->nAcceptedEvents();
      m_MD_initialSumW        = allEventsCBK->sumOfEventWeights();
      //m_MD_initialSumWSquared = allEventsCBK->sumOfEventWeightsSquared();
    }

    if ( m_isDerivation && !DxAODEventsCBK ) {
      ANA_MSG_ERROR( "No CutBookkeeper corresponds to the selected Derivation Framework algorithm name. Check it with your DF experts! Aborting.");
      return EL::StatusCode::FAILURE;
    }

    //m_MD_finalNevents      = ( m_isDerivation ) ? DxAODEventsCBK->nAcceptedEvents() : m_MD_initialNevents;
    //m_MD_finalSumW      = ( m_isDerivation ) ? DxAODEventsCBK->sumOfEventWeights() : m_MD_initialSumW;
    //m_MD_finalSumWSquared   = ( m_isDerivation ) ? DxAODEventsCBK->sumOfEventWeightsSquared() : m_MD_initialSumWSquared;
  
    // the assigned value for m_MD_initialSumW will be saved to output ntuple in the main event loop

  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode compareFELinks :: changeInput (bool /*firstFile*/)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode compareFELinks :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  ANA_MSG_INFO( "Initializing compareFELinks Interface... ");
  ANA_MSG_INFO( "compareFELinks Interface successfully initialized!" );

  // set conversion factor if we want to save pT and energy in GeV
  if( m_convertMeVToGeV ) m_conversionFactor = 0.001;
  else m_conversionFactor = 1;

  // can't use declareProperty in xAH so just comment these lines out
  // declareProperty("ElectronNeutralFEDecorKey", m_electronNeutralFEReadDecorKey = "Electrons.neutralFELinks");
  // declareProperty("ElectronChargedFEDecorKey", m_electronChargedFEReadDecorKey = "Electrons.chargedFELinks");
  // declareProperty("PhotonNeutralFEDecorKey", m_photonNeutralFEReadDecorKey = "Photons.neutralFELinks");
  // declareProperty("MuonNeutralFEDecorKey", m_muonNeutralFEReadDecorKey = "Muons.neutralFELinks");
  // declareProperty("MuonChargedFEDecorKey", m_muonChargedFEReadDecorKey = "Muons.chargedFELinks");
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

  // // decor keys for original <-> global FE links
  // ANA_CHECK(m_neutralOriginalToGlobalFEReadDecorKey.assign("JetETMissNeutralParticleFlowObjects.LinkToNeutralGlobalFlowElement"));
  // ANA_CHECK(m_chargedOriginalToGlobalFEReadDecorKey.assign("JetETMissChargedParticleFlowObjects.LinkToChargedGlobalFlowElement"));
  // ANA_CHECK(m_chargedOriginalToNeutralGlobalFEReadDecorKey.assign("JetETMissChargedParticleFlowObjects.LinksToNeutralGlobalFlowElements"));
  // ANA_CHECK(m_neutralOriginalToGlobalFEReadDecorKey.initialize());
  // ANA_CHECK(m_chargedOriginalToGlobalFEReadDecorKey.initialize());
  // ANA_CHECK(m_chargedOriginalToNeutralGlobalFEReadDecorKey.initialize());

  // declare jet container jey and initialize ... using METJetAssocTool as example
  //declareProperty( "JetContKey", m_jetContKey );
  if( !m_pflowJetContainerName.empty() ) {
    ANA_CHECK( m_jetContKey.assign(m_pflowJetContainerName));
    ANA_CHECK( m_jetContKey.initialize());
  }

  // std::string outputContainerBase = "CHSParticleFlowObjects";
  // m_outAllFEKey = outputContainerBase;
  // ANA_CHECK(m_outAllFEKey.initialize());
  
  // std::string inputContainerBase = "CHS";
  // m_inChargedFEKey = inputContainerBase + "ChargedParticleFlowObjects";
  // m_inNeutralFEKey = inputContainerBase + "NeutralParticleFlowObjects";

  // ANA_CHECK(m_inChargedFEKey.initialize());
  // ANA_CHECK(m_inNeutralFEKey.initialize());

  //std::string inputContainerBaseJet = "JetETMiss";
  std::string inputContainerBaseJet = "Global";
  m_inGlobalChargedFEKey = inputContainerBaseJet + "ChargedParticleFlowObjects";
  m_inGlobalNeutralFEKey = inputContainerBaseJet + "NeutralParticleFlowObjects";

  ANA_CHECK(m_inGlobalChargedFEKey.initialize());
  ANA_CHECK(m_inGlobalNeutralFEKey.initialize());

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode compareFELinks :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.


  // get output file (required for the method for running on grid)
  TFile *m_file = wk()->getOutputFile( "ntuples" );


  /******************************************************************
   * Selection: for Zee/Zmumu/SinglePhoton samples:
   * only save event if:
   *  - at least 2 reco electrons / 2 muons / 1 photon in event
   * and, if truth container is also specified:
   *  - check e/mu/gam can be deltaR and pT matched to truth ones
   * for ttbar, check for min. 1 truth lepton, then look for matching reco lepton
   * for dijet events, just save all electrons / muons / photons
   ******************************************************************/

  // based on what object containers are specified, figure out the sample we're looking at
  if( !m_electronContainerName.empty() && m_muonContainerName.empty() && m_photonContainerName.empty() ) m_isZee = true;
  else if( m_electronContainerName.empty() && !m_muonContainerName.empty() && m_photonContainerName.empty() ) m_isZmumu = true;
  else if( !m_electronContainerName.empty() && !m_muonContainerName.empty() && m_photonContainerName.empty() ) m_isttbar = true;
  else if( m_electronContainerName.empty() && m_muonContainerName.empty() && !m_photonContainerName.empty() ) m_isSinglePhoton = true;
  else if( !m_electronContainerName.empty() && !m_muonContainerName.empty() && !m_photonContainerName.empty() ) m_isDijet = true;

  // Selections if looking at Zee events
  if( m_isZee ) {
    const xAOD::ElectronContainer *inContainer = nullptr;
    ANA_CHECK( evtStore()->retrieve( inContainer, m_electronContainerName ) );

    // check that event contains at least 2 electrons (otherwise, skip event)
    // also store pt, eta and phi of first two electrons (will be used for truth matching below, if specified)
    int electronCount = 0;
    double leadingElectronPt = 0;
    double leadingElectronEta = 0;
    double leadingElectronPhi = 0;
    double subleadingElectronPt = 0;
    double subleadingElectronEta = 0;
    double subleadingElectronPhi = 0;
    for (const xAOD::Electron *electron : *inContainer) {
      electronCount++;
      if( electronCount == 1 ) {
  	leadingElectronPt = electron->pt() * m_conversionFactor;
  	leadingElectronEta = electron->eta();
  	leadingElectronPhi = electron->phi();
      } else if( electronCount > 1 ) {
  	if( electron->pt() * m_conversionFactor > leadingElectronPt ) {
  	  // set previous electron as subleading, and current electron as leading
  	  subleadingElectronPt = leadingElectronPt;
  	  subleadingElectronEta = leadingElectronEta;
  	  subleadingElectronPhi = leadingElectronPhi;
  	  leadingElectronPt = electron->pt() * m_conversionFactor;
  	  leadingElectronEta = electron->eta();
  	  leadingElectronPhi = electron->phi();
  	} else {
  	  // set current electron as subleading
  	  subleadingElectronPt = electron->pt() * m_conversionFactor;
  	  subleadingElectronEta = electron->eta();
  	  subleadingElectronPhi = electron->phi();
  	}

  	// stop iterating once we have two electrons
  	break;
      } // end if electronCount > 1
    } // end loop over electrons

    // if event contains fewer than two electrons: move onto next event
    if( electronCount < 2 ) return EL::StatusCode::SUCCESS;

    // if we have truth container: check that leading electron pair can be matched (otherwise, also skip event)
    if( !m_truthContainerName.empty() ) {
      const xAOD::TruthParticleContainer *truthContainer = nullptr;
      ANA_CHECK( evtStore()->retrieve( truthContainer, m_truthContainerName ) );

      m_truthPDGID = 11; //PDG ID of electron
      int truthCount = 0;
      double leadingTruthPt = 0;
      double leadingTruthEta = 0;
      double leadingTruthPhi = 0;
      double subleadingTruthPt = 0;
      double subleadingTruthEta = 0;
      double subleadingTruthPhi = 0;
      for (const xAOD::TruthParticle* truth : *truthContainer) {
  	// skip if it's not the particle we're looking for
  	if( std::abs( truth->pdgId() ) != m_truthPDGID ) continue;

  	truthCount++;

  	if( truthCount == 1 ) {
  	  leadingTruthPt = truth->pt() * m_conversionFactor;
  	  leadingTruthEta = truth->eta();
  	  leadingTruthPhi = truth->phi();
  	} else if( truthCount > 1 ) {
  	  if( truth->pt() * m_conversionFactor > leadingTruthPt ) {
  	    // set previous truth as subleading, and current truth as leading
  	    subleadingTruthPt = leadingTruthPt;
  	    subleadingTruthEta = leadingTruthEta;
  	    subleadingTruthPhi = leadingTruthPhi;
  	    leadingTruthPt = truth->pt() * m_conversionFactor;
  	    leadingTruthEta = truth->eta();
  	    leadingTruthPhi = truth->phi();
  	  } else {
  	    // set current truth as subleading
  	    subleadingTruthPt = truth->pt() * m_conversionFactor;
  	    subleadingTruthEta = truth->eta();
  	    subleadingTruthPhi = truth->phi();
  	  }

  	  // stop iterating once we have two truths
  	  break;
  	} // end if truthCount > 1
      } // end loop over truth container

      // calculate deltaR and ptDiffRatio for truth vs reco electrons, for both leading and subleading
      const double leadingTruthRecoEtaDiff = leadingElectronEta - leadingTruthEta;
      const double leadingTruthRecoPhiDiff = leadingElectronPhi - leadingTruthPhi;
      const double leadingTruthRecoDeltaR = sqrt( leadingTruthRecoEtaDiff*leadingTruthRecoEtaDiff + leadingTruthRecoPhiDiff*leadingTruthRecoPhiDiff );
      const double leadingTruthRecoPtDiffRatio = std::abs( leadingElectronPt - leadingTruthPt ) / ( leadingElectronPt + leadingTruthPt );

      const double subleadingTruthRecoEtaDiff = subleadingElectronEta - subleadingTruthEta;
      const double subleadingTruthRecoPhiDiff = subleadingElectronPhi - subleadingTruthPhi;
      const double subleadingTruthRecoDeltaR = sqrt( subleadingTruthRecoEtaDiff*subleadingTruthRecoEtaDiff + subleadingTruthRecoPhiDiff*subleadingTruthRecoPhiDiff );
      const double subleadingTruthRecoPtDiffRatio = std::abs( subleadingElectronPt - subleadingTruthPt ) / ( subleadingElectronPt + subleadingTruthPt );

      // skip event if we don't have BOTH leading and subleading matched, each for BOTH deltaR and pT comparisons
      if( (leadingTruthRecoDeltaR > m_deltaRcutoff) || (leadingTruthRecoPtDiffRatio > m_ptDiffCutoff) || (subleadingTruthRecoDeltaR > m_deltaRcutoff) || (subleadingTruthRecoPtDiffRatio > m_ptDiffCutoff) ) {
      	return EL::StatusCode::SUCCESS;
      }
    }
  }

  // single photon
  if( m_isSinglePhoton ) {
    const xAOD::PhotonContainer *inContainer = nullptr;
    ANA_CHECK( evtStore()->retrieve( inContainer, m_photonContainerName ) );

    // check that event contains at least 1 photon (otherwise, skip event)
    // also store pt, eta and phi of first photon (will be used for truth matching below, if specified)
    int photonCount = 0;
    double leadingPhotonPt = 0;
    double leadingPhotonEta = 0;
    double leadingPhotonPhi = 0;
    for (const xAOD::Photon *photon : *inContainer) {
      photonCount++;
      if( photonCount == 1 ) {
  	leadingPhotonPt = photon->pt() * m_conversionFactor;
  	leadingPhotonEta = photon->eta();
  	leadingPhotonPhi = photon->phi();

  	// stop iterating once we have a photon
  	break;
      }
    } // end loop over photons

    // if event contains no photons: move onto next event
    if( photonCount < 1 ) return EL::StatusCode::SUCCESS;

    // if we have truth container: check that leading photon can be matched (otherwise, also skip event)
    if( !m_truthContainerName.empty() ) {
      const xAOD::TruthParticleContainer *truthContainer = nullptr;
      ANA_CHECK( evtStore()->retrieve( truthContainer, m_truthContainerName ) );

      m_truthPDGID = 22;
      int truthCount = 0;
      double leadingTruthPt = 0;
      double leadingTruthEta = 0;
      double leadingTruthPhi = 0;
      for (const xAOD::TruthParticle* truth : *truthContainer) {
  	// skip if it's not the particle we're looking for
  	if( std::abs( truth->pdgId() ) != m_truthPDGID ) continue;

  	truthCount++;

  	if( truthCount == 1 ) {
  	  leadingTruthPt = truth->pt() * m_conversionFactor;
  	  leadingTruthEta = truth->eta();
  	  leadingTruthPhi = truth->phi();

  	  break;
  	}
      } // end loop over truth container

      // calculate deltaR and ptDiffRatio for truth vs reco photon
      const double leadingTruthRecoEtaDiff = leadingPhotonEta - leadingTruthEta;
      const double leadingTruthRecoPhiDiff = leadingPhotonPhi - leadingTruthPhi;
      const double leadingTruthRecoDeltaR = sqrt( leadingTruthRecoEtaDiff*leadingTruthRecoEtaDiff + leadingTruthRecoPhiDiff*leadingTruthRecoPhiDiff );
      const double leadingTruthRecoPtDiffRatio = std::abs( leadingPhotonPt - leadingTruthPt ) / ( leadingPhotonPt + leadingTruthPt );

      // skip event if we don't have a truth--reco match
      if( (leadingTruthRecoDeltaR > m_deltaRcutoff) || (leadingTruthRecoPtDiffRatio > m_ptDiffCutoff) ) {
  	return EL::StatusCode::SUCCESS;
      }
    }
  }

  // Zmumu
  if( m_isZmumu ) {
    const xAOD::MuonContainer *inContainer = nullptr;
    ANA_CHECK( evtStore()->retrieve( inContainer, m_muonContainerName ) );

    // check that event contains at least 2 muons (otherwise, skip event)
    // also store pt, eta and phi of first two muons (will be used for truth matching below, if specified)
    int muonCount = 0;
    double leadingMuonPt = 0;
    double leadingMuonEta = 0;
    double leadingMuonPhi = 0;
    double subleadingMuonPt = 0;
    double subleadingMuonEta = 0;
    double subleadingMuonPhi = 0;
    for (const xAOD::Muon *muon : *inContainer) {
      muonCount++;
      if( muonCount == 1 ) {
  	leadingMuonPt = muon->pt() * m_conversionFactor;
  	leadingMuonEta = muon->eta();
  	leadingMuonPhi = muon->phi();
      } else if( muonCount > 1 ) {
  	if( muon->pt() * m_conversionFactor > leadingMuonPt ) {
  	  // set previous muon as subleading, and current muon as leading
  	  subleadingMuonPt = leadingMuonPt;
  	  subleadingMuonEta = leadingMuonEta;
  	  subleadingMuonPhi = leadingMuonPhi;
  	  leadingMuonPt = muon->pt() * m_conversionFactor;
  	  leadingMuonEta = muon->eta();
  	  leadingMuonPhi = muon->phi();
  	} else {
  	  // set current muon as subleading
  	  subleadingMuonPt = muon->pt() * m_conversionFactor;
  	  subleadingMuonEta = muon->eta();
  	  subleadingMuonPhi = muon->phi();
  	}

  	// stop iterating once we have two muons
  	break;
      } // end if muonCount > 1
    } // end loop over muons

    // if event contains fewer than two muons: move onto next event
    if( muonCount < 2 ) return EL::StatusCode::SUCCESS;

    // if we have truth container: check that leading muon pair can be matched (otherwise, also skip event)
    if( !m_truthContainerName.empty() ) {
      const xAOD::TruthParticleContainer *truthContainer = nullptr;
      ANA_CHECK( evtStore()->retrieve( truthContainer, m_truthContainerName ) );

      m_truthPDGID = 13;
      int truthCount = 0;
      double leadingTruthPt = 0;
      double leadingTruthEta = 0;
      double leadingTruthPhi = 0;
      double subleadingTruthPt = 0;
      double subleadingTruthEta = 0;
      double subleadingTruthPhi = 0;
      for (const xAOD::TruthParticle* truth : *truthContainer) {
  	// skip if it's not the particle we're looking for
  	if( std::abs( truth->pdgId() ) != m_truthPDGID ) continue;

  	truthCount++;

  	if( truthCount == 1 ) {
  	  leadingTruthPt = truth->pt() * m_conversionFactor;
  	  leadingTruthEta = truth->eta();
  	  leadingTruthPhi = truth->phi();
  	} else if( truthCount > 1 ) {
  	  if( truth->pt() * m_conversionFactor > leadingTruthPt ) {
  	    // set previous truth as subleading, and current truth as leading
  	    subleadingTruthPt = leadingTruthPt;
  	    subleadingTruthEta = leadingTruthEta;
  	    subleadingTruthPhi = leadingTruthPhi;
  	    leadingTruthPt = truth->pt() * m_conversionFactor;
  	    leadingTruthEta = truth->eta();
  	    leadingTruthPhi = truth->phi();
  	  } else {
  	    // set current truth as subleading
  	    subleadingTruthPt = truth->pt() * m_conversionFactor;
  	    subleadingTruthEta = truth->eta();
  	    subleadingTruthPhi = truth->phi();
  	  }

  	  // stop iterating once we have two truths
  	  break;
  	} // end if truthCount > 1
      } // end loop over truth container

      // calculate deltaR and ptDiffRatio for truth vs reco muons, for both leading and subleading
      const double leadingTruthRecoEtaDiff = leadingMuonEta - leadingTruthEta;
      const double leadingTruthRecoPhiDiff = leadingMuonPhi - leadingTruthPhi;
      const double leadingTruthRecoDeltaR = sqrt( leadingTruthRecoEtaDiff*leadingTruthRecoEtaDiff + leadingTruthRecoPhiDiff*leadingTruthRecoPhiDiff );
      const double leadingTruthRecoPtDiffRatio = std::abs( leadingMuonPt - leadingTruthPt ) / ( leadingMuonPt + leadingTruthPt );

      const double subleadingTruthRecoEtaDiff = subleadingMuonEta - subleadingTruthEta;
      const double subleadingTruthRecoPhiDiff = subleadingMuonPhi - subleadingTruthPhi;
      const double subleadingTruthRecoDeltaR = sqrt( subleadingTruthRecoEtaDiff*subleadingTruthRecoEtaDiff + subleadingTruthRecoPhiDiff*subleadingTruthRecoPhiDiff );
      const double subleadingTruthRecoPtDiffRatio = std::abs( subleadingMuonPt - subleadingTruthPt ) / ( subleadingMuonPt + subleadingTruthPt );

      // skip event if we don't have BOTH leading and subleading matched, each for BOTH deltaR and pT comparisons
      if( (leadingTruthRecoDeltaR > m_deltaRcutoff) || (leadingTruthRecoPtDiffRatio > m_ptDiffCutoff) || (subleadingTruthRecoDeltaR > m_deltaRcutoff) || (subleadingTruthRecoPtDiffRatio > m_ptDiffCutoff) ) {
  	return EL::StatusCode::SUCCESS;
      }
    }
  }

  // ttbar
  //if( m_isttbar ) {
  if( false ) { //temp: dunno how to save TruthBosonsWithDecayParticles using CA method so just turn off this for now...
    const xAOD::TruthParticleContainer *truthContainer = nullptr;
    ANA_CHECK( evtStore()->retrieve( truthContainer, m_truthContainerName ) );

    // want events with at least one lepton from top decay
    bool foundTruthLepton = false;
    std::vector<int> leptonID;
    std::vector<float> truthLeptonPt;
    std::vector<float> truthLeptonEta;
    std::vector<float> truthLeptonPhi;
    std::vector<float> truthLeptonE;
    for( const xAOD::TruthParticle *truth : *truthContainer ) {
      if( std::abs( truth->pdgId() ) != 11 && std::abs( truth->pdgId() ) != 13 ) continue; //skip if not electron or muon

      foundTruthLepton = true;
      leptonID.push_back( truth->pdgId() );
      truthLeptonPt.push_back( truth->pt() * m_conversionFactor );
      truthLeptonEta.push_back( truth->eta() );
      truthLeptonPhi.push_back( truth->phi() );
      truthLeptonE.push_back( truth->e() * m_conversionFactor );

      if( leptonID.size() == 2 ) break;
    } //end loop over truth particles

    // skip event if both tops decayed hadronically ... though this shouldn't be the case for nonallhad ttbar
    if( !foundTruthLepton ) return EL::StatusCode::SUCCESS;

    // if found a lepton(s), check that reco containers contain minimum required number of same lepton
    int numElectronsExpected = 0;
    int numMuonsExpected = 0;
    int numElectronsFound = 0;
    int numMuonsFound = 0;
    std::vector<float> recoLeptonPt;
    std::vector<float> recoLeptonEta;
    std::vector<float> recoLeptonPhi;
    std::vector<float> recoLeptonE;
    for( int lepton = 0; static_cast<std::vector<int>::size_type>( lepton ) < leptonID.size(); lepton++ ) { //note: the static cast will do funny things if "lepton" is negative! (zero is ok)
      if( std::abs( leptonID.at( lepton ) ) == 11 ) { //check for reco electron
	numElectronsExpected++;

	const xAOD::ElectronContainer *electrons = nullptr;
	ANA_CHECK( evtStore()->retrieve( electrons, m_electronContainerName ) );

	int electronCount = 0;
	for( const xAOD::Electron *electron : *electrons ) {
	  if( electronCount == numElectronsFound ) {
	    recoLeptonPt.push_back( electron->pt() * m_conversionFactor );
	    recoLeptonEta.push_back( electron->eta() );
	    recoLeptonPhi.push_back( electron->phi() );
	    recoLeptonE.push_back( electron->e() * m_conversionFactor );
	    numElectronsFound++;
	    break;
	  }

	  electronCount++;
	}
      } else if( std::abs( leptonID.at( lepton ) ) == 13 ) { //check for reco muon
	numMuonsExpected++;

	const xAOD::MuonContainer *muons = nullptr;
	ANA_CHECK( evtStore()->retrieve( muons, m_muonContainerName ) );

	int muonCount = 0;
	for( const xAOD::Muon *muon : *muons ) {
	  if( muonCount == numMuonsFound ) {
	    recoLeptonPt.push_back( muon->pt() * m_conversionFactor );
	    recoLeptonEta.push_back( muon->eta() );
	    recoLeptonPhi.push_back( muon->phi() );
	    recoLeptonE.push_back( muon->e() * m_conversionFactor );
	    numMuonsFound++;
	    break;
	  }

	  muonCount++;
	}
      }
    } //break out of loop over truth leptons

    if( numElectronsExpected != numElectronsFound || numMuonsExpected != numMuonsFound )
      return EL::StatusCode::SUCCESS;

    // check that electron or muon passes truth-matching requirements based on DeltaR and pTdiffRatio
    for( int lepton = 0; static_cast<std::vector<int>::size_type>( lepton ) < leptonID.size(); lepton++ ) {
      const float etaDiff = truthLeptonEta.at( lepton ) - recoLeptonEta.at( lepton );
      const float phiDiff = truthLeptonPhi.at( lepton ) - recoLeptonPhi.at( lepton );
      const float truthRecoDeltaR = sqrt( etaDiff*etaDiff + phiDiff*phiDiff );
      const float truthRecoPtDiffRatio = std::abs( truthLeptonPt.at( lepton ) - recoLeptonPt.at( lepton ) ) / (truthLeptonPt.at( lepton ) + recoLeptonPt.at( lepton ) );

      if( truthRecoDeltaR > m_deltaRcutoff || truthRecoPtDiffRatio > m_ptDiffCutoff )
	return EL::StatusCode::SUCCESS;
    }

    // pre-emptively fill truth lepton tree in output ntuple since we already did the work to find them above
    m_truthTree->SetDirectory( m_file->GetDirectory( m_outFileName.c_str() ) );
    for( int lepton = 0; static_cast<std::vector<int>::size_type>( lepton ) < leptonID.size(); lepton++ ) {
      m_truthID.push_back( leptonID.at( lepton ) );
      m_pt.push_back( truthLeptonPt.at( lepton ) );
      m_eta.push_back( truthLeptonEta.at( lepton ) );
      m_phi.push_back( truthLeptonPhi.at( lepton ) );
      m_e.push_back( truthLeptonE.at( lepton ) );
    }

    m_truthTree->Fill();
    m_truthID.clear();
    m_pt.clear();
    m_eta.clear();
    m_phi.clear();
    m_e.clear();
  }



  /********************************************************
   * If we passed preliminary checks, start filling trees
   ********************************************************/

  // fetch primary vertex information for this event for later reference (needed to find number of jet tracks)
  size_t vtxIdx = 0;
  m_nv = 0; //reset for each event
  m_npv = 0; //reset for each event
  m_npuv = 0; //reset for each event
  const xAOD::Vertex *priVtx = nullptr;
  const xAOD::VertexContainer* vertices = nullptr;
  if(evtStore()->retrieve(vertices, "PrimaryVertices").isSuccess()) {
    for(const xAOD::Vertex* vtx : *vertices) {
      if(vtx->vertexType() == xAOD::VxType::PriVtx) {
	if( !priVtx ) priVtx = vtx;
	m_npv++;
	//break;
      }

      if(vtx->vertexType() == xAOD::VxType::PileUp) {
	m_npuv++; //count number of pile up vertices in event
      }

      //count total number of vertices in the "PrimaryVertices" container?
      m_nv++;
    }
  } else {
    ANA_MSG_WARNING( "Failed to retrieve primary vertex container" );
  }
  vtxIdx = priVtx->index();

  // Grab run and event number, fill info tree
  const xAOD::EventInfo *ei = nullptr;
  ANA_CHECK( evtStore()->retrieve( ei , "EventInfo" ) );
  m_runNumber = ei->runNumber();
  m_evtNumber = ei->eventNumber();
  m_mu = ei->averageInteractionsPerCrossing();
  m_mcChannelNumber = ei->mcChannelNumber();
  m_mcEventWeights = ei->mcEventWeights().at(0);
  m_infoTree->Fill();

  /// PFLOW JETS
  if( !m_pflowJetContainerName.empty() ) {
    // tell the tree to go to the file
    m_pflowJetTree->SetDirectory( m_file->GetDirectory( m_outFileName.c_str() ) );

    // // get charged and neutral PFO (FE) containers and combine them into one "CHSParticleFlowObjects" container
    // SG::ReadHandle<xAOD::FlowElementContainer> inNeutralFEHandle = makeHandle(m_inNeutralFEKey);
    // SG::ReadHandle<xAOD::FlowElementContainer> inChargedFEHandle = makeHandle(m_inChargedFEKey);

    SG::ReadHandle<xAOD::FlowElementContainer> inGlobalNeutralFEHandle = makeHandle(m_inGlobalNeutralFEKey);
    SG::ReadHandle<xAOD::FlowElementContainer> inGlobalChargedFEHandle = makeHandle(m_inGlobalChargedFEKey);

    // ConstDataVector<xAOD::FlowElementContainer> *m_FlowElements = new ConstDataVector<xAOD::FlowElementContainer>(SG::VIEW_ELEMENTS);

    // (*m_FlowElements).assign((*inNeutralFEHandle).begin(), (*inNeutralFEHandle).end());
    // (*m_FlowElements).insert((*m_FlowElements).end(),
    // 			     (*inChargedFEHandle).begin(),
    // 			     (*inChargedFEHandle).end());

    // ANA_CHECK(evtStore()->record(m_FlowElements, "CHSParticleFlowObjects"));

    // prepare PFlow jet container
    SG::ReadHandle<xAOD::JetContainer> inContainer(m_jetContKey);

    // old(?) way of retrieving jet container
    // const xAOD::JetContainer *inContainer = nullptr;
    // ANA_CHECK( evtStore()->retrieve( inContainer, m_pflowJetContainerName ) );

    // also retrieve truth particle container -- needed for comparing calpfo_... values later
    const xAOD::TruthParticleContainer *truthParticles = nullptr;
    ANA_CHECK( evtStore()->retrieve( truthParticles, "TruthParticles" ) );

    // placeholder if no scale factor information
    static const std::vector<float> junk(1,-999);

    // loop over jets
    for( const xAOD::Jet *jet : *inContainer ) {

      // for each jet, prepare vector to store PFO index of each constituent
      std::vector<int> chargedPFOindex;
      std::vector<int> neutralPFOindex;
      std::vector<double> chargedPFOenergy;
      std::vector<double> neutralPFOenergy;

      // vectors to store truth info from calpfo decoration
      std::vector<std::vector<int>> calpfo_NLeadingTruthParticlePdgId;
      std::vector<std::vector<int>> calpfo_NLeadingTruthParticleBarcode;
      std::vector<std::vector<double>> calpfo_NLeadingTruthParticleEnergy;

      // loop over jet constituents
      ANA_MSG_VERBOSE( "Number of jet constituents: " << jet->numConstituents() );
      for( size_t consti = 0; consti < jet->numConstituents(); consti++) {

      	//Old comment: Points to the combined CHSParticleFlowObjects container
	//Updated comment: Should now point to CHSGParticleFlowObjects
      	const xAOD::FlowElement *fe = static_cast<const xAOD::FlowElement*>(jet->rawConstituent(consti));

      	//Need to figure out now if this is a charged or neutral constituent (use index for this) (note: doing this b/c no aux data for fe->isCharged())
      	//Check if constituent in CHSChargedParticleFlowObjects with index of fe is the same as the fe we're looking at right now
      	//if(inChargedFEHandle->size() > fe->index() && inChargedFEHandle->at(fe->index()) == fe){
	if( fe->isCharged() ){
      	  // Get the original constituent (before four-vector modifications, etc) from the JetETMissChargedParticleFlowObjects container (using index)
      	  // This object can be then used to retrieve additional moments like the links to other objects
      	  const xAOD::FlowElement *fe_original = inGlobalChargedFEHandle->at(fe->index());
	  //const xAOD::FlowElement *fe_original = fe; //lol probs not needed

	  chargedPFOindex.push_back(fe_original->index());
	  chargedPFOenergy.push_back(fe_original->e() * m_conversionFactor);

	  // // not sure if this is necessary
	  // neutralPFOindex.push_back(-10);
	  // neutralPFOenergy.push_back(-10);

      	  ANA_MSG_VERBOSE( "Charged Flow Element index: " << fe_original->index() );
      	}
      	//Same but for neutral FEs
      	//Check if constituent in CHSNeutralParticleFlowObjects with index of fe is the same as the fe we're looking at right now
      	//else if(inNeutralFEHandle->size() > fe->index() && inNeutralFEHandle->at(fe->index())== fe){
	else if( !fe->isCharged() ){
      	  // Get the original constituent (before four-vector modifications, etc) from the JetETMissNeutralParticleFlowObjects container (using index) 
      	  // This object can be then used to retrieve additional moments like the links to other objects
      	  const xAOD::FlowElement *fe_original = inGlobalNeutralFEHandle->at(fe->index());
	  //const xAOD::FlowElement *fe_original = fe;

	  neutralPFOindex.push_back(fe_original->index());
	  neutralPFOenergy.push_back(fe_original->e() * m_conversionFactor);

      	  ANA_MSG_VERBOSE( "Neutral Flow Element index: " << fe_original->index() );

	  // for neutral PFOs, also save truth information on what particles deposited how much energy
	  // note: this is a vector because it's possible for multiple particles to have contributed to this PFO (the 3 largest contributions are saved in this decoration)
	  SG::AuxElement::ConstAccessor< std::vector<std::pair<unsigned int,double>> > calpfoVec("calpfo_NLeadingTruthParticleBarcodeEnergyPairs");
	  std::vector<std::pair<unsigned int,double>> barcodeEnergyPair = calpfoVec(*fe_original);
	  std::vector<Int_t> truthIDs;
	  std::vector<Int_t> truthBarcodes;
	  std::vector<Double_t> truthEnergies; //note: this is the amount of energy deposited in the nPFO by the truth particle

      	  ANA_MSG_VERBOSE( "Got the calpfo vector for this NPFO! Here are its details:" );

	  for( Size_t truthContrib = 0; truthContrib < barcodeEnergyPair.size(); truthContrib++ ) {
	    //loop over truth particles to find corresponding barcode
	    bool foundMatchingBarcode = false;
	    for( const xAOD::TruthParticle *truthParticle : *truthParticles ) {
	      if( barcodeEnergyPair.at(truthContrib).first != truthParticle->barcode() ) continue;
	      foundMatchingBarcode = true;
	      truthIDs.push_back( truthParticle->pdgId() );
	      truthBarcodes.push_back( truthParticle->barcode() );
	      truthEnergies.push_back( barcodeEnergyPair.at(truthContrib).second * m_conversionFactor );
	      //truthEnergies.push_back( truthParticle->e() * m_conversionFactor );

	      ANA_MSG_VERBOSE( "[jet constituent NPFO]   PDGID: " << truthParticle->pdgId() << "   Barcode: " << truthParticle->barcode() << "   energy deposited in NPFO: " << barcodeEnergyPair.at(truthContrib).second * m_conversionFactor << " GeV" );
	      break; //found matching truth particle, move onto the next calpfo contribution
	    }
	    if( !foundMatchingBarcode ) ANA_MSG_INFO( "Huh..? Didn't find any matching barcodes... This is unexpected." );
	  } //end loop over calpfo vector
	  calpfo_NLeadingTruthParticlePdgId.push_back(truthIDs);
	  calpfo_NLeadingTruthParticleBarcode.push_back(truthBarcodes);
	  calpfo_NLeadingTruthParticleEnergy.push_back(truthEnergies);
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

      // fetch b-tag and scale factor information
      SG::AuxElement::ConstAccessor< char > isTag("BTag_DL1r_FixedCutBEff_77"); // int: 0 = not b-tagged, 1 = tagged
      SG::AuxElement::ConstAccessor< std::vector<float> > sf("BTag_SF_DL1r_FixedCutBEff_77"); // scale factors

      m_isTag.push_back( isTag.isAvailable(*jet) ? isTag(*jet) : -1 );
      m_sf.push_back( sf.isAvailable(*jet) ? sf(*jet) : junk);

      // test that tagging information is extracted properly
      if( isTag.isAvailable(*jet) && msgLvl(MSG::VERBOSE) ) std::cout << "Jet b-tag flag: " << (int)isTag(*jet) << "\n";

      // find number of tracks associated with jet
      static const SG::AuxElement::ConstAccessor< std::vector<int> > acc("NumTrkPt500");
      int numTracks = acc(*jet).at(vtxIdx);
      m_nTrk.push_back( numTracks );

      // print number of tracks
      ANA_MSG_VERBOSE( "Number of jet tracks: " << numTracks );

      // fetch JVT cut info
      static const SG::AuxElement::ConstAccessor< char > passJVT("JetJVT_Passed_Tight");
      m_passJVT.push_back( passJVT.isAvailable(*jet) ? passJVT(*jet) : -1 );

      // check JVT pass/fail values
      ANA_MSG_VERBOSE( "Pass JVT cut: " << (int)passJVT(*jet) );

      // get constituent-scale quantities
      xAOD::JetFourMom_t constit4vec = jet->getAttribute<xAOD::JetFourMom_t>("JetConstitScaleMomentum");
      m_pt_constit.push_back( constit4vec.Pt() * m_conversionFactor );
      m_eta_constit.push_back( constit4vec.Eta() );
      m_phi_constit.push_back( constit4vec.Phi() );
      m_e_constit.push_back( constit4vec.E() * m_conversionFactor );

      // for each jet, loop over all other jets, check if there are any with pt > 7 GeV within deltaR < 0.6
      bool foundClosebyJet = false;
      for( const xAOD::Jet *jet2 : *inContainer ) {
	// obviously, don't compare jet with itself
	if( jet == jet2 )
	  continue;

	// // this method uses constituent-level jet pT for the comparison
	// // only check against jets with pT > 7 GeV
	// xAOD::JetFourMom_t constit4vec2 = jet2->getAttribute<xAOD::JetFourMom_t>("JetConstitScaleMomentum");
	// if( constit4vec2.Pt() < 7000 )
	//   continue;

	// // now check deltaR between jets
	// double jetPhiDiff = constit4vec.Phi() - constit4vec2.Phi();
	// double jetEtaDiff = constit4vec.Eta() - constit4vec2.Eta();
	// double jetDeltaR = sqrt( jetPhiDiff*jetPhiDiff + jetEtaDiff*jetEtaDiff );

	// this method uses calibrated jet pT for the comparison
	// only check against jets with pT > 7 GeV
	if( jet2->pt() < 7000 )
	  continue;

	// now check deltaR between jets
	double jetPhiDiff = jet->phi() - jet2->phi();
	double jetEtaDiff = jet->eta() - jet2->eta();
	double jetDeltaR = sqrt( jetPhiDiff*jetPhiDiff + jetEtaDiff*jetEtaDiff );

	if( jetDeltaR < 0.6 ) {
	  foundClosebyJet = true;
	  break;
	}
      } //end inner loop over jets

      if( !foundClosebyJet ) m_isIsoJetDR0p6.push_back( 1 );
      else m_isIsoJetDR0p6.push_back( 0 );

      // get detector eta
      m_detEta.push_back( jet->getAttribute<float>("DetectorEta") );

    } //end loop over jets

    // fill tree
    m_pflowJetTree->Fill();

    // clear vectors before going to next event
    m_pt.clear();
    m_eta.clear();
    m_phi.clear();
    m_e.clear();
    m_isTag.clear();
    m_sf.clear();
    m_nTrk.clear();
    m_passJVT.clear();
    m_pt_constit.clear();
    m_eta_constit.clear();
    m_phi_constit.clear();
    m_e_constit.clear();
    m_isIsoJetDR0p6.clear();
    m_detEta.clear();
    m_neutralPFOindex.clear();
    m_chargedPFOindex.clear();
    m_neutralPFOenergy.clear();
    m_chargedPFOenergy.clear();
    m_calpfo_NLeadingTruthParticlePdgId.clear();
    m_calpfo_NLeadingTruthParticleBarcode.clear();
    m_calpfo_NLeadingTruthParticleEnergy.clear();
  }

  /// TOPO JETS
  if( !m_topoJetContainerName.empty() ) {
    // tell the tree to go to the file
    m_topoJetTree->SetDirectory( m_file->GetDirectory( m_outFileName.c_str() ) );

    const xAOD::JetContainer *inContainer = nullptr;
    ANA_CHECK( evtStore()->retrieve( inContainer, m_topoJetContainerName ) );

    // placeholder if no scale factor information
    static const std::vector<float> junk(1,-999);

    // fetch primary vertex information for this event for later reference (needed to find number of jet tracks)
    size_t vtxIdx = 0;
    const xAOD::Vertex *priVtx = nullptr;
    const xAOD::VertexContainer* vertices = nullptr;
    if(evtStore()->retrieve(vertices, "PrimaryVertices").isSuccess()) {
      for(const xAOD::Vertex* vtx : *vertices) {
	if(vtx->vertexType() == xAOD::VxType::PriVtx) {
	  priVtx = vtx;
	  break;
	}
      }
    } else {
      ANA_MSG_WARNING( "Failed to retrieve primary vertex container" );
    }
    vtxIdx = priVtx->index();

    // loop over jets
    for( const xAOD::Jet *jet : *inContainer ) {
      m_pt.push_back( jet->pt() * m_conversionFactor );
      m_eta.push_back( jet->eta() );
      m_phi.push_back( jet->phi() );
      m_e.push_back( jet->e() * m_conversionFactor );

      // fetch b-tag and scale factor information
      SG::AuxElement::ConstAccessor< char > isTag("BTag_DL1r_FixedCutBEff_77"); // int: 0 = not b-tagged, 1 = tagged
      SG::AuxElement::ConstAccessor< std::vector<float> > sf("BTag_SF_DL1r_FixedCutBEff_77"); // scale factors

      m_isTag.push_back( isTag.isAvailable(*jet) ? isTag(*jet) : -1 );
      m_sf.push_back( sf.isAvailable(*jet) ? sf(*jet) : junk);

      // test that tagging information is extracted properly
      if( isTag.isAvailable(*jet) && msgLvl(MSG::VERBOSE) ) std::cout << "Jet b-tag flag: " << (int)isTag(*jet) << "\n";

      // find number of tracks associated with jet
      static const SG::AuxElement::ConstAccessor< std::vector<int> > acc("NumTrkPt500");
      int numTracks = acc(*jet).at(vtxIdx);
      m_nTrk.push_back( numTracks );

      // print number of tracks
      ANA_MSG_VERBOSE( "Number of jet tracks: " << numTracks );

      // fetch JVT cut info
      static const SG::AuxElement::ConstAccessor< char > passJVT("JetJVT_Passed_Tight");
      m_passJVT.push_back( passJVT.isAvailable(*jet) ? passJVT(*jet) : -1 );

      // check JVT pass/fail values
      ANA_MSG_VERBOSE( "Pass JVT cut: " << (int)passJVT(*jet) );

      // get detector eta
      m_detEta.push_back( jet->getAttribute<float>("DetectorEta") );
    }

    // fill tree
    m_topoJetTree->Fill();

    // clear vectors before going to next event
    m_pt.clear();
    m_eta.clear();
    m_phi.clear();
    m_e.clear();
    m_isTag.clear();
    m_sf.clear();
    m_nTrk.clear();
    m_passJVT.clear();
    m_detEta.clear();
  }

  /// ELECTRONS
  if( !m_electronContainerName.empty() ) {
    // tell the tree to go to the file
    m_electronTree->SetDirectory( m_file->GetDirectory( m_outFileName.c_str() ) );

    const xAOD::ElectronContainer *inContainer = nullptr;
    ANA_CHECK( evtStore()->retrieve( inContainer, m_electronContainerName ) );

    // also retrieve truth particle container and define neutral PFO container -- need this to access and use calpfo vector
    // note that these correspond to the Global containers
    SG::ReadHandle<xAOD::FlowElementContainer> inGlobalNeutralFEHandle = makeHandle(m_inGlobalNeutralFEKey);
    SG::ReadHandle<xAOD::FlowElementContainer> inGlobalChargedFEHandle = makeHandle(m_inGlobalChargedFEKey);
    const xAOD::TruthParticleContainer *truthParticles = nullptr;
    ANA_CHECK( evtStore()->retrieve( truthParticles, "TruthParticles" ) );

    // note: each electron will point to a VECTOR of FELinks
    // here, the outer vector is for each electron in the event; nested vector is for each FE linked to corresponding electron
    std::vector<std::vector<ElementLink<xAOD::FlowElementContainer>>> electronNFELinks;
    std::vector<std::vector<ElementLink<xAOD::FlowElementContainer>>> electronCFELinks;

    Int_t countElectron = 0;
    for( const xAOD::Electron *electron : *inContainer ) {
      // accessor for checking whether electron passed identification selections
      SG::AuxElement::ConstAccessor< char > passSel("passSel");

      // prepare vector to store PFO indices and energies for the current electron
      std::vector<int> neutralPFOindex;
      std::vector<int> chargedPFOindex;
      std::vector<double> neutralPFOenergy;
      std::vector<double> chargedPFOenergy;

      // OLD COMMENT: fetch neutral FE links: electron->JetETMissNeutral and JetETMissNeutral->GlobalNeutral
      // MODELV2: decoration points straight to Global container
      SG::ReadDecorHandle<xAOD::ElectronContainer, std::vector<ElementLink<xAOD::FlowElementContainer>>> neutralFEReadDecorHandle(m_electronNeutralFEReadDecorKey);
      //SG::ReadDecorHandle<xAOD::FlowElementContainer, ElementLink<xAOD::FlowElementContainer>> neutralGlobalFEReadDecorHandle(m_neutralOriginalToGlobalFEReadDecorKey);
      electronNFELinks.push_back( neutralFEReadDecorHandle(*electron) ); // note that this corresponds to FE in the JetETMiss container

      // OLD COMMENT: fetch charged FE links: electron->JetETMissCharged, JetETMissCharged->GlobalCharged, and JetETMissCharged->GlobalNeutral (for cFEs that create new nFEs in Global)
      // MODELV2: decoration points straight to Global container
      SG::ReadDecorHandle<xAOD::ElectronContainer, std::vector<ElementLink<xAOD::FlowElementContainer>>> chargedFEReadDecorHandle(m_electronChargedFEReadDecorKey);
      //SG::ReadDecorHandle<xAOD::FlowElementContainer, ElementLink<xAOD::FlowElementContainer>> chargedGlobalFEReadDecorHandle(m_chargedOriginalToGlobalFEReadDecorKey);
      //SG::ReadDecorHandle<xAOD::FlowElementContainer, std::vector<ElementLink<xAOD::FlowElementContainer>>> chargedOriginalToNeutralGlobalFEReadDecorHandle(m_chargedOriginalToNeutralGlobalFEReadDecorKey);
      electronCFELinks.push_back( chargedFEReadDecorHandle(*electron) ); // note that this corresponds to FE in the JetETMiss container

      // vectors to store calpfo info
      // outer vector is for each neutral FE linked to the current electron; inner vector is for each entry in the calpfo vector for the corresponding neutral FE
      std::vector<std::vector<int>> calpfo_NLeadingTruthParticlePdgId;
      std::vector<std::vector<int>> calpfo_NLeadingTruthParticleBarcode;
      std::vector<std::vector<double>> calpfo_NLeadingTruthParticleEnergy;

      // loop over each neutral FE linked to the electron
      for( ElementLink<xAOD::FlowElementContainer> feLink : electronNFELinks.at(countElectron) ) {
	if( feLink.isValid() ) {
	  // OLD COMMENT: this corresponds to neutral FE in the JetETMiss container
	  // MODELV2: this now points to Global
	  //const xAOD::FlowElement *electronNeutralOriginalFlowElement = *feLink;
	  const xAOD::FlowElement *electronNeutralGlobalFlowElement = *feLink;

	  // // fetch link from neutral JetETMiss->Global
	  // ElementLink<xAOD::FlowElementContainer> feLinkToGlobal = neutralGlobalFEReadDecorHandle(*electronNeutralOriginalFlowElement);
	  // if( !feLinkToGlobal.isValid() ) ANA_MSG_INFO( "FE Link to neutral Global is NOT valid... this message should NEVER get printed." );
	  // else ANA_MSG_VERBOSE( "FE Link to neutral Global is valid. As it should be (always)." );

	  // this now corresponds to neutral FE in the Global container
	  // MODELV2: this is no longer needed
	  //const xAOD::FlowElement *electronNeutralGlobalFlowElement = *feLinkToGlobal;

	  // get the index and energy according to FE in the Global container
	  neutralPFOindex.push_back(electronNeutralGlobalFlowElement->index());
	  neutralPFOenergy.push_back(electronNeutralGlobalFlowElement->e() * m_conversionFactor);

	  // printouts for debug
	  ANA_MSG_VERBOSE( "Electron number " << countElectron << " is linked to a neutral FE with index: " << electronNeutralGlobalFlowElement->index() << " with energy: " << electronNeutralGlobalFlowElement->e() * m_conversionFactor << " GeV." );
	  if( electronNeutralGlobalFlowElement->isCharged() )
	    ANA_MSG_INFO( "Umm... expecting a neutral FE but apparently it is charged?? Lol this shouldn't ever get printed." );

	  // save calibration hit information for neutrals: what particles deposited energy (and how much)
	  // note: this is a vector because it's possible for multiple particles to have contributed to this PFO (default: 3 largest contributions are in the calpfo decoration)
	  // get the corresponding neutral PFO from the FlowElement (should be the same index)
	  ANA_MSG_VERBOSE( "Fetching calpfo vector for the neutral FE linked to this electron..." );
	  SG::AuxElement::ConstAccessor< std::vector<std::pair<unsigned int,double>> > calpfoVec("calpfo_NLeadingTruthParticleBarcodeEnergyPairs");
	  const xAOD::FlowElement *electronNPFO = inGlobalNeutralFEHandle->at(electronNeutralGlobalFlowElement->index()); // fetch the info from the Global container
	  std::vector<std::pair<unsigned int,double>> barcodeEnergyPair = calpfoVec(*electronNPFO);

	  std::vector<Int_t> truthIDs;
	  std::vector<Int_t> truthBarcodes;
	  std::vector<Double_t> truthEnergies; //note: this is the amount of energy DEPOSITED in the neutral FE by the truth particle (not the energy of the particle itself)
	  ANA_MSG_VERBOSE( "Got the calpfo vector! Here are its details:" );
	  for( Size_t truthContrib = 0; truthContrib < barcodeEnergyPair.size(); truthContrib++ ) {
	    //loop over truth particles to find corresponding barcode
	    bool foundMatchingBarcode = false;
	    for( const xAOD::TruthParticle *truthParticle : *truthParticles ) {
	      if( barcodeEnergyPair.at(truthContrib).first != truthParticle->barcode() ) continue;
	      foundMatchingBarcode = true;
	      truthIDs.push_back( truthParticle->pdgId() );
	      truthBarcodes.push_back( truthParticle->barcode() );
	      truthEnergies.push_back( barcodeEnergyPair.at(truthContrib).second * m_conversionFactor );

	      ANA_MSG_VERBOSE( "[electron-linked neutral FE]   PDGID: " << truthParticle->pdgId() << "   Barcode: " << truthParticle->barcode() << "   energy deposited: " << barcodeEnergyPair.at(truthContrib).second * m_conversionFactor << " GeV" );
	      break; //found matching truth particle, move onto the next calpfo contribution
	    }
	    if( !foundMatchingBarcode ) ANA_MSG_INFO( "Huh..? Didn't find any matching barcodes... This is unexpected." );
	  } //end loop over calpfo vector
	  calpfo_NLeadingTruthParticlePdgId.push_back(truthIDs);
	  calpfo_NLeadingTruthParticleBarcode.push_back(truthBarcodes);
	  calpfo_NLeadingTruthParticleEnergy.push_back(truthEnergies);
	}
      }

      for( ElementLink<xAOD::FlowElementContainer> feLink : electronCFELinks.at(countElectron) ) {
	if( feLink.isValid() ) {
	  // OLD COMMENT: this corresponds to charged FE in the JetETMiss container
	  // MODELV2: this now points to Global
	  //const xAOD::FlowElement *electronChargedOriginalFlowElement = *feLink;
	  const xAOD::FlowElement *electronChargedGlobalFlowElement = *feLink;

	  // MODELV2: add the block below
	  // get the index and energy according to FE in the charged Global container
	  chargedPFOindex.push_back(electronChargedGlobalFlowElement->index());
	  chargedPFOenergy.push_back(electronChargedGlobalFlowElement->e() * m_conversionFactor);

	  // printouts
	  ANA_MSG_VERBOSE( "Electron number " << countElectron << ", charged Global FE link has index: " << electronChargedGlobalFlowElement->index() << " and energy: " << electronChargedGlobalFlowElement->e() * m_conversionFactor << " GeV." );

	  // // OLD COMMENT: fetch link from charged JetETMiss -> [charged/neutral] Global
	  // // MODELV2: no longer needed
	  // ElementLink<xAOD::FlowElementContainer> feLinkToGlobal = chargedGlobalFEReadDecorHandle(*electronChargedOriginalFlowElement);
	  // if( feLinkToGlobal.isValid() ) {
	  //   ANA_MSG_VERBOSE( "FE Link from charged JetETMiss to charged Global is valid." );

	  //   // this now corresponds to charged FE in the charged Global container
	  //   const xAOD::FlowElement *electronChargedGlobalFlowElement = *feLinkToGlobal;

	  //   // get the index and energy according to FE in the charged Global container
	  //   chargedPFOindex.push_back(electronChargedGlobalFlowElement->index());
	  //   chargedPFOenergy.push_back(electronChargedGlobalFlowElement->e() * m_conversionFactor);

	  //   // printouts
	  //   ANA_MSG_VERBOSE( "Electron number " << countElectron << ", charged Global FE link has index: " << electronChargedGlobalFlowElement->index() );
	  //   ANA_MSG_VERBOSE( "Compare the index and energy of CHARGED original vs global FE linked to electron... Original: index " << electronChargedOriginalFlowElement->index() << " with energy " << electronChargedOriginalFlowElement->e() * m_conversionFactor << " GeV ... Global: index " << electronChargedGlobalFlowElement->index() << " with energy " << electronChargedGlobalFlowElement->e() * m_conversionFactor << " GeV." );

	  // } else {
	  //   // if element link from charged JetETMiss to charged Global is NOT valid, check if there are links from charged JetETMiss to NEUTRAL Global instead
	  //   // loop over all the neutral Global FEs that the charged JetETMiss FE is linked to
	  //   std::vector<ElementLink<xAOD::FlowElementContainer>> feLinkChargedOriginalToNeutralGlobal = chargedOriginalToNeutralGlobalFEReadDecorHandle(*electronChargedOriginalFlowElement);
	  //   for( ElementLink<xAOD::FlowElementContainer> feLinkToNeutralGlobal : feLinkChargedOriginalToNeutralGlobal ) {
	  //     if( feLinkToNeutralGlobal.isValid() ) { 
	  // 	ANA_MSG_VERBOSE( "FE Link to charged Global is NOT valid... link points to FEs in neutral Global instead." );

	  // 	// this now corresponds to neutral FE in the Global container
	  // 	const xAOD::FlowElement *electronNeutralGlobalFlowElement = *feLinkToNeutralGlobal;

	  // 	// get the index and energy according to FE in the Global container
	  // 	neutralPFOindex.push_back(electronNeutralGlobalFlowElement->index());
	  // 	neutralPFOenergy.push_back(electronNeutralGlobalFlowElement->e() * m_conversionFactor);

	  // 	// save calibration hit information for neutrals: what particles deposited energy (and how much)
	  // 	// note: this is a vector because it's possible for multiple particles to have contributed to this PFO (default: 3 largest contributions are in the calpfo decoration)
	  // 	// get the corresponding neutral PFO from the FlowElement (should be the same index)
	  // 	ANA_MSG_VERBOSE( "Fetching calpfo vector for the neutral FE linked to this electron..." );
	  // 	SG::AuxElement::ConstAccessor< std::vector<std::pair<unsigned int,double>> > calpfoVec("calpfo_NLeadingTruthParticleBarcodeEnergyPairs");
	  // 	const xAOD::FlowElement *electronNPFO = inGlobalNeutralFEHandle->at(electronNeutralGlobalFlowElement->index()); // fetch the info from the Global container
	  // 	std::vector<std::pair<unsigned int,double>> barcodeEnergyPair = calpfoVec(*electronNPFO);

	  // 	std::vector<Int_t> truthIDs;
	  // 	std::vector<Int_t> truthBarcodes;
	  // 	std::vector<Double_t> truthEnergies; //note: this is the amount of energy DEPOSITED in the neutral FE by the truth particle (not the energy of the particle itself)
	  // 	ANA_MSG_VERBOSE( "Got the calpfo vector! Here are its details:" );
	  // 	for( Size_t truthContrib = 0; truthContrib < barcodeEnergyPair.size(); truthContrib++ ) {
	  // 	  //loop over truth particles to find corresponding barcode
	  // 	  bool foundMatchingBarcode = false;
	  // 	  for( const xAOD::TruthParticle *truthParticle : *truthParticles ) {
	  // 	    if( barcodeEnergyPair.at(truthContrib).first != truthParticle->barcode() ) continue;
	  // 	    foundMatchingBarcode = true;
	  // 	    truthIDs.push_back( truthParticle->pdgId() );
	  // 	    truthBarcodes.push_back( truthParticle->barcode() );
	  // 	    truthEnergies.push_back( barcodeEnergyPair.at(truthContrib).second * m_conversionFactor );

	  // 	    ANA_MSG_VERBOSE( "[electron-linked neutral FE]   PDGID: " << truthParticle->pdgId() << "   Barcode: " << truthParticle->barcode() << "   energy deposited: " << barcodeEnergyPair.at(truthContrib).second * m_conversionFactor << " GeV" );
	  // 	    break; //found matching truth particle, move onto the next calpfo contribution
	  // 	  }
	  // 	  if( !foundMatchingBarcode ) ANA_MSG_INFO( "Huh..? Didn't find any matching barcodes... This is unexpected." );
	  // 	} //end loop over calpfo vector
	  // 	calpfo_NLeadingTruthParticlePdgId.push_back(truthIDs);
	  // 	calpfo_NLeadingTruthParticleBarcode.push_back(truthBarcodes);
	  // 	calpfo_NLeadingTruthParticleEnergy.push_back(truthEnergies);
	  //     }

	  //   } // end loop over linked neutral Global FEs
	  // }
	}
      }

      countElectron++;

      m_pt.push_back( electron->pt() * m_conversionFactor );
      m_eta.push_back( electron->eta() );
      m_phi.push_back( electron->phi() );
      m_e.push_back( electron->e() * m_conversionFactor );
      m_passSel.push_back( passSel.isAvailable(*electron) ? passSel(*electron) : -1 );
      m_neutralPFOindex.push_back( neutralPFOindex );
      m_chargedPFOindex.push_back( chargedPFOindex );
      m_neutralPFOenergy.push_back( neutralPFOenergy );
      m_chargedPFOenergy.push_back( chargedPFOenergy );
      m_calpfo_NLeadingTruthParticlePdgId.push_back( calpfo_NLeadingTruthParticlePdgId );
      m_calpfo_NLeadingTruthParticleBarcode.push_back( calpfo_NLeadingTruthParticleBarcode );
      m_calpfo_NLeadingTruthParticleEnergy.push_back( calpfo_NLeadingTruthParticleEnergy );

      // // stop once we have 2 electrons
      // if( m_pt.size() == 2 ) break;
    }

    // fill tree
    m_electronTree->Fill();

    // clear vectors before going to next event
    m_pt.clear();
    m_eta.clear();
    m_phi.clear();
    m_e.clear();
    m_passSel.clear();
    m_neutralPFOindex.clear();
    m_chargedPFOindex.clear();
    m_neutralPFOenergy.clear();
    m_chargedPFOenergy.clear();
    m_calpfo_NLeadingTruthParticlePdgId.clear();
    m_calpfo_NLeadingTruthParticleBarcode.clear();
    m_calpfo_NLeadingTruthParticleEnergy.clear();
  }

  /// PHOTONS
  if( !m_photonContainerName.empty() ) {
    // tell the tree to go to the file
    m_photonTree->SetDirectory( m_file->GetDirectory( m_outFileName.c_str() ) );

    const xAOD::PhotonContainer *inContainer = nullptr;
    ANA_CHECK( evtStore()->retrieve( inContainer, m_photonContainerName ) );

    std::vector<std::vector<ElementLink<xAOD::FlowElementContainer>>> photonNFELinks; //EACH photon will have a VECTOR of FELinks

    Int_t countPhoton = 0;
    for( const xAOD::Photon *photon : *inContainer ) {
      // accessor for checking whether photon passed identification selections
      SG::AuxElement::ConstAccessor< char > passSel("passSel");

      // prepare vector to store PFO indices and energies
      std::vector<int> neutralPFOindex;
      std::vector<double> neutralPFOenergy;

      SG::ReadDecorHandle<xAOD::PhotonContainer, std::vector<ElementLink<xAOD::FlowElementContainer>>> neutralFEReadDecorHandle(m_photonNeutralFEReadDecorKey);
      photonNFELinks.push_back( neutralFEReadDecorHandle(*photon) );

      // printouts for verbose message level
      for( ElementLink<xAOD::FlowElementContainer> feLink : photonNFELinks.at(countPhoton) ) {
	if( feLink.isValid() ) {
	  const xAOD::FlowElement *photonNFE = *feLink;
	  neutralPFOindex.push_back(photonNFE->index());
	  neutralPFOenergy.push_back(photonNFE->e() * m_conversionFactor);
	  ANA_MSG_VERBOSE( "Photon number " << countPhoton << ", link has index: " << photonNFE->index() );
	}
      }

      countPhoton++;

      m_pt.push_back( photon->pt() * m_conversionFactor );
      m_eta.push_back( photon->eta() );
      m_phi.push_back( photon->phi() );
      m_e.push_back( photon->e() * m_conversionFactor );
      m_passSel.push_back( passSel.isAvailable(*photon) ? passSel(*photon) : -1 );
      m_neutralPFOindex.push_back( neutralPFOindex );
      m_neutralPFOenergy.push_back( neutralPFOenergy );

      // // stop after 1 photon
      // if( m_pt.size() == 1 ) break;
    }

    // fill tree
    m_photonTree->Fill();

    // clear vectors before going to next event
    m_pt.clear();
    m_eta.clear();
    m_phi.clear();
    m_e.clear();
    m_passSel.clear();
    m_neutralPFOindex.clear();
    m_neutralPFOenergy.clear();
  }

  /// MUONS
  if( !m_muonContainerName.empty() ) {
    // tell the tree to go to the file
    m_muonTree->SetDirectory( m_file->GetDirectory( m_outFileName.c_str() ) );

    const xAOD::MuonContainer *inContainer = nullptr;
    ANA_CHECK( evtStore()->retrieve( inContainer, m_muonContainerName ) );

    std::vector<std::vector<ElementLink<xAOD::FlowElementContainer>>> muonNFELinks; //EACH muon will have a VECTOR of FELinks
    std::vector<std::vector<ElementLink<xAOD::FlowElementContainer>>> muonCFELinks;

    Int_t countMuon = 0;
    for( const xAOD::Muon *muon : *inContainer ) {
      // accessor for checking whether muon passed identification selections
      SG::AuxElement::ConstAccessor< char > passSel("passSel");

      // prepare vectors to store PFO indices and energies
      std::vector<int> neutralPFOindex;
      std::vector<int> chargedPFOindex;
      std::vector<double> neutralPFOenergy;
      std::vector<double> chargedPFOenergy;

      // links to neutral FE
      SG::ReadDecorHandle<xAOD::MuonContainer, std::vector<ElementLink<xAOD::FlowElementContainer>>> neutralFEReadDecorHandle(m_muonNeutralFEReadDecorKey);
      muonNFELinks.push_back( neutralFEReadDecorHandle(*muon) );

      // links to charged FE
      SG::ReadDecorHandle<xAOD::MuonContainer, std::vector<ElementLink<xAOD::FlowElementContainer>>> chargedFEReadDecorHandle(m_muonChargedFEReadDecorKey);
      muonCFELinks.push_back( chargedFEReadDecorHandle(*muon) );

      // printouts for verbose message level
      for( ElementLink<xAOD::FlowElementContainer> feLink : muonNFELinks.at(countMuon) ) {
	if( feLink.isValid() ) {
	  const xAOD::FlowElement *muonNFE = *feLink;
	  neutralPFOindex.push_back(muonNFE->index());
	  neutralPFOenergy.push_back(muonNFE->e() * m_conversionFactor);
	  ANA_MSG_VERBOSE( "Muon number " << countMuon << ", neutral FE link has index: " << muonNFE->index() );
	}
      }

      for( ElementLink<xAOD::FlowElementContainer> feLink : muonCFELinks.at(countMuon) ) {
	if( feLink.isValid() ) {
	  const xAOD::FlowElement *muonCFE = *feLink;
	  chargedPFOindex.push_back(muonCFE->index());
	  chargedPFOenergy.push_back(muonCFE->e() * m_conversionFactor);
	  ANA_MSG_VERBOSE( "Muon number " << countMuon << ", charged FE link has index: " << muonCFE->index() );
	}
      }

      countMuon++;

      m_pt.push_back( muon->pt() * m_conversionFactor );
      m_eta.push_back( muon->eta() );
      m_phi.push_back( muon->phi() );
      m_e.push_back( muon->e() * m_conversionFactor );
      m_passSel.push_back( passSel.isAvailable(*muon) ? passSel(*muon) : -1 );
      m_neutralPFOindex.push_back( neutralPFOindex );
      m_chargedPFOindex.push_back( chargedPFOindex );
      m_neutralPFOenergy.push_back( neutralPFOenergy );
      m_chargedPFOenergy.push_back( chargedPFOenergy );

      // // stop once we have 2 muons
      // if( m_pt.size() == 2 ) break;
    }

    // fill tree
    m_muonTree->Fill();

    // clear vectors before going to next event
    m_pt.clear();
    m_eta.clear();
    m_phi.clear();
    m_e.clear();
    m_passSel.clear();
    m_neutralPFOindex.clear();
    m_chargedPFOindex.clear();
    m_neutralPFOenergy.clear();
    m_chargedPFOenergy.clear();
  }

  /// TRUTH PARTICLES for Z+jets samples; expect 2 relevant truth electrons / muons
  if( !m_truthContainerName.empty() && (m_isZee || m_isZmumu) ) {
    // set PDG ID of truth lepton
    if( m_isZee ) m_truthPDGID = 11;
    else m_truthPDGID = 13;

    // tell the tree to go to the file
    m_truthTree->SetDirectory( m_file->GetDirectory( m_outFileName.c_str() ) );

    const xAOD::TruthParticleContainer *inContainer = nullptr;
    ANA_CHECK( evtStore()->retrieve( inContainer, m_truthContainerName ) );

    // need to be careful: these aren't sorted by highest pT first
    // make sure we store leading particle as the 0th element in the vector
    int truthCount = 0;
    double leadingTruthPt = 0;
    double leadingTruthEta = 0;
    double leadingTruthPhi = 0;
    double leadingTruthE = 0;
    int leadingTruthID = 0;
    int leadingTruthBarcode = 0;
    double subleadingTruthPt = 0;
    double subleadingTruthEta = 0;
    double subleadingTruthPhi = 0;
    double subleadingTruthE = 0;
    int subleadingTruthID = 0;
    int subleadingTruthBarcode = 0;
    for( const xAOD::TruthParticle* truth : *inContainer ) {
      // skip if it's not the particle we're looking for
      if( std::abs( truth->pdgId() ) != m_truthPDGID ) continue;

  	truthCount++;

  	if( truthCount == 1 ) {
  	  leadingTruthPt = truth->pt() * m_conversionFactor;
  	  leadingTruthEta = truth->eta();
  	  leadingTruthPhi = truth->phi();
  	  leadingTruthE = truth->e() * m_conversionFactor;
  	  leadingTruthID = truth->pdgId();
	  leadingTruthBarcode = truth->barcode();
  	} else if( truthCount > 1 ) {
  	  if( truth->pt() * m_conversionFactor > leadingTruthPt ) {
  	    // set previous truth as subleading, and current truth as leading
  	    subleadingTruthPt = leadingTruthPt;
  	    subleadingTruthEta = leadingTruthEta;
  	    subleadingTruthPhi = leadingTruthPhi;
  	    subleadingTruthE = leadingTruthE;
  	    subleadingTruthID = leadingTruthID;
	    subleadingTruthBarcode = leadingTruthBarcode;
  	    leadingTruthPt = truth->pt() * m_conversionFactor;
  	    leadingTruthEta = truth->eta();
  	    leadingTruthPhi = truth->phi();
  	    leadingTruthE = truth->e() * m_conversionFactor;
  	    leadingTruthID = truth->pdgId();
	    leadingTruthBarcode = truth->barcode();
  	  } else {
  	    // set current truth as subleading
  	    subleadingTruthPt = truth->pt() * m_conversionFactor;
  	    subleadingTruthEta = truth->eta();
  	    subleadingTruthPhi = truth->phi();
  	    subleadingTruthE = truth->e() * m_conversionFactor;
  	    subleadingTruthID = truth->pdgId();
	    subleadingTruthBarcode = truth->barcode();
  	  }

  	  // stop iterating once we have two truths
  	  break;
  	} // end if truthCount > 1
    } // end loop over truth container

    m_pt.push_back( leadingTruthPt );
    m_pt.push_back( subleadingTruthPt );
    m_eta.push_back( leadingTruthEta );
    m_eta.push_back( subleadingTruthEta );
    m_phi.push_back( leadingTruthPhi );
    m_phi.push_back( subleadingTruthPhi );
    m_e.push_back( leadingTruthE );
    m_e.push_back( subleadingTruthE );
    m_truthID.push_back( leadingTruthID );
    m_truthID.push_back( subleadingTruthID );
    m_truthBarcode.push_back( leadingTruthBarcode );
    m_truthBarcode.push_back( subleadingTruthBarcode );

    // fill tree
    m_truthTree->Fill();

    // clear vectors before going to next event
    m_pt.clear();
    m_eta.clear();
    m_phi.clear();
    m_e.clear();
    m_truthID.clear();
    m_truthBarcode.clear();
  }

  /// TRUTH PARTICLES for SinglePhoton samples; expect 1 truth photon per event
  if( !m_truthContainerName.empty() && m_isSinglePhoton ) {
    // tell the tree to go to the file
    m_truthTree->SetDirectory( m_file->GetDirectory( m_outFileName.c_str() ) );

    const xAOD::TruthParticleContainer *inContainer = nullptr;
    ANA_CHECK( evtStore()->retrieve( inContainer, m_truthContainerName ) );

    m_truthPDGID = 22;
    double leadingTruthPt = 0;
    double leadingTruthEta = 0;
    double leadingTruthPhi = 0;
    double leadingTruthE = 0;
    int leadingTruthID = 0;
    for( const xAOD::TruthParticle* truth : *inContainer ) {
      // skip if it's not the particle we're looking for
      if( std::abs( truth->pdgId() ) != m_truthPDGID ) continue;

      leadingTruthPt = truth->pt() * m_conversionFactor;
      leadingTruthEta = truth->eta();
      leadingTruthPhi = truth->phi();
      leadingTruthE = truth->e() * m_conversionFactor;
      leadingTruthID = truth->pdgId();

      // stop iterating once we have a photon
      break;
    } // end loop over truth container

    m_pt.push_back( leadingTruthPt );
    m_eta.push_back( leadingTruthEta );
    m_phi.push_back( leadingTruthPhi );
    m_e.push_back( leadingTruthE );
    m_truthID.push_back( leadingTruthID );

    // fill tree
    m_truthTree->Fill();

    // clear vectors before going to next event
    m_pt.clear();
    m_eta.clear();
    m_phi.clear();
    m_e.clear();
    m_truthID.clear();
  }

  /// Truth jets
  if( !m_truthJetContainerName.empty() ) {
    // tell the tree to go to the file
    m_truthJetTree->SetDirectory( m_file->GetDirectory( m_outFileName.c_str() ) );

    const xAOD::JetContainer *inContainer = nullptr;
    ANA_CHECK( evtStore()->retrieve( inContainer, m_truthJetContainerName ) );

    for( const xAOD::Jet *jet : *inContainer ) {
      m_pt.push_back( jet->pt() * m_conversionFactor );
      m_eta.push_back( jet->eta() );
      m_phi.push_back( jet->phi() );
      m_e.push_back( jet->e() * m_conversionFactor );

      // check if there are any truth jets with pT > 7 GeV within DeltaR < 1 of the current truth jet
      bool foundClosebyJet = false;
      for( const xAOD::Jet *jet2 : *inContainer ) {
	// obviously, don't compare jet with itself
	if( jet == jet2 )
	  continue;

	// only check against jets with pT > 7 GeV
	if( jet2->pt() < 7000 )
	  continue;

	// now check deltaR between jets
	double jetPhiDiff = jet->phi() - jet2->phi();
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

  /// Truth WZ jets
  if( !m_truthWZJetContainerName.empty() ) {
    // tell the tree to go to the file
    m_truthWZJetTree->SetDirectory( m_file->GetDirectory( m_outFileName.c_str() ) );

    const xAOD::JetContainer *inContainer = nullptr;
    ANA_CHECK( evtStore()->retrieve( inContainer, m_truthWZJetContainerName ) );

    for( const xAOD::Jet *jet : *inContainer ) {
      m_pt.push_back( jet->pt() * m_conversionFactor );
      m_eta.push_back( jet->eta() );
      m_phi.push_back( jet->phi() );
      m_e.push_back( jet->e() * m_conversionFactor );

      // check if there are any truth jets with pT > 7 GeV within DeltaR < 1 of the current truth jet
      bool foundClosebyJet = false;
      for( const xAOD::Jet *jet2 : *inContainer ) {
	// obviously, don't compare jet with itself
	if( jet == jet2 )
	  continue;

	// only check against jets with pT > 7 GeV
	if( jet2->pt() < 7000 )
	  continue;

	// now check deltaR between jets
	double jetPhiDiff = jet->phi() - jet2->phi();
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

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode compareFELinks :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.

  ANA_MSG_DEBUG("Calling postExecute");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode compareFELinks :: finalize ()
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



EL::StatusCode compareFELinks :: histFinalize ()
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

  // // The below is only for the "regular" ROOT method of writing to a file (does not work for grid runs)
  // // important: need to cd() to output file otherwise nothing will be written
  // m_file->cd();

  // // write any trees
  // if( !m_inJetContainerName.empty() ) m_jetTree->Write();
  // if( !m_inElectronContainerName.empty() ) m_electronTree->Write();
  // if( !m_inPhotonContainerName.empty() ) m_photonTree->Write();
  // if( !m_inMuonContainerName.empty() ) m_muonTree->Write();
  // //  if( !m_truthContainerName.empty() ) m_truthTree->Write();

  // // done! close file and then we can exit
  // m_file->Close();

  return EL::StatusCode::SUCCESS;
}
