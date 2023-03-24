/*******************************************************************************
 *
 * This algorithm produces ntuples for the purpose of overlap removal studies.
 * Depending on the object containers that are specified in the configuration,
 * the event topology is assumed, as follows:
 *   - only electrons .............. Z->ee
 *   - only muons .................. Z->mm
 *   - only photons ................ single photon
 *   - electrons and muons ......... ttbar
 *   - electrons, muons, photons ... dijet
 *
 * In all cases, a jet container is expected to be specified.
 *
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
  m_infoTree->Branch("initialSumW", &m_MD_initialSumW, "initialSumW/D"); //weighted number of events

  // PFlow jets
  if( !m_pflowJetContainerName.empty() ) {
    m_pflowJetTree = new TTree(m_pflowJetTreeName.c_str(),"pflowJetTree");
    m_pflowJetTree->Branch("pt", &m_pt);
    m_pflowJetTree->Branch("eta", &m_eta);
    m_pflowJetTree->Branch("phi", &m_phi);
    m_pflowJetTree->Branch("e", &m_e);
    m_pflowJetTree->Branch("isTag", &m_isTag);
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

  // EMTopo jets
  if( !m_topoJetContainerName.empty() ) {
    m_topoJetTree = new TTree(m_topoJetTreeName.c_str(),"topoJetTree");
    m_topoJetTree->Branch("pt", &m_pt);
    m_topoJetTree->Branch("eta", &m_eta);
    m_topoJetTree->Branch("phi", &m_phi);
    m_topoJetTree->Branch("e", &m_e);
    m_topoJetTree->Branch("isTag", &m_isTag);
    m_topoJetTree->Branch("nTrk", &m_nTrk);
    m_topoJetTree->Branch("passJVT", &m_passJVT);
    m_topoJetTree->Branch("detEta", &m_detEta);
  }

  // electrons
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

  // photons
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

  // muons
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
  std::string chargedContainerName = "ChargedParticleFlowObjects";
  std::string neutralContainerName = "NeutralParticleFlowObjects";
  m_inGlobalChargedFEKey = inputContainerBaseGlobal + chargedContainerName;
  m_inGlobalNeutralFEKey = inputContainerBaseGlobal + neutralContainerName;
  m_inCHSGChargedFEKey = inputContainerBaseCHSG + chargedContainerName;
  m_inCHSGNeutralFEKey = inputContainerBaseCHSG + neutralContainerName;
  ANA_CHECK(m_inGlobalChargedFEKey.initialize());
  ANA_CHECK(m_inGlobalNeutralFEKey.initialize());
  ANA_CHECK(m_inCHSGChargedFEKey.initialize());
  ANA_CHECK(m_inCHSGNeutralFEKey.initialize());

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
   *
   *  Selection for Zee/Zmumu/SinglePhoton samples:
   *  Only save event if:
   *   - at least 2 reconstructed electrons / 2 muons / 1 photon are present in event
   *  If a truth object container is also specified:
   *   - check that reconstructed and truth electrons/muons/photons are deltaR and pT-matched
   *
   *  For ttbar events: check for at least 1 truth lepton(s), then look for matching reco lepton(s)
   *
   *  For dijet events: save all electrons / muons / photons
   *
   ******************************************************************/

  // based on what object container strings are non-empty, figure out the sample we're looking at, then apply selections accordingly
  if( !m_electronContainerName.empty() && m_muonContainerName.empty() && m_photonContainerName.empty() ) m_isZee = true;
  else if( m_electronContainerName.empty() && !m_muonContainerName.empty() && m_photonContainerName.empty() ) m_isZmumu = true;
  else if( !m_electronContainerName.empty() && !m_muonContainerName.empty() && m_photonContainerName.empty() ) m_isttbar = true;
  else if( m_electronContainerName.empty() && m_muonContainerName.empty() && !m_photonContainerName.empty() ) m_isSinglePhoton = true;
  else if( !m_electronContainerName.empty() && !m_muonContainerName.empty() && !m_photonContainerName.empty() ) m_isDijet = true;

  // PRESELECTION FOR ZEE EVENTS
  if( m_isZee ) {
    SG::ReadHandle<xAOD::ElectronContainer> inContainer(m_eleContKey);

    // check that event contains at least 2 electrons, otherwise skip event
    // pt, eta and phi of first two electrons will be used for truth matching (if truth container is available)
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

    if( electronCount < 2 ) return EL::StatusCode::SUCCESS;

    if( !m_truthContainerName.empty() ) {
      SG::ReadHandle<xAOD::TruthParticleContainer> truthContainer(m_truthContKey);

      m_truthPDGID = 11; //PDG ID of electron
      int truthCount = 0;
      double leadingTruthPt = 0;
      double leadingTruthEta = 0;
      double leadingTruthPhi = 0;
      double subleadingTruthPt = 0;
      double subleadingTruthEta = 0;
      double subleadingTruthPhi = 0;
      for (const xAOD::TruthParticle* truth : *truthContainer) {
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

      // skip event if either of the reconstructed leading or subleading electrons don't match the truth
      const double leadingTruthRecoEtaDiff = leadingElectronEta - leadingTruthEta;
      const double leadingTruthRecoPhiDiff = leadingElectronPhi - leadingTruthPhi;
      const double leadingTruthRecoDeltaR = sqrt( leadingTruthRecoEtaDiff*leadingTruthRecoEtaDiff + leadingTruthRecoPhiDiff*leadingTruthRecoPhiDiff );
      const double leadingTruthRecoPtDiffRatio = std::abs( leadingElectronPt - leadingTruthPt ) / ( leadingElectronPt + leadingTruthPt );

      const double subleadingTruthRecoEtaDiff = subleadingElectronEta - subleadingTruthEta;
      const double subleadingTruthRecoPhiDiff = subleadingElectronPhi - subleadingTruthPhi;
      const double subleadingTruthRecoDeltaR = sqrt( subleadingTruthRecoEtaDiff*subleadingTruthRecoEtaDiff + subleadingTruthRecoPhiDiff*subleadingTruthRecoPhiDiff );
      const double subleadingTruthRecoPtDiffRatio = std::abs( subleadingElectronPt - subleadingTruthPt ) / ( subleadingElectronPt + subleadingTruthPt );

      if( (leadingTruthRecoDeltaR > m_deltaRcutoff)
	  || (leadingTruthRecoPtDiffRatio > m_ptDiffCutoff)
	  || (subleadingTruthRecoDeltaR > m_deltaRcutoff)
	  || (subleadingTruthRecoPtDiffRatio > m_ptDiffCutoff)
	  ) {
      	return EL::StatusCode::SUCCESS;
      }
    }
  }

  // PRESELECTION FOR SINGLE PHOTON EVENTS
  if( m_isSinglePhoton ) {
    SG::ReadHandle<xAOD::PhotonContainer> inContainer(m_phoContKey);

    // check that event contains at least 1 photon, otherwise skip event
    // pt, eta and phi of first photon will be used for truth matching (if truth container is available)
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

    if( photonCount < 1 ) return EL::StatusCode::SUCCESS;

    if( !m_truthContainerName.empty() ) {
      SG::ReadHandle<xAOD::TruthParticleContainer> truthContainer(m_truthContKey);

      m_truthPDGID = 22; //PDG ID of photon
      int truthCount = 0;
      double leadingTruthPt = 0;
      double leadingTruthEta = 0;
      double leadingTruthPhi = 0;
      for (const xAOD::TruthParticle* truth : *truthContainer) {
  	if( std::abs( truth->pdgId() ) != m_truthPDGID ) continue;

  	truthCount++;
  	if( truthCount == 1 ) {
  	  leadingTruthPt = truth->pt() * m_conversionFactor;
  	  leadingTruthEta = truth->eta();
  	  leadingTruthPhi = truth->phi();

  	  break;
  	}
      } // end loop over truth container

      // skip event if reconstructed and truth photon does not match
      const double leadingTruthRecoEtaDiff = leadingPhotonEta - leadingTruthEta;
      const double leadingTruthRecoPhiDiff = leadingPhotonPhi - leadingTruthPhi;
      const double leadingTruthRecoDeltaR = sqrt( leadingTruthRecoEtaDiff*leadingTruthRecoEtaDiff + leadingTruthRecoPhiDiff*leadingTruthRecoPhiDiff );
      const double leadingTruthRecoPtDiffRatio = std::abs( leadingPhotonPt - leadingTruthPt ) / ( leadingPhotonPt + leadingTruthPt );

      if( (leadingTruthRecoDeltaR > m_deltaRcutoff) || (leadingTruthRecoPtDiffRatio > m_ptDiffCutoff) ) {
  	return EL::StatusCode::SUCCESS;
      }
    }
  }

  // PRESELECTION FOR ZMUMU EVENTS
  if( m_isZmumu ) {
    SG::ReadHandle<xAOD::MuonContainer> inContainer(m_muContKey);

    // check that event contains at least 2 muons, otherwise skip event
    // pt, eta and phi of first two muons will be used for truth matching below (if truth container is available)
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

    if( muonCount < 2 ) return EL::StatusCode::SUCCESS;

    if( !m_truthContainerName.empty() ) {
      SG::ReadHandle<xAOD::TruthParticleContainer> truthContainer(m_truthContKey);

      m_truthPDGID = 13; //PDG ID of muon
      int truthCount = 0;
      double leadingTruthPt = 0;
      double leadingTruthEta = 0;
      double leadingTruthPhi = 0;
      double subleadingTruthPt = 0;
      double subleadingTruthEta = 0;
      double subleadingTruthPhi = 0;
      for (const xAOD::TruthParticle* truth : *truthContainer) {
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

      // skip event if either of the reconstructed leading or subleading muons don't match the truth
      const double leadingTruthRecoEtaDiff = leadingMuonEta - leadingTruthEta;
      const double leadingTruthRecoPhiDiff = leadingMuonPhi - leadingTruthPhi;
      const double leadingTruthRecoDeltaR = sqrt( leadingTruthRecoEtaDiff*leadingTruthRecoEtaDiff + leadingTruthRecoPhiDiff*leadingTruthRecoPhiDiff );
      const double leadingTruthRecoPtDiffRatio = std::abs( leadingMuonPt - leadingTruthPt ) / ( leadingMuonPt + leadingTruthPt );

      const double subleadingTruthRecoEtaDiff = subleadingMuonEta - subleadingTruthEta;
      const double subleadingTruthRecoPhiDiff = subleadingMuonPhi - subleadingTruthPhi;
      const double subleadingTruthRecoDeltaR = sqrt( subleadingTruthRecoEtaDiff*subleadingTruthRecoEtaDiff + subleadingTruthRecoPhiDiff*subleadingTruthRecoPhiDiff );
      const double subleadingTruthRecoPtDiffRatio = std::abs( subleadingMuonPt - subleadingTruthPt ) / ( subleadingMuonPt + subleadingTruthPt );

      if( (leadingTruthRecoDeltaR > m_deltaRcutoff)
	  || (leadingTruthRecoPtDiffRatio > m_ptDiffCutoff)
	  || (subleadingTruthRecoDeltaR > m_deltaRcutoff)
	  || (subleadingTruthRecoPtDiffRatio > m_ptDiffCutoff)
	  ) {
  	return EL::StatusCode::SUCCESS;
      }
    }
  }

  // PRESELECTION FOR TTBAR EVENTS
  int numElectronsExpected = 0;
  int numMuonsExpected = 0;
  if( m_isttbar ) {
    SG::ReadHandle<xAOD::TruthParticleContainer> truthContainer(m_truthContKey);

    // check for at least one truth lepton from top decay in the event
    bool foundTruthLepton = false;
    std::vector<int> leptonID;
    std::vector<float> truthLeptonPt;
    std::vector<float> truthLeptonEta;
    std::vector<float> truthLeptonPhi;
    std::vector<float> truthLeptonE;

    for( const xAOD::TruthParticle *truth : *truthContainer ) {
      // skip any particles that aren't electrons or muons
      if( std::abs( truth->pdgId() ) != 11 && std::abs( truth->pdgId() ) != 13 ) continue;

      foundTruthLepton = true;
      leptonID.push_back( truth->pdgId() );
      truthLeptonPt.push_back( truth->pt() * m_conversionFactor );
      truthLeptonEta.push_back( truth->eta() );
      truthLeptonPhi.push_back( truth->phi() );
      truthLeptonE.push_back( truth->e() * m_conversionFactor );

      if( leptonID.size() == 2 ) break;
    } //end loop over truth particles

    if( !foundTruthLepton ) return EL::StatusCode::SUCCESS;

    // if a lepton(s) is found, check that corresponding reconstructed containers have at least the same number of leptons
    int numElectronsFound = 0;
    int numMuonsFound = 0;
    std::vector<float> recoLeptonPt;
    std::vector<float> recoLeptonEta;
    std::vector<float> recoLeptonPhi;
    std::vector<float> recoLeptonE;
    for( int lepton = 0; static_cast<std::vector<int>::size_type>( lepton ) < leptonID.size(); lepton++ ) {
      if( std::abs( leptonID.at( lepton ) ) == 11 ) {
	numElectronsExpected++;

	// fetch the corresponding electron from the reconstructed electron container
	SG::ReadHandle<xAOD::ElectronContainer> electrons(m_eleContKey);

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
      } else if( std::abs( leptonID.at( lepton ) ) == 13 ) {
	numMuonsExpected++;

	// fetch the corresponding muon from the reconstructed muon container
	SG::ReadHandle<xAOD::MuonContainer> muons(m_muContKey);

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

    // skip event if any of the reconstructed lepton(s) don't match the truth
    for( int lepton = 0; static_cast<std::vector<int>::size_type>( lepton ) < leptonID.size(); lepton++ ) {
      const float etaDiff = truthLeptonEta.at( lepton ) - recoLeptonEta.at( lepton );
      const float phiDiff = truthLeptonPhi.at( lepton ) - recoLeptonPhi.at( lepton );
      const float truthRecoDeltaR = sqrt( etaDiff*etaDiff + phiDiff*phiDiff );
      const float truthRecoPtDiffRatio = std::abs( truthLeptonPt.at( lepton ) - recoLeptonPt.at( lepton ) ) / ( truthLeptonPt.at( lepton ) + recoLeptonPt.at( lepton ) );

      if( truthRecoDeltaR > m_deltaRcutoff || truthRecoPtDiffRatio > m_ptDiffCutoff )
	return EL::StatusCode::SUCCESS;
      else {
	ANA_MSG_VERBOSE( "Match found! Reco vs truth lepton:" );
	ANA_MSG_VERBOSE( "\tLepton number " << lepton << ": reco pT = " << recoLeptonPt.at( lepton ) << " GeV, eta = " << recoLeptonEta.at( lepton ) << " and phi = " << recoLeptonPhi.at( lepton ) );
	ANA_MSG_VERBOSE( "\tLepton number " << lepton << ": truth pT = " << truthLeptonPt.at( lepton ) << " GeV, eta = " << truthLeptonEta.at( lepton ) << " and phi = " << truthLeptonPhi.at( lepton ) );
      }
    }

    // fill the truth lepton tree in the output ntuple now, since we already did the work to find them
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
   *  END OF PRESELECTION; START FILLING TREES
   ********************************************************/

  // fetch primary vertex information for this event for later reference (needed to find number of jet tracks)
  size_t vtxIdx = 0;
  m_nv = 0;
  m_npv = 0;
  m_npuv = 0;
  const xAOD::Vertex *priVtx = nullptr;
  const xAOD::VertexContainer* vertices = nullptr;
  if(evtStore()->retrieve(vertices, "PrimaryVertices").isSuccess()) {
    for(const xAOD::Vertex* vtx : *vertices) {
      m_nv++;

      if(vtx->vertexType() == xAOD::VxType::PriVtx) {
	if( !priVtx ) priVtx = vtx;
	m_npv++;
      }

      if(vtx->vertexType() == xAOD::VxType::PileUp) {
	m_npuv++;
      }
    }
  } else {
    ANA_MSG_WARNING( "Failed to retrieve primary vertex container" );
  }
  vtxIdx = priVtx->index();

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
	  SG::AuxElement::ConstAccessor< std::vector<std::pair<unsigned int,double>> > calpfoVec("calpfo_NLeadingTruthParticleBarcodeEnergyPairs");
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

      // b-tag and scale factor
      SG::AuxElement::ConstAccessor< char > isTag("BTag_DL1dv01_FixedCutBEff_77"); // actually an int: 0 = not b-tagged, 1 = tagged
      SG::AuxElement::ConstAccessor< std::vector<float> > sf("BTag_SF_DL1dv01_FixedCutBEff_77");
      m_isTag.push_back( isTag.isAvailable(*jet) ? isTag(*jet) : -1 );
      m_sf.push_back( sf.isAvailable(*jet) ? sf(*jet) : junk);
      if( isTag.isAvailable(*jet) && msgLvl(MSG::VERBOSE) ) std::cout << "Jet b-tag flag: " << (int)isTag(*jet) << "\n";

      // number of tracks associated with jet
      static const SG::AuxElement::ConstAccessor< std::vector<int> > acc("NumTrkPt500");
      int numTracks = acc(*jet).at(vtxIdx);
      m_nTrk.push_back( numTracks );
      ANA_MSG_VERBOSE( "Number of jet tracks: " << numTracks );

      // JVT
      static const SG::AuxElement::ConstAccessor< char > passJVT("JetJVT_Passed_Tight");
      m_passJVT.push_back( passJVT.isAvailable(*jet) ? passJVT(*jet) : -1 );
      ANA_MSG_VERBOSE( "Pass JVT cut: " << (int)passJVT(*jet) );

      // constituent-scale quantities
      xAOD::JetFourMom_t constit4vec = jet->getAttribute<xAOD::JetFourMom_t>("JetConstitScaleMomentum");
      m_pt_constit.push_back( constit4vec.Pt() * m_conversionFactor );
      m_eta_constit.push_back( constit4vec.Eta() );
      m_phi_constit.push_back( constit4vec.Phi() );
      m_e_constit.push_back( constit4vec.E() * m_conversionFactor );

      // isolation check: see if there are any other jets with pt > 7 GeV within deltaR < 0.6
      bool foundClosebyJet = false;
      for( const xAOD::Jet *jet2 : *inContainer ) {
	if( jet2 == jet || jet2->pt() < 7000 ) {
	  continue;
	}

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

      // detector eta
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

      // b-tag and scale factor -- actually nothing for EMTopo jets available...
      SG::AuxElement::ConstAccessor< char > isTag("BTag_DL1dv01_FixedCutBEff_77"); // actually an int: 0 = not b-tagged, 1 = tagged
      SG::AuxElement::ConstAccessor< std::vector<float> > sf("BTag_SF_DL1dv01_FixedCutBEff_77");
      m_isTag.push_back( isTag.isAvailable(*jet) ? isTag(*jet) : -1 );
      m_sf.push_back( sf.isAvailable(*jet) ? sf(*jet) : junk);
      if( isTag.isAvailable(*jet) && msgLvl(MSG::VERBOSE) ) std::cout << "Jet b-tag flag: " << (int)isTag(*jet) << "\n";

      // number of tracks associated with jet
      static const SG::AuxElement::ConstAccessor< std::vector<int> > acc("NumTrkPt500");
      int numTracks = acc(*jet).at(vtxIdx);
      m_nTrk.push_back( numTracks );
      ANA_MSG_VERBOSE( "Number of jet tracks: " << numTracks );

      // JVT
      static const SG::AuxElement::ConstAccessor< char > passJVT("JetJVT_Passed_Tight");
      m_passJVT.push_back( passJVT.isAvailable(*jet) ? passJVT(*jet) : -1 );
      ANA_MSG_VERBOSE( "Pass JVT cut: " << (int)passJVT(*jet) );

      // detector eta
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
    SG::ReadHandle<xAOD::ElectronContainer> inContainer(m_eleContKey);
    SG::ReadHandle<xAOD::TruthParticleContainer> truthParticles(m_truthParticlesContKey);

    Int_t countElectron = 0;
    for( const xAOD::Electron *electron : *inContainer ) {
      // identification selection check
      SG::AuxElement::ConstAccessor< char > passSel("passSel");

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
	  SG::AuxElement::ConstAccessor< std::vector<std::pair<unsigned int,double>> > calpfoVec("calpfo_NLeadingTruthParticleBarcodeEnergyPairs");
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
      m_passSel.push_back( passSel.isAvailable(*electron) ? passSel(*electron) : -1 );
      m_neutralPFOindex.push_back( neutralPFOindex );
      m_chargedPFOindex.push_back( chargedPFOindex );
      m_neutralPFOenergy.push_back( neutralPFOenergy );
      m_chargedPFOenergy.push_back( chargedPFOenergy );
      m_calpfo_NLeadingTruthParticlePdgId.push_back( calpfo_NLeadingTruthParticlePdgId );
      m_calpfo_NLeadingTruthParticleBarcode.push_back( calpfo_NLeadingTruthParticleBarcode );
      m_calpfo_NLeadingTruthParticleEnergy.push_back( calpfo_NLeadingTruthParticleEnergy );

      countElectron++;
      if( m_doTruthObjectMatch && ( (m_isZee && countElectron == 2) || (m_isttbar && countElectron == numElectronsExpected) ) ) {
	break;
      }
    } //end loop over electrons

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
    SG::ReadHandle<xAOD::PhotonContainer> inContainer(m_phoContKey);

    Int_t countPhoton = 0;
    for( const xAOD::Photon *photon : *inContainer ) {
      // identification selection check
      SG::AuxElement::ConstAccessor< char > passSel("passSel");

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
      m_passSel.push_back( passSel.isAvailable(*photon) ? passSel(*photon) : -1 );
      m_neutralPFOindex.push_back( neutralPFOindex );
      m_neutralPFOenergy.push_back( neutralPFOenergy );

      countPhoton++;
      if( m_doTruthObjectMatch && m_isSinglePhoton && countPhoton == 1 ) {
	break;
      }
    } //end loop over photons

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
    SG::ReadHandle<xAOD::MuonContainer> inContainer(m_muContKey);

    Int_t countMuon = 0;
    for( const xAOD::Muon *muon : *inContainer ) {
      // identification selection check
      SG::AuxElement::ConstAccessor< char > passSel("passSel");

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
      m_passSel.push_back( passSel.isAvailable(*muon) ? passSel(*muon) : -1 );
      m_neutralPFOindex.push_back( neutralPFOindex );
      m_chargedPFOindex.push_back( chargedPFOindex );
      m_neutralPFOenergy.push_back( neutralPFOenergy );
      m_chargedPFOenergy.push_back( chargedPFOenergy );

      countMuon++;
      if( m_doTruthObjectMatch && ( (m_isZmumu && countMuon == 2) || (m_isttbar && countMuon == numMuonsExpected) ) ) {
	break;
      }
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

  // TRUTH PARTICLES (for Z+jets samples; expect 2 truth electrons/muons)
  // side note: for ttbar topology, we already saved truth particles during the preselection step
  if( !m_truthContainerName.empty() && (m_isZee || m_isZmumu) ) {
    if( m_isZee ) m_truthPDGID = 11;
    else m_truthPDGID = 13;

    SG::ReadHandle<xAOD::TruthParticleContainer> inContainer(m_truthContKey);

    // save pT-ordered truth particles
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

  // TRUTH PARTICLES (for SinglePhoton samples; expect 1 truth photon)
  if( !m_truthContainerName.empty() && m_isSinglePhoton ) {
    SG::ReadHandle<xAOD::TruthParticleContainer> inContainer(m_truthContKey);

    m_truthPDGID = 22;
    double leadingTruthPt = 0;
    double leadingTruthEta = 0;
    double leadingTruthPhi = 0;
    double leadingTruthE = 0;
    int leadingTruthID = 0;
    for( const xAOD::TruthParticle* truth : *inContainer ) {
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

  return EL::StatusCode::SUCCESS;
}
