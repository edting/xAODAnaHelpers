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

  // initialise decoration keys for links between objects and FlowElements
  ANA_CHECK(m_electronNeutralFEReadDecorKey.assign("Electrons.neutralGlobalFELinks"));
  ANA_CHECK(m_electronChargedFEReadDecorKey.assign("Electrons.chargedGlobalFELinks"));
  ANA_CHECK(m_electronNeutralFEReadDecorKey.initialize());
  ANA_CHECK(m_electronChargedFEReadDecorKey.initialize());

  // initialise containers
  ANA_CHECK(m_jetContKey.assign("AntiKt4EMPFlowJets"));
  ANA_CHECK(m_eleContKey.assign("Electrons"));
  ANA_CHECK(m_truthParticleContKey.assign("TruthParticles"));
  ANA_CHECK(m_jetContKey.initialize());
  ANA_CHECK(m_eleContKey.initialize());
  ANA_CHECK(m_truthParticleContKey.initialize());

  // initialise jet constituent container(s)
  std::string inputContainerBaseGlobal = "Global";
  std::string inputContainerBaseCHSG = "CHSG";
  std::string inputContainerBaseCSSKG = "CSSKG";
  std::string chargedContainerName = "ChargedParticleFlowObjects";
  std::string neutralContainerName = "NeutralParticleFlowObjects";
  m_inGlobalChargedFEKey = inputContainerBaseGlobal + chargedContainerName;
  m_inGlobalNeutralFEKey = inputContainerBaseGlobal + neutralContainerName;
  m_inCHSGChargedFEKey = inputContainerBaseCHSG + chargedContainerName;
  m_inCHSGNeutralFEKey = inputContainerBaseCHSG + neutralContainerName;
  m_inCSSKGChargedFEKey = inputContainerBaseCSSKG + chargedContainerName;
  m_inCSSKGNeutralFEKey = inputContainerBaseCSSKG + neutralContainerName;
  ANA_CHECK(m_inGlobalChargedFEKey.initialize());
  ANA_CHECK(m_inGlobalNeutralFEKey.initialize());
  ANA_CHECK(m_inCHSGChargedFEKey.initialize());
  ANA_CHECK(m_inCHSGNeutralFEKey.initialize());
  ANA_CHECK(m_inCSSKGChargedFEKey.initialize());
  ANA_CHECK(m_inCSSKGNeutralFEKey.initialize());

  // more decorations
  ANA_CHECK(m_neutralFECellsRemovedReadDecorKey.assign("GlobalNeutralParticleFlowObjects.cellsRemovedFromNeutralFE_barcodeEnergyPair"));
  ANA_CHECK(m_neutralFECellsRemovedReadDecorKey.initialize());

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode mySimpleAlg :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.


  SG::ReadHandle<xAOD::FlowElementContainer> inGlobalNeutralFEHandle = makeHandle(m_inGlobalNeutralFEKey);
  SG::ReadHandle<xAOD::FlowElementContainer> inGlobalChargedFEHandle = makeHandle(m_inGlobalChargedFEKey);
  SG::ReadHandle<xAOD::FlowElementContainer> inCHSGNeutralFEHandle = makeHandle(m_inCHSGNeutralFEKey);
  SG::ReadHandle<xAOD::FlowElementContainer> inCHSGChargedFEHandle = makeHandle(m_inCHSGChargedFEKey);
  SG::ReadHandle<xAOD::FlowElementContainer> inCSSKGNeutralFEHandle = makeHandle(m_inCSSKGNeutralFEKey);
  SG::ReadHandle<xAOD::FlowElementContainer> inCSSKGChargedFEHandle = makeHandle(m_inCSSKGChargedFEKey);
  SG::ReadHandle<xAOD::ElectronContainer> inContainer(m_eleContKey);
  SG::ReadHandle<xAOD::TruthParticleContainer> truthParticles(m_truthParticleContKey);
  SG::ReadHandle<xAOD::JetContainer> ufoJets("AntiKt10UFOCSSKJets");
  SG::ReadHandle<xAOD::JetContainer> pflowJets("AntiKt4EMPFlowJets");

  // this block tests fetching UFO jet constituents and the otherObjects() and chargedObjects() of those constituents
  if( false ) {
    for( const xAOD::Jet *jet : *ufoJets ) {
      for( size_t consti = 0; consti < jet->numConstituents(); consti++) {
	const xAOD::FlowElement *constit = static_cast<const xAOD::FlowElement*>(jet->rawConstituent(consti)); //assumption: this returns constituents in UFOCSSK
	//ANA_MSG_INFO( "Trying to access UFO constituent energy... " << constit->e() << " Did it work? If so, rawConstituent is valid." );

	if( constit->isCharged() ) {
	  std::vector<const xAOD::IParticle*> chargedObjects = constit->chargedObjects(); //assumption: this returns charged objects in CSSKGChargedParticleFlowObjects
	  Int_t chargedObjSize = chargedObjects.size();
	  //ANA_MSG_INFO( "Size of UFO jet constituent chargedObjects() vector: " << chargedObjSize );
	
	  if( chargedObjSize > 0 && chargedObjects.at(0) ) ANA_MSG_INFO( "Valid charged object..." );
	  else ANA_MSG_INFO( "Charged object is a nullptr." );
	}

	if( !constit->isCharged() ) {
	  std::vector<const xAOD::IParticle*> otherObjects = constit->otherObjects(); //assumption: this returns other objects in CSSKGNeutralParticleFlowObjects
	  Int_t otherObjSize = otherObjects.size();
	  //ANA_MSG_INFO( "Size of UFO jet constituent otherObjects() vector: " << otherObjSize );
	
	  if( otherObjSize > 0 && otherObjects.at(0) ) ANA_MSG_INFO( "Valid other object..." );
	  else ANA_MSG_INFO( "Other object is a nullptr." );
	}
      }
    }
  }

  // this block tests the isCharged() accessor for FEs
  if( false ) {
    for( const xAOD::FlowElement *fe : *inCHSGNeutralFEHandle ) {
      const xAOD::FlowElement *fe_orig = static_cast<const xAOD::FlowElement*>(getOriginalObject(*fe));
      if( fe_orig->isCharged() ) ANA_MSG_INFO( "This FE is charged!" );
      else ANA_MSG_INFO( "This FE is neutral!" );
    }
  }

  // this block checks the size of the calpfo vector
  if( true ) {
    for( const xAOD::FlowElement *fe : *inGlobalNeutralFEHandle ) {
      SG::AuxElement::ConstAccessor< std::vector<std::pair<unsigned int,double>> > calpfoVec("calpfo_100LeadingTruthParticleBarcodeEnergyPairs");
      std::vector<std::pair<unsigned int,double>> barcodeEnergyPair = calpfoVec(*fe);
      ANA_MSG_INFO( "Got the calpfo vector! It has size: " << barcodeEnergyPair.size() );
    }
  }

  // this block tests the calibration hits of neutral FEs linked to electrons
  if( false ) {
    Int_t countElectron = 0;
    for( const xAOD::Electron *electron : *inContainer ) {

      // decorations
      SG::ReadDecorHandle<xAOD::ElectronContainer, std::vector<ElementLink<xAOD::FlowElementContainer>>> neutralFEReadDecorHandle(m_electronNeutralFEReadDecorKey);
      SG::ReadDecorHandle<xAOD::ElectronContainer, std::vector<ElementLink<xAOD::FlowElementContainer>>> chargedFEReadDecorHandle(m_electronChargedFEReadDecorKey);
      std::vector<ElementLink<xAOD::FlowElementContainer>> electronNFELinks = neutralFEReadDecorHandle(*electron);
      std::vector<ElementLink<xAOD::FlowElementContainer>> electronCFELinks = chargedFEReadDecorHandle(*electron);

      // vectors to store calibration hit information for each neutral FE
      // by default, (up to) the three largest contributions is saved to the calpfo vector
      std::vector<std::vector<int>> calpfo_NLeadingTruthParticlePdgId;
      std::vector<std::vector<int>> calpfo_NLeadingTruthParticleBarcode;
      std::vector<std::vector<double>> calpfo_NLeadingTruthParticleEnergy;

      // // also prepare vectors for CaloCalTopoCluster calibration hits
      // std::vector<std::vector<int>> calclus_NLeadingTruthParticlePdgId;
      // std::vector<std::vector<int>> calclus_NLeadingTruthParticleBarcode;
      // std::vector<std::vector<double>> calclus_NLeadingTruthParticleEnergy;

      // loop over each neutral FE linked to the electron
      for( ElementLink<xAOD::FlowElementContainer> feLink : electronNFELinks ) {
	if( feLink.isValid() ) {
	  const xAOD::FlowElement *electronNeutralGlobalFlowElement = *feLink;

	  ANA_MSG_INFO( "Electron number " << countElectron << " is linked to a neutral FE with index: " << electronNeutralGlobalFlowElement->index() << " with energy: " << electronNeutralGlobalFlowElement->e() * 0.001 << " GeV." );

	  // cells removed from the neutral PFO
	  SG::ReadDecorHandle<xAOD::FlowElementContainer, std::vector<std::pair<unsigned int, double>>> neutralFECellsRemovedReadDecorHandle(m_neutralFECellsRemovedReadDecorKey);
	  std::vector<std::pair<unsigned int, double>> caloCellsRemoved = neutralFECellsRemovedReadDecorHandle(*electronNeutralGlobalFlowElement);
	  ANA_MSG_INFO( "Size of the cellsRemovedFromNeutralFE decoration: " << caloCellsRemoved.size() );

	  // calibration hit information
	  ANA_MSG_INFO( "Fetching calpfo vector for the neutral FE linked to this electron..." );
	  SG::AuxElement::ConstAccessor< std::vector<std::pair<unsigned int,double>> > calpfoVec("calpfo_100LeadingTruthParticleBarcodeEnergyPairs");
	  std::vector<std::pair<unsigned int,double>> barcodeEnergyPair = calpfoVec(*electronNeutralGlobalFlowElement);
	  ANA_MSG_INFO( "Got the calpfo vector! Here are its details:" );

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
	      truthEnergies.push_back( barcodeEnergyPair.at(truthContrib).second * 0.001 );

	      ANA_MSG_INFO( "[electron-linked neutral FE]   PDGID: " << truthParticle->pdgId() << "   Barcode: " << truthParticle->barcode() << "   energy deposited: " << barcodeEnergyPair.at(truthContrib).second * 0.001 << " GeV" );
	      break;
	    }
	    if( !foundMatchingBarcode ) ANA_MSG_INFO( "Huh..? Didn't find any matching barcodes... This is unexpected." );
	  } //end loop over calpfo vector
	  calpfo_NLeadingTruthParticlePdgId.push_back(truthIDs);
	  calpfo_NLeadingTruthParticleBarcode.push_back(truthBarcodes);
	  calpfo_NLeadingTruthParticleEnergy.push_back(truthEnergies);

	  // // save calibration hit information for the topo-cluster corresponding to this neutral FE (if decoration is available)
	  // if( m_calclusIsAvail ) {
	  //   xAOD::CaloCluster *linkedCluster = (xAOD::CaloCluster*) electronNeutralGlobalFlowElement->otherObjects().at(0);

	  //   ANA_MSG_INFO( "Fetching calclus vector for the topo-cluster linked to the neutral FE that is linked to this electron..." );
	  //   SG::AuxElement::ConstAccessor< std::vector<std::pair<unsigned int,double>> > calpfoVec("calclus_NLeadingTruthParticleBarcodeEnergyPairs");
	  //   barcodeEnergyPair = calpfoVec(*linkedCluster);
	  //   ANA_MSG_INFO( "Got the calclus vector! Here are its details:" );

	  //   // find truth particle with matching barcode and check its PDG ID
	  //   truthIDs.clear();
	  //   truthBarcodes.clear();
	  //   truthEnergies.clear();
	  //   for( Size_t truthContrib = 0; truthContrib < barcodeEnergyPair.size(); truthContrib++ ) {
	  //     bool foundMatchingBarcode = false;
	  //     for( const xAOD::TruthParticle *truthParticle : *truthParticles ) {
	  //       if( barcodeEnergyPair.at(truthContrib).first != truthParticle->barcode() ) continue;
	  //       foundMatchingBarcode = true;
	  //       truthIDs.push_back( truthParticle->pdgId() );
	  //       truthBarcodes.push_back( truthParticle->barcode() );
	  //       truthEnergies.push_back( barcodeEnergyPair.at(truthContrib).second * 0.001 );

	  //       ANA_MSG_INFO( "[topo-cluster linked to the electron-linked neutral FE]   PDGID: " << truthParticle->pdgId() << "   Barcode: " << truthParticle->barcode() << "   energy deposited: " << barcodeEnergyPair.at(truthContrib).second * 0.001 << " GeV" );
	  //       break;
	  //     }
	  //     if( !foundMatchingBarcode ) ANA_MSG_INFO( "Huh..? Didn't find any matching barcodes... This is unexpected." );
	  //   } //end loop over calpfo vector
	  //   calclus_NLeadingTruthParticlePdgId.push_back(truthIDs);
	  //   calclus_NLeadingTruthParticleBarcode.push_back(truthBarcodes);
	  //   calclus_NLeadingTruthParticleEnergy.push_back(truthEnergies);
	  // }
	}
      } //end loop over linked neutral FEs
    } //end loop over electrons
  }

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
