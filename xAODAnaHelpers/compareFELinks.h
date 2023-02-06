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

#ifndef xAODAnaHelpers_compareFELinks_H
#define xAODAnaHelpers_compareFELinks_H

// algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"

#include "AsgTools/ToolHandle.h"
#include "AsgDataHandles/ReadDecorHandle.h"
#include "AsgDataHandles/WriteDecorHandleKey.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODPFlow/FlowElement.h"
#include "xAODPFlow/FlowElementContainer.h"
#include "xAODPFlow/FlowElementAuxContainer.h"
#include "xAODPFlow/PFO.h"
#include "xAODPFlow/PFOContainer.h"
#include "xAODPFlow/PFOAuxContainer.h"

// ROOT includes
#include "TTree.h"

class compareFELinks : public xAH::Algorithm
{
public:
  // The variables in this public scope can be changed by the user in the xAH config file

  // name of output ntuple root file
  std::string m_outFileName = "compareFELinks.root";

  // names of output trees
  std::string m_infoTreeName = "info";
  std::string m_pflowJetTreeName = "pflowJets";
  std::string m_topoJetTreeName = "topoJets";
  std::string m_electronTreeName = "electrons";
  std::string m_photonTreeName = "photons";
  std::string m_muonTreeName = "muons";
  std::string m_truthTreeName = "truthParticles";
  std::string m_truthJetTreeName = "truthJets";
  std::string m_truthWZJetTreeName = "truthWZJets";

  // name of input xAOD containers (to be specified in configuration)
  std::string m_pflowJetContainerName = "";
  std::string m_topoJetContainerName = "";
  std::string m_electronContainerName = "";
  std::string m_photonContainerName = "";
  std::string m_muonContainerName = "";

  // name of truth containers (to be specified in configuration)
  std::string m_truthContainerName = "";
  std::string m_truthJetContainerName = "";
  std::string m_truthWZJetContainerName = "";

  // PDG ID of truth objects
  int m_truthPDGID = 0;

  // parameters for truth--reco object matching
  double m_deltaRcutoff = 10000.;
  double m_ptDiffCutoff = 10000.;

  // convert pT and energy to GeV or leave it in the default MeV
  bool m_convertMeVToGeV = true;

  // Metadata
  std::string m_derivationName = "";
  bool m_useMetaData = true;

private:

  // boolean flags to be set depending on non-empty input container strings
  bool m_isZee = false;
  bool m_isZmumu = false;
  bool m_isttbar = false;
  bool m_isSinglePhoton = false;
  bool m_isDijet = false;

  TTree *m_infoTree = nullptr;
  TTree *m_pflowJetTree = nullptr;
  TTree *m_topoJetTree = nullptr;
  TTree *m_electronTree = nullptr;
  TTree *m_photonTree = nullptr;
  TTree *m_muonTree = nullptr;
  TTree *m_truthTree = nullptr;
  TTree *m_truthJetTree = nullptr;
  TTree *m_truthWZJetTree = nullptr;

  int m_runNumber;
  int m_evtNumber;
  int m_mcChannelNumber;
  int m_nv; //number of vertices in the "PrimaryVertices" container
  int m_npv; //number of primary vertices in event
  int m_npuv; //number of pileup vertices in event
  int m_mu; //mean number of interactions per bunch crossing
  double m_MD_initialSumW; //weighted number of events
  double m_mcEventWeights;
  std::vector<double> m_pt;
  std::vector<double> m_eta;
  std::vector<double> m_phi;
  std::vector<double> m_e;
  std::vector<int> m_truthID;
  std::vector<int> m_truthBarcode;

  // identification selection check
  std::vector<int> m_passSel;

  // jet constituent-level quantities
  std::vector<double> m_pt_constit;
  std::vector<double> m_eta_constit;
  std::vector<double> m_phi_constit;
  std::vector<double> m_e_constit;

  // jet detector eta
  std::vector<float> m_detEta;

  // jet b-tag and scale factors
  std::vector<int> m_isTag;
  std::vector<std::vector<float>> m_sf;

  // number of tracks associated with jet
  std::vector<int> m_nTrk;

  // jet JVT check
  std::vector<int> m_passJVT;

  // jet isolation check
  std::vector<int> m_isIsoJetDR0p6; //for reco jets
  std::vector<int> m_isIsoJetDR1p0; //for truth jets

  // charged and neutral FEs
  std::vector<std::vector<int>> m_neutralPFOindex;
  std::vector<std::vector<int>> m_chargedPFOindex;
  std::vector<std::vector<double>> m_neutralPFOenergy;
  std::vector<std::vector<double>> m_chargedPFOenergy;

  // calibration hit for neutral FEs
  std::vector<std::vector<std::vector<int>>> m_calpfo_NLeadingTruthParticlePdgId;
  std::vector<std::vector<std::vector<int>>> m_calpfo_NLeadingTruthParticleBarcode;
  std::vector<std::vector<std::vector<double>>> m_calpfo_NLeadingTruthParticleEnergy;

  // unit conversion factor (MeV -> GeV)
  double m_conversionFactor = 1.0;

  // FE links for objects
  SG::ReadDecorHandleKey<xAOD::ElectronContainer> m_electronNeutralFEReadDecorKey;
  SG::ReadDecorHandleKey<xAOD::ElectronContainer> m_electronChargedFEReadDecorKey;
  SG::ReadDecorHandleKey<xAOD::PhotonContainer> m_photonNeutralFEReadDecorKey;
  SG::ReadDecorHandleKey<xAOD::MuonContainer> m_muonNeutralFEReadDecorKey;
  SG::ReadDecorHandleKey<xAOD::MuonContainer> m_muonChargedFEReadDecorKey;

  // Global PFlow container keys
  SG::ReadHandleKey<xAOD::FlowElementContainer> m_inGlobalChargedFEKey;
  SG::ReadHandleKey<xAOD::FlowElementContainer> m_inGlobalNeutralFEKey;

  // jet container key
  SG::ReadHandleKey<xAOD::JetContainer> m_jetContKey;

public:

  // this is a standard constructor
  compareFELinks ();

  // destructor
  virtual ~compareFELinks();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  /// @cond
  // this is needed to distribute the algorithm to the workers
  ClassDef(compareFELinks, 1);
  /// @endcond

};

#endif
