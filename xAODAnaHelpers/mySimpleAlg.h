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

#ifndef xAODAnaHelpers_mySimpleAlg_H
#define xAODAnaHelpers_mySimpleAlg_H

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
#include "xAODCaloEvent/CaloClusterContainer.h"

// ROOT includes
#include "TTree.h"

class mySimpleAlg : public xAH::Algorithm
{
public:
  // The variables in this public scope can be changed by the user in the xAH config file

    // name of output ntuple root file
  std::string m_outFileName = "simpleNtuple.root";

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

  // name of input xAOD containers
  std::string m_pflowJetContainerName = "AntiKt4EMPFlowJets";
  std::string m_topoJetContainerName = "";//"AntiKt4EMTopoJets";
  std::string m_electronContainerName = "";//"Electrons";
  std::string m_photonContainerName = "";//"Photons";
  std::string m_muonContainerName = "";//"Muons";

  // name of truth containers (to be specified in configuration)
  std::string m_truthContainerName = "";//"TruthParticles";
  std::string m_truthJetContainerName = "";//"AntiKt4TruthJets";
  std::string m_truthWZJetContainerName = "AntiKt4TruthWZJets";

  // PDG ID of truth objects
  int m_truthPDGID = 0;

  // parameters for truth--reco object matching
  double m_deltaRcutoff = 10000.;
  double m_ptDiffCutoff = 10000.;

  // toggle TRUE makes it so only truth-matched objects are saved
  bool m_doTruthObjectMatch = false;

  // whether topocluster calibration hits is available in input or not
  bool m_calclusIsAvail = false;

  // convert pT and energy to GeV or leave it in the default MeV
  bool m_convertMeVToGeV = true;

private:

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
  double m_mcEventWeights;
  std::vector<double> m_pt;
  std::vector<double> m_eta;
  std::vector<double> m_phi;
  std::vector<double> m_e;
  std::vector<int> m_truthID;
  std::vector<int> m_truthBarcode;

  // number of tracks associated with jet
  std::vector<int> m_nTrk;

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

  // calibration hit for CaloCaloTopoClusters
  std::vector<std::vector<std::vector<int>>> m_calclus_NLeadingTruthParticlePdgId;
  std::vector<std::vector<std::vector<int>>> m_calclus_NLeadingTruthParticleBarcode;
  std::vector<std::vector<std::vector<double>>> m_calclus_NLeadingTruthParticleEnergy;

  // unit conversion factor (MeV -> GeV)
  double m_conversionFactor = 0.001;

  // FE links for objects
  SG::ReadDecorHandleKey<xAOD::ElectronContainer> m_electronNeutralFEReadDecorKey;
  SG::ReadDecorHandleKey<xAOD::ElectronContainer> m_electronChargedFEReadDecorKey;
  SG::ReadDecorHandleKey<xAOD::PhotonContainer> m_photonNeutralFEReadDecorKey;
  SG::ReadDecorHandleKey<xAOD::MuonContainer> m_muonNeutralFEReadDecorKey;
  SG::ReadDecorHandleKey<xAOD::MuonContainer> m_muonChargedFEReadDecorKey;

  // Global PFlow container keys
  SG::ReadHandleKey<xAOD::FlowElementContainer> m_inGlobalChargedFEKey;
  SG::ReadHandleKey<xAOD::FlowElementContainer> m_inGlobalNeutralFEKey;
  SG::ReadHandleKey<xAOD::FlowElementContainer> m_inCHSGChargedFEKey;
  SG::ReadHandleKey<xAOD::FlowElementContainer> m_inCHSGNeutralFEKey;

  // container keys
  SG::ReadHandleKey<xAOD::JetContainer> m_jetContKey;
  SG::ReadHandleKey<xAOD::JetContainer> m_jetTopoContKey;
  SG::ReadHandleKey<xAOD::JetContainer> m_jetTruthContKey;
  SG::ReadHandleKey<xAOD::JetContainer> m_jetTruthWZContKey;
  SG::ReadHandleKey<xAOD::ElectronContainer> m_eleContKey;
  SG::ReadHandleKey<xAOD::PhotonContainer> m_phoContKey;
  SG::ReadHandleKey<xAOD::MuonContainer> m_muContKey;
  SG::ReadHandleKey<xAOD::TruthParticleContainer> m_truthContKey;
  SG::ReadHandleKey<xAOD::TruthParticleContainer> m_truthParticlesContKey;

public:

  // this is a standard constructor
  mySimpleAlg ();

  // destructor
  virtual ~mySimpleAlg();

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
  ClassDef(mySimpleAlg, 1);
  /// @endcond

};

#endif
