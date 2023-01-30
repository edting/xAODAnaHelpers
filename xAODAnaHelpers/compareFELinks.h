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

#ifndef xAODAnaHelpers_compareFELinks_H
#define xAODAnaHelpers_compareFELinks_H

// algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"

#include "AsgTools/ToolHandle.h"
#include "AsgDataHandles/ReadDecorHandle.h" //needed to get FlowElementLinks from electron/photon/muon
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
//#include "TFile.h"
#include "TTree.h"

class compareFELinks : public xAH::Algorithm
{
public:
  // The variables in this public scope can be changed by the user in the xAH config file

  // name of output root file produced by algorithm
  // note that this is RECREATED if the algorithm is called more than once
  // so any existing file will be OVERWRITTEN unless a different name is specified in the job configuration
  std::string m_outFileName = "compareFELinks.root";

  // name of trees to write to output file
  std::string m_infoTreeName = "info";
  std::string m_pflowJetTreeName = "pflowJets";
  std::string m_topoJetTreeName = "topoJets";
  std::string m_electronTreeName = "electrons";
  std::string m_photonTreeName = "photons";
  std::string m_muonTreeName = "muons";
  std::string m_truthTreeName = "truthParticles";
  std::string m_truthJetTreeName = "truthJets";
  std::string m_truthWZJetTreeName = "truthWZJets";

  // name of container(s) containing calibrated objects
  std::string m_pflowJetContainerName = "";
  std::string m_topoJetContainerName = "";
  std::string m_electronContainerName = "";
  std::string m_photonContainerName = "";
  std::string m_muonContainerName = "";

  // name of truth containers for objects / jets / jets without WZ contributions
  std::string m_truthContainerName = "";
  std::string m_truthJetContainerName = "";
  std::string m_truthWZJetContainerName = "";

  // PDG ID of truth objects
  int m_truthPDGID = 0;

  // parameters for truth--reco matching
  double m_deltaRcutoff = 10000.;
  double m_ptDiffCutoff = 10000.;

  // convert pT and energy to GeV or leave it in the default MeV
  bool m_convertMeVToGeV = true;

  // Metadata
  std::string m_derivationName = "";
  bool m_useMetaData = true;

private:

  // one of these will be set to true depending on which object container names are passed to the algorithm
  // only electron container name specified -> Zee
  //      muon -> Zmumu
  //      electron + muon -> ttbar
  //      photon -> SinglePhoton
  //      electron + muon + photon -> dijet
  bool m_isZee = false;
  bool m_isZmumu = false;
  bool m_isttbar = false;
  bool m_isSinglePhoton = false;
  bool m_isDijet = false;

  //  TFile *m_file = nullptr;
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
  int m_nv; //number of vertices inthe "PrimaryVertices" container
  int m_npv; //number of primary vertices in event
  int m_npuv; //number of pileup vertices in event
  int m_mu; //mean number of pp interactions per bunch crossing
  double m_MD_initialSumW; //weighted number of events
  double m_mcEventWeights;
  std::vector<double> m_pt;
  std::vector<double> m_eta;
  std::vector<double> m_phi;
  std::vector<double> m_e;
  std::vector<int> m_truthID;
  std::vector<int> m_truthBarcode;

  // decoration for whether object passes its corresponding identification criteria
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

  // decoration for whether or not jet passes JVT cut
  std::vector<int> m_passJVT;

  // check whether a jet has any other jets within some deltaR of it
  std::vector<int> m_isIsoJetDR0p6; //for reco jets; true if no other jet with pT > 7 GeV in deltaR < 0.6
  std::vector<int> m_isIsoJetDR1p0; //for truth jets; deltaR < 1

  // for each electron/muon/photon, vector containing PFO indices and energies (each object can be associated with multiple PFOs)
  // for jets, entries will correspond to the PFO index of each jet constituent (each constituent is associated with 1 PFO)
  std::vector<std::vector<int>> m_neutralPFOindex;
  std::vector<std::vector<int>> m_chargedPFOindex;
  std::vector<std::vector<double>> m_neutralPFOenergy;
  std::vector<std::vector<double>> m_chargedPFOenergy;

  // for each neutral PFO of each jet, use calpfo_NLeadingTruthParticleBarcodeEnergyPairs to determine PDGID of truth particle that deposited energy (and also save that energy value)
  std::vector<std::vector<std::vector<int>>> m_calpfo_NLeadingTruthParticlePdgId; //the inner nested vector<int> is because there can be multiple particles contributing to one PFO
  std::vector<std::vector<std::vector<int>>> m_calpfo_NLeadingTruthParticleBarcode;
  std::vector<std::vector<std::vector<double>>> m_calpfo_NLeadingTruthParticleEnergy;

  // controls whether we convert pT and energy from MeV to GeV, depending on m_convertMeVToGeV
  double m_conversionFactor = 1.0;

  // flow element link decoration keys -- for links from physics objects to FEs
  SG::ReadDecorHandleKey<xAOD::ElectronContainer> m_electronNeutralFEReadDecorKey;
  SG::ReadDecorHandleKey<xAOD::ElectronContainer> m_electronChargedFEReadDecorKey;
  SG::ReadDecorHandleKey<xAOD::PhotonContainer> m_photonNeutralFEReadDecorKey;
  SG::ReadDecorHandleKey<xAOD::MuonContainer> m_muonNeutralFEReadDecorKey;
  SG::ReadDecorHandleKey<xAOD::MuonContainer> m_muonChargedFEReadDecorKey;

  /* // decor keys for links between JetETMiss and Global containers */
  /* SG::ReadDecorHandleKey<xAOD::FlowElementContainer> m_neutralOriginalToGlobalFEReadDecorKey; */
  /* SG::ReadDecorHandleKey<xAOD::FlowElementContainer> m_chargedOriginalToGlobalFEReadDecorKey; */
  /* SG::ReadDecorHandleKey<xAOD::FlowElementContainer> m_chargedOriginalToNeutralGlobalFEReadDecorKey; */

  // flow element container key
  SG::ReadHandleKey<xAOD::FlowElementContainer> m_inGlobalChargedFEKey;//{this, "InGlobalChargedFEKey", "", "ReadHandleKey for modified Charged FlowElements"};
  SG::ReadHandleKey<xAOD::FlowElementContainer> m_inGlobalNeutralFEKey;//{this, "InGlobalNeutralFEKey", "", "ReadHandleKey for modified Neutral FlowElements"};

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
