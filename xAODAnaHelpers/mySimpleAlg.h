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

private:

  // FE links for objects
  SG::ReadDecorHandleKey<xAOD::ElectronContainer> m_electronNeutralFEReadDecorKey;
  SG::ReadDecorHandleKey<xAOD::ElectronContainer> m_electronChargedFEReadDecorKey;

  // Global PFlow container keys
  SG::ReadHandleKey<xAOD::FlowElementContainer> m_inGlobalChargedFEKey;
  SG::ReadHandleKey<xAOD::FlowElementContainer> m_inGlobalNeutralFEKey;
  SG::ReadHandleKey<xAOD::FlowElementContainer> m_inCHSGChargedFEKey;
  SG::ReadHandleKey<xAOD::FlowElementContainer> m_inCHSGNeutralFEKey;
  SG::ReadHandleKey<xAOD::FlowElementContainer> m_inCSSKGChargedFEKey;
  SG::ReadHandleKey<xAOD::FlowElementContainer> m_inCSSKGNeutralFEKey;

  // more decor keys
  SG::ReadDecorHandleKey<xAOD::FlowElementContainer> m_neutralFECellsRemovedReadDecorKey;

  // container keys
  SG::ReadHandleKey<xAOD::JetContainer> m_jetContKey;
  SG::ReadHandleKey<xAOD::ElectronContainer> m_eleContKey;
  SG::ReadHandleKey<xAOD::TruthParticleContainer> m_truthParticleContKey;

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
