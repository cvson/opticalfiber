//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: SILICAOpticalPhysics.cc 102356 2017-01-23 16:22:42Z gcosmo $
//
/// \file optical/wls/src/SILICAOpticalPhysics.cc
/// \brief Implementation of the SILICAOpticalPhysics class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"

#include "SILICAOpticalPhysics.hh"

SILICAOpticalPhysics::SILICAOpticalPhysics(G4bool toggle)
    : G4VPhysicsConstructor("Optical")
{
  fSILICAProcess                = NULL;
  fScintProcess              = NULL;
  fCerenkovProcess           = NULL;
  fBoundaryProcess           = NULL;
  fAbsorptionProcess         = NULL;
  fRayleighScattering        = NULL;
  fMieHGScatteringProcess    = NULL;

  fAbsorptionOn              = toggle;
}

SILICAOpticalPhysics::~SILICAOpticalPhysics() { }

#include "G4OpticalPhoton.hh"

void SILICAOpticalPhysics::ConstructParticle()
{
  G4OpticalPhoton::OpticalPhotonDefinition();
}

#include "G4ProcessManager.hh"

void SILICAOpticalPhysics::ConstructProcess()
{
    G4cout << "SILICAOpticalPhysics:: Add Optical Physics Processes"
           << G4endl;

  fSILICAProcess = new G4OpWLS();

  fScintProcess = new G4Scintillation();
  fScintProcess->SetScintillationYieldFactor(1.);
  fScintProcess->SetTrackSecondariesFirst(true);

  fCerenkovProcess = new G4Cerenkov();
  fCerenkovProcess->SetMaxNumPhotonsPerStep(300);
  fCerenkovProcess->SetTrackSecondariesFirst(true);

  fAbsorptionProcess      = new G4OpAbsorption();
  fRayleighScattering     = new G4OpRayleigh();
  fMieHGScatteringProcess = new G4OpMieHG();
  fBoundaryProcess        = new G4OpBoundaryProcess();

  G4ProcessManager* pManager =
                G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

  if (!pManager) {
     std::ostringstream o;
     o << "Optical Photon without a Process Manager";
     G4Exception("SILICAOpticalPhysics::ConstructProcess()","",
                  FatalException,o.str().c_str());
  }

  if (fAbsorptionOn) pManager->AddDiscreteProcess(fAbsorptionProcess);

  //pManager->AddDiscreteProcess(fRayleighScattering);
  //pManager->AddDiscreteProcess(fMieHGScatteringProcess);

  pManager->AddDiscreteProcess(fBoundaryProcess);
    //this is for WLS fiber
  //fSILICAProcess->UseTimeProfile("delta");
  fSILICAProcess->UseTimeProfile("exponential");

  pManager->AddDiscreteProcess(fSILICAProcess);

  fScintProcess->SetScintillationYieldFactor(1.);
  fScintProcess->SetScintillationExcitationRatio(0.0);
  fScintProcess->SetTrackSecondariesFirst(true);

  // Use Birks Correction in the Scintillation process

  G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
  fScintProcess->AddSaturation(emSaturation);

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while ( (*particleIterator)() ){

    G4ParticleDefinition* particle = particleIterator->value();
    G4String particleName = particle->GetParticleName();

    pManager = particle->GetProcessManager();
    if (!pManager) {
       std::ostringstream o;
       o << "Particle " << particleName << "without a Process Manager";
       G4Exception("SILICAOpticalPhysics::ConstructProcess()","",
                    FatalException,o.str().c_str());
    }

    if(fCerenkovProcess->IsApplicable(*particle)){
      pManager->AddProcess(fCerenkovProcess);
      pManager->SetProcessOrdering(fCerenkovProcess,idxPostStep);
    }
    if(fScintProcess->IsApplicable(*particle)){
      pManager->AddProcess(fScintProcess);
      pManager->SetProcessOrderingToLast(fScintProcess,idxAtRest);
      pManager->SetProcessOrderingToLast(fScintProcess,idxPostStep);
    }

  }
}

void SILICAOpticalPhysics::SetNbOfPhotonsCerenkov(G4int maxNumber)
{
  fCerenkovProcess->SetMaxNumPhotonsPerStep(maxNumber);
}
