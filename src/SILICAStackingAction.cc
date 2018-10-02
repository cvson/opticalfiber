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
// $Id: SILICAStackingAction.cc 69561 2013-05-08 12:25:56Z gcosmo $
//
/// \file optical/wls/src/SILICAStackingAction.cc
/// \brief Implementation of the SILICAStackingAction class
//
//
#include "SILICAStackingAction.hh"

#include "G4RunManager.hh"

#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SILICAStackingAction::SILICAStackingAction(SILICAHistoManager* histo) : fHistoManager(histo),fPhotonCounter(0) { }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SILICAStackingAction::~SILICAStackingAction() { }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
      SILICAStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  G4ParticleDefinition* particleType = aTrack->GetDefinition();

  // keep primary particle
  if (aTrack->GetParentID() == 0) return fUrgent;

  if (particleType == G4OpticalPhoton::OpticalPhotonDefinition()) {
     // keep optical photon
     fPhotonCounter++;
     return fUrgent;
  } else {
     // discard all other secondaries
     // return fKill;
  }
  return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICAStackingAction::NewStage() {
   G4cout << "Number of optical photons produces in this event : "
          << fPhotonCounter << G4endl;
    fHistoManager->FillHisto(4, fPhotonCounter);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICAStackingAction::PrepareNewEvent() { fPhotonCounter = 0; }
