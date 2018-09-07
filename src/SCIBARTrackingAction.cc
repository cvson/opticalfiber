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
// $Id: SCIBARTrackingAction.cc 90240 2015-05-21 09:08:13Z gcosmo $
//
/// \file optical/wls/src/SCIBARTrackingAction.cc
/// \brief Implementation of the SCIBARTrackingAction class
//
//
#include "globals.hh"
#include "G4RunManager.hh"

#include "SCIBARTrajectory.hh"

#include "SCIBARUserTrackInformation.hh"

#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4TrackingManager.hh"

#include "SCIBARTrackingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  //Let this be up to the user via vis.mac
  //  fpTrackingManager->SetStoreTrajectory(true);

  //Use custom trajectory class
  fpTrackingManager->SetTrajectory(new SCIBARTrajectory(aTrack));

  SCIBARUserTrackInformation* trackInformation = new SCIBARUserTrackInformation();

  if (aTrack->GetMomentumDirection().z()>0.0) {
     trackInformation->AddStatusFlag(right);
  } else {
     trackInformation->AddStatusFlag(left);
  }

  G4String PVName = aTrack->GetVolume()->GetName();

  if (PVName == "SCIBARFiber" || PVName == "Clad1" || PVName == "Clad2")
     trackInformation->AddStatusFlag(InsideOfFiber);

  fpTrackingManager->SetUserTrackInformation(trackInformation);
}
