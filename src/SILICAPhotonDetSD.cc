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
// $Id: SILICAPhotonDetSD.cc 69561 2013-05-08 12:25:56Z gcosmo $
//
/// \file optical/wls/src/SILICAPhotonDetSD.cc
/// \brief Implementation of the SILICAPhotonDetSD class
//
//
#include "SILICAPhotonDetSD.hh"
#include "SILICAPhotonDetHit.hh"
#include "SILICAUserTrackInformation.hh"

#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SILICAPhotonDetSD::SILICAPhotonDetSD(G4String name)
  : G4VSensitiveDetector(name), fPhotonDetHitCollection(0)
{
  collectionName.insert("PhotonDetHitCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SILICAPhotonDetSD::~SILICAPhotonDetSD() { }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICAPhotonDetSD::Initialize(G4HCofThisEvent* HCE)
{
  fPhotonDetHitCollection =
       new SILICAPhotonDetHitsCollection(SensitiveDetectorName,collectionName[0]);
  //Store collection with event and keep ID
  static G4int HCID = -1;
  if (HCID<0) HCID = GetCollectionID(0);
  HCE->AddHitsCollection( HCID, fPhotonDetHitCollection );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SILICAPhotonDetSD::ProcessHits(G4Step* , G4TouchableHistory* )
{
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SILICAPhotonDetSD::ProcessHits_constStep(const G4Step* aStep,
                                             G4TouchableHistory* )
//Generates a hit and uses the postStepPoint; PostStepPoint because the hit
//is generated manually when the photon hits the detector
{
  if (aStep == NULL) return false;
    //Retrieve the energy deposited from the step
    G4double energyDeposit = aStep -> GetTotalEnergyDeposit();//newadd3
    
  G4Track* theTrack = aStep->GetTrack();

  // Need to know if this is an optical photon
  if(theTrack->GetDefinition()
     != G4OpticalPhoton::OpticalPhotonDefinition()) return false;

  // Find out information regarding the hit
  G4StepPoint* thePostPoint = aStep->GetPostStepPoint();
 
  SILICAUserTrackInformation* trackInformation
      = (SILICAUserTrackInformation*)theTrack->GetUserInformation();
 
  G4TouchableHistory* theTouchable
      = (G4TouchableHistory*)(thePostPoint->GetTouchable());
 
  G4ThreeVector photonExit   = trackInformation -> GetExitPosition();
  G4ThreeVector photonArrive = thePostPoint -> GetPosition();
  G4double      arrivalTime  = theTrack -> GetGlobalTime();

  // Convert the global coordinate for arriving photons into
  // the local coordinate of the detector
  photonArrive = theTouchable->GetHistory()->
                                GetTopTransform().TransformPoint(photonArrive);

  // Creating the hit and add it to the collection
  /*fPhotonDetHitCollection->
            insert(new SILICAPhotonDetHit(photonExit, photonArrive, arrivalTime));*/
    fPhotonDetHitCollection->
    insert(new SILICAPhotonDetHit(photonExit, photonArrive, arrivalTime,energyDeposit));

  return true;
}
