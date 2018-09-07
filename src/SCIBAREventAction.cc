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
// $Id: SCIBAREventAction.cc 90240 2015-05-21 09:08:13Z gcosmo $
//
/// \file optical/wls/src/SCIBAREventAction.cc
/// \brief Implementation of the SCIBAREventAction class
//
//
#include "SCIBAREventAction.hh"

#include "SCIBARRunAction.hh"

#include "SCIBAREventActionMessenger.hh"

#include "SCIBARPhotonDetHit.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "SCIBARTrajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"

#include "Randomize.hh"

// Purpose: Accumulates statistics regarding hits
//          in the PhotonDet detector

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SCIBAREventAction::SCIBAREventAction(SCIBARRunAction* runaction)
 : fRunAction(runaction), fVerboseLevel(0)
{
  fMPPCCollID = 0;

  fEventMessenger = new SCIBAREventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SCIBAREventAction::~SCIBAREventAction()
{
  delete fEventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBAREventAction::BeginOfEventAction(const G4Event* evt)
{
 G4int evtNb = evt->GetEventID();

 if(fVerboseLevel>0)
    G4cout << "<<< Event  " << evtNb << " started." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Threading.hh"

void SCIBAREventAction::EndOfEventAction(const G4Event* evt)
{
  if (fVerboseLevel>0)
     G4cout << "<<< Event  " << evt->GetEventID() << " ended." << G4endl;
 
  if (fRunAction->GetRndmFreq() == 2)
    {
     std::ostringstream os;
     os<<"endOfEvent_"<<G4Threading::G4GetThreadId()<<".rndm";
     G4Random::saveEngineStatus(os.str().c_str());
    }

  // Get Hits from the detector if any
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colName = "PhotonDetHitCollection";
  fMPPCCollID = SDman->GetCollectionID(colName);

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  SCIBARPhotonDetHitsCollection* mppcHC = 0;

  // Get the hit collections
  if (HCE)
  {
     if (fMPPCCollID>=0) mppcHC = 
                        (SCIBARPhotonDetHitsCollection*)(HCE->GetHC(fMPPCCollID));
  }

  // Get hit information about photons that reached the detector in this event
  if (mppcHC)
  {
	//newadd
	G4int n_hit = mppcHC->entries();
	G4cout << "<<< Event  " << evt->GetEventID() << " having hit."<<n_hit << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int SCIBAREventAction::GetEventNo()
{
  return fpEventManager->GetConstCurrentEvent()->GetEventID();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBAREventAction::SetEventVerbose(G4int level)
{
  fVerboseLevel = level;
}
