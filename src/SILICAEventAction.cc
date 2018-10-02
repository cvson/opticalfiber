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
// $Id: SILICAEventAction.cc 90240 2015-05-21 09:08:13Z gcosmo $
//
/// \file optical/wls/src/SILICAEventAction.cc
/// \brief Implementation of the SILICAEventAction class
//
//
#include "SILICAEventAction.hh"
#include "SILICAStackingAction.hh"
#include "SILICARunAction.hh"
#include "SILICAHistoManager.hh"
#include "SILICAEventActionMessenger.hh"

#include "SILICAPhotonDetHit.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "SILICATrajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"

#include "Randomize.hh"

// Purpose: Accumulates statistics regarding hits
//          in the PhotonDet detector

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SILICAEventAction::SILICAEventAction(SILICARunAction* runaction, SILICAHistoManager* histo)
:G4UserEventAction(),
fRunAction(runaction),fHistoManager(histo),fVerboseLevel(0),fMPPCCollID(0),
 fEnergyDep(0.),fTrackLDep(0.), fEnergyDepHit(0.),
fprix(0.), fpriy(0.), fpriz(0.), fpripz(0.), fpritheta(0.),
fnphotonpass2end(0), fnphotonfail2end(0),
fnphotonscint(0), fnphotonscheren(0), fenscint(0.0), fencheren(0.0)
{
 //fMPPCCollID = 0;
 fEventMessenger = new SILICAEventActionMessenger(this); }



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SILICAEventAction::~SILICAEventAction()
{
  delete fEventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICAEventAction::BeginOfEventAction(const G4Event* evt)
{
 G4int evtNb = evt->GetEventID();

 if(fVerboseLevel>0)
    G4cout << "<<< Event  " << evtNb << " started." << G4endl;
 // initialisation per event
 fEnergyDep = 0.;//newadd
 fTrackLDep = 0.;//newadd
    //remember to reset
    fnphotonpass2end = 0;
    fnphotonfail2end = 0;
    fnphotonscint = 0;
    fnphotonscheren = 0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Threading.hh"

void SILICAEventAction::EndOfEventAction(const G4Event* evt)
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
  SILICAPhotonDetHitsCollection* mppcHC = 0;

  // Get the hit collections
  if (HCE)
  {
     if (fMPPCCollID>=0) mppcHC = 
                        (SILICAPhotonDetHitsCollection*)(HCE->GetHC(fMPPCCollID));
  }

  // Get hit information about photons that reached the detector in this event
    G4int n_hit;
  if (mppcHC)
  {
	//newadd
    n_hit = mppcHC->entries();
	G4cout << "<<< Event  " << evt->GetEventID() << " having hit."<<n_hit << G4endl;
  }
  else n_hit = 0;
    //get number of photon electron?
    G4int n_photondet;
    
    if (n_hit==0) {
        n_photondet=0;
    } else {
    for(int l=0; l<n_hit; l++){
        //to count number of photon?
        //n_photondet = n_hit;//TODO
        SILICAPhotonDetHit* hitth = (*mppcHC)[l];
        fEnergyDepHit += hitth->GetEnergyDeposit();//newadd3
    }//end for l
    }//end else n_hit
	//NEW ADD
	  //accumulates statistic
  //
  fRunAction->fillPerEvent(fEnergyDep, fTrackLDep);
  
  //fill histograms
  //
  fHistoManager->FillHisto(0, fEnergyDep);
  fHistoManager->FillHisto(1, fTrackLDep);
  fHistoManager->FillHisto(2, fEnergyDep*1.0/fTrackLDep);
    fHistoManager->FillHisto(3, n_hit);
    //this class is not work? why?
    //move this to StackAction class
    /*SILICAStackingAction * pfstachaction = new SILICAStackingAction();
    G4int n_opticalphoton = pfstachaction->GetNumberofOpticalPhoton();
    fHistoManager->FillHisto(4, n_opticalphoton);*/
    fHistoManager->FillHisto(5, fEnergyDepHit);//newadd3
  //fill ntuple
  //
  fHistoManager->FillNtuple(fEnergyDep, fTrackLDep);
 fHistoManager->FillPrimTrack(fprix, fpriy, fpriz, fpripz, fpritheta,fnphotonpass2end,fnphotonfail2end,fnphotonscint,fnphotonscheren,fenscint,fencheren);//newadd
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int SILICAEventAction::GetEventNo()
{
  return fpEventManager->GetConstCurrentEvent()->GetEventID();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICAEventAction::SetEventVerbose(G4int level)
{
  fVerboseLevel = level;
}
