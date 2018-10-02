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
// $Id: SILICARunAction.cc 70603 2013-06-03 11:23:16Z gcosmo $
//
/// \file optical/wls/src/SILICARunAction.cc
/// \brief Implementation of the SILICARunAction class
//
//
#include "SILICARunAction.hh"
#include "SILICAHistoManager.hh"//newadd
#include "SILICARunActionMessenger.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh" //newadd

#include "Randomize.hh"

#include "SILICADetectorConstruction.hh"
#include "SILICASteppingAction.hh"

#include <ctime>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*SILICARunAction::SILICARunAction()
  : fSaveRndm(0), fAutoSeed(false)
{
  fRunMessenger = new SILICARunActionMessenger(this);
}*/
SILICARunAction::SILICARunAction(SILICAHistoManager* histo)
: G4UserRunAction(),
  fHistoManager(histo),
  fSumEDep(0.), fSum2EDep(0.),
  fSumLDep(0.), fSum2LDep(0.), fSaveRndm(1), fAutoSeed(false)
{
fRunMessenger = new SILICARunActionMessenger(this);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SILICARunAction::~SILICARunAction()
{
  delete fRunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICARunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  G4RunManager::GetRunManager()->SetRandomNumberStoreDir("random/");

  if (fAutoSeed) {
     // automatic (time-based) random seeds for each run
     G4cout << "*******************" << G4endl;
     G4cout << "*** AUTOSEED ON ***" << G4endl;
     G4cout << "*******************" << G4endl;
     long seeds[2];
     time_t systime = time(NULL);
     seeds[0] = (long) systime;
     seeds[1] = (long) (systime*G4UniformRand());
     G4Random::setTheSeeds(seeds);
     G4Random::showEngineStatus();
  } else {
     G4Random::showEngineStatus();
  }

  if (fSaveRndm > 0) G4Random::saveEngineStatus("BeginOfRun.rndm");
//newadd
 //
  fSumEDep = fSum2EDep = 0.;
  fSumLDep = fSum2LDep = 0.;
  
  //histograms
  //
  fHistoManager->Book();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICARunAction::fillPerEvent(G4double EDep,G4double LDep)
{ 
  //accumulate statistic
  //
  fSumEDep += EDep;  fSum2EDep += EDep*EDep;
  
  fSumLDep += LDep;  fSum2LDep += LDep*LDep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICARunAction::EndOfRunAction(const G4Run* aRun )
{
  if (fSaveRndm == 1)
  {
     G4Random::showEngineStatus();
     G4Random::saveEngineStatus("endOfRun.rndm");
//newadd
 G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  
  //compute statistics: mean and rms
  //
  fSumEDep /= NbOfEvents; fSum2EDep /= NbOfEvents;
  G4double rmsEDep = fSum2EDep - fSumEDep*fSumEDep;
  if (rmsEDep >0.) rmsEDep = std::sqrt(rmsEDep); else rmsEDep = 0.;
  

  
  fSumLDep /= NbOfEvents; fSum2LDep /= NbOfEvents;
  G4double rmsLDep = fSum2LDep - fSumLDep*fSumLDep;
  if (rmsLDep >0.) rmsLDep = std::sqrt(rmsLDep); else rmsLDep = 0.;
  
 
  
  //print
  //
  G4cout
     << "\n--------------------End of Run------------------------------\n"
     << "\n mean Energy in scintillator : " << G4BestUnit(fSumEDep,"Energy")
     << " +- "                          << G4BestUnit(rmsEDep,"Energy")
     << G4endl;
     
  G4cout
     << "\n mean trackLength in scintillator : " << G4BestUnit(fSumLDep,"Length")
     << " +- "                               << G4BestUnit(rmsLDep,"Length")  
     << "\n------------------------------------------------------------\n"
     << G4endl;
     
  //save histograms
  //
  fHistoManager->PrintStatistic();
  fHistoManager->Save();   
  }

}
