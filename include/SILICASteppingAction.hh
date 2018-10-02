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
// $Id: SILICASteppingAction.hh 69561 2013-05-08 12:25:56Z gcosmo $
//
/// \file optical/wls/include/SILICASteppingAction.hh
/// \brief Definition of the SILICASteppingAction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef SILICASteppingAction_h
#define SILICASteppingAction_h 1

#include "G4UserSteppingAction.hh"

class SILICADetectorConstruction;
class SILICAEventAction;//newadd
class SILICASteppingActionMessenger;

class G4Track;
class G4StepPoint;

class G4OpBoundaryProcess;

class SILICASteppingAction : public G4UserSteppingAction
{
  public:

    //SILICASteppingAction(SILICADetectorConstruction*);
    SILICASteppingAction(SILICADetectorConstruction*, SILICAEventAction*);
    virtual ~SILICASteppingAction();

    virtual void UserSteppingAction(const G4Step*);
 
    // Set the bounce limit, 0 for no limit
    void  SetBounceLimit(G4int);
 
    G4int GetNumberOfBounces();
    G4int GetNumberOfClad1Bounces();
    G4int GetNumberOfClad2Bounces();
    G4int GetNumberOfSILICABounces();
    //unless you need to store for each step
    G4int GetNumberOfPhotonsPass2end(){return fCounterEnd;};//newadd2
    G4int GetNumberOfPhotonsFail2end(){return fCounterMid;};//newadd2
    G4int GetNumberOfPhotonsScintillation(){return fScintillationCounter;};
    G4int GetNumberOfPhotonsCherenkov(){return fCerenkovCounter;};
    // return number of successful events and reset the counter
    G4int ResetSuccessCounter();
 
  private:

    // Artificially kill the photon after it has bounced more than this number
    G4int fBounceLimit;
    // number of photons that reach the end
    G4int fCounterEnd;
    // number of photons that didn't make it to the end
    G4int fCounterMid;
    // total number of bounces that a photon been through
    G4int fCounterBounce;
    // number of bounces that a photon been through within the fibre
    G4int fCounterSILICABounce;
    // number of bounces that a photon been through from Cladding 1 to 2
    G4int fCounterClad1Bounce;
    // number of bounces that a photon been through from Cladding 2 to World
    G4int fCounterClad2Bounce;
    
    G4int fScintillationCounter;//newadd4
    G4int fCerenkovCounter;//newadd4

    // initial gamma of the photon
    G4double fInitGamma;
    // initial theta of the photon
    G4double fInitTheta;

    G4OpBoundaryProcess* fOpProcess;

    // maximum number of save states
    static G4int fMaxRndmSave;
 
    SILICADetectorConstruction* fDetector;
    SILICAEventAction* fEventAction;//newadd

    SILICASteppingActionMessenger* fSteppingMessenger;

    inline void ResetCounters()
    { 
      fCounterBounce = fCounterSILICABounce =
      fCounterClad1Bounce = fCounterClad2Bounce = 0;
        fScintillationCounter = fCerenkovCounter = 0;//newadd4
      fInitGamma = fInitTheta = -1;
    }

    // save the random status into a sub-directory
    // Pre: subDir must be empty or ended with "/"
    inline void saveRandomStatus(G4String subDir);

};

#endif
