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
// $Id: SILICASteppingAction.cc 75292 2013-10-30 09:25:15Z gcosmo $
//
/// \file optical/wls/src/SILICASteppingAction.cc
/// \brief Implementation of the SILICASteppingAction class
//
//
#include "G4Run.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"

#include "SILICASteppingAction.hh"
#include "SILICADetectorConstruction.hh"
#include "SILICAEventAction.hh"
#include "SILICASteppingActionMessenger.hh"
#include "SILICAPhotonDetSD.hh"

#include "G4ParticleTypes.hh"

#include "SILICAUserTrackInformation.hh"

#include "G4ProcessManager.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"

#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include <sstream>

// Purpose: Save relevant information into User Track Information

static const G4ThreeVector ZHat = G4ThreeVector(0.0,0.0,1.0);

G4int SILICASteppingAction::fMaxRndmSave = 10000;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SILICASteppingAction::SILICASteppingAction(SILICADetectorConstruction* detector,  SILICAEventAction* event)
: fDetector(detector),fEventAction(event)
{
    fSteppingMessenger = new SILICASteppingActionMessenger(this);
    
    fCounterEnd = 0;
    fCounterMid = 0;
    fBounceLimit = 1000000;
    
    fOpProcess = NULL;
    
    ResetCounters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SILICASteppingAction::~SILICASteppingAction()
{
    delete fSteppingMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  SILICASteppingAction::SetBounceLimit(G4int i)   {fBounceLimit = i;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int SILICASteppingAction::GetNumberOfBounces()      {return fCounterBounce;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int SILICASteppingAction::GetNumberOfClad1Bounces() {return fCounterClad1Bounce;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int SILICASteppingAction::GetNumberOfClad2Bounces() {return fCounterClad2Bounce;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int SILICASteppingAction::GetNumberOfSILICABounces()   {return fCounterSILICABounce;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int SILICASteppingAction::ResetSuccessCounter()     {
    G4int temp = fCounterEnd; fCounterEnd = 0; return temp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void SILICASteppingAction::saveRandomStatus(G4String subDir) 
// save the random status into a sub-directory
// Pre: subDir must be empty or ended with "/"
{
    
    // don't save if the maximum amount has been reached
    if (SILICASteppingAction::fMaxRndmSave == 0) return;
    
    G4RunManager* theRunManager = G4RunManager::GetRunManager();
    G4String randomNumberStatusDir = theRunManager->GetRandomNumberStoreDir();
    
    G4String fileIn  = randomNumberStatusDir + "currentEvent.rndm";
    
    std::ostringstream os;
    
    os << "run" << theRunManager->GetCurrentRun()->GetRunID() << "evt"
    << theRunManager->GetCurrentEvent()->GetEventID() << ".rndm" << '\0';
    
    G4String fileOut = randomNumberStatusDir + subDir + os.str();
    
    G4String copCmd = "/control/shell cp "+fileIn+" "+fileOut;
    G4UImanager::GetUIpointer()->ApplyCommand(copCmd);
    
    SILICASteppingAction::fMaxRndmSave--;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICASteppingAction::UserSteppingAction(const G4Step* theStep)
{
    G4StepPoint* thePrePoint  = theStep->GetPreStepPoint();
    G4StepPoint* thePostPoint = theStep->GetPostStepPoint();
    
    G4VPhysicalVolume* thePrePV  = thePrePoint->GetPhysicalVolume();
    G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();
    
    G4String thePrePVname  = " ";
    G4String thePostPVname = " ";
    
    if (thePostPV) {
        thePrePVname  = thePrePV->GetName();
        thePostPVname = thePostPV->GetName();
    }
    ////////newadd
    // collect energy and track length step by step
    G4double edep = 0;
    
    //should make for inside Scintillator only
    G4double stepl = 0.;
    if (theStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0. && (thePrePVname=="FiberCore" || thePrePVname=="FiberCladding")  ){
        stepl = theStep->GetStepLength();
        edep = theStep->GetTotalEnergyDeposit();
       fEventAction->AddDep(edep,stepl);
    }
    ////////////
    G4Track* theTrack = theStep->GetTrack();
    SILICAUserTrackInformation* trackInformation
    = (SILICAUserTrackInformation*)theTrack->GetUserInformation();
    
    //newadd4 for counting scintillator photon and cherenkov photon
    //this loop over all step
    G4RunManager* theRunManager = G4RunManager::GetRunManager();
    G4String ParticleName = theTrack->GetDynamicParticle()->
    GetParticleDefinition()->GetParticleName();
    //if (ParticleName == "opticalphoton") return;
    
    const std::vector<const G4Track*>* secondaries =
    theStep->GetSecondaryInCurrentStep();
    G4double enscintth = 0;
    G4double encherenth = 0;
    if (ParticleName != "opticalphoton" && secondaries->size()>0) {
        for(unsigned int i=0; i<secondaries->size(); ++i) {
            if (secondaries->at(i)->GetParentID()>0) {
                if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()
                   == G4OpticalPhoton::OpticalPhotonDefinition()){
                    if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
                        == "Scintillation"){fScintillationCounter++;
                        enscintth = secondaries->at(i)->GetKineticEnergy();
                    }
                    if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
                        == "Cerenkov"){fCerenkovCounter++;
                        encherenth = secondaries->at(i)->GetKineticEnergy();
                    }
                }
            }
        }
         G4cout << "event "<<theRunManager->GetCurrentEvent()->GetEventID()<<" Number of scintillation photons produces "<< fScintillationCounter <<" cherenkov "<< fCerenkovCounter << G4endl;
        fEventAction->AddNoPhotonType(fScintillationCounter, fCerenkovCounter);//newadd6
        fEventAction->AddPhotonEnergy(enscintth, encherenth);//newadd6
    }
    //end newadd4
    //fEventAction->AddPriTrackNophotonScint(fScintillationCounter);//newadd5
    //fEventAction->AddPriTrackNophotonCheren(fCerenkovCounter);//newadd5
    

    
    //Recording data for start
    //to retrive primary track information
    
    if (theTrack->GetParentID()==0) {
        //This is a primary track
        if ( theTrack->GetCurrentStepNumber() == 1 ) {
            //newadd
            G4double xval  = theTrack->GetVertexPosition().x();
            fEventAction->AddPriTrackx(xval);
            G4double yval  = theTrack->GetVertexPosition().y();
            fEventAction->AddPriTracky(yval);
            G4double zval  = theTrack->GetVertexPosition().z();
            fEventAction->AddPriTrackz(zval);
            G4double pzval = theTrack->GetVertexMomentumDirection().z();
            fEventAction->AddPriTrackpz(pzval);
            G4double fInitThetaval =
            theTrack->GetVertexMomentumDirection().angle(ZHat);
            fEventAction->AddPriTracktheta(fInitThetaval);
            
        }
    }
    
    // Retrieve the status of the photon
    G4OpBoundaryProcessStatus theStatus = Undefined;
    
    //is this to look optical photon only?
    G4ProcessManager* OpManager =
    G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
    
    if (OpManager) {
        G4int MAXofPostStepLoops =
        OpManager->GetPostStepProcessVector()->entries();
        G4ProcessVector* fPostStepDoItVector =
        OpManager->GetPostStepProcessVector(typeDoIt);
        
        for ( G4int i=0; i<MAXofPostStepLoops; i++) {
            G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
            fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
            if (fOpProcess) { theStatus = fOpProcess->GetStatus(); break;}
        }
    }
    
    // Find the skewness of the ray at first change of boundary
    //turn off isStatus(InsideOfFiber) doesn't change the tracking
    if ( fInitGamma == -1 &&
        (theStatus == TotalInternalReflection
         || theStatus == FresnelReflection
         || theStatus == FresnelRefraction)
        && trackInformation->isStatus(InsideOfFiber) ) {
        
        G4double px = theTrack->GetVertexMomentumDirection().x();
        G4double py = theTrack->GetVertexMomentumDirection().y();
        G4double x  = theTrack->GetPosition().x();
        G4double y  = theTrack->GetPosition().y();
        
        fInitGamma = x * px + y * py;
        
        fInitGamma =
        fInitGamma / std::sqrt(px*px + py*py) / std::sqrt(x*x + y*y);
        
        fInitGamma = std::acos(fInitGamma*rad);
        
        if ( fInitGamma / deg > 90.0)  { fInitGamma = 180 * deg - fInitGamma;}
    }
    // Record Photons that missed the photon detector but escaped from readout
    if ( !thePostPV && trackInformation->isStatus(EscapedFromReadOut) ) {
        //     UpdateHistogramSuccess(thePostPoint,theTrack);
        ResetCounters();
        
        return;
    }
    
    // Assumed photons are originated at the fiber OR
    // the fiber is the first material the photon hits
    switch (theStatus) {
            //boundary for photons from scintillation
            
            // Exiting the fiber
        case FresnelRefraction:
        case SameMaterial:
            
            G4bool isFiber;
            /*isFiber = thePostPVname == "SILICAFiber"
            || thePostPVname == "Clad1"
            || thePostPVname == "Clad2" ||thePostPVname == "Hole" || thePostPVname =="Scintillator" ;*/
            
            isFiber = thePostPVname == "FiberCladding"
            || thePostPVname == "FiberCore";
            if ( isFiber ) {
                
                if (trackInformation->isStatus(OutsideOfFiber))
                    trackInformation->AddStatusFlag(InsideOfFiber);
                
                // Set the Exit flag when the photon refracted out of the fiber
            } else if (trackInformation->isStatus(InsideOfFiber)) {
                
                // EscapedFromReadOut if the z position is the same as fiber's end
                if (theTrack->GetPosition().z() == fDetector->GetSILICAFiberEnd())
                {
                    trackInformation->AddStatusFlag(EscapedFromReadOut);
                    fCounterEnd++;
                }
                else // Escaped from side
                {
                    trackInformation->AddStatusFlag(EscapedFromSide);
                    trackInformation->SetExitPosition(theTrack->GetPosition());
                    
                    //              UpdateHistogramEscape(thePostPoint,theTrack);
                    
                    fCounterMid++;
                    ResetCounters();
                }
                
                trackInformation->AddStatusFlag(OutsideOfFiber);
                trackInformation->SetExitPosition(theTrack->GetPosition());
                
            }
            
            return;
            
            // Internal Reflections
        case TotalInternalReflection:
            
            // Kill the track if it's number of bounces exceeded the limit
            if (fBounceLimit > 0 && fCounterBounce >= fBounceLimit)
            {
                theTrack->SetTrackStatus(fStopAndKill);
                trackInformation->AddStatusFlag(murderee);
                ResetCounters();
                G4cout << "\n Bounce Limit Exceeded" << G4endl;
                return;
            }
            
        case FresnelReflection:
            
            fCounterBounce++;
            
            if ( thePrePVname == "FiberCore") fCounterSILICABounce++;
            
            else if ( thePrePVname == "FiberCladding") fCounterClad1Bounce++;
            
            
            // Determine if the photon has reflected off the read-out end
            if (theTrack->GetPosition().z() == fDetector->GetSILICAFiberEnd())
            {
                if (!trackInformation->isStatus(ReflectedAtReadOut) &&
                    trackInformation->isStatus(InsideOfFiber))
                {
                    trackInformation->AddStatusFlag(ReflectedAtReadOut);
                    
                    if (fDetector->IsPerfectFiber() &&
                        theStatus == TotalInternalReflection)
                    {
                        theTrack->SetTrackStatus(fStopAndKill);
                        trackInformation->AddStatusFlag(murderee);
                        //                UpdateHistogramReflect(thePostPoint,theTrack);
                        ResetCounters();
                        return;
                    }
                }
            }
            return;
            
            // Reflection of the mirror
        case LambertianReflection:
        case LobeReflection:
        case SpikeReflection:
            
            // Check if it hits the mirror
            if ( thePostPVname == "Mirror" )
                trackInformation->AddStatusFlag(ReflectedAtMirror);
            
            return;
            
            // Detected by a detector
        case Detection:
            
            // Check if the photon hits the detector and process the hit if it does
            if ( thePostPVname == "PhotonDet" ) {
                
                G4SDManager* SDman = G4SDManager::GetSDMpointer();
                G4String SDname="SILICA/PhotonDet";
                SILICAPhotonDetSD* mppcSD =
                (SILICAPhotonDetSD*)SDman->FindSensitiveDetector(SDname);
                
                if (mppcSD) mppcSD->ProcessHits_constStep(theStep,NULL);
                
                // Record Photons that escaped at the end
                //          if (trackInformation->isStatus(EscapedFromReadOut))
                //                              UpdateHistogramSuccess(thePostPoint,theTrack);
                
                // Stop Tracking when it hits the detector's surface
                ResetCounters();
                theTrack->SetTrackStatus(fStopAndKill);
                
                return;
            }
            
            break;
            
        default: break;
            
    }//end of switch status
    //G4cout << "event "<<theRunManager->GetCurrentEvent()->GetEventID()<<" Number of photons pass2end "<< fCounterEnd <<" fail2end "<< fCounterMid << G4endl;
    //add information of photon pass or fail to reach sensitive detector
    //fEventAction->AddPriTrackNophotonpass2end(fCounterEnd);//newadd2
    //fEventAction->AddPriTrackNophotonfail2end(fCounterMid);//newadd2
    if (theStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.){
        fEventAction->AddNoPhotonPass2end(fCounterEnd);//newadd6
        fEventAction->AddNoPhotonFail2end(fCounterMid);//newadd6
    }
    // Check for absorbed photons
    if (theTrack->GetTrackStatus() != fAlive  &&
        trackInformation->isStatus(InsideOfFiber))
    {
        //     UpdateHistogramAbsorb(thePostPoint,theTrack);
        ResetCounters();
        return;
    }
    
 
}
