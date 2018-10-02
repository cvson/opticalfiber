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
// $Id: SILICAActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file optical/wls/src/SILICAActionInitialization.cc
/// \brief Implementation of the SILICAActionInitialization class

#include "SILICAActionInitialization.hh"
#include "SILICADetectorConstruction.hh"
#include "SILICAHistoManager.hh"//newadd
#include "SILICAPrimaryGeneratorAction.hh"

#include "SILICARunAction.hh"
#include "SILICAEventAction.hh"
#include "SILICATrackingAction.hh"
#include "SILICASteppingAction.hh"
#include "SILICAStackingAction.hh"
#include "SILICASteppingVerbose.hh"
#include "G4GeneralParticleSource.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SILICAActionInitialization::SILICAActionInitialization(SILICADetectorConstruction* det)
 : G4VUserActionInitialization(), fDetector(det)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SILICAActionInitialization::~SILICAActionInitialization()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICAActionInitialization::BuildForMaster() const
{
  // Histo manager
  SILICAHistoManager*  histo = new SILICAHistoManager();//newadd
  //SetUserAction(new SILICARunAction());
   SetUserAction(new SILICARunAction(histo));//newadd
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICAActionInitialization::Build() const
{
   // Histo manager
  SILICAHistoManager*  histo = new SILICAHistoManager();//newadd
  SetUserAction(new SILICAPrimaryGeneratorAction(fDetector));

  //SILICARunAction* runAction = new SILICARunAction();
  //SILICAEventAction* eventAction = new SILICAEventAction(runAction);
 SILICARunAction* runAction = new SILICARunAction(histo);
 SILICAEventAction* eventAction = new SILICAEventAction(runAction,histo);

  SetUserAction(runAction);
  SetUserAction(eventAction);
  SetUserAction(new SILICATrackingAction());
  //SetUserAction(new SILICASteppingAction(fDetector));
  SetUserAction(new SILICASteppingAction(fDetector,eventAction));//newadd
 // SetUserAction(new SILICAStackingAction());
    SetUserAction(new SILICAStackingAction(histo));
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VSteppingVerbose* SILICAActionInitialization::InitializeSteppingVerbose() const
{
  return new SILICASteppingVerbose();
}
