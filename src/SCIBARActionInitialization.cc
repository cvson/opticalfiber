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
// $Id: SCIBARActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file optical/wls/src/SCIBARActionInitialization.cc
/// \brief Implementation of the SCIBARActionInitialization class

#include "SCIBARActionInitialization.hh"
#include "SCIBARDetectorConstruction.hh"

#include "SCIBARPrimaryGeneratorAction.hh"

#include "SCIBARRunAction.hh"
#include "SCIBAREventAction.hh"
#include "SCIBARTrackingAction.hh"
#include "SCIBARSteppingAction.hh"
#include "SCIBARStackingAction.hh"
#include "SCIBARSteppingVerbose.hh"
#include "G4GeneralParticleSource.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SCIBARActionInitialization::SCIBARActionInitialization(SCIBARDetectorConstruction* det)
 : G4VUserActionInitialization(), fDetector(det)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SCIBARActionInitialization::~SCIBARActionInitialization()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARActionInitialization::BuildForMaster() const
{
  SetUserAction(new SCIBARRunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARActionInitialization::Build() const
{
  SetUserAction(new SCIBARPrimaryGeneratorAction(fDetector));

  SCIBARRunAction* runAction = new SCIBARRunAction();
  SCIBAREventAction* eventAction = new SCIBAREventAction(runAction);

  SetUserAction(runAction);
  SetUserAction(eventAction);
  SetUserAction(new SCIBARTrackingAction());
  SetUserAction(new SCIBARSteppingAction(fDetector));
  SetUserAction(new SCIBARStackingAction());
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VSteppingVerbose* SCIBARActionInitialization::InitializeSteppingVerbose() const
{
  return new SCIBARSteppingVerbose();
}
