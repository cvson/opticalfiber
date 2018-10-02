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
// $Id: SILICARunActionMessenger.hh 69561 2013-05-08 12:25:56Z gcosmo $
//
/// \file optical/wls/include/SILICARunActionMessenger.hh
/// \brief Definition of the SILICARunActionMessenger class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef SILICARunActionMessenger_h
#define SILICARunActionMessenger_h 1

#include "globals.hh"

#include "G4UImessenger.hh"

class SILICARunAction;

class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;

class SILICARunActionMessenger : public G4UImessenger
{
  public:

    SILICARunActionMessenger(SILICARunAction* );
    virtual ~SILICARunActionMessenger();

    virtual void SetNewValue(G4UIcommand* ,G4String );

  private:

    SILICARunAction*              fRunAction;

    G4UIdirectory*             fRndmDir;
    G4UIcmdWithAnInteger*      fRndmSaveCmd;
    G4UIcmdWithAString*        fRndmReadCmd;
    G4UIcmdWithABool*          fSetAutoSeedCmd;

};

#endif
