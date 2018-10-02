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
// $Id: SILICADetectorConstruction.cc 84718 2014-10-20 07:40:45Z gcosmo $
//
/// \file optical/wls/src/SILICADetectorConstruction.cc
/// \brief Implementation of the SILICADetectorConstruction class
//
//
#include "G4ios.hh"
#include "globals.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4RunManager.hh"

#include "SILICADetectorConstruction.hh"
#include "SILICADetectorMessenger.hh"
#include "SILICAMaterials.hh"
#include "SILICAPhotonDetSD.hh"

#include "G4UserLimits.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
//to include color
#include "G4VisAttributes.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SILICADetectorConstruction::SILICADetectorConstruction()
: fMaterials(NULL), fLogicHole(NULL), fLogicWorld(NULL),
fPhysiWorld(NULL), fPhysiHole(NULL)
{
    fDetectorMessenger = new SILICADetectorMessenger(this);
    
    fNumOfCladLayers = 0;//0 default
    
    fSurfaceRoughness = 1;
    
    fMirrorToggle = true;
    fMirrorPolish = 1.;
    fMirrorReflectivity = 1.;
    
    fMPPCPolish = 1.;
    fMPPCReflectivity = 0.;
    
    fExtrusionPolish = 1.;
    fExtrusionReflectivity = 1.;
    
    fXYRatio = 1.0;
    
    fSILICAfiberZ     = 0.65*m;
    fSILICAfiberRY  = 0.5*mm;//0.5*mm;
    fSILICAfiberOrigin = 0.0;
    
    fMPPCShape = "Circle";
    fMPPCHalfL = fSILICAfiberRY;
    fMPPCDist  = 0.00*mm;
    fMPPCTheta = 0.0*deg;
    fMPPCZ     = 0.05*mm;
    
    fClrfiberZ  = fMPPCZ + 10.*nm;
    fMirrorZ    = 0.1*mm;
    
    fBarLength        = 0.6*m;//New change half of length, INGRID type 0.2
    fBarBase          = 9.6*mm;
    fBarWidth	    = 25*mm;//width 5cm
    fBarThick	    = 4.91*mm;//thickness 1cm
    fHoleRadius       = 0.9*mm;
    fHoleLength       = fSILICAfiberZ;//fBarLength;
    fCoatingThickness = 0.25*mm;
    fCoatingRadius    = 1.875*mm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SILICADetectorConstruction::~SILICADetectorConstruction()
{
    if (fDetectorMessenger) delete fDetectorMessenger;
    if (fMaterials)         delete fMaterials;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* SILICADetectorConstruction::Construct()
{
    if (fPhysiWorld) {
        G4GeometryManager::GetInstance()->OpenGeometry();
        G4PhysicalVolumeStore::GetInstance()->Clean();
        G4LogicalVolumeStore::GetInstance()->Clean();
        G4SolidStore::GetInstance()->Clean();
        G4LogicalSkinSurface::CleanSurfaceTable();
        G4LogicalBorderSurface::CleanSurfaceTable();
    }
    
    fMaterials = SILICAMaterials::GetInstance();
    
    UpdateGeometryParameters();
    
    return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* SILICADetectorConstruction::ConstructDetector()
{
    //--------------------------------------------------
    // World
    //--------------------------------------------------
    
    G4VSolid* solidWorld =
    new G4Box("World", fWorldSizeX, fWorldSizeY, fWorldSizeZ);
    
    fLogicWorld = new G4LogicalVolume(solidWorld,
                                      FindMaterial("G4_AIR"),
                                      "World");
    
    fPhysiWorld = new G4PVPlacement(0,
                                    G4ThreeVector(),
                                    fLogicWorld,
                                    "World",
                                    0,
                                    false,
                                    0);
    
    
    //--------------------------------------------------
    // Fiber
    //--------------------------------------------------
    //Optical properties for the vacuum and aluminum
    // Boundary Surface Properties
    /*G4OpticalSurface* opSurface = NULL;
    if (fSurfaceRoughness < 1.)
        opSurface = new G4OpticalSurface("RoughSurface",          // Surface Name
                                         glisur,                  // SetModel
                                         ground,                  // SetFinish
                                         dielectric_dielectric,   // SetType
                                         fSurfaceRoughness);      // SetPolish
    */
    double fiberLength = fSILICAfiberZ*2.;
    G4Tubs *cladding = new G4Tubs("Cladding",0.8*mm*0.5,0.83*mm*0.5,fiberLength*0.5,0.*deg, 360.*deg);
    
    G4Tubs *core = new G4Tubs("Core",0.0,0.8*mm*0.5,fiberLength*0.5,0.*deg, 360.*deg);
    G4LogicalVolume* logicClad = new G4LogicalVolume( cladding, FindMaterial("TECS"), "FiberCladding", 0,0,0);
    G4VisAttributes* vscint_intVisAtt1 = new G4VisAttributes(G4Color(0.,0.0,1.));//blue
    vscint_intVisAtt1->SetForceSolid(true);
    logicClad->SetVisAttributes(vscint_intVisAtt1);
    new G4PVPlacement(0,
                      G4ThreeVector(0.0,0.0,0.0),
                      logicClad,
                      "FiberCladding",
                      fLogicWorld,
                      false,
                      0);
    
    G4LogicalVolume* logicCore = new G4LogicalVolume( core, FindMaterial("Silica"), "FiberCore", 0,0,0);
    G4VisAttributes* vscint_intVisAtt2 = new G4VisAttributes(G4Color(1.,1.,1.));//white
    vscint_intVisAtt2->SetForceSolid(true);
    logicCore->SetVisAttributes(vscint_intVisAtt2);
    new G4PVPlacement(0,
                      G4ThreeVector(0.0,0.0,0.0),
                      logicCore,
                      "FiberCore",
                      logicClad,
                      false,
                      0);
    
    //finish the fiber
    // Clear Fiber (Coupling Layer)
    G4VSolid* solidCouple = new G4Box("Couple",fCoupleRX,fCoupleRY,fCoupleZ);
    
    G4LogicalVolume*   logicCouple = new G4LogicalVolume(solidCouple,
                                                         FindMaterial("G4_AIR"),
                                                         "Couple");
    
    new G4PVPlacement(0,
                      G4ThreeVector(0.0,0.0,fCoupleOrigin),
                      logicCouple,
                      "Couple",
                      fLogicWorld,
                      false,
                      0);

    //--------------------------------------------------
    // A logical layer in front of PhotonDet
    //--------------------------------------------------
    
    // Purpose: Preventing direct dielectric to metal contact
    
    // Check for valid placement of PhotonDet
    if (fMPPCTheta > std::atan(fMPPCDist / fMPPCHalfL)) {
        
        fMPPCTheta = 0;
        fMPPCOriginX  = std::sin(fMPPCTheta) * (fMPPCDist + fClrfiberZ);
        fMPPCOriginZ  = -fCoupleZ+std::cos(fMPPCTheta)*(fMPPCDist+fClrfiberZ);
        G4cerr << "Invalid alignment.  Alignment Reset to 0" << G4endl;
    }
    
    // Clear Fiber (Coupling Layer)
    G4VSolid* solidClrfiber;
    
    if ( fMPPCShape == "Square" )
        solidClrfiber =
        new G4Box("ClearFiber",fClrfiberHalfL,fClrfiberHalfL,fClrfiberZ);
    else
        solidClrfiber =
        new G4Tubs("ClearFiber",0.,fClrfiberHalfL,fClrfiberZ,0.0*rad,twopi*rad);
    
    G4LogicalVolume*   logicClrfiber =
    new G4LogicalVolume(solidClrfiber,
                        FindMaterial("G4_AIR"),
                        "ClearFiber");
    
    new G4PVPlacement(new G4RotationMatrix(CLHEP::HepRotationY(-fMPPCTheta)),
                      G4ThreeVector(fMPPCOriginX,0.0,fMPPCOriginZ),
                      logicClrfiber,
                      "ClearFiber",
                      logicCouple,
                      false,
                      0);
    
    //should clear the fiber at the end
    //--------------------------------------------------
    // PhotonDet (Sensitive Detector)
    //--------------------------------------------------
    
    // Physical Construction
    G4VSolid* solidPhotonDet;
    
    if ( fMPPCShape == "Square" )
        solidPhotonDet = new G4Box("PhotonDet",fMPPCHalfL,fMPPCHalfL,fMPPCZ);
    else
        solidPhotonDet =
        new G4Tubs("PhotonDet",0.,fMPPCHalfL,fMPPCZ,0.0*rad,twopi*rad);
    
    G4LogicalVolume*   logicPhotonDet =
    new G4LogicalVolume(solidPhotonDet,
                        FindMaterial("G4_Al"),
                        "PhotonDet_LV");
    G4VisAttributes* mppcVisAtt = new G4VisAttributes(G4Color(0.7,0.,0.7)); // magenta
    mppcVisAtt->SetForceSolid(true);
    logicPhotonDet->SetVisAttributes(mppcVisAtt);
    
    new G4PVPlacement(0,
                      G4ThreeVector(0.0,0.0,0.0),
                      logicPhotonDet,
                      "PhotonDet",
                      logicClrfiber,
                      false,
                      0);
    
    // PhotonDet Surface Properties
    G4OpticalSurface* photonDetSurface = new G4OpticalSurface("PhotonDetSurface",
                                                              glisur,
                                                              ground,
                                                              dielectric_metal,
                                                              fMPPCPolish);
    
    G4MaterialPropertiesTable* photonDetSurfaceProperty =
    new G4MaterialPropertiesTable();
    
    G4double p_mppc[] = {2.00*eV, 3.47*eV};
    const G4int nbins = sizeof(p_mppc)/sizeof(G4double);
    G4double refl_mppc[] = {fMPPCReflectivity,fMPPCReflectivity};
    assert(sizeof(refl_mppc) == sizeof(p_mppc));
    G4double effi_mppc[] = {1, 1};
    assert(sizeof(effi_mppc) == sizeof(p_mppc));
    
    photonDetSurfaceProperty->AddProperty("REFLECTIVITY",p_mppc,refl_mppc,nbins);
    photonDetSurfaceProperty->AddProperty("EFFICIENCY",p_mppc,effi_mppc,nbins);
    
    photonDetSurface->SetMaterialPropertiesTable(photonDetSurfaceProperty);
    
    new G4LogicalSkinSurface("PhotonDetSurface",logicPhotonDet,photonDetSurface);
    
    //--------------------------------------------------
    // End of Construction
    //--------------------------------------------------
    
    return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::ConstructSDandField()
{
    if (!fmppcSD.Get()) {
        G4String mppcSDName = "SILICA/PhotonDet";
        SILICAPhotonDetSD* mppcSD = new SILICAPhotonDetSD(mppcSDName);
        fmppcSD.Put(mppcSD);
    }
    SetSensitiveDetector("PhotonDet_LV", fmppcSD.Get(), true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::UpdateGeometryParameters()
{
    fSILICAfiberRX  = fXYRatio * fSILICAfiberRY;
    
    fClad1RX = fSILICAfiberRX + 0.03*fSILICAfiberRX;
    fClad1RY = fSILICAfiberRY + 0.03*fSILICAfiberRY;
    fClad1Z  = fSILICAfiberZ;
    
    fClad2RX = fClad1RX + 0.03*fSILICAfiberRX;
    fClad2RY = fClad1RY + 0.03*fSILICAfiberRY;
    fClad2Z  = fSILICAfiberZ;
    
    fWorldSizeX = fClad2RX   + fMPPCDist + fMPPCHalfL + 5.*cm;//1.*cm;
    fWorldSizeY = fClad2RY   + fMPPCDist + fMPPCHalfL + 1.*cm;
    fWorldSizeZ = fSILICAfiberZ + fMPPCDist + fMPPCHalfL + 5.*cm;//1.0*cm
    
    fCoupleRX = fWorldSizeX;
    fCoupleRY = fWorldSizeY;
    fCoupleZ  = (fWorldSizeZ - fSILICAfiberZ) / 2;
    
    fClrfiberHalfL = fMPPCHalfL;
    
    fMirrorRmax = fClad2RY;
    
    fCoupleOrigin = fSILICAfiberOrigin + fSILICAfiberZ + fCoupleZ;
    fMirrorOrigin = fSILICAfiberOrigin - fSILICAfiberZ - fMirrorZ;
    fMPPCOriginX  = std::sin(fMPPCTheta) * (fMPPCDist + fClrfiberZ);
    fMPPCOriginZ  = -fCoupleZ + std::cos(fMPPCTheta) * (fMPPCDist + fClrfiberZ);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RotationMatrix
SILICADetectorConstruction::StringToRotationMatrix(G4String rotation)
{
    // We apply successive rotations OF THE OBJECT around the FIXED
    // axes of the parent's local coordinates; rotations are applied
    // left-to-right (rotation="r1,r2,r3" => r1 then r2 then r3).
    
    G4RotationMatrix rot;
    
    unsigned int place = 0;
    
    while (place < rotation.size()) {
        
        G4double angle;
        char* p;
        
        const G4String tmpstring=rotation.substr(place+1);
        
        angle = strtod(tmpstring.c_str(),&p) * deg;
        
        if (!p || (*p != (char)',' && *p != (char)'\0')) {
            G4cerr << "Invalid rotation specification: " <<
            rotation.c_str() << G4endl;
            return rot;
        }
        
        G4RotationMatrix thisRotation;
        
        switch(rotation.substr(place,1).c_str()[0]) {
            case 'X': case 'x':
                thisRotation = G4RotationMatrix(CLHEP::HepRotationX(angle));
                break;
            case 'Y': case 'y':
                thisRotation = G4RotationMatrix(CLHEP::HepRotationY(angle));
                break;
            case 'Z': case 'z':
                thisRotation = G4RotationMatrix(CLHEP::HepRotationZ(angle));
                break;
            default:
                G4cerr << " Invalid rotation specification: "
                << rotation << G4endl;
                return rot;
        }
        
        rot = thisRotation * rot;
        place = rotation.find(',',place);
        if (place > rotation.size()) break;
        ++place;
    }
    
    return rot;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetPhotonDetGeometry (G4String shape)
// Set the Geometry of the PhotonDet detector
// Pre:  shape must be either "Circle" and "Square"
{
    if (shape == "Circle" || shape == "Square" ) fMPPCShape = shape;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetNumberOfCladding(G4int num)
// Set the number of claddings
// Pre: 0 <= num <= 2
{
    fNumOfCladLayers = num;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetSILICALength (G4double length)
// Set the TOTAL length of the SILICA fiber
{
    fSILICAfiberZ = length;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetSILICARadius (G4double radius)
// Set the Y radius of SILICA fiber
{
    fSILICAfiberRY = radius;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetClad1Radius (G4double radius)
// Set the Y radius of Cladding 1
{
    fClad1RY = radius;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetClad2Radius (G4double radius)
// Set the Y radius of Cladding 2
{
    fClad2RY = radius;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetPhotonDetHalfLength(G4double halfL)
// Set the half length of the PhotonDet detector
// The half length will be the radius if PhotonDet is circular
{
    fMPPCHalfL = halfL;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetGap (G4double gap)
// Set the distance between fiber end and PhotonDet
{ 
    fMPPCDist = gap;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetPhotonDetAlignment(G4double theta)
// Set the Aligment of PhotonDet with respect to the z axis
// If theta is 0 deg, then the detector is perfectly aligned
// PhotonDet will be deviated by theta from z axis
// facing towards the center of the fiber
{
    fMPPCTheta = theta;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetSurfaceRoughness(G4double roughness)
// Set the Surface Roughness between Cladding 1 and SILICA fiber
// Pre: 0 < roughness <= 1
{
    fSurfaceRoughness = roughness;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetMirrorPolish(G4double polish)
// Set the Polish of the mirror, polish of 1 is a perfect mirror surface
// Pre: 0 < polish <= 1
{
    fMirrorPolish = polish;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetMirrorReflectivity(G4double reflectivity)
// Set the Reflectivity of the mirror, reflectivity of 1 is a perfect mirror
// Pre: 0 < reflectivity <= 1
{
    fMirrorReflectivity = reflectivity;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetPhotonDetPolish(G4double polish)
// Set the Polish of the PhotonDet, polish of 1 is a perfect mirror surface
// Pre: 0 < polish <= 1
{
    fMPPCPolish = polish;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetPhotonDetReflectivity(G4double reflectivity)
// Set the Reflectivity of the PhotonDet, reflectivity of 1 is a perfect mirror
// Pre: 0 < reflectivity <= 1
{
    fMPPCReflectivity = reflectivity;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetMirror(G4bool flag)
// Toggle to place the mirror or not at one end (-z end) of the fiber
// True means place the mirror, false means otherwise
{
    fMirrorToggle = flag;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetXYRatio(G4double r)
// Set the ratio of the x and y radius of the ellipse (x/y)
// a ratio of 1 would produce a circle
{
    fXYRatio = r;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetBarLength (G4double length)
// Set the length of the scintillator bar
{
    fBarLength = length;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetBarBase (G4double side)
// Set the side of the scintillator bar
{
    fBarBase = side;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void SILICADetectorConstruction::SetBarWidth (G4double side)
// Set the side of the scintillator bar
{
    fBarWidth = side;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void SILICADetectorConstruction::SetBarThick (G4double side)
// Set the side of the scintillator bar
{
    fBarThick = side;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetHoleRadius (G4double radius)
// Set the radius of the fiber hole
{
    fHoleRadius = radius;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetCoatingThickness (G4double thick)
// Set thickness of the coating on the bars
{
    fCoatingThickness = thick;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICADetectorConstruction::SetCoatingRadius (G4double radius)
// Set inner radius of the corner bar coating
{
    fCoatingRadius = radius;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SILICADetectorConstruction::GetSILICAFiberLength() { return fSILICAfiberZ; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SILICADetectorConstruction::GetBarLength() { return fBarLength; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SILICADetectorConstruction::GetBarBase() { return fBarBase; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SILICADetectorConstruction::GetBarWidth() { return fBarWidth; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SILICADetectorConstruction::GetBarThick() { return fBarThick; }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SILICADetectorConstruction::GetHoleRadius() { return fHoleRadius; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SILICADetectorConstruction::GetHoleLength() { return fHoleLength; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SILICADetectorConstruction::GetFiberRadius() { return GetSILICAFiberRMax(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SILICADetectorConstruction::GetCoatingThickness()
{ return fCoatingThickness; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SILICADetectorConstruction::GetCoatingRadius() { return fCoatingRadius; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SILICADetectorConstruction::GetSILICAFiberEnd()
{
    return fSILICAfiberOrigin + fSILICAfiberZ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SILICADetectorConstruction::GetSILICAFiberRMax()
{
    if (fNumOfCladLayers == 2) return fClad2RY;
    if (fNumOfCladLayers == 1) return fClad1RY;
    return fSILICAfiberRY;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SILICADetectorConstruction::GetSurfaceRoughness()
{
    return fSurfaceRoughness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Return True if the fiber construction is ideal
G4bool SILICADetectorConstruction::IsPerfectFiber()
{
    return     fSurfaceRoughness == 1. && fXYRatio == 1.
    && (!fMirrorToggle    ||
        (fMirrorPolish    == 1. && fMirrorReflectivity == 1.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* SILICADetectorConstruction::FindMaterial(G4String name) {
    G4Material* material = G4Material::GetMaterial(name,true);
    return material;
}
