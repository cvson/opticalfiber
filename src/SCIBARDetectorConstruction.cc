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
// $Id: SCIBARDetectorConstruction.cc 84718 2014-10-20 07:40:45Z gcosmo $
//
/// \file optical/wls/src/SCIBARDetectorConstruction.cc
/// \brief Implementation of the SCIBARDetectorConstruction class
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

#include "SCIBARDetectorConstruction.hh"
#include "SCIBARDetectorMessenger.hh"
#include "SCIBARMaterials.hh"
#include "SCIBARPhotonDetSD.hh"

#include "G4UserLimits.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
//to include color
#include "G4VisAttributes.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SCIBARDetectorConstruction::SCIBARDetectorConstruction()
  : fMaterials(NULL), fLogicHole(NULL), fLogicWorld(NULL),
    fPhysiWorld(NULL), fPhysiHole(NULL)
{
  fDetectorMessenger = new SCIBARDetectorMessenger(this);

  fNumOfCladLayers = 0;
 
  fSurfaceRoughness = 1;
 
  fMirrorToggle = true;
  fMirrorPolish = 1.;
  fMirrorReflectivity = 1.;

  fMPPCPolish = 1.;
  fMPPCReflectivity = 0.;

  fExtrusionPolish = 1.;
  fExtrusionReflectivity = 1.;
 
  fXYRatio = 1.0;

  fSCIBARfiberZ     = 0.65*m;
  fSCIBARfiberRY  = 0.5*mm;
  fSCIBARfiberOrigin = 0.0;
 
  fMPPCShape = "Circle";
  fMPPCHalfL = fSCIBARfiberRY;
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
  fHoleLength       = fSCIBARfiberZ;//fBarLength;
  fCoatingThickness = 0.25*mm;
  fCoatingRadius    = 1.875*mm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SCIBARDetectorConstruction::~SCIBARDetectorConstruction()
{
  if (fDetectorMessenger) delete fDetectorMessenger;
  if (fMaterials)         delete fMaterials;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* SCIBARDetectorConstruction::Construct()
{
  if (fPhysiWorld) {
     G4GeometryManager::GetInstance()->OpenGeometry();
     G4PhysicalVolumeStore::GetInstance()->Clean();
     G4LogicalVolumeStore::GetInstance()->Clean();
     G4SolidStore::GetInstance()->Clean();
     G4LogicalSkinSurface::CleanSurfaceTable();
     G4LogicalBorderSurface::CleanSurfaceTable();
  }

  fMaterials = SCIBARMaterials::GetInstance();

  UpdateGeometryParameters();

  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* SCIBARDetectorConstruction::ConstructDetector()
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
  //  Extrusion
  //--------------------------------------------------
  //newadd
  G4VSolid* solidExtrusion =
        new G4Box("Extrusion",GetBarWidth(),GetBarThick(),GetBarLength());
	//new G4Box("Extrusion",GetBarBase()/2,GetBarBase()/2,GetBarLength()/2);

  G4LogicalVolume* logicExtrusion =
                      new G4LogicalVolume(solidExtrusion,
                                          FindMaterial("Coating"),
                                          "Extrusion");
   G4VisAttributes* vref_intVisAtt = new G4VisAttributes(G4Color(1.0,0.0,0.0));//red
  vref_intVisAtt->SetForceSolid(true);
  //logicExtrusion->SetVisAttributes(vref_intVisAtt);

  G4OpticalSurface* TiO2Surface = new G4OpticalSurface("TiO2Surface",
                                                       glisur,
                                                       ground,
                                                       dielectric_metal,
                                                       fExtrusionPolish);

  G4MaterialPropertiesTable* TiO2SurfaceProperty =
                                             new G4MaterialPropertiesTable();

  G4double p_TiO2[] = {2.00*eV, 3.47*eV};
  const G4int nbins = sizeof(p_TiO2)/sizeof(G4double);

  G4double refl_TiO2[] = {fExtrusionReflectivity,fExtrusionReflectivity};
  assert(sizeof(refl_TiO2) == sizeof(p_TiO2));
  G4double effi_TiO2[] = {0, 0};
  assert(sizeof(effi_TiO2) == sizeof(p_TiO2));

  TiO2SurfaceProperty -> AddProperty("REFLECTIVITY",p_TiO2,refl_TiO2,nbins);
  TiO2SurfaceProperty -> AddProperty("EFFICIENCY",p_TiO2,effi_TiO2,nbins);

  TiO2Surface -> SetMaterialPropertiesTable(TiO2SurfaceProperty);

  new G4PVPlacement(0,
                    G4ThreeVector(),
                    logicExtrusion,
                    "Extrusion",
                    fLogicWorld,
                    false,
                    0);

  new G4LogicalSkinSurface("TiO2Surface",logicExtrusion,TiO2Surface);

  //--------------------------------------------------
  // Scintillator
  //--------------------------------------------------
   //newadd
  /*G4VSolid* solidScintillator = new G4Box("Scintillator",
                                GetBarBase()/2-GetCoatingThickness()
                                           -GetCoatingRadius(),
                                GetBarBase()/2-GetCoatingThickness()
                                           -GetCoatingRadius(),
                                GetBarLength()/2);*/
  G4VSolid* solidScintillator = new G4Box("Scintillator",
                                GetBarWidth()-GetCoatingThickness()
                                           -GetCoatingRadius(),
                                GetBarThick()-GetCoatingThickness()
                                           -GetCoatingRadius(),
                                GetBarLength());

  G4LogicalVolume* logicScintillator =
                             new G4LogicalVolume(solidScintillator,
                                                 FindMaterial("Polystyrene"),
                                                 "Scintillator");
   // G4VisAttributes* vscint_intVisAtt3 = new G4VisAttributes(G4Color(1.,0.0,1.));//magenta
  G4VisAttributes* vscint_intVisAtt3 = new G4VisAttributes(G4Color(1.,1.0,1.));//white  
//vscint_intVisAtt->SetForceWireframe(true);
  vscint_intVisAtt3->SetForceSolid(true);
  logicScintillator->SetVisAttributes(vscint_intVisAtt3);

  new G4PVPlacement(0,
                    G4ThreeVector(),
                    logicScintillator,
                    "Scintillator",
                    logicExtrusion,
                    false,
                    0);

  if (GetCoatingRadius() > 0.*mm) {
     /*G4VSolid* solidScintside = new G4Box("SideOfBar",
                                GetBarBase()/2-GetCoatingThickness()
                                           -GetCoatingRadius(),
                                GetCoatingRadius()/2,
                                GetBarLength()/2);*/
     //newadd
    
     G4VSolid* solidScintside = new G4Box("SideOfBar",
                                GetBarWidth()-GetCoatingThickness()
                                           -GetCoatingRadius(),
                                GetCoatingRadius()/2,
                                GetBarLength());
	//newadd, not width is not same as thickness, can't just rotate
     G4VSolid* solidScintsideS = new G4Box("SideOfBarS",
				GetCoatingRadius()/2,
                                GetBarThick()-GetCoatingThickness()
                                           -GetCoatingRadius(),
                                GetBarLength());

     
     G4VSolid* solidScintcrnr = new G4Tubs("CrnrOfBar",
                                 0.0*cm,
                                 GetCoatingRadius(),
                                 GetBarLength(),
                                 0.*deg,
                                 90.*deg);
     G4LogicalVolume* logicScintSide =
                             new G4LogicalVolume(solidScintside,
                                                 FindMaterial("Polystyrene"),
                                                 "SideOfBar");
// G4VisAttributes* vscint_intVisAtt1 = new G4VisAttributes(G4Color(0.,1.0,0.));//green
G4VisAttributes* vscint_intVisAtt1 = new G4VisAttributes(G4Color(0.5,0.5,0.5));//grey
  vscint_intVisAtt1->SetForceSolid(true);
  logicScintSide->SetVisAttributes(vscint_intVisAtt1);


 G4LogicalVolume* logicScintSideS =
                             new G4LogicalVolume(solidScintsideS,
                                                 FindMaterial("Polystyrene"),
                                                 "SideOfBarS");
 G4VisAttributes* vscint_intVisAtt10 = new G4VisAttributes(G4Color(0.5,0.5,0.5));//grey
  vscint_intVisAtt10->SetForceSolid(true);
  logicScintSideS->SetVisAttributes(vscint_intVisAtt10);


     G4LogicalVolume* logicScintCrnr =
                             new G4LogicalVolume(solidScintcrnr,
                                                 FindMaterial("Polystyrene"),
                                                 "CrnrOfBar");

   //G4VisAttributes* vscint_intVisAtt2 = new G4VisAttributes(G4Color(0.,0.0,1.));//blue
  G4VisAttributes* vscint_intVisAtt2 = new G4VisAttributes(G4Color(0.5,0.5,0.5));//grey
  //vscint_intVisAtt->SetForceWireframe(true);
  vscint_intVisAtt2->SetForceSolid(true);
  logicScintCrnr->SetVisAttributes(vscint_intVisAtt2);

     //G4double x = GetBarBase()/2-GetCoatingThickness()-GetCoatingRadius()/2;
     //G4double y = GetBarBase()/2-GetCoatingThickness()-GetCoatingRadius()/2;
	//newadd
	G4double x = GetBarWidth()-GetCoatingThickness()-GetCoatingRadius()/2;
        G4double y = GetBarThick()-GetCoatingThickness()-GetCoatingRadius()/2;

     new G4PVPlacement(0,
                       G4ThreeVector(0,-y,0),
                       logicScintSide,
                       "SideOfBar",
                       logicExtrusion,
                       false,
                       0);
     new G4PVPlacement(0,
                       G4ThreeVector(0, y,0),
                       logicScintSide,
                       "SideOfBar",
                       logicExtrusion,
                       false,
                       1);

    new G4PVPlacement(0,
                       G4ThreeVector(-x,0,0),
                       logicScintSideS,
                       "SideOfBarS",
                       logicExtrusion,
                       false,
                       0);
     new G4PVPlacement(0,
                       G4ThreeVector(x, 0,0),
                       logicScintSideS,
                       "SideOfBarS",
                       logicExtrusion,
                       false,
                       1);

     G4RotationMatrix* g4rot = new G4RotationMatrix();
     *g4rot = StringToRotationMatrix("Z90");
     *g4rot = g4rot->inverse();
     if (*g4rot == G4RotationMatrix()) g4rot = NULL;

     /*new G4PVPlacement(g4rot,
                       G4ThreeVector(x,0,0),
                       logicScintSide,
                       "SideOfBar",
                       logicExtrusion,
                       false,
                       2);
     new G4PVPlacement(g4rot,
                       G4ThreeVector(-x,0,0),
                       logicScintSide,
                       "SideOfBar",
                       logicExtrusion,
                       false,
                       3);
*/
     //x = GetBarBase()/2-GetCoatingThickness()-GetCoatingRadius();
     //y = GetBarBase()/2-GetCoatingThickness()-GetCoatingRadius();
	//newadd
	x = GetBarWidth()-GetCoatingThickness()-GetCoatingRadius();
        y = GetBarThick()-GetCoatingThickness()-GetCoatingRadius();

     new G4PVPlacement(0,
                       G4ThreeVector(x,y,0),
                       logicScintCrnr,
                       "CrnrOfBar",
                       logicExtrusion,
                       false,
                       0);

     new G4PVPlacement(g4rot,
                       G4ThreeVector(-x,y,0),
                       logicScintCrnr,
                       "CrnrOfBar",
                       logicExtrusion,
                       false,
                       1);

     g4rot = new G4RotationMatrix();
     *g4rot = StringToRotationMatrix("Z180");
     *g4rot = g4rot->inverse();
     if (*g4rot == G4RotationMatrix()) g4rot = NULL;

     new G4PVPlacement(g4rot,
                       G4ThreeVector(-x,-y,0),
                       logicScintCrnr,
                       "CrnrOfBar",
                       logicExtrusion,
                       false,
                       2);

     g4rot = new G4RotationMatrix();
     *g4rot = StringToRotationMatrix("Z270");
     *g4rot = g4rot->inverse();
     if (*g4rot == G4RotationMatrix()) g4rot = NULL;

     new G4PVPlacement(g4rot,
                       G4ThreeVector(x,-y,0),
                       logicScintCrnr,
                       "CrnrOfBar",
                       logicExtrusion,
                       false,
                       3);

  }

  //make a hole
  if (GetFiberRadius()<GetHoleRadius()) {

        G4VSolid* solidHole = new G4Tubs("Hole",
                                         0.0*cm,
                                         GetHoleRadius(),
                                         GetHoleLength(),
                                         0.*deg,
                                         360.*deg);
        fLogicHole = new G4LogicalVolume(solidHole,
                                         FindMaterial("G4_AIR"),
                                         "Hole");

        fPhysiHole = new G4PVPlacement(0,
                                       G4ThreeVector(),
                                       fLogicHole,
                                       "Hole",
                                       logicScintillator,
                                       false,
                                       0);
  }

  //--------------------------------------------------
  // Fiber
  //--------------------------------------------------

  ConstructFiber();

  //--------------------------------------------------
  // End of Construction
  //--------------------------------------------------

  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::ConstructFiber()
{
  if (!(fLogicHole) || !(fPhysiHole) ) {
     std::ostringstream o;
     o << "The Fiber Hole has not been constructed";
     G4Exception("SCIBARDetectorConstruction::ConstructFiber","",
                  FatalException,o.str().c_str());
  }

  // Pointers to the most recently constructed volume
  G4LogicalVolume* logicPlacement = fLogicHole;
  G4VPhysicalVolume* physiPlacement = fPhysiHole;

  //--------------------------------------------------
  // Fiber Construction
  //-------------------------------------------------- 

  // Boundary Surface Properties
  G4OpticalSurface* opSurface = NULL;
 
  if (fSurfaceRoughness < 1.)
     opSurface = new G4OpticalSurface("RoughSurface",          // Surface Name
                                      glisur,                  // SetModel
                                      ground,                  // SetFinish
                                      dielectric_dielectric,   // SetType
                                      fSurfaceRoughness);      // SetPolish

  G4LogicalVolume   *logicClad1, *logicClad2;
  G4VPhysicalVolume *physiClad1, *physiClad2;

  // Determine the number of cladding layers to be built
  switch ( fNumOfCladLayers ) {
 
    case 2:

     //--------------------------------------------------
     // Cladding 2
     //--------------------------------------------------

     G4VSolid* solidClad2;
 
     if (fXYRatio == 1.)
       solidClad2 = new G4Tubs("Clad2",0.,fClad2RX,fClad2Z,0.0*rad,twopi*rad);
     else
       solidClad2 = new G4EllipticalTube("Clad2",fClad2RX,fClad2RY,fClad2Z);

     logicClad2  = new G4LogicalVolume(solidClad2,
                                       FindMaterial("FPethylene"),
                                       "Clad2");

     physiClad2 = new G4PVPlacement(0,
                                    G4ThreeVector(0.0,0.0,fSCIBARfiberOrigin),
                                    logicClad2,
                                    "Clad2",
                                    logicPlacement,
                                    false,
                                    0);

     // Place the rough surface only if needed
     if (opSurface) {
       new G4LogicalBorderSurface("surfaceClad2Out",
                                  physiClad2,
                                  physiPlacement,
                                  opSurface);
       new G4LogicalBorderSurface("surfaceClad2In",
                                  physiPlacement,
                                  physiClad2,
                                  opSurface);
     }

     logicPlacement = logicClad2;
     physiPlacement = physiClad2;

    case 1:

     //--------------------------------------------------
     // Cladding 1
     //--------------------------------------------------

     G4VSolid* solidClad1;

     if (fXYRatio == 1.)
       solidClad1 = new G4Tubs("Clad1",0.,fClad1RX,fClad1Z,0.0*rad,twopi*rad);
     else
       solidClad1 = new G4EllipticalTube("Clad1",fClad1RX,fClad1RY,fClad1Z);

     logicClad1 = new G4LogicalVolume(solidClad1,
                                      FindMaterial("Pethylene"),
                                      "Clad1");

     physiClad1 = new G4PVPlacement(0,
                                    G4ThreeVector(0.0,0.0,fSCIBARfiberOrigin),
                                    logicClad1,
                                    "Clad1",
                                    logicPlacement,
                                    false,
                                    0);

     // Place the rough surface only if needed
     if (opSurface) {
       new G4LogicalBorderSurface("surfaceClad1Out",
                                  physiClad1,
                                  physiPlacement,
                                  opSurface);
       new G4LogicalBorderSurface("surfaceClad1In",
                                  physiPlacement,
                                  physiClad1,
                                  opSurface);
     }

     logicPlacement = logicClad1;
     physiPlacement = physiClad1;

    default:

     //--------------------------------------------------
     // SCIBAR Fiber
     //--------------------------------------------------

     G4VSolid* solidSCIBARfiber;

     if (fXYRatio == 1.)
       solidSCIBARfiber =
           new G4Tubs("SCIBARFiber",0.,fSCIBARfiberRX,fSCIBARfiberZ,0.0*rad,twopi*rad);
     else
       solidSCIBARfiber =
           new G4EllipticalTube("SCIBARFiber",fSCIBARfiberRX,fSCIBARfiberRY,fSCIBARfiberZ);

     G4LogicalVolume*   logicSCIBARfiber =
                                     new G4LogicalVolume(solidSCIBARfiber,
                                                         FindMaterial("PMMA"),
                                                         "SCIBARFiber");

     logicSCIBARfiber->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,10*ms));

	      G4VisAttributes* fiberVisAtt = new G4VisAttributes(G4Color(0.,1.,0.)); // green
  fiberVisAtt->SetForceSolid(true);
  logicSCIBARfiber->SetVisAttributes(fiberVisAtt);

     G4VPhysicalVolume* physiSCIBARfiber = new G4PVPlacement(0,
                                       G4ThreeVector(0.0,0.0,fSCIBARfiberOrigin),
                                       logicSCIBARfiber,
                                       "SCIBARFiber",
                                       logicPlacement,
                                       false,
                                       0);

     // Place the rough surface only if needed
     if (opSurface) {
       new G4LogicalBorderSurface("surfaceSCIBAROut",
                                  physiSCIBARfiber,
                                  physiPlacement,
                                  opSurface);
       new G4LogicalBorderSurface("surfaceSCIBARIn",
                                  physiPlacement,
                                  physiSCIBARfiber,
                                  opSurface);
     }
  }

  //--------------------------------------------------
  // Mirror for reflection at one of the end
  //--------------------------------------------------

  // Place the mirror only if the user wants the mirror
  if (fMirrorToggle) {  

     G4VSolid* solidMirror = new G4Box("Mirror",
                                       fMirrorRmax,
                                       fMirrorRmax,
                                       fMirrorZ);
 
     G4LogicalVolume* logicMirror = new G4LogicalVolume(solidMirror,
                                                        FindMaterial("G4_Al"),
                                                        "Mirror");

     G4OpticalSurface* mirrorSurface = new G4OpticalSurface("MirrorSurface",
                                                             glisur,
                                                             ground,
                                                             dielectric_metal,
                                                             fMirrorPolish);

     G4MaterialPropertiesTable* mirrorSurfaceProperty =
                                              new G4MaterialPropertiesTable();

     G4double p_mirror[] = {2.00*eV, 3.47*eV};
     const G4int nbins = sizeof(p_mirror)/sizeof(G4double);
     G4double refl_mirror[] = {fMirrorReflectivity,fMirrorReflectivity};
     assert(sizeof(refl_mirror) == sizeof(p_mirror));
     G4double effi_mirror[] = {0, 0};
     assert(sizeof(effi_mirror) == sizeof(effi_mirror));

     mirrorSurfaceProperty->
                       AddProperty("REFLECTIVITY",p_mirror,refl_mirror,nbins);
     mirrorSurfaceProperty->
                       AddProperty("EFFICIENCY",p_mirror,effi_mirror,nbins);

     mirrorSurface -> SetMaterialPropertiesTable(mirrorSurfaceProperty);

     new G4PVPlacement(0,
                       G4ThreeVector(0.0,0.0,fMirrorOrigin),
                       logicMirror,
                       "Mirror",
                       fLogicWorld,
                       false,
                       0);

     new G4LogicalSkinSurface("MirrorSurface",logicMirror,mirrorSurface);
  }

  //--------------------------------------------------
  // Coupling at the read-out end
  //--------------------------------------------------

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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::ConstructSDandField()
{
  if (!fmppcSD.Get()) {
     G4String mppcSDName = "SCIBAR/PhotonDet";
     SCIBARPhotonDetSD* mppcSD = new SCIBARPhotonDetSD(mppcSDName);
     fmppcSD.Put(mppcSD);
  }
  SetSensitiveDetector("PhotonDet_LV", fmppcSD.Get(), true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::UpdateGeometryParameters()
{
  fSCIBARfiberRX  = fXYRatio * fSCIBARfiberRY;

  fClad1RX = fSCIBARfiberRX + 0.03*fSCIBARfiberRX;
  fClad1RY = fSCIBARfiberRY + 0.03*fSCIBARfiberRY;
  fClad1Z  = fSCIBARfiberZ;

  fClad2RX = fClad1RX + 0.03*fSCIBARfiberRX;
  fClad2RY = fClad1RY + 0.03*fSCIBARfiberRY;
  fClad2Z  = fSCIBARfiberZ;

  fWorldSizeX = fClad2RX   + fMPPCDist + fMPPCHalfL + 5.*cm;//1.*cm;
  fWorldSizeY = fClad2RY   + fMPPCDist + fMPPCHalfL + 1.*cm;
  fWorldSizeZ = fSCIBARfiberZ + fMPPCDist + fMPPCHalfL + 5.*cm;//1.0*cm
 
  fCoupleRX = fWorldSizeX;
  fCoupleRY = fWorldSizeY;
  fCoupleZ  = (fWorldSizeZ - fSCIBARfiberZ) / 2;
 
  fClrfiberHalfL = fMPPCHalfL;
 
  fMirrorRmax = fClad2RY;
 
  fCoupleOrigin = fSCIBARfiberOrigin + fSCIBARfiberZ + fCoupleZ;
  fMirrorOrigin = fSCIBARfiberOrigin - fSCIBARfiberZ - fMirrorZ;
  fMPPCOriginX  = std::sin(fMPPCTheta) * (fMPPCDist + fClrfiberZ);
  fMPPCOriginZ  = -fCoupleZ + std::cos(fMPPCTheta) * (fMPPCDist + fClrfiberZ);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RotationMatrix
            SCIBARDetectorConstruction::StringToRotationMatrix(G4String rotation)
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

void SCIBARDetectorConstruction::SetPhotonDetGeometry (G4String shape)
// Set the Geometry of the PhotonDet detector
// Pre:  shape must be either "Circle" and "Square"
{
  if (shape == "Circle" || shape == "Square" ) fMPPCShape = shape;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetNumberOfCladding(G4int num)
// Set the number of claddings
// Pre: 0 <= num <= 2
{
  fNumOfCladLayers = num;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetSCIBARLength (G4double length)
// Set the TOTAL length of the SCIBAR fiber
{
  fSCIBARfiberZ = length;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetSCIBARRadius (G4double radius)
// Set the Y radius of SCIBAR fiber
{
  fSCIBARfiberRY = radius;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetClad1Radius (G4double radius)
// Set the Y radius of Cladding 1
{
  fClad1RY = radius;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetClad2Radius (G4double radius)
// Set the Y radius of Cladding 2
{
  fClad2RY = radius;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetPhotonDetHalfLength(G4double halfL)
// Set the half length of the PhotonDet detector
// The half length will be the radius if PhotonDet is circular
{
  fMPPCHalfL = halfL;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetGap (G4double gap)
// Set the distance between fiber end and PhotonDet
{ 
  fMPPCDist = gap;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetPhotonDetAlignment(G4double theta)
// Set the Aligment of PhotonDet with respect to the z axis
// If theta is 0 deg, then the detector is perfectly aligned
// PhotonDet will be deviated by theta from z axis
// facing towards the center of the fiber
{
  fMPPCTheta = theta;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetSurfaceRoughness(G4double roughness)
// Set the Surface Roughness between Cladding 1 and SCIBAR fiber
// Pre: 0 < roughness <= 1
{
  fSurfaceRoughness = roughness;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetMirrorPolish(G4double polish)
// Set the Polish of the mirror, polish of 1 is a perfect mirror surface
// Pre: 0 < polish <= 1
{
  fMirrorPolish = polish;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetMirrorReflectivity(G4double reflectivity)
// Set the Reflectivity of the mirror, reflectivity of 1 is a perfect mirror
// Pre: 0 < reflectivity <= 1
{
  fMirrorReflectivity = reflectivity;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetPhotonDetPolish(G4double polish)
// Set the Polish of the PhotonDet, polish of 1 is a perfect mirror surface
// Pre: 0 < polish <= 1
{
  fMPPCPolish = polish;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetPhotonDetReflectivity(G4double reflectivity)
// Set the Reflectivity of the PhotonDet, reflectivity of 1 is a perfect mirror
// Pre: 0 < reflectivity <= 1
{
  fMPPCReflectivity = reflectivity;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetMirror(G4bool flag)
// Toggle to place the mirror or not at one end (-z end) of the fiber
// True means place the mirror, false means otherwise
{
  fMirrorToggle = flag;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetXYRatio(G4double r)
// Set the ratio of the x and y radius of the ellipse (x/y)
// a ratio of 1 would produce a circle
{
  fXYRatio = r;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetBarLength (G4double length)
// Set the length of the scintillator bar
{
  fBarLength = length;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetBarBase (G4double side)
// Set the side of the scintillator bar
{
  fBarBase = side;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void SCIBARDetectorConstruction::SetBarWidth (G4double side)
// Set the side of the scintillator bar
{
  fBarWidth = side;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void SCIBARDetectorConstruction::SetBarThick (G4double side)
// Set the side of the scintillator bar
{
  fBarThick = side;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetHoleRadius (G4double radius)
// Set the radius of the fiber hole
{
  fHoleRadius = radius;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetCoatingThickness (G4double thick)
// Set thickness of the coating on the bars
{
  fCoatingThickness = thick;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARDetectorConstruction::SetCoatingRadius (G4double radius)
// Set inner radius of the corner bar coating
{
  fCoatingRadius = radius;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SCIBARDetectorConstruction::GetSCIBARFiberLength() { return fSCIBARfiberZ; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SCIBARDetectorConstruction::GetBarLength() { return fBarLength; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SCIBARDetectorConstruction::GetBarBase() { return fBarBase; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SCIBARDetectorConstruction::GetBarWidth() { return fBarWidth; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SCIBARDetectorConstruction::GetBarThick() { return fBarThick; }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SCIBARDetectorConstruction::GetHoleRadius() { return fHoleRadius; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SCIBARDetectorConstruction::GetHoleLength() { return fHoleLength; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SCIBARDetectorConstruction::GetFiberRadius() { return GetSCIBARFiberRMax(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SCIBARDetectorConstruction::GetCoatingThickness()
                                                   { return fCoatingThickness; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SCIBARDetectorConstruction::GetCoatingRadius() { return fCoatingRadius; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SCIBARDetectorConstruction::GetSCIBARFiberEnd()
{
  return fSCIBARfiberOrigin + fSCIBARfiberZ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SCIBARDetectorConstruction::GetSCIBARFiberRMax()
{
  if (fNumOfCladLayers == 2) return fClad2RY;
  if (fNumOfCladLayers == 1) return fClad1RY;
  return fSCIBARfiberRY;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SCIBARDetectorConstruction::GetSurfaceRoughness()
{
  return fSurfaceRoughness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Return True if the fiber construction is ideal
G4bool SCIBARDetectorConstruction::IsPerfectFiber()
{
  return     fSurfaceRoughness == 1. && fXYRatio == 1.
             && (!fMirrorToggle    ||
             (fMirrorPolish    == 1. && fMirrorReflectivity == 1.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* SCIBARDetectorConstruction::FindMaterial(G4String name) {
    G4Material* material = G4Material::GetMaterial(name,true);
    return material;
}
