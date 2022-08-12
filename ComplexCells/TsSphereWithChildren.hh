// Component for TsSphereWithChildren

#ifndef TsSphereWithChildren_hh
#define TsSphereWithChildren_hh

#include "TsVGeometryComponent.hh"

class TsSphereWithChildren : public TsVGeometryComponent
{
public:
	TsSphereWithChildren(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
						 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	virtual ~TsSphereWithChildren() = default;

	virtual G4VPhysicalVolume* Construct();

protected:
	std::vector<G4VPhysicalVolume*> RandomlyPlaceSolid(G4VSolid* solid, G4int n, G4double radialOffset, G4VPhysicalVolume* parent = nullptr, G4bool independentLogicals = false, G4String name = "");

	std::vector<G4VPhysicalVolume*> ConstructSphericalChildren(G4String name, G4double radialOffset = 0);
	std::vector<G4VPhysicalVolume*> ConstructSphericalChildren(G4String name, G4int n, G4double radius, G4double radialOffset, G4VPhysicalVolume* parent = nullptr, G4bool independentLogicals = false);
	std::vector<G4VPhysicalVolume*> ConstructEllipsoidalChildren(G4String name, G4double radialOffset = 0);
	std::vector<G4VPhysicalVolume*> ConstructEllipsoidalChildren(G4String name, G4int n, G4double semiAxisA, G4double semiAxisB, G4double semiAxisC, G4double radialOffset, G4VPhysicalVolume* parent = nullptr, G4bool independentLogicals = false);

protected:
	G4double fRadius;
};

#endif
