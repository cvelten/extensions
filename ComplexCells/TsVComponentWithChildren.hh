// Extra Class for TsVComponentWithChildren

#ifndef TsSphereWithChildren_hh
#define TsSphereWithChildren_hh

#include "TsVGeometryComponent.hh"

class TsVComponentWithChildren : public TsVGeometryComponent
{
public:
	TsVComponentWithChildren(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
							 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	virtual ~TsVComponentWithChildren() = default;

	virtual G4VPhysicalVolume* Construct() = 0;

public:
	G4ThreeVector GetPointWithinVolume(G4VSolid* solid, G4double minDistanceFromSurface = 0) const;

	virtual void DuplicateSource();

	static G4String ConstructParameterName(const char* component, const char* parmName);
	static std::string StringReplace(std::string str, const std::string& from, const std::string& to);

protected:
	std::vector<G4VPhysicalVolume*> RandomlyPlaceSolid(G4VSolid* solid, G4int n = 1, G4double radialOffset = 0, G4VPhysicalVolume* parent = nullptr, G4bool independentLogicals = false, G4String name = "");

	std::vector<G4VPhysicalVolume*> ConstructSphericalChildren(G4String name, G4double radialOffset = 0);
	std::vector<G4VPhysicalVolume*> ConstructSphericalChildren(G4String name, G4int n, G4double radius, G4double radialOffset, G4VPhysicalVolume* parent = nullptr, G4bool independentLogicals = false);
	std::vector<G4VPhysicalVolume*> ConstructEllipsoidalChildren(G4String name, G4double radialOffset = 0);
	std::vector<G4VPhysicalVolume*> ConstructEllipsoidalChildren(G4String name, G4int n, G4double semiAxisA, G4double semiAxisB, G4double semiAxisC, G4double radialOffset, G4VPhysicalVolume* parent = nullptr, G4bool independentLogicals = false);

protected:
	const G4int fMaxNumAttempts = 10000;
};

#endif
