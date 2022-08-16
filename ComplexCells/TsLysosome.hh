// Component for TsLysosome

#ifndef TsLysosome_hh
#define TsLysosome_hh

#include "TsVComponentWithChildren.hh"

class TsLysosome : public TsVComponentWithChildren
{
public:
	TsLysosome(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			   TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsLysosome() = default;

	G4VPhysicalVolume* Construct() override;

	SurfaceType GetSurfaceID(G4String surfaceName) override;
	G4bool IsOnBoundary(G4ThreeVector localpos, G4VSolid* solid, SurfaceType surfaceID) override;
	G4double GetAreaOfSelectedSurface(G4VSolid* solid, SurfaceType surfaceID, G4int i, G4int j, G4int k) override;

private:
	G4double fSemiAxisA, fSemiAxisB, fSemiAxisC;
};

#endif
