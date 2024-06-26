// Component for TsComplexCell

#ifndef TsComplexCell_hh
#define TsComplexCell_hh

#include "VComponentWithChildren.hh"

class TsComplexCell : public VComponentWithChildren
{
public:
	TsComplexCell(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsComplexCell() = default;

	G4VPhysicalVolume* Construct() override;
	void ConstructNucleus();
	void ConstructLysosomes();
	void ConstructMitochondria();

	SurfaceType GetSurfaceID(G4String surfaceName) override;
	G4bool IsOnBoundary(G4ThreeVector localpos, G4VSolid* solid, SurfaceType surfaceID) override;
	// G4double GetAreaOfSelectedSurface(G4VSolid* solid, SurfaceType surfaceID, G4int i, G4int j, G4int k) override;

public:
private:
	G4bool fUseParameterSystem;
	G4double fRadius, fNucleusRadius;
	std::vector<G4VPhysicalVolume*> fLysosomePhysicals, fMitochondriaPhysicals;

	G4int fLysosomesN, fMitochondriaN;

	const G4double fNanomaterialRadiusDefault = 25 * CLHEP::nm;
};

#endif
