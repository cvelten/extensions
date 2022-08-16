// Component for TsComplexSphericalCell

#ifndef TsComplexSphericalCell_hh
#define TsComplexSphericalCell_hh

#include "TsVComponentWithChildren.hh"

class TsComplexSphericalCell : public TsVComponentWithChildren
{
public:
	TsComplexSphericalCell(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
						   TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsComplexSphericalCell() = default;

	G4VPhysicalVolume* Construct() override;

protected:
	virtual void ConstructNucleus();
	virtual void ConstructMitochondria();
	virtual void ConstructLysosomes();

protected:
	G4double fCellRadius, fNucleusRadius;
};

#endif
