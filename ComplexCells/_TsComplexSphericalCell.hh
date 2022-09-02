// Component for __TsComplexSphericalCell

#ifndef _TsComplexSphericalCell_hh
#define _TsComplexSphericalCell_hh

#include "TsVComponentWithChildren.hh"

class _TsComplexSphericalCell : public TsVComponentWithChildren
{
public:
	_TsComplexSphericalCell(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
							TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~_TsComplexSphericalCell() = default;

	G4VPhysicalVolume* Construct() override;

protected:
	virtual void ConstructNucleus();
	virtual void ConstructMitochondria();
	virtual void ConstructLysosomes();

protected:
	G4double fCellRadius, fNucleusRadius;
};

#endif
