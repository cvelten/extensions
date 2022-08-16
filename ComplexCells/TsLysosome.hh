// Component for TsLysosome

#ifndef TsLysosome_hh
#define TsLysosome_hh

#include "TsSphereWithChildren.hh"

class TsLysosome : public TsSphereWithChildren
{
public:
	TsLysosome(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			   TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsLysosome() = default;

	G4VPhysicalVolume* Construct() override;

protected:
	virtual void ConstructLysosome();

protected:
	G4double fSemiAxisA, fSemiAxisB, fSemiAxisC;
};

#endif
