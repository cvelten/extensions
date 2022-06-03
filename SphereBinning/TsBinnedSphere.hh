// Component for TsBinnedSphere

#ifndef TsBinnedSphere_hh
#define TsBinnedSphere_hh

#include "TsVGeometryComponent.hh"

class TsBinnedSphere : public TsVGeometryComponent
{
public:
	TsBinnedSphere(TsParameterManager *pM, TsExtensionManager *eM, TsMaterialManager *mM, TsGeometryManager *gM,
				   TsVGeometryComponent *parentComponent, G4VPhysicalVolume *parentVolume, G4String &name);
	~TsBinnedSphere();

	G4VPhysicalVolume *Construct();

	// Methods used by scorers
	G4int GetIndex(G4Step *aStep);
	G4int GetIndex(G4int iR, G4int iPhi, G4int iTheta);
	G4int GetBin(G4int index, G4int iBin);
	SurfaceType GetSurfaceID(G4String surfaceName);
	G4bool IsOnBoundary(G4ThreeVector localpos, G4VSolid *solid, SurfaceType surfaceID);
	G4double GetAreaOfSelectedSurface(G4VSolid *solid, SurfaceType surfaceID, G4int i, G4int j, G4int k);
	void PreCalculateThetaRatios();

	// Methods used by parameterization
	G4Material *ComputeMaterial(const G4int repNo, G4VPhysicalVolume *pvol, const G4VTouchable *parent);
	void ComputeTransformation(const G4int, G4VPhysicalVolume *pvol) const;
	void ComputeDimensions(G4Sphere &, const G4int, const G4VPhysicalVolume *) const;
	void ComputeDimensions(G4Tubs &, const G4int, const G4VPhysicalVolume *) const {};

	// Precalculated ratios to account for area theta divisions
	std::vector<G4double> fThetaAreaRatios;

	static void CreateDefaults(TsParameterManager *pM, G4String &childName, G4String &parentName);

private:
	G4double fTotalRMin;
	G4double fTotalRMax;
	G4double fTotalSTheta;

	std::vector<G4RotationMatrix *> fPhiDivisionRotations;

	RadialBinning fRadialBinning;
	std::vector<G4double> fRBinValues;
	G4double fDeltaR;

	G4bool fConstructParameterized;
	G4bool fUseOldParameterizationForSphere;

	const std::vector<G4double> &GetRBinValues();
};

#endif
