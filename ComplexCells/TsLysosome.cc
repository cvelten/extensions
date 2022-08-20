// Component for TsLysosome

#include "TsLysosome.hh"

#include "TsParameterManager.hh"

#include "G4Ellipsoid.hh"
#include "G4LogicalVolume.hh"
#include "G4Orb.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
// #include "Randomize.hh"

#define um3 (um * um * um)
#define nm3 (nm * nm * nm)

TsLysosome::TsLysosome(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM, TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name)
	: TsVComponentWithChildren(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
	fIsDividable = false;
	fCanCalculateSurfaceArea = true;
}

G4VPhysicalVolume* TsLysosome::Construct()
{
	BeginConstruction();

	//
	// Lysosome

	fSemiAxisA = fPm->GetDoubleParameter(GetFullParmName("SemiAxisA"), "Length");
	fSemiAxisB = fPm->GetDoubleParameter(GetFullParmName("SemiAxisB"), "Length");
	fSemiAxisC = fPm->GetDoubleParameter(GetFullParmName("SemiAxisC"), "Length");

	G4Ellipsoid* solid = new G4Ellipsoid(fName, fSemiAxisA, fSemiAxisB, fSemiAxisC);
	fEnvelopeLog = CreateLogicalVolume(solid);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

	//
	// Enzyme Complexes (with nanomaterial)

	auto nChildren = fPm->GetIntegerParameter(GetFullParmName("Complexes", "N"));
	auto rChildren = fPm->ParameterExists(GetFullParmName("Complexes/Radius")) ? fPm->GetDoubleParameter(GetFullParmName("Complexes", "Radius"), "Length") : 0.1 * micrometer;

	// If no nanomaterial densities are given, they are all -0-
	auto nanomaterialDensities = std::vector<G4double>(nChildren, 0.0);
	if (fPm->ParameterExists(GetFullParmName("Complexes", "NanomaterialDensity")))
	{
		nanomaterialDensities = std::vector<G4double>(nChildren, fPm->GetUnitlessParameter(GetFullParmName("Complexes", "NanomaterialDensity")));
	}
	if (fPm->ParameterExists(GetFullParmName("Complexes", "NanomaterialDensities")))
	{
		if (fVerbosity > 0 && fPm->ParameterExists(GetFullParmName("Complexes", "NanomaterialDensity")))
			G4cerr << "Nanomaterial density was initially set through " << GetFullParmName("Complexes", "NanomaterialDensity")
				   << "but " << GetFullParmName("Complexes", "NanomaterialDensities") << " is also provided, so we'll use that instead!"
				   << G4endl;

		auto arr = fPm->GetUnitlessVector(GetFullParmName("Complexes", "NanomaterialDensities"));
		nanomaterialDensities = std::vector<G4double>(arr, arr + nChildren);

		if (!std::all_of(nanomaterialDensities.cbegin(), nanomaterialDensities.cend(), [](G4double n) { return n >= 0 && n <= 1; }))
			OutOfRange(GetFullParmName("Complexes", "NanomaterialDensities"), "must have values within [0, 1]");

		// Add color & drawing style if not specified
		if (!fPm->ParameterExists(GetFullParmName("Complexes/Nanomaterial", "Color")))
			fPm->AddParameter("s:" + GetFullParmName("Complexes/Nanomaterial", "Color"), "\"green\"");
		if (!fPm->ParameterExists(GetFullParmName("Complexes/Nanomaterial", "DrawingStyle")))
			fPm->AddParameter("s:" + GetFullParmName("Complexes/Nanomaterial", "DrawingStyle"), "\"solid\"");
	}

	std::vector<G4VPhysicalVolume*> lysosome_children = ConstructSphericalChildren("Complexes", nChildren, rChildren, 0, fEnvelopePhys, true);
	for (size_t i = 0; i < lysosome_children.size(); ++i)
	{
		if (std::ceil(nanomaterialDensities[i] * fNanomaterialMaxN) < 1) continue;

		auto nanomaterial_radius = fNanomaterialRadiusPerComplexRadius * rChildren;
		auto nanomaterial_volume = 4. / 3. * CLHEP::pi * std::pow(nanomaterial_radius, 3);
		auto nanomaterial_n = (G4int)std::ceil(nanomaterialDensities[i] * fNanomaterialMaxN);
		auto child_volume = 4. / 3. * CLHEP::pi * std::pow(rChildren, 3);

		G4cerr << "n = " << nanomaterial_n << " / V = " << nanomaterial_n * nanomaterial_volume / um3 << " um3 / VComplex = " << child_volume / um3 << " um3 = " << nanomaterial_n * nanomaterial_volume / child_volume * 1000 << " mM" << G4endl;

		ConstructSphericalChildren("Complexes/Nanomaterial", nanomaterial_n, nanomaterial_radius, 0, lysosome_children[i], true);
	}

	if (fVerbosity > 0)
	{
		G4cout << "Nanomaterial volumes in component '" << fName << "' (V = " << GetCubicVolume() / um3 << " um3):" << G4endl;
		G4double totalNanomaterialVolume = 0;
		for (size_t i = 0; i < lysosome_children.size(); ++i)
		{
			auto nanomaterial_radius = fNanomaterialRadiusPerComplexRadius * rChildren;
			auto nanomaterial_volume = 4. / 3. * CLHEP::pi * std::pow(nanomaterial_radius, 3);
			auto nanomaterial_n = (G4int)std::ceil(nanomaterialDensities[i] * fNanomaterialMaxN);

			G4cout << "> " << lysosome_children[i]->GetName() << " (V=" << lysosome_children[i]->GetLogicalVolume()->GetSolid()->GetCubicVolume() / um3 << " um3, n="
				   << nanomaterial_n << "): " << nanomaterial_volume * nanomaterial_n / um3 << " um3"
				   << G4endl;

			totalNanomaterialVolume += nanomaterial_volume * nanomaterial_n;
		}
		G4cout << "Total nanomaterial volume: " << totalNanomaterialVolume / um3 << " um3" << G4endl;
	}

	//
	// Additional children

	if (fPm->ParameterExists(GetFullParmName("Spheres", "N")))
		ConstructSphericalChildren("Spheres");

	if (fPm->ParameterExists(GetFullParmName("Ellipsoids", "N")))
		ConstructSphericalChildren("Ellipsoids");

	InstantiateChildren(fEnvelopePhys);

	return fEnvelopePhys;
}

TsVGeometryComponent::SurfaceType TsLysosome::GetSurfaceID(G4String surfaceName)
{
	SurfaceType surfaceID;
	G4String surfaceNameLower = surfaceName;
	surfaceNameLower.toLower();
	if (surfaceNameLower == "outercurvedsurface")
		surfaceID = OuterCurvedSurface;
	else if (surfaceNameLower == "anysurface")
		surfaceID = AnySurface;
	else {
		surfaceID = None;
		G4cerr << "Topas is exiting due to a serious error in scoring." << G4endl;
		G4cerr << "Scorer name: " << GetName() << " has unknown surface name: " << surfaceName << G4endl;
		fPm->AbortSession(1);
	}
	return surfaceID;
}

G4bool TsLysosome::IsOnBoundary(G4ThreeVector localpos, G4VSolid* solid, SurfaceType surfaceID)
{
	switch (surfaceID)
	{
	case AnySurface:
	case OuterCurvedSurface:
		return ((G4Ellipsoid*)solid)->Inside(localpos) == kSurface;
	default:
		G4cerr << "Topas is exiting due to a serious error." << G4endl;
		G4cerr << "TsLysosome::IsOnBoundary called for unknown surface of component: " << fName << G4endl;
		fPm->AbortSession(1);
		return false;
	}
}

G4double TsLysosome::GetAreaOfSelectedSurface(G4VSolid* solid, SurfaceType surfaceID, G4int, G4int, G4int copyNo)
{
	switch (surfaceID) {
	case AnySurface:
	case OuterCurvedSurface:
		return ((G4Ellipsoid*)solid)->GetSurfaceArea();
	default:
		G4cerr << "Topas is exiting due to a serious error." << G4endl;
		G4cerr << "TsLysosome::GetAreaOfSelectedSurface called for unknown surface of component: " << fName << G4endl;
		fPm->AbortSession(1);
		return 0.;
	}
}