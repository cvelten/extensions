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

TsLysosome::TsLysosome(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM, TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name)
	: TsSphereWithChildren(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
	fIsDividable = false;
	fCanCalculateSurfaceArea = false;
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
	if (fPm->ParameterExists(GetFullParmName("Complexes", "NanomaterialDensities")))
	{
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
		if (nanomaterialDensities[i] < 1e-6) continue;

		const auto nanomaterial_max = 500;

		// auto child_volume = (*it)->GetLogicalVolume()->GetSolid()->GetCubicVolume();
		// auto nanomaterial_radius = std::cbrt(nanomaterial_volume / (4. / 3. * CLHEP::pi));
		auto nanomaterial_radius = 0.05 * rChildren;
		auto nanomaterial_volume = 4. / 3. * CLHEP::pi * std::pow(nanomaterial_radius, 3);
		auto nanomaterial_n = (G4int)std::ceil(nanomaterialDensities[i] * nanomaterial_max);
		auto child_volume = 4. / 3. * CLHEP::pi * std::pow(rChildren, 3);

		const G4double um3 = (micrometer * micrometer * micrometer);
		G4cerr << "n = " << nanomaterial_n << " / V = " << nanomaterial_n * nanomaterial_volume / um3 << " um3 / VComplex = " << child_volume / um3 << " um3 = " << nanomaterial_n * nanomaterial_volume / child_volume * 1000 << " mM" << G4endl;

		ConstructSphericalChildren("Complexes/Nanomaterial", nanomaterial_n, nanomaterial_radius, 0, lysosome_children[i]);
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

void TsLysosome::ConstructLysosome()
{
}
