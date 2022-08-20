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
#define micromolar (1E-3 * mole / m3)

TsLysosome::TsLysosome(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM, TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name)
	: TsVComponentWithChildren(pM, eM, mM, gM, parentComponent, parentVolume, name),
	  fNanomaterialRadius(fNanomaterialRadiusDefault)
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
	fNanomaterialNumbers = std::vector<G4int>(nChildren, 0.0);
	if (fPm->ParameterExists(GetFullParmName("Complexes/Nanomaterial", "Number")) || fPm->ParameterExists(GetFullParmName("Complexes/Nanomaterial", "Numbers")))
	{
		if (fPm->ParameterExists(GetFullParmName("Complexes/Nanomaterial", "Number")))
		{
			fNanomaterialNumbers = std::vector<G4int>(nChildren, fPm->GetIntegerParameter(GetFullParmName("Complexes/Nanomaterial", "Number")));
		}
		if (fPm->ParameterExists(GetFullParmName("Complexes/Nanomaterial", "Numbers")))
		{
			if (fVerbosity > 0 && fPm->ParameterExists(GetFullParmName("Complexes/Nanomaterial", "Number")))
				G4cerr << "Nanomaterial density was initially set through " << GetFullParmName("Complexes/Nanomaterial", "Number")
					   << "but " << GetFullParmName("Complexes/Nanomaterial", "Numbers") << " is also provided, so we'll use that instead!"
					   << G4endl;

			// Numbers
			auto arr = fPm->GetIntegerVector(GetFullParmName("Complexes/Nanomaterial", "Numbers"));
			fNanomaterialNumbers = std::vector<G4int>(arr, arr + nChildren);

			auto max_num = std::floor(std::pow(rChildren / fNanomaterialRadius, 3));
			if (!std::all_of(fNanomaterialNumbers.cbegin(), fNanomaterialNumbers.cend(), [&max_num](G4double n) { return n >= 0 && n <= max_num; }))
				OutOfRange(GetFullParmName("Complexes", "NanomaterialNumbers"), G4String("must have values within [0, ") + std::to_string(max_num));
		}
		// Radius
		if (fPm->ParameterExists(GetFullParmName("Complexes/Nanomaterial", "Radius")))
			fNanomaterialRadius = fPm->GetDoubleParameter(GetFullParmName("Complexes/Nanomaterial", "Radius"), "Length");
		else
		{
			G4cerr << GetFullParmName("Complexes/Nanomaterial", "Radius") << " was not given, so we'll assume " << fNanomaterialRadiusDefault / nm << G4endl;
			fNanomaterialRadius = fNanomaterialRadiusDefault;
		}
		// Add color & drawing style if not specified
		if (!fPm->ParameterExists(GetFullParmName("Complexes/Nanomaterial", "Color")))
			fPm->AddParameter("s:" + GetFullParmName("Complexes/Nanomaterial", "Color"), "\"green\"");
		if (!fPm->ParameterExists(GetFullParmName("Complexes/Nanomaterial", "DrawingStyle")))
			fPm->AddParameter("s:" + GetFullParmName("Complexes/Nanomaterial", "DrawingStyle"), "\"solid\"");
	}

	std::vector<G4VPhysicalVolume*> lysosome_children = ConstructSphericalChildren("Complexes", nChildren, rChildren, 0, fEnvelopePhys, true);
	for (size_t i = 0; i < lysosome_children.size(); ++i)
	{
		if (fNanomaterialNumbers[i] < 1) continue;
		ConstructSphericalChildren("Complexes/Nanomaterial", fNanomaterialNumbers[i], fNanomaterialRadius, 0, lysosome_children[i], true);
	}

	if (fVerbosity > 0)
	{
		G4cout << "Nanomaterial volumes in component '" << fName << "' (V = " << GetCubicVolume() / um3 << " um3):" << G4endl;
		G4int totalNanomaterialNumber = 0;
		G4double totalNanomaterialVolume = 0;
		for (size_t i = 0; i < lysosome_children.size(); ++i)
		{
			auto nanomaterial_volume = 4. / 3. * CLHEP::pi * std::pow(fNanomaterialRadius, 3);
			auto child_volume = lysosome_children[i]->GetLogicalVolume()->GetSolid()->GetCubicVolume();
			auto child_molarity = (fNanomaterialNumbers[i] / Avogadro) / child_volume;

			totalNanomaterialNumber += fNanomaterialNumbers[i];
			totalNanomaterialVolume += nanomaterial_volume * fNanomaterialNumbers[i];

			G4cout << "> " << lysosome_children[i]->GetName() << " (V=" << child_volume / um3 << " um3, n="
				   << fNanomaterialNumbers[i] << ", c=" << child_molarity / micromolar << " uM): "
				   << nanomaterial_volume * fNanomaterialNumbers[i] / um3 << " um3"
				   << G4endl;
		}
		G4cout << "Total nanomaterial volume:   " << totalNanomaterialVolume / um3 << " um3" << G4endl;
		G4cout << "Total nanomaterial molarity: " << (totalNanomaterialNumber / Avogadro) / GetCubicVolume() / micromolar << " uM"
			   << G4endl;
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