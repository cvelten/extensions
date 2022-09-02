// Component for TsComplexCell

#include "TsComplexCell.hh"

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

TsComplexCell::TsComplexCell(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM, TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name)
	: TsVComponentWithChildren(pM, eM, mM, gM, parentComponent, parentVolume, name),
	  fLysosomesN(0), fMitochondriaN(0)
{
	fIsDividable = false;
	// fCanCalculateSurfaceArea = false;
}

G4VPhysicalVolume* TsComplexCell::Construct()
{
	BeginConstruction();

	// Cell
	fRadius = fPm->GetDoubleParameter(GetFullParmName("Radius"), "Length");

	G4Orb* cellSolid = new G4Orb(fName, fRadius);
	fEnvelopeLog = CreateLogicalVolume(cellSolid);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

	// Optional
	ConstructNucleus();

	// Lysosomes
	ConstructLysosomes();

	// Mitochondria
	ConstructMitochondria();

	//
	// Additional children

	if (fPm->ParameterExists(GetFullParmName("Spheres", "N")))
		ConstructSphericalChildren("Spheres");

	if (fPm->ParameterExists(GetFullParmName("Ellipsoids", "N")))
		ConstructSphericalChildren("Ellipsoids");

	InstantiateChildren(fEnvelopePhys);

	return fEnvelopePhys;
}

void TsComplexCell::ConstructNucleus()
{
	if (!fPm->ParameterExists(GetFullParmName("Nucleus/Radius")))
		return;

	auto fNucleusRadius = fPm->GetDoubleParameter(GetFullParmName("Nucleus/Radius"), "Length");

	if (!fPm->ParameterExists(GetFullParmName("Nucleus", "Material")))
		fPm->CloneParameter(GetFullParmName("Material"), GetFullParmName("Nucleus", "Material"));
	if (!fPm->ParameterExists(GetFullParmName("Nucleus", "DrawingStyle")))
		fPm->AddParameter("s:" + GetFullParmName("Nucleus", "DrawingStyle"), "\"wireframe\"");
	if (!fPm->ParameterExists(GetFullParmName("Nucleus", "Color")))
		fPm->AddParameter("s:" + GetFullParmName("Nucleus", "Color"), "\"red\"");

	G4RotationMatrix* nucRot = new G4RotationMatrix();
	nucRot->rotateX(0);
	nucRot->rotateY(0);

	auto transNucX = 0 * um;
	auto transNucY = 0 * um;
	auto transNucZ = 0 * um;
	if (fPm->ParameterExists(GetFullParmName("Nucleus/TransX"))) {
		transNucX = fPm->GetDoubleParameter(GetFullParmName("Nucleus/TransX"), "Length");
		if (transNucX > fRadius - fNucleusRadius)
			OutOfRange(GetFullParmName("Nucleus/TransX"), "sets the nucleus outside of the cell");
	}
	if (fPm->ParameterExists(GetFullParmName("Nucleus/TransY"))) {
		transNucY = fPm->GetDoubleParameter(GetFullParmName("Nucleus/TransY"), "Length");
		if (transNucY > fRadius - fNucleusRadius)
			OutOfRange(GetFullParmName("Nucleus/TransY"), "sets the nucleus outside of the cell");
	}
	if (fPm->ParameterExists(GetFullParmName("Nucleus/TransZ"))) {
		transNucZ = fPm->GetDoubleParameter(GetFullParmName("Nucleus/TransZ"), "Length");
		if (transNucZ > fRadius - fNucleusRadius)
			OutOfRange(GetFullParmName("Nucleus/TransZ"), "sets the nucleus outside of the cell");
	}

	G4ThreeVector* nucPos = new G4ThreeVector(transNucX, transNucY, transNucZ);

	auto name = "Nucleus";
	G4Orb* nucSolid = new G4Orb(name, fNucleusRadius);
	G4LogicalVolume* nucLogical = CreateLogicalVolume(name, nucSolid);
	G4VPhysicalVolume* nucleusPhysical = CreatePhysicalVolume(name, nucLogical, nucRot, nucPos, fEnvelopePhys);

	if (nucleusPhysical->CheckOverlaps())
	{
		G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
		G4cerr << "Nucleus overlaps with the cell." << G4endl;
		fPm->AbortSession(1);
	}
}

void TsComplexCell::ConstructLysosomes()
{
	if (fPm->ParameterExists(GetFullParmName("Lysosome", "N")))
		fLysosomesN = fPm->GetIntegerParameter(GetFullParmName("Lysosome", "N"));
	if (fLysosomesN == 0)
		return;

	auto semiAxisA = fPm->GetDoubleParameter(GetFullParmName("Lysosome", "SemiAxisA"), "Length");
	auto semiAxisB = fPm->GetDoubleParameter(GetFullParmName("Lysosome", "SemiAxisB"), "Length");
	auto semiAxisC = fPm->GetDoubleParameter(GetFullParmName("Lysosome", "SemiAxisC"), "Length");

	std::vector<G4String> lysosomeParams;
	fPm->GetParameterNamesStartingWith(GetFullParmName("Lysosome"), &lysosomeParams);

	//
	// Enzyme Complexes (with nanomaterial)

	auto nComplexes = fPm->GetIntegerParameter(GetFullParmName("Lysosome/Complexes", "N"));
	auto rComplexes = fPm->ParameterExists(GetFullParmName("Lysosome/Complexes/Radius")) ? fPm->GetDoubleParameter(GetFullParmName("Lysosome/Complexes", "Radius"), "Length") : 0.1 * micrometer;

	if (!fPm->ParameterExists(GetFullParmName("Lysosome/Complexes", "Color")))
		fPm->AddParameter("s:" + GetFullParmName("Lysosome/Complexes", "Color"), "\"gray\"");

	// Nanomaterial

	G4double nanomaterialRadius;
	std::vector<G4int> nanomaterialNumbers;

	// If no nanomaterial densities are given, they are all -0-
	nanomaterialNumbers = std::vector<G4int>(nComplexes, 0.0);
	if (fPm->ParameterExists(GetFullParmName("Lysosome/Complexes/Nanomaterial", "Number")) || fPm->ParameterExists(GetFullParmName("Lysosome/Complexes/Nanomaterial", "Numbers")))
	{
		// Radius
		if (fPm->ParameterExists(GetFullParmName("Lysosome/Complexes/Nanomaterial", "Radius")))
			nanomaterialRadius = fPm->GetDoubleParameter(GetFullParmName("Lysosome/Complexes/Nanomaterial", "Radius"), "Length");
		else
		{
			G4cerr << GetFullParmName("Lysosome/Complexes/Nanomaterial", "Radius")
				   << " was not given, so we'll assume "
				   << fNanomaterialRadiusDefault / nm << " nm"
				   << G4endl;
			nanomaterialRadius = fNanomaterialRadiusDefault;
		}

		if (fPm->ParameterExists(GetFullParmName("Lysosome/Complexes/Nanomaterial", "Number")))
		{
			nanomaterialNumbers = std::vector<G4int>(nComplexes, fPm->GetIntegerParameter(GetFullParmName("Lysosome/Complexes/Nanomaterial", "Number")));
		}
		if (fPm->ParameterExists(GetFullParmName("Lysosome/Complexes/Nanomaterial", "Numbers")))
		{
			if (fVerbosity > 0 && fPm->ParameterExists(GetFullParmName("Lysosome/Complexes/Nanomaterial", "Number")))
				G4cerr << "Nanomaterial density was initially set through " << GetFullParmName("Lysosome/Complexes/Nanomaterial", "Number")
					   << "but " << GetFullParmName("Lysosome/Complexes/Nanomaterial", "Numbers") << " is also provided, so we'll use that instead!"
					   << G4endl;

			// Numbers
			auto arr = fPm->GetIntegerVector(GetFullParmName("Lysosome/Complexes/Nanomaterial", "Numbers"));
			nanomaterialNumbers = std::vector<G4int>(arr, arr + nComplexes);

			auto max_num = std::floor(std::pow(rComplexes / nanomaterialRadius, 3));
			if (!std::all_of(nanomaterialNumbers.cbegin(), nanomaterialNumbers.cend(), [&max_num](G4double n) { return n >= 0 && n <= max_num; }))
				OutOfRange(GetFullParmName("Lysosome/Complexes", "NanomaterialNumbers"), G4String("must have values within [0, ") + std::to_string(max_num));
		}

		// Add color & drawing style if not specified
		if (!fPm->ParameterExists(GetFullParmName("Lysosome/Complexes/Nanomaterial", "Color")))
			fPm->AddParameter("s:" + GetFullParmName("Lysosome/Complexes/Nanomaterial", "Color"), "\"green\"");
		if (!fPm->ParameterExists(GetFullParmName("Lysosome/Complexes/Nanomaterial", "DrawingStyle")))
			fPm->AddParameter("s:" + GetFullParmName("Lysosome/Complexes/Nanomaterial", "DrawingStyle"), "\"solid\"");
	}

	G4Ellipsoid* solid = new G4Ellipsoid("Lysosome", semiAxisA, semiAxisB, semiAxisC);
	std::vector<G4VPhysicalVolume*> lysosomes = RandomlyPlaceSolid(solid, fLysosomesN, 0, fEnvelopePhys, true);
	for (auto it_lyso = lysosomes.begin(); it_lyso != lysosomes.end(); ++it_lyso)
	{
		std::vector<G4VPhysicalVolume*> lysosome_children = ConstructSphericalChildren("Lysosome/Complexes", nComplexes, rComplexes, 0, *it_lyso, true);
		for (size_t j = 0; j < lysosome_children.size(); ++j)
		{
			if (nanomaterialNumbers[j] < 1) continue;
			ConstructSphericalChildren("Lysosome/Complexes/Nanomaterial", nanomaterialNumbers[j], nanomaterialRadius, 0, lysosome_children[j], true);
		}
	}
}

void TsComplexCell::ConstructMitochondria()
{
	if (!fPm->ParameterExists(GetFullParmName("Mitochondria", "N")))
		return;

	G4int fMitochondriaN = fPm->GetIntegerParameter(GetFullParmName("Mitochondria", "N"));

	if (!fPm->ParameterExists(GetFullParmName("Mitochondria", "Material")))
		fPm->CloneParameter(GetFullParmName("Material"), GetFullParmName("Mitochondria", "Material"));
	if (!fPm->ParameterExists(GetFullParmName("Mitochondria", "DrawingStyle")))
		fPm->AddParameter("s:" + GetFullParmName("Mitochondria", "DrawingStyle"), "\"wireframe\"");
	if (!fPm->ParameterExists(GetFullParmName("Mitochondria", "Color")))
		fPm->AddParameter("s:" + GetFullParmName("Mitochondria", "Color"), "\"green\"");

	// defaults, may be overriden by fPm [TBD]
	G4double semiAxisA = 0.4 * micrometer;
	G4double semiAxisB = 0.15 * micrometer;
	G4double semiAxisC = 0.15 * micrometer;

	std::vector<G4VPhysicalVolume*> mitochondria = ConstructEllipsoidalChildren("Mitochondria", fMitochondriaN, semiAxisA, semiAxisB, semiAxisC, fNucleusRadius);
}

TsVGeometryComponent::SurfaceType TsComplexCell::GetSurfaceID(G4String surfaceName)
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

G4bool TsComplexCell::IsOnBoundary(G4ThreeVector localpos, G4VSolid* solid, SurfaceType surfaceID)
{
	switch (surfaceID)
	{
	case AnySurface:
	case OuterCurvedSurface:
		return ((G4Orb*)solid)->Inside(localpos) == kSurface;
	default:
		G4cerr << "Topas is exiting due to a serious error." << G4endl;
		G4cerr << "TsComplexCell::IsOnBoundary called for unknown surface of component: " << fName << G4endl;
		fPm->AbortSession(1);
		return false;
	}
}
