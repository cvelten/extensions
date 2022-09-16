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
	  fUseParameterSystem(false), fLysosomesN(0), fMitochondriaN(0)
{
	fIsDividable = false;
	// fCanCalculateSurfaceArea = false;
}

G4VPhysicalVolume* TsComplexCell::Construct()
{
	BeginConstruction();

	if (fPm->ParameterExists(GetFullParmName("UseParameterSystem")))
		fUseParameterSystem = fPm->GetBooleanParameter(GetFullParmName("UseParameterSystem"));

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

	// Remove if using parameter system
	if (fUseParameterSystem)
	{
		for (auto it = fLysosomePhysicals.begin(); it != fLysosomePhysicals.end(); ++it)
			fEnvelopePhys->GetLogicalVolume()->RemoveDaughter(*it);
		fLysosomePhysicals.clear();
		for (auto it = fMitochondriaPhysicals.begin(); it != fMitochondriaPhysicals.end(); ++it)
			fEnvelopePhys->GetLogicalVolume()->RemoveDaughter(*it);
		fMitochondriaPhysicals.clear();
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
	fLysosomePhysicals = RandomlyPlaceSolid(solid, fLysosomesN, 0, fEnvelopePhys, true);

	if (fUseParameterSystem)
	{
		for (auto it = fLysosomePhysicals.begin(); it != fLysosomePhysicals.end(); ++it)
		{
			fPm->AddParameter("s:" + ConstructParameterName((*it)->GetName(), "Type"), "\"TsLysosome\"");
			fPm->AddParameter("s:" + ConstructParameterName((*it)->GetName(), "Parent"), "\"" + GetName() + "\"");
			fPm->CloneParameter(GetFullParmName("Lysosome/Material"), ConstructParameterName((*it)->GetName(), "Material"));
			fPm->CloneParameter(GetFullParmName("Lysosome/SemiAxisA"), ConstructParameterName((*it)->GetName(), "SemiAxisA"));
			fPm->CloneParameter(GetFullParmName("Lysosome/SemiAxisB"), ConstructParameterName((*it)->GetName(), "SemiAxisB"));
			fPm->CloneParameter(GetFullParmName("Lysosome/SemiAxisC"), ConstructParameterName((*it)->GetName(), "SemiAxisC"));
			fPm->CloneParameter(GetFullParmName("Lysosome/Color"), ConstructParameterName((*it)->GetName(), "Color"));
			fPm->CloneParameter(GetFullParmName("Lysosome/Complexes/N"), ConstructParameterName((*it)->GetName(), "Complexes/N"));
			fPm->CloneParameter(GetFullParmName("Lysosome/Complexes/Radius"), ConstructParameterName((*it)->GetName(), "Complexes/Radius"));
			fPm->CloneParameter(GetFullParmName("Lysosome/Complexes/Material"), ConstructParameterName((*it)->GetName(), "Complexes/Material"));
			fPm->CloneParameter(GetFullParmName("Lysosome/Complexes/Nanomaterial/Number"), ConstructParameterName((*it)->GetName(), "Complexes/Nanomaterial/Number"));
			fPm->CloneParameter(GetFullParmName("Lysosome/Complexes/Nanomaterial/Radius"), ConstructParameterName((*it)->GetName(), "Complexes/Nanomaterial/Radius"));
			fPm->CloneParameter(GetFullParmName("Lysosome/Complexes/Nanomaterial/Material"), ConstructParameterName((*it)->GetName(), "Complexes/Nanomaterial/Material"));
			fPm->CloneParameter(GetFullParmName("Lysosome/Complexes/Nanomaterial/Invisible"), ConstructParameterName((*it)->GetName(), "Complexes/Nanomaterial/Invisible"));

			auto trans = (*it)->GetTranslation();
			auto rot = (*it)->GetRotation();
			auto psi = std::acos(rot->yy());
			auto phi = std::acos(rot->xx());
			fPm->AddParameter("d:" + ConstructParameterName((*it)->GetName(), "TransX"), std::to_string(trans.x() / um) + " um");
			fPm->AddParameter("d:" + ConstructParameterName((*it)->GetName(), "TransY"), std::to_string(trans.y() / um) + " um");
			fPm->AddParameter("d:" + ConstructParameterName((*it)->GetName(), "TransZ"), std::to_string(trans.z() / um) + " um");
			fPm->AddParameter("d:" + ConstructParameterName((*it)->GetName(), "RotX"), std::to_string(psi / radian) + " radian");
			fPm->AddParameter("d:" + ConstructParameterName((*it)->GetName(), "RotY"), std::to_string(phi / radian) + " radian");
		}

		if (fPm->ParameterExists(GetFullParmName("Lysosome/DuplicateScorerForVolumes")))
		{
			auto scorerName = fPm->GetStringParameter(GetFullParmName("Lysosome/DuplicateScorerForVolumes"));
			std::vector<G4String> parameterNames;
			fPm->GetParameterNamesStartingWith("Sc/" + scorerName + "/", &parameterNames);
			for (auto it = fLysosomePhysicals.begin(); it != fLysosomePhysicals.end(); ++it)
			{
				for (std::string parName : parameterNames)
				{
					G4String nameWithoutSlash = StringReplace((*it)->GetName(), "/", "_");
					G4String newParName = StringReplace(parName, scorerName, scorerName + "_" + nameWithoutSlash);

					if (newParName.find("/OutputFile") != std::string::npos)
					{
						G4String outputFile = fPm->GetStringParameter(parName) + "_" + nameWithoutSlash;
						fPm->AddParameter("s:Sc/" + scorerName + "_" + nameWithoutSlash + "/OutputFile", '"' + outputFile + '"');
						continue;
					}

					fPm->CloneParameter(parName, newParName);
				}
			}
		}
	}
	else
	{
		for (auto it = fLysosomePhysicals.begin(); it != fLysosomePhysicals.end(); ++it)
		{
			G4cerr << (*it)->GetName() << ": " << (*it)->GetTranslation() << G4endl;

			std::vector<G4VPhysicalVolume*> lysosome_children = ConstructSphericalChildren("Lysosome/Complexes", nComplexes, rComplexes, 0, *it, true);
			for (size_t j = 0; j < lysosome_children.size(); ++j)
			{
				if (nanomaterialNumbers[j] < 1) continue;
				ConstructSphericalChildren("Lysosome/Complexes/Nanomaterial", nanomaterialNumbers[j], nanomaterialRadius, 0, lysosome_children[j], true);
			}
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
	{
		if (fPm->ParameterExists(GetFullParmName("Mitochondria/SemiAxisA")))
			semiAxisA = fPm->GetDoubleParameter(GetFullParmName("Mitochondria/SemiAxisA"), "Length");
		else
			fPm->AddParameter("d:" + GetFullParmName("Mitochondria/SemiAxisA"), std::to_string(semiAxisA / um) + " um");
		if (fPm->ParameterExists(GetFullParmName("Mitochondria/SemiAxisB")))
			semiAxisB = fPm->GetDoubleParameter(GetFullParmName("Mitochondria/SemiAxisB"), "Length");
		else
			fPm->AddParameter("d:" + GetFullParmName("Mitochondria/SemiAxisB"), std::to_string(semiAxisB / um) + " um");
		if (fPm->ParameterExists(GetFullParmName("Mitochondria/SemiAxisC")))
			semiAxisC = fPm->GetDoubleParameter(GetFullParmName("Mitochondria/SemiAxisC"), "Length");
		else
			fPm->AddParameter("d:" + GetFullParmName("Mitochondria/SemiAxisC"), std::to_string(semiAxisC / um) + " um");
	}

	fMitochondriaPhysicals = ConstructEllipsoidalChildren("Mitochondria", fMitochondriaN, semiAxisA, semiAxisB, semiAxisC, fNucleusRadius);

	if (fUseParameterSystem)
	{
		for (auto it = fMitochondriaPhysicals.begin(); it != fMitochondriaPhysicals.end(); ++it)
		{
			fPm->AddParameter("s:" + ConstructParameterName((*it)->GetName(), "Type"), "\"G4Ellipsoid\"");
			fPm->AddParameter("s:" + ConstructParameterName((*it)->GetName(), "Parent"), "\"" + GetName() + "\"");
			fPm->CloneParameter(GetFullParmName("Mitochondria/Material"), ConstructParameterName((*it)->GetName(), "Material"));
			fPm->CloneParameter(GetFullParmName("Mitochondria/SemiAxisA"), ConstructParameterName((*it)->GetName(), "HLX"));
			fPm->CloneParameter(GetFullParmName("Mitochondria/SemiAxisB"), ConstructParameterName((*it)->GetName(), "HLY"));
			fPm->CloneParameter(GetFullParmName("Mitochondria/SemiAxisC"), ConstructParameterName((*it)->GetName(), "HLZ"));
			fPm->CloneParameter(GetFullParmName("Mitochondria/Color"), ConstructParameterName((*it)->GetName(), "Color"));

			auto trans = (*it)->GetTranslation();
			auto rot = (*it)->GetRotation();
			auto psi = std::acos(rot->yy());
			auto phi = std::acos(rot->xx());
			fPm->AddParameter("d:" + ConstructParameterName((*it)->GetName(), "TransX"), std::to_string(trans.x() / um) + " um");
			fPm->AddParameter("d:" + ConstructParameterName((*it)->GetName(), "TransY"), std::to_string(trans.y() / um) + " um");
			fPm->AddParameter("d:" + ConstructParameterName((*it)->GetName(), "TransZ"), std::to_string(trans.z() / um) + " um");
			fPm->AddParameter("d:" + ConstructParameterName((*it)->GetName(), "RotX"), std::to_string(psi / radian) + " radian");
			fPm->AddParameter("d:" + ConstructParameterName((*it)->GetName(), "RotY"), std::to_string(phi / radian) + " radian");
		}
	}
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
