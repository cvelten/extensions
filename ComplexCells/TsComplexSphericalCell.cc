// Component for TsComplexSphericalCell

#include "TsComplexSphericalCell.hh"

#include "TsParameterManager.hh"

#include "G4Ellipsoid.hh"
#include "G4LogicalVolume.hh"
#include "G4Orb.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
// #include "Randomize.hh"

TsComplexSphericalCell::TsComplexSphericalCell(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM, TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name)
	: TsVComponentWithChildren(pM, eM, mM, gM, parentComponent, parentVolume, name),
	  fCellRadius(), fNucleusRadius()
{}

G4VPhysicalVolume* TsComplexSphericalCell::Construct()
{
	BeginConstruction();

	//***********************************************************************
	//              Envelope Geometry : spherical cell
	//***********************************************************************

	fCellRadius = fPm->GetDoubleParameter(GetFullParmName("Radius"), "Length");

	G4Orb* cellSolid = new G4Orb(fName, fCellRadius);
	fEnvelopeLog = CreateLogicalVolume(cellSolid);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

	// Optional
	ConstructNucleus();

	// Lysosomes
	ConstructLysosomes();

	// Mitochondria
	ConstructMitochondria();

	// Additional spherical children
	if (fPm->ParameterExists(GetFullParmName("Spheres", "N")))
		ConstructSphericalChildren("Spheres");

	if (fPm->ParameterExists(GetFullParmName("Ellipsoids", "N")))
		ConstructSphericalChildren("Ellipsoids");

	InstantiateChildren(fEnvelopePhys);

	return fEnvelopePhys;
}

void TsComplexSphericalCell::ConstructMitochondria()
{
	if (!fPm->ParameterExists(GetFullParmName("Mitochondria", "N")))
		return;

	G4int n = fPm->GetIntegerParameter(GetFullParmName("Mitochondria", "N"));

	// defaults, may be overriden by fPm
	G4double semiAxisA = 0.4 * micrometer;
	G4double semiAxisB = 0.15 * micrometer;
	G4double semiAxisC = 0.15 * micrometer;

	std::vector<G4VPhysicalVolume*> mitochondria = ConstructEllipsoidalChildren("Mitochondria", n, semiAxisA, semiAxisB, semiAxisC, fNucleusRadius);
}

void TsComplexSphericalCell::ConstructLysosomes()
{
	if (!fPm->ParameterExists(GetFullParmName("Lysosomes", "N")))
		return;

	G4int n = fPm->GetIntegerParameter(GetFullParmName("Lysosomes", "N"));

	// defaults, may be overriden by fPm
	G4double semiAxisA = 0.5 * micrometer;
	G4double semiAxisB = 0.5 * micrometer;
	G4double semiAxisC = 0.3 * micrometer;

	std::vector<G4VPhysicalVolume*> lysosomes = ConstructEllipsoidalChildren("Lysosomes", n, semiAxisA, semiAxisB, semiAxisC, fNucleusRadius, nullptr, false);

	if (fPm->ParameterExists(GetFullParmName("Lysosomes/Complexes", "N")))
	{
		auto nChildren = fPm->GetIntegerParameter(GetFullParmName("Lysosomes/Complexes", "N"));
		auto rChildren = fPm->ParameterExists(GetFullParmName("Lysosomes/Complexes/Radius")) ? fPm->GetDoubleParameter(GetFullParmName("Lysosomes/Complexes", "Radius"), "Length") : 0.1 * micrometer;

		// If no nanomaterial densities are given, they are all -0-
		auto nanomaterialDensities = std::vector<G4double>(nChildren, 0.0);
		if (fPm->ParameterExists(GetFullParmName("Lysosomes/Complexes", "NanomaterialDensities")))
		{
			auto arr = fPm->GetUnitlessVector(GetFullParmName("Lysosomes/Complexes", "NanomaterialDensities"));
			nanomaterialDensities = std::vector<G4double>(arr, arr + nChildren);

			if (!std::all_of(nanomaterialDensities.cbegin(), nanomaterialDensities.cend(), [](G4double n) { return n >= 0 && n <= 1; }))
				OutOfRange(GetFullParmName("Lysosomes/Complexes", "NanomaterialDensities"), "must have values within [0, 1]");

			// Add color & drawing style if not specified
			if (!fPm->ParameterExists(GetFullParmName("Lysosomes/Compelxes/Nanomaterial", "Color")))
				fPm->AddParameter("s:" + GetFullParmName("Lysosomes/Complexes/Nanomaterial", "Color"), "\"green\"");
			if (!fPm->ParameterExists(GetFullParmName("Lysosomes/Compelxes/Nanomaterial", "DrawingStyle")))
				fPm->AddParameter("s:" + GetFullParmName("Lysosomes/Complexes/Nanomaterial", "DrawingStyle"), "\"solid\"");
		}

		std::map<G4VPhysicalVolume*, std::vector<G4VPhysicalVolume*>> lysosome_children;
		for (auto it = lysosomes.begin(); it != lysosomes.end(); ++it)
		{
			// auto radius = 0.1 * micrometer;
			lysosome_children[*it] = ConstructSphericalChildren("Lysosomes/Complexes", nChildren, rChildren, 0, *it, true);

			for (size_t i = 0; i < lysosome_children[*it].size(); ++i)
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

				ConstructSphericalChildren("Lysosomes/Complexes/Nanomaterial", nanomaterial_n, nanomaterial_radius, 0, lysosome_children[*it][i]);
			}
			break;
		}
	}
}

void TsComplexSphericalCell::ConstructNucleus()
{
	if (!fPm->ParameterExists(GetFullParmName("Nucleus/Radius")))
		return;

	auto fNucleusRadius = fPm->GetDoubleParameter(GetFullParmName("Nucleus/Radius"), "Length");

	G4RotationMatrix* nucRot = new G4RotationMatrix();
	nucRot->rotateX(0);
	nucRot->rotateY(0);

	auto transNucX = 0 * um;
	auto transNucY = 0 * um;
	auto transNucZ = 0 * um;
	if (fPm->ParameterExists(GetFullParmName("Nucleus/TransX"))) {
		transNucX = fPm->GetDoubleParameter(GetFullParmName("Nucleus/TransX"), "Length");
		if (transNucX > fNucleusRadius)
			OutOfRange(GetFullParmName("Nucleus/TransX"), "sets the nucleus outside of the cell");
	}
	if (fPm->ParameterExists(GetFullParmName("Nucleus/TransY"))) {
		transNucY = fPm->GetDoubleParameter(GetFullParmName("Nucleus/TransY"), "Length");
		if (transNucY > fNucleusRadius)
			OutOfRange(GetFullParmName("Nucleus/TransY"), "sets the nucleus outside of the cell");
	}
	if (fPm->ParameterExists(GetFullParmName("Nucleus/TransZ"))) {
		transNucZ = fPm->GetDoubleParameter(GetFullParmName("Nucleus/TransZ"), "Length");
		if (transNucZ > fNucleusRadius)
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
