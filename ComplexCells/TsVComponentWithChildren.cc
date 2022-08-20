// Extra Class for TsVComponentWithChildren

#include "TsVComponentWithChildren.hh"

#include "G4VPhysicalVolume.hh"
#include "TsParameterManager.hh"

#include "G4Ellipsoid.hh"
#include "G4LogicalVolume.hh"
#include "G4Orb.hh"
#include "G4PhysicalConstants.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <vector>

TsVComponentWithChildren::TsVComponentWithChildren(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM, TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name)
	: TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{}

std::vector<G4VPhysicalVolume*> TsVComponentWithChildren::ConstructSphericalChildren(G4String name, G4double radialOffset)
{
	if (!fPm->ParameterExists(GetFullParmName(name, "N")))
		return std::vector<G4VPhysicalVolume*>();
	auto n = fPm->GetIntegerParameter(GetFullParmName(name, "N"));
	auto radius = fPm->GetDoubleParameter(GetFullParmName(name, "Radius"), "Length");

	return ConstructSphericalChildren(name, n, radius, radialOffset);
}

std::vector<G4VPhysicalVolume*> TsVComponentWithChildren::ConstructSphericalChildren(G4String name, G4int n, G4double radius, G4double radialOffset, G4VPhysicalVolume* parent, G4bool independentLogicals)
{
	if (fPm->ParameterExists(GetFullParmName(name, "N")))
		n = fPm->GetIntegerParameter(GetFullParmName(name, "N"));
	if (n == 0)
		return std::vector<G4VPhysicalVolume*>();

	if (fPm->ParameterExists(GetFullParmName(name, "Radius")))
		radius = fPm->GetDoubleParameter(GetFullParmName(name, "Radius"), "Length");
	else {
		// Post warning
	}

	if (n < 0) OutOfRange(GetFullParmName(name, "N"), "cannot be negative.");
	if (radius <= 0) OutOfRange(GetFullParmName(name, "Radius"), "cannot be less than or equal to zero.");

	G4Sphere* solid = new G4Sphere(name, 0, radius, 0, 360 * deg, 0, 180 * deg);
	// G4LogicalVolume* logical = CreateLogicalVolume(name, solid);

	return RandomlyPlaceSolid(solid, n, radialOffset, parent, independentLogicals, name);
}

std::vector<G4VPhysicalVolume*> TsVComponentWithChildren::ConstructEllipsoidalChildren(G4String name, G4double radialOffset)
{
	if (!fPm->ParameterExists(GetFullParmName(name, "N")))
		return std::vector<G4VPhysicalVolume*>();
	auto n = fPm->GetIntegerParameter(GetFullParmName(name, "N"));
	auto semiAxisA = fPm->GetDoubleParameter(GetFullParmName(name, "a"), "Length");
	auto semiAxisB = fPm->GetDoubleParameter(GetFullParmName(name, "b"), "Length");
	auto semiAxisC = fPm->GetDoubleParameter(GetFullParmName(name, "c"), "Length");

	return ConstructEllipsoidalChildren(name, n, semiAxisA, semiAxisB, semiAxisC, radialOffset);
}

std::vector<G4VPhysicalVolume*> TsVComponentWithChildren::ConstructEllipsoidalChildren(G4String name, G4int n, G4double semiAxisA, G4double semiAxisB, G4double semiAxisC, G4double radialOffset, G4VPhysicalVolume* parent, G4bool independentLogicals)
{
	if (n <= 0 && fPm->ParameterExists(GetFullParmName(name, "N")))
		n = fPm->GetIntegerParameter(GetFullParmName(name, "N"));
	if (n <= 0)
		return std::vector<G4VPhysicalVolume*>();

	if (fPm->ParameterExists(GetFullParmName(name, "a")))
		semiAxisA = fPm->GetDoubleParameter(GetFullParmName(name, "a"), "Length");
	if (fPm->ParameterExists(GetFullParmName(name, "b")))
		semiAxisB = fPm->GetDoubleParameter(GetFullParmName(name, "b"), "Length");
	if (fPm->ParameterExists(GetFullParmName(name, "c")))
		semiAxisC = fPm->GetDoubleParameter(GetFullParmName(name, "c"), "Length");

	if (semiAxisA <= 0) OutOfRange(GetFullParmName(name, "a"), "cannot be less than or equal to zero.");
	if (semiAxisB <= 0) OutOfRange(GetFullParmName(name, "b"), "cannot be less than or equal to zero.");
	if (semiAxisC <= 0) OutOfRange(GetFullParmName(name, "c"), "cannot be less than or equal to zero.");

	G4Ellipsoid* solid = new G4Ellipsoid(name, semiAxisA, semiAxisB, semiAxisC);
	// G4LogicalVolume* logical = CreateLogicalVolume(name, solid);

	return RandomlyPlaceSolid(solid, n, radialOffset, parent, independentLogicals, name);
}

std::vector<G4VPhysicalVolume*> TsVComponentWithChildren::RandomlyPlaceSolid(G4VSolid* solid, G4int n, G4double radialOffset, G4VPhysicalVolume* parent, G4bool independentLogicals, G4String name)
{
	if (parent == nullptr)
		parent = fEnvelopePhys;
	if (name.empty())
		name = solid->GetName();

	G4LogicalVolume* logical = CreateLogicalVolume(name, solid);

	auto volume = parent->GetLogicalVolume()->GetSolid()->GetCubicVolume();
	auto radius = std::cbrt(volume * 3. / 4. / CLHEP::pi);

	auto children = std::vector<G4VPhysicalVolume*>(n);

	for (int j = 0; j < n; ++j)
	{
		G4VPhysicalVolume* physical = nullptr;
		auto overlapCheck = true;
		while (overlapCheck)
		{
			G4double u = G4UniformRand() * 2 * pi;
			G4double v = std::acos(2 * G4UniformRand() - 1);
			G4double dr = G4UniformRand() * (radius - radialOffset);
			G4double phi = G4UniformRand() * 2 * pi;
			G4double psi = G4UniformRand() * 2 * pi;
			G4double x = 0.0;
			G4double y = 0.0;
			G4double z = 0.0;

			x = (radialOffset + dr) * std::cos(u) * std::sin(v);
			y = (radialOffset + dr) * std::sin(u) * std::sin(v);
			z = (radialOffset + dr) * std::cos(v);

			G4ThreeVector* pos = new G4ThreeVector(x, y, z);
			G4RotationMatrix* rot = new G4RotationMatrix();
			rot->rotateX(psi);
			rot->rotateY(phi);

			physical = CreatePhysicalVolume(name, j, false, logical, rot, pos, parent);

			overlapCheck = physical->CheckOverlaps(1000, 0, false);
			if (overlapCheck)
			{
				parent->GetLogicalVolume()->RemoveDaughter(physical);
				physical = nullptr;
				G4cout << "**** Finding new position for volume " << name << ":" << j << " ****" << G4endl;
			}
		}
		if (physical != nullptr)
			children[j] = physical;
		if (independentLogicals)
			logical = CreateLogicalVolume(name, solid);
	}

	return children;
}
