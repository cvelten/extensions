// Extra Class for VComponentWithChildren

#include "VComponentWithChildren.hh"

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

VComponentWithChildren::VComponentWithChildren(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM, TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name)
	: TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{}

std::vector<G4VPhysicalVolume*> VComponentWithChildren::ConstructSphericalChildren(G4String name, G4double radialOffset)
{
	if (!fPm->ParameterExists(GetFullParmName(name, "N")))
		return std::vector<G4VPhysicalVolume*>();
	auto n = fPm->GetIntegerParameter(GetFullParmName(name, "N"));
	auto radius = fPm->GetDoubleParameter(GetFullParmName(name, "Radius"), "Length");

	return ConstructSphericalChildren(name, n, radius, radialOffset);
}

std::vector<G4VPhysicalVolume*> VComponentWithChildren::ConstructSphericalChildren(G4String name, G4int n, G4double radius, G4double radialOffset, G4VPhysicalVolume* parent, G4bool independentLogicals)
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

std::vector<G4VPhysicalVolume*> VComponentWithChildren::ConstructEllipsoidalChildren(G4String name, G4double radialOffset)
{
	if (!fPm->ParameterExists(GetFullParmName(name, "N")))
		return std::vector<G4VPhysicalVolume*>();
	auto n = fPm->GetIntegerParameter(GetFullParmName(name, "N"));
	auto semiAxisA = fPm->GetDoubleParameter(GetFullParmName(name, "a"), "Length");
	auto semiAxisB = fPm->GetDoubleParameter(GetFullParmName(name, "b"), "Length");
	auto semiAxisC = fPm->GetDoubleParameter(GetFullParmName(name, "c"), "Length");

	return ConstructEllipsoidalChildren(name, n, semiAxisA, semiAxisB, semiAxisC, radialOffset);
}

std::vector<G4VPhysicalVolume*> VComponentWithChildren::ConstructEllipsoidalChildren(G4String name, G4int n, G4double semiAxisA, G4double semiAxisB, G4double semiAxisC, G4double radialOffset, G4VPhysicalVolume* parent, G4bool independentLogicals)
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

std::vector<G4VPhysicalVolume*> VComponentWithChildren::RandomlyPlaceSolid(G4VSolid* solid, G4int n, G4double radialOffset, G4VPhysicalVolume* parent, G4bool independentLogicals, G4String name)
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
		auto nattempts = 0;
		while (overlapCheck)
		{
			++nattempts;

			if (nattempts > fMaxNumAttempts)
			{
				G4cerr << "Topas is exiting due to a serious error in geometry." << G4endl;
				G4cerr << "Component " << GetName() << " has unsuccessfully tried "
					   << nattempts << " times to instantiate '" << name
					   << G4endl;
				fPm->AbortSession(1);
			}

			G4double u = G4UniformRand() * 2 * pi;
			G4double v = std::acos(2 * G4UniformRand() - 1);
			G4double dr = G4UniformRand() * (radius - radialOffset);
			G4double phi = G4UniformRand() * pi;
			G4double psi = G4UniformRand() * pi;
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

G4ThreeVector VComponentWithChildren::GetPointWithinVolume(G4VSolid* solid, G4double minDistanceFromSurface) const
{
	G4ThreeVector surfacePoint = solid->GetPointOnSurface();
	G4ThreeVector surfaceNormal = solid->SurfaceNormal(surfacePoint);
	G4double internalChordLength = solid->DistanceToOut(surfacePoint, -surfaceNormal);

	surfacePoint += minDistanceFromSurface * (-surfaceNormal);

	G4double stepLength = G4UniformRand() * (internalChordLength - 2 * (minDistanceFromSurface));

	return surfacePoint + stepLength * (-surfaceNormal);
}

void VComponentWithChildren::DuplicateSource()
{
	std::vector<G4String> clonedSources, cloningVolumes;

	if (fPm->ParameterExists("DuplicateSourcesForVolumes/Sources")) {
		auto n = fPm->GetVectorLength("DuplicateSourcesForVolumes/Sources");
		auto sources = fPm->GetStringVector("DuplicateSourcesForVolumes/Sources");
		clonedSources = std::vector<G4String>(sources, sources + n);
		delete[] sources;

		n = fPm->GetVectorLength("DuplicateSourcesForVolumes/Volumes");
		auto volumes = fPm->GetStringVector("DuplicateSourcesForVolumes/Volumes");
		cloningVolumes = std::vector<G4String>(volumes, volumes + n);
		delete[] volumes;
	}

	if (fPm->ParameterExists(GetFullParmName("DuplicateSourceForVolume")))
	{
		auto sourceName = fPm->GetStringParameter(GetFullParmName("DuplicateSourceForVolume"));
		if (!std::binary_search(clonedSources.begin(), clonedSources.end(), sourceName))
			clonedSources.push_back(sourceName);
		if (!std::binary_search(cloningVolumes.begin(), cloningVolumes.end(), GetName()))
			cloningVolumes.push_back(GetName());
	}

	if (std::binary_search(cloningVolumes.begin(), cloningVolumes.end(), GetName())) {
		for (auto sourceName : clonedSources) {
			std::vector<G4String> parameterNames;
			fPm->GetParameterNamesStartingWith("So/" + sourceName + "/", &parameterNames);

			for (auto parName : parameterNames) {
				G4String nameWithoutSlash = StringReplace(GetName(), "/", "_");
				G4String newParName = StringReplace(parName, sourceName, sourceName + "_" + nameWithoutSlash);

				if (newParName.find("/Component") != std::string::npos)
				{
					G4String component = fPm->GetStringParameter(parName) + "_" + nameWithoutSlash;
					fPm->AddParameter("s:So/" + sourceName + "_" + nameWithoutSlash + "/Component", '"' + GetName() + '"');
					G4cerr << "s:So/" + sourceName + "_" + nameWithoutSlash + "/Component = \"" << GetName() << '"' << G4endl;
					continue;
				}
				fPm->CloneParameter(parName, newParName);
			}
			G4cout << GetName() << " cloned source parameters for source '" << sourceName << "'." << G4endl;
		}
	}
}

G4String VComponentWithChildren::ConstructParameterName(const char* component, const char* parmName)
{
	// G4String s_component = G4String(component
	G4String fullName = "Ge/" + G4String(component) + "/" + G4String(parmName);
	return fullName;
}

std::string VComponentWithChildren::StringReplace(std::string str, const std::string& from, const std::string& to)
{
	size_t start_pos = str.find(from);
	if (start_pos != std::string::npos)
		str.replace(start_pos, from.length(), to);
	return str;
}