// Component for TsBinnedSphere

#include "TsBinnedSphere.hh"

#include "TsParameterManager.hh"
#include "TsParameterisation.hh"
#include "TsMath.hh"

#include "G4GeometryTolerance.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"

// #include <cfloat>
#include <cmath>

TsBinnedSphere::TsBinnedSphere(TsParameterManager *pM, TsExtensionManager *eM, TsMaterialManager *mM, TsGeometryManager *gM,
							   TsVGeometryComponent *parentComponent, G4VPhysicalVolume *parentVolume, G4String &name)
	: TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name),
	  fTotalRMin(0), fTotalSTheta(0), fRadialBinning(RadialBinning::Equal),
	  fConstructParameterized(false), fUseOldParameterizationForSphere(false)
{
	fIsDividable = true;
	fCanCalculateSurfaceArea = true;
	fHasDifferentVolumePerDivision = true;

	fDivisionNames[0] = "R";
	fDivisionNames[1] = "Phi";
	fDivisionNames[2] = "Theta";
	fDivisionUnits[0] = "cm";
	fDivisionUnits[1] = "deg";
	fDivisionUnits[2] = "deg";
}

TsBinnedSphere::~TsBinnedSphere()
{
}

G4VPhysicalVolume *TsBinnedSphere::Construct()
{
	BeginConstruction();

	fPhiDivisionRotations.clear();
	fThetaAreaRatios.clear();

	if (!fIsCopy)
	{
		for (G4int i = 0; i < 3; i++)
		{
			if (fPm->ParameterExists(GetBinParmName(i)))
			{
				fDivisionCounts[i] = fPm->GetIntegerParameter(GetBinParmName(i));
				if (fDivisionCounts[i] <= 0.)
					OutOfRange(GetBinParmName(i), "must be larger than zero");
			}
		}
	}

	if (fPm->ParameterExists(GetFullParmName("VoxelMaterials")))
	{
		G4int nDivisions = fDivisionCounts[0] * fDivisionCounts[1] * fDivisionCounts[2];
		if (fPm->GetVectorLength(GetFullParmName("VoxelMaterials")) != nDivisions)
		{
			G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
			G4cerr << GetName() << " has " << nDivisions << " voxels," << G4endl;
			G4cerr << "but " << GetFullParmName("VoxelMaterials") << " has length of " << fPm->GetVectorLength(GetFullParmName("VoxelMaterials")) << G4endl;
			fPm->AbortSession(1);
		}
		if (fVerbosity > 0)
		{
			G4cout << "Component " << GetName() << " will be constructed as a parameterized volume" << G4endl;
			G4cout << "since individual voxel materials are set." << G4endl;
			if (fDivisionCounts[1] == 1 && fDivisionCounts[2] == 1)
				G4cout << "If all voxel materials are the same consider just using Ge/.../Material since the component is only binned in R." << G4endl;
		}
		fConstructParameterized = true;
	}

	if (fDivisionCounts[0] > 1 || fDivisionCounts[1] > 1 || fDivisionCounts[2] > 1)
		SetHasDividedCylinderOrSphere();

	if (fPm->ParameterExists(GetFullParmName("RMin")))
		fTotalRMin = fPm->GetDoubleParameter(GetFullParmName("RMin"), "Length");
	else
		fTotalRMin = 0.;

	fTotalRMax = fPm->GetDoubleParameter(GetFullParmName("RMax"), "Length");
	fDeltaR = (fTotalRMax - fTotalRMin) / fDivisionCounts[0];

	G4double totalSPhi;
	if (fPm->ParameterExists(GetFullParmName("SPhi")))
		totalSPhi = fPm->GetDoubleParameter(GetFullParmName("SPhi"), "Angle");
	else
		totalSPhi = 0. * deg;

	G4double totalDPhi;
	if (fPm->ParameterExists(GetFullParmName("DPhi")))
		totalDPhi = fPm->GetDoubleParameter(GetFullParmName("DPhi"), "Angle");
	else
		totalDPhi = 360. * deg;

	if (fPm->ParameterExists(GetFullParmName("STheta")))
		fTotalSTheta = fPm->GetDoubleParameter(GetFullParmName("STheta"), "Angle");
	else
		fTotalSTheta = 0. * deg;

	G4double totalDTheta;
	if (fPm->ParameterExists(GetFullParmName("DTheta")))
		totalDTheta = fPm->GetDoubleParameter(GetFullParmName("DTheta"), "Angle");
	else
		totalDTheta = 180. * deg;

	if (fTotalRMin < 0.)
		OutOfRange(GetFullParmName("RMin"), "can not be negative");
	if (fTotalRMax < fTotalRMin)
		OutOfRange(GetFullParmName("RMax"), "can not be less than RMin");
	if (totalSPhi < 0.)
		OutOfRange(GetFullParmName("SPhi"), "can not be negative");
	if (totalSPhi > 360. * deg)
		OutOfRange(GetFullParmName("SPhi"), "can not be greater than 360 degrees");
	if (totalDPhi < 0.)
		OutOfRange(GetFullParmName("DPhi"), "can not be negative");
	if (totalDPhi > 360. * deg)
		OutOfRange(GetFullParmName("DPhi"), "can not be greater than 360 degrees");
	if (fTotalSTheta < 0.)
		OutOfRange(GetFullParmName("STheta"), "can not be negative");
	if (fTotalSTheta > 180. * deg)
		OutOfRange(GetFullParmName("STheta"), "can not be greater than 180 degrees");
	if (totalDTheta < 0.)
		OutOfRange(GetFullParmName("DTheta"), "can not be negative");
	if (totalDTheta > 180. * deg)
		OutOfRange(GetFullParmName("DTheta"), "can not be greater than 180 degrees");

	fFullWidths[0] = fTotalRMax - fTotalRMin;
	fFullWidths[1] = totalDPhi;
	fFullWidths[2] = totalDTheta;

	if (fPm->ParameterExists(GetFullParmName("RadialBinning")))
	{
		auto rbinning = fPm->GetStringParameter(GetFullParmName("RadialBinning"));
		rbinning.toLower();
		if (rbinning.compare("log") == 0)
			fRadialBinning = RadialBinning::Log;
		else if (rbinning.compare("custom") == 0)
			fRadialBinning = RadialBinning::Custom;
		else if (rbinning.compare("equal") == 0)
			fRadialBinning = RadialBinning::Equal;
		else
		{
			G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
			G4cerr << GetFullParmName("RadialBinning") << " has unknown value of " << rbinning << G4endl;
			G4cerr << "Acceptable values are (case-insensitive): Equal, Log, Custom" << G4endl;
			fPm->AbortSession(1);
		}
		GetRBinValues();
	}

	G4Sphere *envelopeSolid = new G4Sphere(fName, fTotalRMin, fTotalRMax, totalSPhi, totalDPhi, fTotalSTheta, totalDTheta);
	fEnvelopeLog = CreateLogicalVolume(envelopeSolid);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

	G4int nDivisions = fDivisionCounts[0] * fDivisionCounts[1] * fDivisionCounts[2];
	if (nDivisions > 1)
	{
		G4String envelopeMaterialName = GetResolvedMaterialName();

		// Calculate divisions
		// const G4double deltaPhi = totalDPhi / fDivisionCounts[1];
		// const G4double deltaTheta = totalDTheta / fDivisionCounts[2];

		G4String divisionName;

		if (fDivisionCounts[1] * fDivisionCounts[2] == 1 && !fConstructParameterized) // only division in R
		{
			G4cerr << "!! Russian Doll !!" << G4endl;

			G4double rMin, rMax;
			G4VSolid *divisionSolid;
			G4LogicalVolume *divisionLog;
			G4VPhysicalVolume *divisionMother = fEnvelopePhys;

			divisionName = fName + "_R_Division";

			// auto isFullSphere = totalDPhi == 360 * deg && totalDTheta == 180 * deg && fTotalRMin == 0;
			for (auto iRadius = fDivisionCounts[0] - 1; iRadius >= 0; --iRadius)
			{
				if (fRadialBinning == RadialBinning::Equal)
				{
					// rMin = fTotalRMin + (isFullSphere ? 0 : 1) * iRadius * FLT_EPSILON;
					// rMin = fTotalRMin + !isFullSphere * (fDivisionCounts[0] - iRadius - 1) * FLT_EPSILON;
					rMin = fTotalRMin;
					rMax = fTotalRMin + (iRadius + 1) * fDeltaR;
				}
				else
				{
					rMin = fTotalRMin;
					rMax = fRBinValues[iRadius];
				}
				// CreatePhysicalVolume(divisionName, ZDivisionLog, fEnvelopePhys, kZAxis, fDivisionCounts[2], 2. * deltaZ);
				divisionSolid = new G4Sphere(divisionName, rMin, rMax, totalSPhi, totalDPhi, fTotalSTheta, totalDTheta);
				divisionLog = CreateLogicalVolume(divisionName, envelopeMaterialName, divisionSolid);
				divisionMother = CreatePhysicalVolume(divisionName, iRadius, true, divisionLog, nullptr, new G4ThreeVector(), divisionMother);
			}
		}
		else
		{
			G4cerr << "!! Parameterization !!" << G4endl;
			fPm->AbortSession(1);
			/*
			// Use Parameterisation for all three divisions
			divisionName = fName + "_Division";
			G4VSolid *divisionSolid = new G4Sphere(divisionName, fTotalRMin, fTotalRMax, totalSPhi, deltaPhi, fTotalSTheta, deltaTheta);
			G4LogicalVolume *divisionLog = CreateLogicalVolume(divisionName, envelopeMaterialName, divisionSolid);
			CreatePhysicalVolume(divisionName, divisionLog, fEnvelopePhys, kUndefined, nDivisions, new TsParameterisation(this));
			// CreatePhysicalVolume(divisionName, divisionLog, fEnvelopeLog, kUndefined, nDivisions, new TsParameterisation(this));
			divisionLog->SetVisAttributes(GetVisAttributes(""));

			fScoringVolume = divisionLog;
			*/
		}
	}

	// Precalculate rotation matrices for all phi rotations.
	// Be sure to do this after first placement so that these rotations do not appear as fRotations(0).
	for (int iPhi = 0; iPhi < fDivisionCounts[1]; iPhi++)
	{
		G4RotationMatrix *rot = new G4RotationMatrix();
		rot->rotateZ(-iPhi * fFullWidths[1] / fDivisionCounts[1]);
		fPhiDivisionRotations.push_back(rot);

		// Also save into virtual geometry component's overall list of rotation vectors,
		// used when update method needs to delete all rotation vectors.
		RegisterRotation(rot);
	}

	PreCalculateThetaRatios();

	// Instantiate children unless this was the World component.
	// The World's children are constructed in a separate call from the Geometry Manager.
	// This makes the World available to Geant4 when our other constructions ask Geant4 for parallel worlds.
	if (fParentVolume)
		InstantiateChildren();

	return fEnvelopePhys;
}

const std::vector<G4double> &TsBinnedSphere::GetRBinValues()
{
	auto abortSession = false;
	if (fRadialBinning == RadialBinning::Custom)
	{
		auto n = fPm->GetVectorLength(GetFullParmName("RBinValues"));
		auto arr = fPm->GetDoubleVector(GetFullParmName("RBinValues"), "Length");
		fRBinValues = std::vector<G4double>(arr, arr + n);
		fRBinValueUnit = fPm->GetUnitOfParameter(GetFullParmName("RBinValues"));

		if (n != fDivisionCounts[0])
		{
			G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
			G4cerr << GetName() << " has " << fDivisionCounts[0] << " radial bins," << G4endl;
			G4cerr << "but " << GetFullParmName("RBinValues") << " has length of " << n << G4endl;
			abortSession = true;
		}
		if (std::adjacent_find(fRBinValues.cbegin(), fRBinValues.cend(), std::greater_equal<G4double>()) != fRBinValues.cend())
		{
			G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
			G4cerr << GetName() << " has not monotonically increasing values in" << G4endl;
			G4cerr << "vector parameter " << GetFullParmName("RBinValues") << G4endl;
			abortSession = true;
		}
		if (fRBinValues.front() <= fTotalRMin)
		{
			G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
			G4cerr << GetName() << " has " << fRBinValues.front() << " as first custom bin," << G4endl;
			G4cerr << "but " << GetFullParmName("RMin") << " is " << fTotalRMin << ", which is greater or equal!" << G4endl;
			abortSession = true;
		}
		if (fRBinValues.back() != fTotalRMax)
		{
			G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
			G4cerr << GetName() << " has " << fRBinValues.back() << " as last custom bin," << G4endl;
			G4cerr << "but " << GetFullParmName("RMax") << " is " << fTotalRMax << G4endl;
			abortSession = true;
		}
		if (abortSession)
			fPm->AbortSession(1);
	}
	else if (fRadialBinning == RadialBinning::Log)
	{
		G4double first = 0, base = 0;
		if (fPm->ParameterExists(GetFullParmName("RBinFirstValue")))
			first = fPm->GetDoubleParameter(GetFullParmName("RBinFirstValue"), "Length");
		else
		{
			G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
			G4cerr << GetName() << " has defined logarithmic binning," << G4endl;
			G4cerr << "but " << GetFullParmName("RBinFirstValue") << " is not set." << G4endl;
			abortSession = true;
		}
		if (first <= 0 || first <= fTotalRMin)
		{
			G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
			G4cerr << GetName() << " has defined logarithmic binning," << G4endl;
			G4cerr << "but " << GetFullParmName("RBinFirstValue") << " is set to an illogical value of " << first << G4endl;
			abortSession = true;
		}
		base = fPm->ParameterExists(GetFullParmName("RBinLogBase")) ? fPm->GetUnitlessParameter(GetFullParmName("RBinLogBase")) : 10;
		if (base <= 1)
		{
			G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
			G4cerr << GetName() << " has defined logarithmic binning," << G4endl;
			G4cerr << "but " << GetFullParmName("RBinLogBase") << " is set to an illogical value of " << base << G4endl;
			abortSession = true;
		}
		if (abortSession)
			fPm->AbortSession(1);

		fRBinValues = TsMath::LogspacePy(std::log(first) / std::log(base), std::log(fTotalRMax) / std::log(base), fDivisionCounts[0], true, base);
		fRBinValueUnit = fPm->GetUnitOfParameter(GetFullParmName("RBinFirstValue"));
		// fRBinValues.clear();
		// fRBinValues.reserve(fDivisionCounts[0]);
		// std::generate_n(std::back_inserter(fRBinValues), fDivisionCounts[0], Logspace<G4double>(first, 10));
	}
	else // fRadialBinning == RadialBinning:Equal or None or whatever
	{
		fRBinValueUnit = fPm->GetUnitOfParameter(GetFullParmName("RMax"));
	}
	if (fVerbosity > 0 && fRadialBinning != RadialBinning::Equal)
	{
		G4cout << "Component " << GetName() << " is using " << (fRadialBinning == RadialBinning::Log ? "logarithmic" : "custom") << " radial binning." << G4endl;
		G4cout << "Individual (RMin, RMax) in " << fRBinValueUnit << ":" << G4endl;
		G4cout << std::setprecision(14);
		// auto rmin = fTotalRMin;
		for (size_t i = 0; i < fRBinValues.size(); ++i)
		{
			G4cout << " i = " << std::setw(std::floor(std::log10(fRBinValues.size())) + 1) << i << " | ";
			G4cout << std::scientific << "(" << (i == 0 ? fTotalRMin : fRBinValues[i - 1]) / fPm->GetUnitValue(fRBinValueUnit) << ", " << fRBinValues[i] / fPm->GetUnitValue(fRBinValueUnit) << ")" << G4endl;
			// rmin = fRBinValues[i];
		}
		G4cout << G4endl;
		if (fVerbosity > 1 || (fPm->ParameterExists(GetFullParmName("PrintRBinValuesString")) && fPm->GetBooleanParameter(GetFullParmName("PrintRBinValuesString"))))
		{
			G4cout << "dv:" << GetFullParmName("RBinValues") << " = " << fRBinValues.size() << " ";
			for (auto it = fRBinValues.cbegin(); it != fRBinValues.cend(); ++it)
				G4cout << std::scientific << *it / fPm->GetUnitValue(fRBinValueUnit) << " ";
			G4cout << fRBinValueUnit << G4endl;
		}
		G4cout << std::setprecision(6);
	}
	return fRBinValues;
}

G4int TsBinnedSphere::GetIndex(G4Step *aStep)
{
	// G4cerr << "TsBinnedSphere::GetIndex(G4Step)" << G4endl;
	if (fDivisionCounts[0] * fDivisionCounts[1] * fDivisionCounts[2] == 1)
	{
		return 0;
	}
	else
	{
		const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();

		G4int index;
		if (fDivisionCounts[1] * fDivisionCounts[2] == 1 && !(fConstructParameterized || fUseOldParameterizationForSphere)) // only divisions in R
																															// index = touchable->GetReplicaNumber(0);
		{
			index = touchable->GetVolume()->GetCopyNo();
			index = touchable->GetReplicaNumber(0);
			// G4cerr << touchable->GetVolume()->GetName() << G4endl;
			// G4cerr << ">> copyno = " << index << G4endl;
		}
		else
		{
			index = touchable->GetReplicaNumber(0);
			// G4cerr << ">> replica = " << index << G4endl;
		}

		if (index < 0 || index >= fDivisionCounts[0] * fDivisionCounts[1] * fDivisionCounts[2])
		{
			if (fNumberOfDetailedErrorReports < fMaximumNumberOfDetailedErrorReports)
			{
				G4cerr << "\nTopas experienced a potentially serious error in scoring." << G4endl;
				G4cerr << "A step in " << GetNameWithCopyId() << " returned index: " << index << " outside of the valid range of 0 to " << fDivisionCounts[0] * fDivisionCounts[1] * fDivisionCounts[2] - 1 << G4endl;
			}
			return -1;
		}

		return index;
	}
}

G4int TsBinnedSphere::GetIndex(G4int iR, G4int iPhi, G4int iTheta)
{
	G4cerr << "TsBinnedSphere::GetIndex(G4int, G4int, G4int)" << G4endl;
	if (fDivisionCounts[0] * fDivisionCounts[1] * fDivisionCounts[2] == 1)
	{
		return 0;
	}
	else
	{
		return iR * fDivisionCounts[1] * fDivisionCounts[2] + iPhi * fDivisionCounts[2] + iTheta;
	}
}

G4int TsBinnedSphere::GetBin(G4int index, G4int iBin)
{
	G4cerr << "TsBinnedSphere::GetBin(G4int, G4int)" << G4endl;
	G4int binR = int(index / (fDivisionCounts[1] * fDivisionCounts[2]));
	if (iBin == 0)
		return binR;

	G4int binPhi = int((index - binR * fDivisionCounts[1] * fDivisionCounts[2]) / fDivisionCounts[2]);
	if (iBin == 1)
		return binPhi;

	G4int binTheta = index - binR * fDivisionCounts[1] * fDivisionCounts[2] - binPhi * fDivisionCounts[2];
	if (iBin == 2)
		return binTheta;

	return -1;
}

TsVGeometryComponent::SurfaceType TsBinnedSphere::GetSurfaceID(G4String surfaceName)
{
	SurfaceType surfaceID;
	G4String surfaceNameLower = surfaceName;
	surfaceNameLower.toLower();
	if (surfaceNameLower == "outercurvedsurface")
		surfaceID = OuterCurvedSurface;
	else if (surfaceNameLower == "innercurvedsurface")
		surfaceID = InnerCurvedSurface;
	else if (surfaceNameLower == "phiplussurface")
		surfaceID = PhiPlusSurface;
	else if (surfaceNameLower == "phiminussurface")
		surfaceID = PhiMinusSurface;
	else if (surfaceNameLower == "thetaplussurface")
		surfaceID = ThetaPlusSurface;
	else if (surfaceNameLower == "thetaminussurface")
		surfaceID = ThetaMinusSurface;
	else if (surfaceNameLower == "anysurface")
		surfaceID = AnySurface;
	else
	{
		surfaceID = None;
		G4cerr << "Topas is exiting due to a serious error in scoring." << G4endl;
		G4cerr << "Scorer name: " << GetName() << " has unknown surface name: " << surfaceName << G4endl;
		fPm->AbortSession(1);
	}
	return surfaceID;
}

G4bool TsBinnedSphere::IsOnBoundary(G4ThreeVector localpos, G4VSolid *solid, SurfaceType surfaceID)
{
	G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

	switch (surfaceID)
	{

	case AnySurface:
		return true;

	case OuterCurvedSurface:
	{
		G4double localR2 = localpos.mag2();
		G4double radius = ((G4Sphere *)solid)->GetOuterRadius();
		return (localR2 > (radius - kCarTolerance) * (radius - kCarTolerance) && localR2 < (radius + kCarTolerance) * (radius + kCarTolerance));
	}

	case InnerCurvedSurface:
	{
		G4double localR2 = localpos.mag2();
		G4double radius;
		if (fDivisionCounts[0] > 1 && fDivisionCounts[1] * fDivisionCounts[2] == 1 && !(fConstructParameterized || fUseOldParameterizationForSphere))
		{
			if (fRadialBinning == RadialBinning::Equal)
				radius = ((G4Sphere *)solid)->GetOuterRadius() - fDeltaR;
			else
				radius = *(std::lower_bound(fRBinValues.cbegin(), fRBinValues.cend(), ((G4Sphere *)solid)->GetOuterRadius()) - 1);
			G4cerr << ((G4Sphere *)solid)->GetOuterRadius() / mm << " mm has calculated inner radius of " << radius / mm << " mm" << G4endl;
		}
		else
			radius = ((G4Sphere *)solid)->GetInnerRadius();
		return (localR2 > (radius - kCarTolerance) * (radius - kCarTolerance) && localR2 < (radius + kCarTolerance) * (radius + kCarTolerance));
	}

	case PhiPlusSurface:
	{
		G4double localPhi = localpos.phi();
		G4double sPhi = ((G4Sphere *)solid)->GetStartPhiAngle();
		G4double dPhi = ((G4Sphere *)solid)->GetDeltaPhiAngle();
		return (localPhi > (sPhi + dPhi - kCarTolerance) && localPhi < (sPhi + dPhi + kCarTolerance));
	}

	case PhiMinusSurface:
	{
		G4double localPhi = localpos.phi();
		G4double sPhi = ((G4Sphere *)solid)->GetStartPhiAngle();
		return (localPhi > (sPhi - kCarTolerance) && localPhi < (sPhi + kCarTolerance));
	}

	case ThetaPlusSurface:
	{
		G4double localTheta = localpos.theta();
		G4double sTheta = ((G4Sphere *)solid)->GetStartThetaAngle();
		G4double dTheta = ((G4Sphere *)solid)->GetDeltaThetaAngle();
		return (localTheta > (sTheta + dTheta - kCarTolerance) && localTheta < (sTheta + dTheta + kCarTolerance));
	}

	case ThetaMinusSurface:
	{
		G4double localTheta = localpos.theta();
		G4double sTheta = ((G4Sphere *)solid)->GetStartThetaAngle();
		return (localTheta > (sTheta - kCarTolerance) && localTheta < (sTheta + kCarTolerance));
	}

	default:
		G4cerr << "Topas is exiting due to a serious error." << G4endl;
		G4cerr << "TsBinnedSphere::IsOnBoundary called for unknown surface of component: " << fName << G4endl;
		fPm->AbortSession(1);
		return false;
	}
}

G4double TsBinnedSphere::GetAreaOfSelectedSurface(G4VSolid *solid, SurfaceType surfaceID, G4int, G4int, G4int copyNo)
{
	// Third index is the copy number. Other two indices are no used for the Sphere.
	// copyNo = iR * nBinsPhi * nBinsTheta + iPhi * nBinsTheta + iTheta
	G4int iR = int(copyNo / (fDivisionCounts[1] * fDivisionCounts[2]));
	G4int iPhi = int((copyNo - iR * fDivisionCounts[1] * fDivisionCounts[2]) / fDivisionCounts[2]);
	G4int iTheta = copyNo - iR * fDivisionCounts[1] * fDivisionCounts[2] - iPhi * fDivisionCounts[2];

	// G4cout << "In GetAreaOfSelectedSurface for iTheta: " << iTheta << " using fThetaAreaRatios[iTheta]: " << fThetaAreaRatios[iTheta] << G4endl;

	G4double r_inner;
	if (fDivisionCounts[0] > 1 && fDivisionCounts[1] * fDivisionCounts[2] == 1 && !(fConstructParameterized || fUseOldParameterizationForSphere))
	{
		if (fRadialBinning == RadialBinning::Equal)
			r_inner = fTotalRMin + iR * fDeltaR;
		else
			r_inner = iR == 0 ? fTotalRMin : fRBinValues[iR];
	}
	else
		r_inner = ((G4Sphere *)solid)->GetInnerRadius();
	G4double r_outer = ((G4Sphere *)solid)->GetOuterRadius();
	G4double delta_phi = ((G4Sphere *)solid)->GetDeltaPhiAngle();
	G4double delta_theta = ((G4Sphere *)solid)->GetDeltaThetaAngle();
	G4double theta_start = ((G4Sphere *)solid)->GetStartThetaAngle();
	G4double theta_end = theta_start + delta_theta;

	switch (surfaceID)
	{
	case OuterCurvedSurface:
		return fThetaAreaRatios[iTheta] * 2. * delta_phi * r_outer * r_outer;

	case InnerCurvedSurface:
		return fThetaAreaRatios[iTheta] * 2. * delta_phi * r_inner * r_inner;

	case PhiPlusSurface:
	case PhiMinusSurface:
		return 0.5 * delta_theta * (r_outer * r_outer - r_inner * r_inner);

	case ThetaPlusSurface:
		return 0.5 * sin(theta_start) * delta_phi * (r_outer * r_outer - r_inner * r_inner);

	case ThetaMinusSurface:
		return 0.5 * sin(theta_end) * delta_phi * (r_outer * r_outer - r_inner * r_inner);

	case AnySurface:
		return (fThetaAreaRatios[iTheta] * 2. * delta_phi * r_outer * r_outer) + (fThetaAreaRatios[iTheta] * 2. * delta_phi * r_inner * r_inner) + (delta_theta * (r_outer * r_outer - r_inner * r_inner)) + (0.5 * sin(theta_start) * delta_phi * (r_outer * r_outer - r_inner * r_inner)) + (0.5 * sin(theta_end) * delta_phi * (r_outer * r_outer - r_inner * r_inner));

	default:
		G4cerr << "Topas is exiting due to a serious error." << G4endl;
		G4cerr << "TsBinnedSphere::GetAreaOfSelectedSurface called for unknown surface of component: " << fName << G4endl;
		fPm->AbortSession(1);
		return 0.;
	}
}

void TsBinnedSphere::PreCalculateThetaRatios()
{
	// Precalculate ratios of theta bin surface areas.
	// Do so by subtracting area of the spherical cap defined by smaller theta boundary
	// from area of the spherical cap defined by the larger theta boundary.
	// If odd number of bins, the middle bin is special. For this one, we calculate half the value then double it.
	//
	// Area of Spherical Cap = 2 Pi R (1 - Sin Theta)
	// Area of Sphere = 4 Pi R
	// Since we only need ratio of spherical cap over sphere, use:
	// Area Factor = (1/2) (1 - Sin Theta)
	//
	G4double halfWidthTheta = fFullWidths[2] / (2. * fDivisionCounts[2]);

	G4double theta;
	G4double areaFactorOfSmallerCap;
	G4double areaFactorOfLargerCap = .5;
	G4double areaFactorOfBin;

	G4bool isOdd;
	if (fDivisionCounts[2] % 2 == 0)
		isOdd = false;
	else
		isOdd = true;

	G4int startTheta = 1;
	if (isOdd)
		startTheta = 0;

	G4int middleBin = fDivisionCounts[2] / 2;
	fThetaAreaRatios.resize(fDivisionCounts[2], 0);

	for (G4int iTheta = startTheta; iTheta < fDivisionCounts[2] / 2. + .5; iTheta++)
	{
		theta = fTotalSTheta + 2. * iTheta * halfWidthTheta;
		if (isOdd)
			theta += halfWidthTheta;

		areaFactorOfSmallerCap = (1. / 2.) * (1. - sin(theta));
		areaFactorOfBin = areaFactorOfLargerCap - areaFactorOfSmallerCap;

		// If odd number of bins, double size of center bin since it actually crosses the equator.
		if (iTheta == 0)
			areaFactorOfBin *= 2.;

		/*G4cout << "TsBinnedSphere::PreCalculateThetaRatios Filling theta bin: " << iTheta << G4endl;
		G4cout << "Theta: " << theta << G4endl;
		G4cout << "areaFactorOfSmallerCap: " << areaFactorOfSmallerCap << G4endl;
		G4cout << "areaFactorOfLargerCap: " << areaFactorOfLargerCap << G4endl;
		G4cout << "areaFactorOfBin: " << areaFactorOfBin << G4endl;
		*/

		if (iTheta == 0)
		{
			fThetaAreaRatios[middleBin] = areaFactorOfBin;
		}
		else
		{
			fThetaAreaRatios[middleBin - iTheta] = areaFactorOfBin;
			if (isOdd)
				fThetaAreaRatios[middleBin + iTheta] = areaFactorOfBin;
			else
				fThetaAreaRatios[middleBin + iTheta - 1] = areaFactorOfBin;
		}

		areaFactorOfLargerCap = areaFactorOfSmallerCap;
	}
}

G4Material *TsBinnedSphere::ComputeMaterial(const G4int repNo, G4VPhysicalVolume *pvol, const G4VTouchable *parent)
{
	if (parent == 0)
		G4cerr << "TsBinnedSphere::ComputeMaterial called with parent touchable zero for repNo: " << repNo << ", pvol: " << pvol->GetName() << G4endl;

	// G4cout << "TsBinnedSphere::ComputeMaterial called with repNo: " << repNo << ", pvol: " << pvol << ", parent: " << parent << G4endl;
	if (parent == 0 || !fHasVariableMaterial)
		return pvol->GetLogicalVolume()->GetMaterial();

	unsigned int matIndex = (*fCurrentMaterialIndex)[repNo];

	G4Material *material = (*fMaterialList)[matIndex];

	pvol->GetLogicalVolume()->SetVisAttributes(fPm->GetColor(GetDefaultMaterialColor(material->GetName())));

	return material;
}

void TsBinnedSphere::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume *physVol) const
{
	// copyNo = iR * nBinsPhi * nBinsTheta + iPhi * nBinsTheta + iTheta
	G4int iR = int(copyNo / (fDivisionCounts[1] * fDivisionCounts[2]));
	G4int iPhi = int((copyNo - iR * fDivisionCounts[1] * fDivisionCounts[2]) / fDivisionCounts[2]);

	if (fPhiDivisionRotations.size() > 0)
		physVol->SetRotation(fPhiDivisionRotations[iPhi]);
}

void TsBinnedSphere::ComputeDimensions(G4Sphere &sphere, const G4int copyNo, const G4VPhysicalVolume *) const
{
	// copyNo = iR * nBinsPhi * nBinsTheta + iPhi * nBinsTheta + iTheta
	G4int iR = int(copyNo / (fDivisionCounts[1] * fDivisionCounts[2]));
	G4int iPhi = int((copyNo - iR * fDivisionCounts[1] * fDivisionCounts[2]) / fDivisionCounts[2]);
	G4int iTheta = copyNo - iR * fDivisionCounts[1] * fDivisionCounts[2] - iPhi * fDivisionCounts[2];

	if (fRadialBinning == RadialBinning::Equal)
	{
		sphere.SetInsideRadius(fTotalRMin + iR * fFullWidths[0] / fDivisionCounts[0]);
		sphere.SetOuterRadius(fTotalRMin + (iR + 1) * fFullWidths[0] / fDivisionCounts[0]);
	}
	else
	{
		sphere.SetInsideRadius(iR == 0 ? fTotalRMin : fRBinValues[iR - 1]);
		sphere.SetOuterRadius(fRBinValues[iR]);
	}
	sphere.SetStartThetaAngle(fTotalSTheta + iTheta * fFullWidths[2] / fDivisionCounts[2]);
}

void TsBinnedSphere::CreateDefaults(TsParameterManager *pM, G4String &childName, G4String &)
{
	G4String parameterName;
	G4String transValue;

	parameterName = "dc:Ge/" + childName + "/RMax";
	transValue = "80. cm";
	pM->AddParameter(parameterName, transValue);

	parameterName = "dc:Ge/" + childName + "/RMin";
	transValue = "0. cm";
	pM->AddParameter(parameterName, transValue);

	parameterName = "ic:Ge/" + childName + "/RBins";
	transValue = "1";
	pM->AddParameter(parameterName, transValue);

	parameterName = "ic:Ge/" + childName + "/PhiBins";
	transValue = "1";
	pM->AddParameter(parameterName, transValue);

	parameterName = "ic:Ge/" + childName + "/ThetaBins";
	transValue = "1";
	pM->AddParameter(parameterName, transValue);

	parameterName = "dc:Ge/" + childName + "/SPhi";
	transValue = "0. deg";
	pM->AddParameter(parameterName, transValue);

	parameterName = "dc:Ge/" + childName + "/DPhi";
	transValue = "360. deg";
	pM->AddParameter(parameterName, transValue);

	parameterName = "dc:Ge/" + childName + "/STheta";
	transValue = "0. deg";
	pM->AddParameter(parameterName, transValue);

	parameterName = "dc:Ge/" + childName + "/DTheta";
	transValue = "180. deg";
	pM->AddParameter(parameterName, transValue);
}
