// Particle Source for DistributedExtended
//
// ********************************************************************
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * TOPAS collaboration.                                             *
// * Use or redistribution of this code is not permitted without the  *
// * explicit approval of the TOPAS collaboration.                    *
// * Contact: Joseph Perl, perl@slac.stanford.edu                     *
// *                                                                  *
// ********************************************************************
//

#include "TsSourceDistributedExtended.hh"

#include "TsParameterManager.hh"
#include "TsVGeometryComponent.hh"

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisExtent.hh"
#include "Randomize.hh"

TsSourceDistributedExtended::TsSourceDistributedExtended(TsParameterManager* pM, TsSourceManager* psM, G4String sourceName)
	: TsSource(pM, psM, sourceName),
	  fNumberOfSourcePoints(0), fNumberOfSourcePointsPerHistory(1), fPreviousNumberOfSourcePoints(-1), fRedistributePointsOnNewRun(false),
	  fRestrictToVolumes()
{
	ResolveParameters();
}

void TsSourceDistributedExtended::UpdateForNewRun(G4bool rebuiltSomeComponents)
{
	TsSource::UpdateForNewRun(rebuiltSomeComponents);

	if (fRedistributePointsOnNewRun || fNumberOfSourcePoints != fPreviousNumberOfSourcePoints)
		PrepareSampledPoints();
}

void TsSourceDistributedExtended::ResolveParameters()
{
	TsSource::ResolveParameters();

	if (fPm->ParameterExists(GetFullParmName("RedistributePointsOnNewRun")))
		fRedistributePointsOnNewRun = fPm->GetBooleanParameter(GetFullParmName("RedistributePointsOnNewRun"));

	fNumberOfSourcePoints = fPm->GetIntegerParameter(GetFullParmName("NumberOfSourcePoints"));
	if (fNumberOfSourcePoints < 1) {
		G4cerr << "TOPAS is exiting due to a serious error in the distributed source: " << GetName() << G4endl;
		G4cerr << "NumberOfSourcePoints must be greater than zero." << G4endl;
		fPm->AbortSession(1);
	}

	if (fPm->ParameterExists(GetFullParmName("NumberOfSourcePointsPerHistory")))
		fNumberOfSourcePointsPerHistory = fPm->GetIntegerParameter(GetFullParmName("NumberOfSourcePointsPerHistory"));
	if (fNumberOfSourcePointsPerHistory < 1 || fNumberOfSourcePointsPerHistory > fNumberOfSourcePoints) {
		G4cerr << "TOPAS is exiting due to a serious error in the distributed source: " << GetName() << G4endl;
		G4cerr << "NumberOfSourcePointsPerHistory must be greater than zero and smaller than NumberOfSourcePoints." << G4endl;
		fPm->AbortSession(1);
	}

	G4String dist;
	if (fPm->ParameterExists(GetFullParmName("PointDistribution")))
		dist = fPm->GetStringParameter(GetFullParmName("PointDistribution"));
	else
		dist = "flat";

	dist.toLower();
	if (dist == "flat") {
		fPointDistribution = FLAT;
	}
	else if (dist == "gaussian") {
		fPointDistribution = GAUSSIAN;
		fPointDistributionSigma = fPm->GetDoubleParameter(GetFullParmName("PointDistributionSigma"), "Length");
		if (fPointDistributionSigma <= 0.) {
			G4cerr << GetFullParmName("PointDistributionSigma") << " must be greater than zero." << G4endl;
			fPm->AbortSession(1);
		}
	}
	else {
		G4cerr << "Particle source \"" << fSourceName << "\" has unknown PointDistribution \""
			   << fPm->GetStringParameter(GetFullParmName("PointDistribution")) << "\"" << G4endl;
		G4cerr << "Accepted values are Flat and Gaussian." << G4endl;
		fPm->AbortSession(1);
	}

	if (fPm->ParameterExists(GetFullParmName("RestrictToVolumesContaining"))) {
		auto n = fPm->GetVectorLength(GetFullParmName("RestrictToVolumesContaining"));
		auto arr = fPm->GetStringVector(GetFullParmName("RestrictToVolumesContaining"));
		for (auto i = 0; i < n; ++i)
		{
			G4String s = *(arr + i);
			s.toLower();
			fRestrictToVolumes.push_back(s);
		}
		delete[] arr;
	}
}

void TsSourceDistributedExtended::PrepareSampledPoints()
{
	fPreviousNumberOfSourcePoints = fNumberOfSourcePoints;
	fSampledPoints.clear();

	G4int maxNumberOfPointsToSample = 1000000;
	if (fPm->ParameterExists(GetFullParmName("MaxNumberOfPointsToSample")))
		maxNumberOfPointsToSample = fPm->GetIntegerParameter(GetFullParmName("MaxNumberOfPointsToSample"));

	G4bool recursivelyIncludeChildren = false;
	if (fPm->ParameterExists(GetFullParmName("RecursivelyIncludeChildren")))
		recursivelyIncludeChildren = fPm->GetBooleanParameter(GetFullParmName("RecursivelyIncludeChildren"));

	G4VisExtent myExtent = fComponent->GetExtent();
	G4double xMin = myExtent.GetXmin();
	G4double xMax = myExtent.GetXmax();
	G4double yMin = myExtent.GetYmin();
	G4double yMax = myExtent.GetYmax();
	G4double zMin = myExtent.GetZmin();
	G4double zMax = myExtent.GetZmax();

	G4TransportationManager* transportationManager = G4TransportationManager::GetTransportationManager();
	G4Navigator* navigator = transportationManager->GetNavigator(transportationManager->GetParallelWorld(fComponent->GetWorldName()));

	G4double testX;
	G4double testY;
	G4double testZ;

	G4VPhysicalVolume* foundVolume;

	std::vector<G4VPhysicalVolume*> volumes = fComponent->GetAllPhysicalVolumes(recursivelyIncludeChildren);
	if (fRestrictToVolumes.size() > 0) {
		for (auto it = volumes.begin(); it != volumes.end();)
		{
			G4String volumeName = (*it)->GetName();
			volumeName.toLower();

			auto found = false;
			for (auto substr : fRestrictToVolumes) {
				if (volumeName.find(substr) != std::string::npos) {
					found = true;
					break;
				}
			}
			if (!found)
				volumes.erase(it);
			else
				++it;
		}
	}

	if (volumes.size() == 0) {
		G4cerr << "TOPAS is exiting due to a serious error in the distributed source: " << GetName() << G4endl;
		G4cerr << "The number of volumes for findin points must be greater than zero." << G4endl;
		fPm->AbortSession(1);
	}

	for (G4int iPoint = 0; iPoint < fNumberOfSourcePoints; iPoint++)
	{
		G4int counter = 0;
		G4bool foundPoint = false;
		while (!foundPoint) {
			// G4cout << "looking for point: " << counter << G4endl;
			//  Randomly sample a point in a big cube
			if (fPointDistribution == FLAT) {
				testX = G4RandFlat::shoot(xMin, xMax);
				testY = G4RandFlat::shoot(yMin, yMax);
				testZ = G4RandFlat::shoot(zMin, zMax);
			}
			else {
				testX = G4RandGauss::shoot(0., fPointDistributionSigma);
				testY = G4RandGauss::shoot(0., fPointDistributionSigma);
				testZ = G4RandGauss::shoot(0., fPointDistributionSigma);
			}

			// Check whether they are inside any of the component's volumes
			foundVolume = navigator->LocateGlobalPointAndSetup(G4ThreeVector(testX, testY, testZ));
			// G4cout << "foundVolume: " << foundVolume->GetName() << G4endl;

			// Protect against case where extent has extended outside of world
			if (foundVolume) {
				if (foundVolume->GetName() != "World") {
					for (size_t t = 0; !foundPoint && t < volumes.size(); t++)
						foundPoint = foundVolume == volumes[t];
				}

				if (counter++ > maxNumberOfPointsToSample) {
					G4cerr << "TOPAS is exiting due to a serious error in the distributed source: " << GetName() << G4endl;
					G4cerr << "In " << maxNumberOfPointsToSample << " attempts, we have never found a suitable starting position." << G4endl;
					fPm->AbortSession(1);
				}
			}
		}
		fSampledPoints.push_back(new G4Point3D(testX, testY, testZ));
	}
}
