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

#ifndef TsSourceDistributedExtended_hh
#define TsSourceDistributedExtended_hh

#include "TsSource.hh"

#include "G4Point3D.hh"

#include <vector>

class TsSourceDistributedExtended : public TsSource
{
public:
	TsSourceDistributedExtended(TsParameterManager* pM, TsSourceManager* psM, G4String sourceName);
	~TsSourceDistributedExtended() = default;

	void UpdateForNewRun(G4bool rebuiltSomeComponents);

	void ResolveParameters();

	void PrepareSampledPoints();

private:
	G4int fNumberOfSourcePoints;
	G4int fNumberOfSourcePointsPerHistory;
	G4int fPreviousNumberOfSourcePoints;
	G4bool fRedistributePointsOnNewRun;

	std::vector<G4Point3D*> fSampledPoints;

	enum DistributionType
	{
		FLAT,
		GAUSSIAN
	};
	DistributionType fPointDistribution;

	G4double fPointDistributionSigma;

public:
	inline const std::vector<G4Point3D*>& GetSampledPoints() const { return fSampledPoints; }
	inline G4int GetNumberOfSourcePointsPerHistory() const { return fNumberOfSourcePointsPerHistory; }
};

#endif
