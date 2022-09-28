// Particle Generator for DistributedExtended
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

#include "TsGeneratorDistributedExtended.hh"

#include "TsParameterManager.hh"
#include "TsSourceDistributedExtended.hh"

#include "G4RandomDirection.hh"

TsGeneratorDistributedExtended::TsGeneratorDistributedExtended(TsParameterManager* pM, TsGeometryManager* gM, TsGeneratorManager* pgM, G4String sourceName) : TsVGenerator(pM, gM, pgM, sourceName)
{
	ResolveParameters();
}

void TsGeneratorDistributedExtended::GeneratePrimaries(G4Event* anEvent)
{
	if (CurrentSourceHasGeneratedEnough())
		return;

	G4int sampledPoints;
	std::vector<G4Point3D*>::const_iterator iter;
	for (iter = ((TsSourceDistributedExtended*)fPs)->GetSampledPoints().cbegin(), sampledPoints = 0;
		 iter != ((TsSourceDistributedExtended*)fPs)->GetSampledPoints().cend() && sampledPoints < ((TsSourceDistributedExtended*)fPs)->GetNumberOfSourcePointsPerHistory();
		 ++iter, ++sampledPoints)
	{
		TsPrimaryParticle p;

		G4double costheta = G4RandFlat::shoot(-1., 1);
		G4double sintheta = sqrt(1. - costheta * costheta);
		G4double phi = 2. * CLHEP::pi * G4UniformRand();
		G4double sinphi = sin(phi);
		G4double cosphi = cos(phi);
		G4double px = sintheta * cosphi;
		G4double py = sintheta * sinphi;
		G4double pz = costheta;
		G4double mag = std::sqrt((px * px) + (py * py) + (pz * pz));

		p.dCos1 = px / mag;
		p.dCos2 = py / mag;
		p.dCos3 = pz / mag;

		p.posX = (*iter)->x();
		p.posY = (*iter)->y();
		p.posZ = (*iter)->z();

		SetEnergy(p);
		SetParticleType(p);

		p.weight = 1.;
		p.isNewHistory = true;

		GenerateOnePrimary(anEvent, p);
	}
	AddPrimariesToEvent(anEvent);
}
