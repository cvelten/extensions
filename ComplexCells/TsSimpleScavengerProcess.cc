// Extra Class for TsEmDNAChemistry

#include "TsSimpleScavengerProcess.hh"

#include "G4MolecularConfiguration.hh"
#include "G4Molecule.hh"
#include "G4MoleculeFinder.hh"
#include "G4SystemOfUnits.hh"

#ifndef State
#define State(theXInfo) (GetState<ScavengerState>()->theXInfo)
#endif

TsSimpleScavengerProcess::TsSimpleScavengerProcess(const G4String& aName, G4ProcessType type)
	: G4VITRestDiscreteProcess(aName, type),
	  fIsInitialized(false),
	  fMolecularConfiguration(nullptr), fProductsMolecularConfigurations(), fScavengingCapacity(0), fHasProducts(false)
{
	enableAtRestDoIt = true;
	enableAlongStepDoIt = false;
	enablePostStepDoIt = true;

	SetProcessSubType(60);
	G4VITProcess::SetInstantiateProcessState(false);

	fProposesTimeStep = true;

	pParticleChange = &fParticleChange;
}

TsSimpleScavengerProcess::ScavengerState::ScavengerState()
	: G4ProcessState(), fPreviousTimeAtPreStepPoint(-1)
{}

void TsSimpleScavengerProcess::BuildPhysicsTable(const G4ParticleDefinition&)
{
	fIsInitialized = true;
}

void TsSimpleScavengerProcess::StartTracking(G4Track* track)
{
	G4VProcess::StartTracking(track);
	G4VITProcess::fpState.reset(new ScavengerState());
	G4VITProcess::StartTracking(track);
}

void TsSimpleScavengerProcess::SetReaction(const G4MolecularConfiguration* molConf,
										   double scavengingCapacity)
{
	if (fIsInitialized) {
		G4ExceptionDescription exceptionDescription;
		exceptionDescription << "TsSimpleScavengerProcess was already initialised. ";
		exceptionDescription << "You cannot set a reaction after initialisation.";
		G4Exception("TsSimpleScavengerProcess::SetReaction", "TsSimpleScavengerProcess001",
					FatalErrorInArgument, exceptionDescription);
	}

	fMolecularConfiguration = molConf;
	fProductsMolecularConfigurations.clear();
	fScavengingCapacity = scavengingCapacity;
	fHasProducts = false;
}

void TsSimpleScavengerProcess::SetReaction(const G4MolecularConfiguration* molConf,
										   const std::vector<G4MolecularConfiguration*>& products,
										   double scavengingCapacity)
{
	if (fIsInitialized) {
		G4ExceptionDescription exceptionDescription;
		exceptionDescription << "TsSimpleScavengerProcess was already initialised. ";
		exceptionDescription << "You cannot set a reaction after initialisation.";
		G4Exception("TsSimpleScavengerProcess::SetReaction", "TsSimpleScavengerProcess001",
					FatalErrorInArgument, exceptionDescription);
	}

	G4cerr << "TsSimpleScavengerProcess currently has NOT implemented scavenging with chemical products - it will just kill the track!" << G4endl;

	fMolecularConfiguration = molConf;
	fProductsMolecularConfigurations = std::vector<G4MolecularConfiguration*>(products);
	fScavengingCapacity = scavengingCapacity;
	fHasProducts = true;
}

G4double TsSimpleScavengerProcess::PostStepGetPhysicalInteractionLength(const G4Track& track, G4double, G4ForceCondition* pForceCond)
{
	return AtRestGetPhysicalInteractionLength(track, pForceCond);
}

G4double TsSimpleScavengerProcess::AtRestGetPhysicalInteractionLength(const G4Track& track, G4ForceCondition* pForceCond)
{
	G4Molecule* mol = GetMolecule(track);
	if (mol == nullptr || mol->GetMolecularConfiguration() != fMolecularConfiguration)
		return DBL_MAX;

	// condition is set to "Not Forced"
	*pForceCond = NotForced;

	G4double previousTimeStep = 0;
	if (State(fPreviousTimeAtPreStepPoint) != -1) {
		previousTimeStep = track.GetLocalTime() - State(fPreviousTimeAtPreStepPoint);
	}
	fpState->currentInteractionLength = 1 / fScavengingCapacity;
	if ((previousTimeStep < 0.0) || (fpState->theNumberOfInteractionLengthLeft <= 0.0)) {
		ResetNumberOfInteractionLengthLeft();
	}
	else if (previousTimeStep > 0.0) {
		SubtractNumberOfInteractionLengthLeft(previousTimeStep);
	}

	fpState->currentInteractionLength = GetMeanLifeTime(track, pForceCond);

	State(fPreviousTimeAtPreStepPoint) = track.GetLocalTime();
	return -(fpState->theNumberOfInteractionLengthLeft * fpState->currentInteractionLength);
}

G4double TsSimpleScavengerProcess::GetMeanFreePath(const G4Track& track, G4double, G4ForceCondition* pForceCond)
{
	return GetMeanLifeTime(track, pForceCond);
}

G4double TsSimpleScavengerProcess::GetMeanLifeTime(const G4Track& track, G4ForceCondition* pForceCond)
{
	G4double previousTimeStep = 0;
	if (State(fPreviousTimeAtPreStepPoint) != -1) {
		previousTimeStep = track.GetLocalTime() - State(fPreviousTimeAtPreStepPoint);
	}

	// Probability to have scavenged since last time step DeltaT, after having made it to T, i.e.
	// P[~S(T+dT)|~S(T)] = P[~S(T)|~S(T+dT)] * P[~S(T+dT)] / P[~S(T)]
	// P[~S(T)|~S(T+dT)] = 1 by def.
	// P[~S(T)] = 1 - W(T) = exp(-k_obs * T); W(T) = 1 - exp(-k_obs * T)
	// P[~S(T+dT)|~S(T)] = exp(-k_obs * dT)
	G4double probNotScavenged = std::exp(-fScavengingCapacity * previousTimeStep);

	G4double rnd = G4UniformRand();
	if (rnd >= probNotScavenged) {
		// Compete for PostStepDoIt, but also set condition to forced
		// This will allow other (more dominant processes to work first, before scavenging may be applied)
		auto mol = GetMolecule(track);
		if (verboseLevel > 0) {
			G4cout << " -- limiting step through scavenging: " << mol->GetName()
				   << " at T(global) = " << track.GetGlobalTime() / ns << " ns"
				   << " and T(local) = " << track.GetLocalTime() / ns << " ns"
				   << " after DeltaT = " << previousTimeStep / ns << " ns"
				   << " with P[~S(T+dT)|~S(T)] = " << probNotScavenged
				   << " and RND = " << rnd
				   << " Nint: " << fpState->theNumberOfInteractionLengthLeft << " * " << fpState->currentInteractionLength * ns << " ns^-1"
				   << G4endl;
		}
		return 0;
	}

	return DBL_MAX;
}

G4VParticleChange* TsSimpleScavengerProcess::PostStepDoIt(const G4Track& track, const G4Step& aStep)
{
	return AtRestDoIt(track, aStep);
}

G4VParticleChange* TsSimpleScavengerProcess::AtRestDoIt(const G4Track& track, const G4Step&)
{
	if (verboseLevel > 0) {
		G4cout << " -- TsSimpleScavengerProcess::PostStepDoIt: " << track.GetTrackID() << " (" << GetMolecule(track)->GetName() << ")" << G4endl;
	}

	/*
		const auto pMoleculeA = GetMolecule(trackA)->GetMolecularConfiguration();
		const auto pMoleculeB = GetMolecule(trackB)->GetMolecularConfiguration();

		const auto pReactionData = fMolReactionTable->GetReactionData(pMoleculeA, pMoleculeB);

		const G4int nbProducts = pReactionData->GetNbProducts();

		if (nbProducts)
		{
			const G4double D1 = pMoleculeA->GetDiffusionCoefficient();
			const G4double D2 = pMoleculeB->GetDiffusionCoefficient();
			const G4double sqrD1 = D1 == 0. ? 0. : std::sqrt(D1);
			const G4double sqrD2 = D2 == 0. ? 0. : std::sqrt(D2);
			const G4double inv_numerator = 1./(sqrD1 + sqrD2);
			const G4ThreeVector reactionSite = sqrD2 * inv_numerator * trackA.GetPosition()
											 + sqrD1 * inv_numerator * trackB.GetPosition();

			for (G4int j = 0; j < nbProducts; ++j)
			{
				auto pProduct = new G4Molecule(pReactionData->GetProduct(j));
				auto pProductTrack = pProduct->BuildTrack(trackA.GetGlobalTime(), reactionSite);

				pProductTrack->SetTrackStatus(fAlive);

				G4ITTrackHolder::Instance()->Push(pProductTrack);

				pChanges->AddSecondary(pProductTrack);
				G4MoleculeFinder::Instance()->Push(pProductTrack);
			}
		}
	*/
	if (fHasProducts) {
		const auto molecule = GetMolecule(track)->GetMolecularConfiguration();
		const G4double D = molecule->GetDiffusionCoefficient();
		for (auto it = fProductsMolecularConfigurations.cbegin(); it != fProductsMolecularConfigurations.cend(); ++it) {
			auto product = new G4Molecule(*it);
			auto productTrack = product->BuildTrack(track.GetGlobalTime(), track.GetPosition());

			productTrack->SetTrackStatus(fAlive);
			G4ITTrackHolder::Instance()->Push(productTrack);
			G4MoleculeFinder::Instance()->Push(productTrack);
		}
	}

	fParticleChange.Initialize(track);
	fParticleChange.ProposeTrackStatus(fStopAndKill);
	State(fPreviousTimeAtPreStepPoint) = -1;
	return &fParticleChange;
}
