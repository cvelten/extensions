#ifndef CvSimpleScavengerProcess_hh
#define CvSimpleScavengerProcess_hh

#include "G4VITRestDiscreteProcess.hh"

class G4MolecularConfiguration;

class CvSimpleScavengerProcess : public G4VITRestDiscreteProcess
{
public:
	CvSimpleScavengerProcess(const G4String& aName = "SimpleScavengerProcess", G4ProcessType type = fDecay);
	~CvSimpleScavengerProcess() override = default;

public:
	void StartTracking(G4Track*) override;

	virtual void SetReaction(const G4MolecularConfiguration*, G4double);
	virtual void SetReaction(const G4MolecularConfiguration*, const std::vector<G4MolecularConfiguration*>&, G4double);

	void BuildPhysicsTable(const G4ParticleDefinition&) override;

	G4double PostStepGetPhysicalInteractionLength(const G4Track&, G4double, G4ForceCondition*) override;
	G4double AtRestGetPhysicalInteractionLength(const G4Track&, G4ForceCondition*) override;

	G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;
	G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&) override;

protected:
	G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*) override;
	G4double GetMeanLifeTime(const G4Track&, G4ForceCondition*) override;

private:
	struct ScavengerState : public G4ProcessState
	{
		ScavengerState();
		~ScavengerState() override = default;
		G4double fPreviousTimeAtPreStepPoint;
	};

	G4VParticleChange fParticleChange;

	G4bool fIsInitialized;

	const G4MolecularConfiguration* fMolecularConfiguration;
	std::vector<G4MolecularConfiguration*> fProductsMolecularConfigurations;
	G4double fScavengingCapacity;
	G4bool fHasProducts;
};

#endif
