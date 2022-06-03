// Extra Class for TsBinnedSphere

#include "TsMath.hh"

std::vector<G4double> TsMath::LogspacePy(G4double start, G4double stop, G4int num, G4bool endpoint, G4double base)
{
    G4double realStart = pow(base, start);
    G4double realBase = endpoint ? pow(base, (stop - start) / (num - 1)) : pow(base, (stop - start) / num);

    std::vector<G4double> retval;
    retval.reserve(num);
    std::generate_n(std::back_inserter(retval), num, Logspace<G4double>(realStart, realBase));
    return retval;
}