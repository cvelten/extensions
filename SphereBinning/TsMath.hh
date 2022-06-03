// Extra Class for TsBinnedSphere

#ifndef TsMath_hh
#define TsMath_hh

#include "G4Types.hh"

#include <vector>

class TsMath
{
    // Logspace
    // From: https://stackoverflow.com/questions/21429294/is-there-something-like-numpy-logspace-in-c
public:
    template <typename T>
    class Logspace
    {
    private:
        T curValue, base;

    public:
        Logspace(T first, T base) : curValue(first), base(base) {}

        T operator()()
        {
            T retval = curValue;
            curValue *= base;
            return retval;
        }
    };
    static std::vector<G4double> LogspacePy(G4double start, G4double stop, G4int num = 50, G4bool endpoint = true, G4double base = 10);
};

#endif