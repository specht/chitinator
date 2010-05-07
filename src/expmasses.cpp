/*
Copyright (c) 2007-2008 Michael Specht

This file is part of Chitinator.

Chitinator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SimQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Chitinator.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <QtCore>
#include "RefPtr.h"
#include "Digestion.h"
#include "Enzyme.h"
#include "Polymer.h"


int main(int ai_ArgumentCount, char** ac_Arguments__)
{
    for (int j = 1; j <= 1000; ++j)
    {
        QList<double> lk_Masses;
        for (int i = 1; i < j; ++i)
        {
            for (int k = 0; k <= i; ++k)
            {
                for (int li_Charge = 1; li_Charge <= 4; ++li_Charge)
                {
                    int li_DCount = k;
                    int li_ACount = i - li_DCount;
                    double ld_Mass = MASS_A * li_ACount + MASS_D * li_DCount - MASS_WATER * (li_ACount + li_DCount - 1);
                    ld_Mass = (ld_Mass + (li_Charge * MASS_HYDROGEN)) / li_Charge;
                    lk_Masses.push_back(ld_Mass);
                }
            }
        }
        qSort(lk_Masses);
        QList<double> lk_MassDelta;
        double ld_MinDelta = 1e20;
        for (int i = 0; i < lk_Masses.size() - 1; ++i)
        {
            double ld_MassDelta = lk_Masses[i + 1] - lk_Masses[i];
            ld_MinDelta = std::min<double>(ld_MinDelta, ld_MassDelta);
        }
        printf("%d;%1.20lf.\n", j, ld_MinDelta);
    }
    exit(1);
    printf("This is the Chitinator.\n");
    printf("Creating polymers... ");
    r_PolymerParameters lr_PolymerParameters(1000.0, 10.0, 0.64);
    QList<RefPtr<k_Polymer> > lk_Polymers;
    for (int i = 0; i < 20000; ++i)
    {
        RefPtr<k_Polymer> lk_pPolymer(new k_Polymer(lr_PolymerParameters));
        lk_Polymers.push_back(lk_pPolymer);
    }
    printf("done.\n");
    r_EnzymeParameters lr_EnzymeParameters(2, 3, 0, 1);
    k_Enzyme lk_Enzyme(lr_EnzymeParameters);
    k_Digestion lk_Digestion;
    QList<RefPtr<k_Polymer> > lk_Digest = lk_Digestion.digest(lk_Polymers, lk_Enzyme);
    printf("Digested %d polymers into %d products.\n", lk_Polymers.size(), lk_Digest.size());
    int li_MinLength = lk_Digest.first()->mi_Length;
    int li_MaxLength = li_MinLength;
    for (int i = 1; i < lk_Digest.size(); ++i)
    {
        li_MinLength = std::min<int>(li_MinLength, lk_Digest[i]->mi_Length);
        li_MaxLength = std::max<int>(li_MinLength, lk_Digest[i]->mi_Length);
    }
    printf("min length: %d, max length: %d.\n", li_MinLength, li_MaxLength);
}
