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

#include "Digestion.h"


k_Digestion::k_Digestion()
{
}


k_Digestion::~k_Digestion()
{
}


QList<QSharedPointer<k_Polymer> > k_Digestion::digest(QList<QSharedPointer<k_Polymer> > ak_Polymers, k_Enzyme ak_Enzyme)
{
    QList<QSharedPointer<k_Polymer> > lk_Digest;
    printf("Digesting... ");
    for (int i = 0; i < ak_Polymers.size(); ++i)
    {
        QSharedPointer<k_Polymer> lk_pPolymer = ak_Polymers[i];
        for (int k = 1; k < lk_pPolymer->getLength(); ++k)
        {
            if (ak_Enzyme.canCutPolymerAt(lk_pPolymer.data(), k))
            {
                QList<QSharedPointer<k_Polymer> > lk_Parts = lk_pPolymer->cleaveAt(ak_Enzyme.mr_EnzymeParameters.mi_Cut + k);
                lk_Digest.push_back(lk_Parts[0]);
                lk_pPolymer = lk_Parts[1];
            }
        }
        lk_Digest.push_back(lk_pPolymer);
    }
    printf("done.\n");
    return lk_Digest;
}
