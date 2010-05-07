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

#include "Polymer.h"


k_Polymer::k_Polymer(r_PolymerParameters ar_PolymerParameters)
{
    mi_Length = (int)ar_PolymerParameters.md_LengthMean;
    mi_PolymerOffset = 0;
    int li_BufferLength = (mi_Length >> 5) + ((mi_Length & 31) ? 1 : 0);
    mui_pPolymer = RefPtr<unsigned int>(new unsigned int[li_BufferLength]);
    memset(mui_pPolymer.get_Pointer(), 0, li_BufferLength * 4);
    for (int i = 0; i < mi_Length; ++i)
    {
        if (((double)rand() / RAND_MAX) < ar_PolymerParameters.md_DegreeOfAcetylation)
            this->setMonomer(i, A);
        else
            this->setMonomer(i, D);
    }
}


k_Polymer::k_Polymer(k_Polymer* ak_Other_, int ai_Offset, int ai_Length)
    // don't actually copy values, but increase reference count
    // of source polymer array
    : mui_pPolymer(ak_Other_->mui_pPolymer)
    , mi_Length(ai_Length)
    , mi_PolymerOffset(ak_Other_->mi_PolymerOffset + ai_Offset)
{
}


k_Polymer::~k_Polymer()
{
}


int k_Polymer::getLength() const
{
    return mi_Length;
}


int k_Polymer::getDCount()
{
    /*
    int li_DCount = 0;
    unsigned int* lui_Block_ = mui_pPolymer.get_Pointer();
    unsigned int* lui_End_ = lui_Block_ + (mi_Length >> 5);
    unsigned lui_Block;
    // count bits in full unsigned ints
    while (lui_Block_ != lui_End_)
    {
        lui_Block = *lui_Block_;
        // Counting bits set, Brian Kernighan's way 
        // modified from: http://www-graphics.stanford.edu/~seander/bithacks.html
        for (; lui_Block; ++li_DCount)
            lui_Block &= lui_Block - 1;
        ++lui_Block_;
    }
    // count remaining bits
    lui_Block = *lui_Block_;
    // mask out excess bits
    lui_Block &= ((1 << (mi_Length & 31)) - 1);
    for (; lui_Block; ++li_DCount)
        lui_Block &= lui_Block - 1;
    */
    // TODO: speed this up but beware of mi_PolymerOffset!
    int li_DCount = 0;
    for (int i = 0; i < mi_Length; ++i)
        if (this->getMonomer(i) == D)
            ++li_DCount;
        
    return li_DCount;
}


double k_Polymer::getMass()
{
    int li_DCount = this->getDCount();
    return MASS_A * (mi_Length - li_DCount) + MASS_D * li_DCount - MASS_WATER * (mi_Length - 1);
}


unsigned int k_Polymer::slice(int ai_Position, int ai_Length)
{
    // TODO: speed this up
    unsigned int lui_Result = 0;
    for (int i = 0; i < ai_Length; ++i)
    {
        lui_Result <<= 1;
        lui_Result |= this->getMonomer(ai_Position + i);
    }
    return lui_Result;
}


QList<RefPtr<k_Polymer> > k_Polymer::cleaveAt(int ai_Position)
{
    // TODO: speed this up
    QList<RefPtr<k_Polymer> > lk_Result;
    lk_Result.push_back(RefPtr<k_Polymer>(new k_Polymer(this, 0, ai_Position)));
    lk_Result.push_back(RefPtr<k_Polymer>(new k_Polymer(this, ai_Position, mi_Length - ai_Position)));
    return lk_Result;
}


void k_Polymer::setMonomer(int ai_Position, int ai_Flag)
{
    int li_Position = ai_Position + mi_PolymerOffset;
    unsigned int lui_UnsignedInt = mui_pPolymer.get_Pointer()[li_Position >> 5];
    if (ai_Flag == 0)
        lui_UnsignedInt &= ~(1 << (li_Position & 31));
    else
        lui_UnsignedInt |= (1 << (li_Position & 31));
    mui_pPolymer.get_Pointer()[li_Position >> 5] = lui_UnsignedInt;
}


int k_Polymer::getMonomer(int ai_Position)
{
    int li_Position = ai_Position + mi_PolymerOffset;
    return (mui_pPolymer.get_Pointer()[li_Position >> 5] & (1 << (li_Position & 31))) == 0 ? 0 : 1;
}


QString k_Polymer::toString()
{
    QString ls_Result;
    for (int i = 0; i < mi_Length; ++i)
        ls_Result += this->getMonomer(i) == A ? "A" : "D";
    return ls_Result;
}
