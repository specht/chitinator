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

#pragma once
#include <QtCore>
#include "RefPtr.h"
#include "Polymer.h"


struct r_EnzymeParameters
{
    r_EnzymeParameters(int ai_Length, unsigned int aui_Mask, unsigned int aui_Pattern, int ai_Cut)
        : mi_Length(ai_Length)
        , mui_Mask(aui_Mask)
        , mui_Pattern(aui_Pattern)
        , mi_Cut(ai_Cut)
    {
    };

    int mi_Length;
    unsigned int mui_Mask;
    unsigned int mui_Pattern;
    int mi_Cut;
};


class k_Enzyme
{
public:
    k_Enzyme(r_EnzymeParameters ar_EnzymeParameters);
    virtual ~k_Enzyme();
    
    virtual bool canCutPolymerAt(k_Polymer* ak_Polymer_, int ai_Position);

    r_EnzymeParameters mr_EnzymeParameters;
};
