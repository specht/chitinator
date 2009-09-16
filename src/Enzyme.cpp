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

#include "Enzyme.h"


k_Enzyme::k_Enzyme(r_EnzymeParameters ar_EnzymeParameters)
	: mr_EnzymeParameters(ar_EnzymeParameters)
{
}


k_Enzyme::~k_Enzyme()
{
}


bool k_Enzyme::canCutPolymerAt(k_Polymer* ak_Polymer_, int ai_Position)
{
	// return false if cut is too far to the right
	if (ai_Position + mr_EnzymeParameters.mi_Length > ak_Polymer_->getLength())
		return false;
	// extract subslice from polymer
	unsigned int lui_Slice = ak_Polymer_->slice(ai_Position, mr_EnzymeParameters.mi_Length);
	// apply enzyme mask
	lui_Slice &= mr_EnzymeParameters.mui_Mask;
	// check enzyme pattern
	return lui_Slice == mr_EnzymeParameters.mui_Pattern;
}
