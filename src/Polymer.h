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


#define A 0
#define D 1
#define MASS_A 203.195
#define MASS_D 161.158
#define MASS_WATER 18.01057
#define MASS_HYDROGEN 1.007947


struct r_PolymerParameters
{
	r_PolymerParameters(double ad_LengthMean, double ad_LengthSd, double ad_DegreeOfAcetylation)
		: md_LengthMean(ad_LengthMean)
		, md_LengthSd(ad_LengthSd)
		, md_DegreeOfAcetylation(ad_DegreeOfAcetylation)
	{
	};
	
	double md_LengthMean;
	double md_LengthSd;
	double md_DegreeOfAcetylation;
};


class k_Polymer
{
public:
	k_Polymer(r_PolymerParameters ar_PolymerParameters);
	k_Polymer(k_Polymer* ak_Other_, int ai_Offset, int ai_Length);
	virtual ~k_Polymer();
	
	int getLength() const;
	double getMass();
	int getDCount();
	unsigned int slice(int ai_Position, int ai_Length);
	QList<RefPtr<k_Polymer> > cleaveAt(int ai_Position);
	
	void setMonomer(int ai_Position, int ai_Flag);
	int getMonomer(int ai_Position);
	
	QString toString();

	int mi_Length;
	int mi_PolymerOffset;
	RefPtr<unsigned int> mui_pPolymer;
};
