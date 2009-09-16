/*
Copyright (c) 2007-2008 Michael Specht

This file is part of SimQuant.

SimQuant is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SimQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SimQuant.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#include <QtCore>
#include <ptb/ScanIterator.h>


#define HYDROGEN 1.007947
#define WATER 18.01057
#define NEUTRON 1.002
#define MASS_D 179.172
#define MASS_A 221.210


class k_ChitoScanner: public k_ScanIterator
{
public:
	k_ChitoScanner(r_ScanType::Enumeration ae_ScanType,
				 QList<tk_IntPair> ak_MsLevels,
				 int ai_MaxIsotopeCount, int ai_MinCharge, int ai_MaxCharge, 
				 double ad_MinSnr, double ad_MassAccuracy);
	virtual ~k_ChitoScanner();
	
	// quantify takes a list of spectra files and a hash of (peptide => protein) entries
	virtual void scan(QStringList ak_SpectraFiles);
	virtual void handleScan(r_Scan& ar_Scan);
	virtual void progressFunction(QString as_ScanId, bool ab_InterestingScan);
	
protected:
	int mi_MaxIsotopeCount;
	int mi_MinCharge;
	int mi_MaxCharge;
	double md_MinSnr;
	double md_MassAccuracy;
	QString ms_CurrentSpot;
};
