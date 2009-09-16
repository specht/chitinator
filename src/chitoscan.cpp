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

#include <QtCore>
#include <stdio.h>
#include "ChitoScanner.h"
#include "RefPtr.h"


void printUsageAndExit()
{
	printf("Usage: chitoscan [options] [spectra files]\n");
	printf("Spectra files may be mzData, mzXML or mzML, optionally compressed (.gz|.bz2|.zip).\n");
	printf("Options:\n");
	printf("  --maxIsotopeCount [int] (default: 3)\n");
	printf("  --minCharge [int] (default: 2)\n");
	printf("  --maxCharge [int] (default: 3)\n");
	printf("  --minSnr [float] (default: 2.0)\n");
	printf("  --massAccuracy (ppm) [float] (default: 5.0)\n");
	printf("      This mass accuracy is used to check for the presence of peaks.\n");
	exit(1);
}


bool stringToBool(QString& as_String)
{
	if (as_String == "yes" || as_String == "true" || as_String == "on" || as_String == "enable" || as_String == "enabled")
		return true;
	else if (as_String == "no" || as_String == "false" || as_String == "off" || as_String == "disable" || as_String == "disabled")
		return false;
	else
	{
		printf("Error: unknown boolean value '%s'.\n", as_String.toStdString().c_str());
		exit(1);
	}
};


int main(int ai_ArgumentCount, char** ac_Arguments__)
{
	QStringList lk_Arguments;
	for (int i = 1; i < ai_ArgumentCount; ++i)
		lk_Arguments << ac_Arguments__[i];
		
	r_ScanType::Enumeration le_ScanType = r_ScanType::All;
	int li_MaxIsotopeCount = 3;
	int li_MinCharge = 1;
	int li_MaxCharge = 3;
	double ld_MinSnr = 2.0;
	double ld_MassAccuracy = 5.0;
	
	QStringList lk_SpectraFiles;
	foreach (QString ls_Path, lk_Arguments)
		lk_SpectraFiles << ls_Path;
	
	k_ChitoScanner lk_ChitoScanner(r_ScanType::All, QList<tk_IntPair>() << tk_IntPair(1, 1),
		li_MaxIsotopeCount, li_MinCharge, li_MaxCharge, ld_MinSnr, ld_MassAccuracy);
		
	lk_ChitoScanner.scan(lk_SpectraFiles);
}
