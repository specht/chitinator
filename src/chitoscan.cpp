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
	printf("  --maxDP [int] (default: 10)\n");
	printf("  --minCharge [int] (default: 1)\n");
	printf("  --maxCharge [int] (default: 3)\n");
	printf("  --maxIsotopeCount [int] (default: 3)\n");
	printf("  --minSnr [float] (default: 2.0)\n");
	printf("  --precursorMassAccuracy (ppm) [float] (default: 5.0)\n");
	printf("  --productMassAccuracy (ppm) [float] (default: 700.0)\n");
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
	
	if (lk_Arguments.empty())
		printUsageAndExit();
		
	r_ScanType::Enumeration le_ScanType = r_ScanType::All;
	int li_MaxDP = 10;
	int li_MaxIsotopeCount = 3;
	int li_MinCharge = 1;
	int li_MaxCharge = 3;
	double ld_MinSnr = 2.0;
	double ld_PrecursorMassAccuracy = 5.0;
	double ld_ProductMassAccuracy = 700.0;
	
	// consume options
	int li_Index;
	
	li_Index = lk_Arguments.indexOf("--maxDP");
	if (li_Index > -1)
	{
		li_MaxDP = QVariant(lk_Arguments[li_Index + 1]).toInt();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--maxIsotopeCount");
	if (li_Index > -1)
	{
		li_MaxIsotopeCount = QVariant(lk_Arguments[li_Index + 1]).toInt();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--minCharge");
	if (li_Index > -1)
	{
		li_MinCharge = QVariant(lk_Arguments[li_Index + 1]).toInt();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--maxCharge");
	if (li_Index > -1)
	{
		li_MaxCharge = QVariant(lk_Arguments[li_Index + 1]).toInt();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--minSnr");
	if (li_Index > -1)
	{
		ld_MinSnr = QVariant(lk_Arguments[li_Index + 1]).toDouble();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--precursorMassAccuracy");
	if (li_Index > -1)
	{
		ld_PrecursorMassAccuracy = QVariant(lk_Arguments[li_Index + 1]).toDouble();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--productMassAccuracy");
	if (li_Index > -1)
	{
		ld_ProductMassAccuracy = QVariant(lk_Arguments[li_Index + 1]).toDouble();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	QStringList lk_SpectraFiles;
	foreach (QString ls_Path, lk_Arguments)
		lk_SpectraFiles << ls_Path;
	
	k_ChitoScanner lk_ChitoScanner(r_ScanType::All, 
		QList<tk_IntPair>() << tk_IntPair(1, 10000), li_MaxIsotopeCount, 
		li_MinCharge, li_MaxCharge, ld_MinSnr, ld_PrecursorMassAccuracy, 
		ld_ProductMassAccuracy);
		
	lk_ChitoScanner.scan(lk_SpectraFiles);
}
