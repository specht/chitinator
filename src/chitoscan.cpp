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
	printf("  --massAccuracy (ppm) [float] (default: 5.0)\n");
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
	double ld_MassAccuracy = 5.0;
	
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
	
	li_Index = lk_Arguments.indexOf("--massAccuracy");
	if (li_Index > -1)
	{
		ld_MassAccuracy = QVariant(lk_Arguments[li_Index + 1]).toDouble();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	if (lk_Arguments.size() == 1 && !QFile::exists(lk_Arguments.first()))
	{
		QMap<double, QString> lk_Targets;
		for (int li_DP = 1; li_DP < li_MaxDP; ++li_DP)
		{
			for (int li_DA = 0; li_DA <= li_DP; ++li_DA)
			{
				for (int li_Charge = li_MinCharge; li_Charge <= li_MaxCharge; ++li_Charge)
				{
					for (int li_Isotope = 0; li_Isotope < li_MaxIsotopeCount; ++li_Isotope)
					{
						int li_A = li_DA;
						int li_D = li_DP - li_DA;
						double ld_Mz = (MASS_A * li_A + MASS_D * li_D + MASS_WATER + MASS_HYDROGEN * li_Charge + MASS_NEUTRON * li_Isotope) / li_Charge;
						lk_Targets[ld_Mz] = QString("A%1D%2+%3 (%4+)").arg(li_A).arg(li_D).arg(li_Isotope).arg(li_Charge);
					}
				}
			}
		}
		double ld_QueryMz = lk_Arguments.first().toFloat();
		double ld_MinError = 1e20;
		QMap<double, QString>::const_iterator lk_BestIter;
		QMap<double, QString>::const_iterator lk_Iter = lk_Targets.constBegin();
		for (; lk_Iter != lk_Targets.constEnd(); ++lk_Iter)
		{
			double ld_Mz = lk_Iter.key();
			double ld_Error = fabs(ld_Mz - ld_QueryMz) / ld_Mz * 1000000.0;
			if (ld_Error < ld_MinError)
			{
				ld_MinError = ld_Error;
				lk_BestIter = lk_Iter;
			}
		}
		while (lk_BestIter != lk_Targets.constBegin())
		{
			--lk_BestIter;
			double ld_Mz = lk_BestIter.key();
			double ld_Error = fabs(ld_Mz - ld_QueryMz) / ld_Mz * 1000000.0;
			if (ld_Error > ld_MassAccuracy)
				break;
		}
		bool lb_BreakNext = false;
		while (true)
		{
			double ld_Mz = lk_BestIter.key();
			double ld_Error = fabs(ld_Mz - ld_QueryMz) / ld_Mz * 1000000.0;
			printf("%9.4f %s (%1.2f ppm)\n", lk_BestIter.key(), lk_BestIter.value().toStdString().c_str(), ld_Error);
			if (lb_BreakNext)
				break;
			++lk_BestIter;
			ld_Mz = lk_BestIter.key();
			ld_Error = fabs(ld_Mz - ld_QueryMz) / ld_Mz * 1000000.0;
			if (ld_Error > ld_MassAccuracy)
				lb_BreakNext = true;
		}
		exit(0);
	}
	
	QStringList lk_SpectraFiles;
	foreach (QString ls_Path, lk_Arguments)
		lk_SpectraFiles << ls_Path;
	
	k_ChitoScanner lk_ChitoScanner(r_ScanType::All, QList<tk_IntPair>() << tk_IntPair(1, 1),
		li_MaxIsotopeCount, li_MinCharge, li_MaxCharge, ld_MinSnr, ld_MassAccuracy);
		
	lk_ChitoScanner.scan(lk_SpectraFiles);
}
