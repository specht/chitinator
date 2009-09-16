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

#include "ChitoScanner.h"
#include <ptb/ScanIterator.h>
#include <QtCore>
#include <math.h> 
#include <limits>


k_ChitoScanner::k_ChitoScanner(r_ScanType::Enumeration ae_ScanType,
						   QList<tk_IntPair> ak_MsLevels,
						   int ai_MaxIsotopeCount, int ai_MinCharge, int ai_MaxCharge, 
						   double ad_MinSnr, double ad_MassAccuracy)
	: k_ScanIterator(ae_ScanType, ak_MsLevels)
	, mi_MaxIsotopeCount(ai_MaxIsotopeCount)
	, mi_MinCharge(ai_MinCharge)
	, mi_MaxCharge(ai_MaxCharge)
	, md_MinSnr(ad_MinSnr)
	, md_MassAccuracy(ad_MassAccuracy)
{
	// test parameter sanity
	if (mi_MinCharge < 1 || mi_MinCharge > 100)
	{
		printf("Error: minimum charge (%d) out of range.\n", mi_MinCharge);
		exit(1);
	}
	if (mi_MaxCharge < 1 || mi_MaxCharge > 100)
	{
		printf("Error: maximum charge (%d) out of range.\n", mi_MaxCharge);
		exit(1);
	}
	if (mi_MaxCharge < mi_MinCharge)
	{
		printf("Error: minimum charge (%d) bigger than maximum charge(%d).\n", mi_MinCharge, mi_MaxCharge);
		exit(1);
	}
	if (md_MassAccuracy < 0.0)
	{
		printf("Error: mass accuracy must not be less than zero.\n");
		exit(1);
	}
}


k_ChitoScanner::~k_ChitoScanner()
{
}


void k_ChitoScanner::scan(QStringList ak_SpectraFiles)
{
	// parse all bands
	foreach (QString ls_Path, ak_SpectraFiles)
	{
		ms_CurrentSpot = QFileInfo(ls_Path).baseName();
		
		// parse spot
		this->parseFile(ls_Path);
		
		printf(" done.\n");
	}
}


void k_ChitoScanner::handleScan(r_Scan& ar_Scan)
{
/*	if (QVariant(ar_Scan.ms_Id).toInt() != 3894)
		return;*/
	
	if (ar_Scan.mr_Spectrum.mi_PeaksCount == 0)
	{
		printf("Warning: Empty spectrum (scan #%s @ %1.2f minutes)!\n", ar_Scan.ms_Id.toStdString().c_str(), ar_Scan.md_RetentionTime);
		return;
	}
	
	// find all peaks in this spectrum
	QList<r_Peak> lk_AllPeaks = k_ScanIterator::findAllPeaks(ar_Scan.mr_Spectrum);
	
	// filter by SNR
	QList<r_Peak> lk_TempPeaks;
	lk_TempPeaks.clear();
	foreach (r_Peak lr_Peak, lk_AllPeaks)
		if (lr_Peak.md_Snr >= md_MinSnr)
			lk_TempPeaks << lr_Peak;
	lk_AllPeaks = lk_TempPeaks;
	
	// filter out lower 5%
	double ld_MaxIntensity = 0.0;
	foreach (r_Peak lr_Peak, lk_AllPeaks)
		ld_MaxIntensity = std::max<double>(ld_MaxIntensity, lr_Peak.md_PeakIntensity);
	lk_TempPeaks.clear();
	foreach (r_Peak lr_Peak, lk_AllPeaks)
		if (lr_Peak.md_PeakIntensity >= ld_MaxIntensity * 0.05)
			lk_TempPeaks << lr_Peak;
	lk_AllPeaks = lk_TempPeaks;
	
	printf("MS%d: %4d ", ar_Scan.mi_MsLevel, lk_AllPeaks.size());
	foreach (r_Peak lr_Peak, lk_AllPeaks)
		printf("%7.2f ", lr_Peak.md_PeakMz);
	printf("\n");
}


void k_ChitoScanner::progressFunction(QString as_ScanId, bool)
{
// 	printf("\r%s: scan #%s...", ms_CurrentSpot.toStdString().c_str(), as_ScanId.toStdString().c_str());
}
