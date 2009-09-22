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
						   double ad_MinSnr, double ad_PrecursorMassAccuracy,
						   double ad_ProductMassAccuracy)
	: k_ScanIterator(ae_ScanType, ak_MsLevels)
	, mi_MaxIsotopeCount(ai_MaxIsotopeCount)
	, mi_MinCharge(ai_MinCharge)
	, mi_MaxCharge(ai_MaxCharge)
	, md_MinSnr(ad_MinSnr)
	, md_PrecursorMassAccuracy(ad_PrecursorMassAccuracy)
	, md_ProductMassAccuracy(ad_ProductMassAccuracy)
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
	if (md_PrecursorMassAccuracy < 0.0)
	{
		printf("Error: mass accuracy must not be less than zero.\n");
		exit(1);
	}
	if (md_ProductMassAccuracy< 0.0)
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
		
		mk_MS1MassAccuracies.clear();
		mi_MS2PeakHitCount = 0;
		mi_MS2PeakHitCountWithIsotope = 0;
		mk_DAAmounts.clear();
		mk_DAAmountsUnweighted.clear();
		mk_DPs.clear();
	
		// parse spot
		this->parseFile(ls_Path);
		
		printf(" done.\n");
		
		qSort(mk_DPs.begin(), mk_DPs.end());
		qSort(mk_MS1MassAccuracies.begin(), mk_MS1MassAccuracies.end());
		
		printf("MS1 mass accuracies (%d entries):\n", mk_MS1MassAccuracies.size());
			for (int i = 8; i <= 10; ++i)
			printf("%3i%% %9.4f / ", i * 10, mk_MS1MassAccuracies[(mk_MS1MassAccuracies.size() - 1) * i / 10]);
		printf("\n\n");
		
		double ld_ProductDA = 0.0;
		double ld_AmountSum = 0.0;
		foreach (tk_IntPair lk_Pair, mk_DAAmounts.keys())
		{
			double ld_Amount = mk_DAAmounts[lk_Pair];
			ld_AmountSum += ld_Amount;
			double ld_DA = (double)lk_Pair.first / (double)(lk_Pair.first + lk_Pair.second);
			ld_ProductDA += ld_DA * ld_Amount;
		}
		ld_ProductDA /= ld_AmountSum;
		printf("Product DA (weighted)  : %1.2f%%\n", ld_ProductDA * 100.0);
		ld_ProductDA = 0.0;
		ld_AmountSum = 0.0;
		foreach (tk_IntPair lk_Pair, mk_DAAmountsUnweighted.keys())
		{
			double ld_Amount = mk_DAAmountsUnweighted[lk_Pair];
			ld_AmountSum += ld_Amount;
			double ld_DA = (double)lk_Pair.first / (double)(lk_Pair.first + lk_Pair.second);
			ld_ProductDA += ld_DA * ld_Amount;
		}
		ld_ProductDA /= ld_AmountSum;
		printf("Product DA (unweighted): %1.2f%%\n\n", ld_ProductDA * 100.0);
		
		printf("Product DP median : %d\n", (int)(mk_DPs[mk_DPs.size() / 2]));
		double ld_DPMean = 0.0;
		double ld_DPSD = 0.0;
		calculateMeanAndStandardDeviation(mk_DPs, &ld_DPMean, &ld_DPSD);
		printf("Product DP mean/SD: %1.2f / %1.2f\n\n", ld_DPMean, ld_DPSD);
		
		printf("MS2 peak hit with isotope: %1.2f%%.\n\n", (double)mi_MS2PeakHitCountWithIsotope / mi_MS2PeakHitCount * 100.0);
	}
}


void k_ChitoScanner::handleScan(r_Scan& ar_Scan)
{
/*	if (QVariant(ar_Scan.ms_Id).toInt() != 165)
		return;*/
	
	if (ar_Scan.mr_Spectrum.mi_PeaksCount == 0)
	{
		printf("Warning: Empty spectrum (scan #%s @ %1.2f minutes)!\n", ar_Scan.ms_Id.toStdString().c_str(), ar_Scan.md_RetentionTime);
		return;
	}
	
	// find all peaks in this spectrum
	QList<r_Peak> lk_AllPeaks = k_ScanIterator::findAllPeaks(ar_Scan.mr_Spectrum);
	
	// filter out lower 5%
	double ld_MaxIntensity = 0.0;
	foreach (r_Peak lr_Peak, lk_AllPeaks)
		ld_MaxIntensity = std::max<double>(ld_MaxIntensity, lr_Peak.md_PeakIntensity);
	QList<r_Peak> lk_TempPeaks;
	lk_TempPeaks.clear();
	foreach (r_Peak lr_Peak, lk_AllPeaks)
		if (lr_Peak.md_PeakIntensity >= ld_MaxIntensity * 0.05)
			lk_TempPeaks << lr_Peak;
	lk_AllPeaks = lk_TempPeaks;
	
	// filter by SNR
	lk_TempPeaks.clear();
	foreach (r_Peak lr_Peak, lk_AllPeaks)
		if (lr_Peak.md_Snr >= md_MinSnr)
			lk_TempPeaks << lr_Peak;
	lk_AllPeaks = lk_TempPeaks;

	if (ar_Scan.mi_MsLevel == 1)
	{
		foreach (r_Peak lr_Peak, lk_AllPeaks)
		{
			r_OligoHit lr_Hit = 
				matchOligomer(lr_Peak.md_PeakMz, 
							   md_PrecursorMassAccuracy, 
							   gk_Undefined, gk_Undefined, gk_Undefined, 
							   gk_Undefined, tk_IntPair(1, 1));
			if (lr_Hit.mb_IsGood)
			{
				mk_MS1MassAccuracies << lr_Hit.md_MassAccuracy;
				tk_IntPair lk_Pair(lr_Hit.mi_A, lr_Hit.mi_D);
				
				if (!mk_DAAmounts.contains(lk_Pair))
					mk_DAAmounts[lk_Pair] = 0.0;
				mk_DAAmounts[lk_Pair] += lr_Peak.md_PeakIntensity;
				
				if (!mk_DAAmountsUnweighted.contains(lk_Pair))
					mk_DAAmountsUnweighted[lk_Pair] = 0.0;
				mk_DAAmountsUnweighted[lk_Pair] += 1.0;
				
				mk_DPs << (lr_Hit.mi_A + lr_Hit.mi_D);
			}
// 			printf("%1.4f %s\n", lr_Peak.md_PeakMz, lr_Hit.toString().toStdString().c_str());
		}
	}
	
	if (ar_Scan.mi_MsLevel == 2)
	{
		r_Precursor lr_Precursor = ar_Scan.mk_Precursors.first();
		r_OligoHit lr_PrecursorHit = 
			matchOligomer(lr_Precursor.md_Mz, md_PrecursorMassAccuracy,
						   gk_Undefined, gk_Undefined,
						   tk_IntPair(lr_Precursor.mi_ChargeState,
									   lr_Precursor.mi_ChargeState));
		if (lr_PrecursorHit.mb_IsGood)
		{
			printf("precursor is %1.2f %d+ %s\n", lr_Precursor.md_Mz, 
					lr_Precursor.mi_ChargeState, 
					lr_PrecursorHit.toString().toStdString().c_str());
			foreach (r_Peak lr_Peak, lk_AllPeaks)
			{
				r_OligoHit lr_Hit = matchOligomer(lr_Peak.md_PeakMz, 
					md_ProductMassAccuracy, tk_IntPair(0, 0), 
					gk_Undefined, tk_IntPair(1, 1), tk_IntPair(0, 1), 
					gk_Undefined, true);
				QString ls_Fragment;
				if (lr_Hit.mb_IsGood)
				{
					++mi_MS2PeakHitCount;
					if (lr_Hit.mi_Isotope > 0)
						++mi_MS2PeakHitCountWithIsotope;
					ls_Fragment = "[";
					for (int i = 0; i < lr_Hit.mi_A; ++i)
						ls_Fragment += "A";
					for (int i = 0; i < lr_Hit.mi_D; ++i)
						ls_Fragment += "D";
					ls_Fragment += "]";
					if (lr_Hit.mi_Label == 1)
						ls_Fragment += "*";
	 				printf("%9.4f / %9.4f  %-12s %s\n", lr_Peak.md_PeakMz, lr_Peak.md_PeakIntensity, ls_Fragment.toStdString().c_str(), lr_Hit.toString().toStdString().c_str());
				}
			}
	// 		printf("\n");
		}
	}
// 	printf("%s: scan #%s: ", ms_CurrentSpot.toStdString().c_str(), ar_Scan.ms_Id.toStdString().c_str());	
// 	printf("MS%d: %4d peaks ", ar_Scan.mi_MsLevel, lk_AllPeaks.size());
}


void k_ChitoScanner::progressFunction(QString as_ScanId, bool)
{
	printf("\r%s: scan #%s...", ms_CurrentSpot.toStdString().c_str(), as_ScanId.toStdString().c_str());
}


double k_ChitoScanner::calculateOligomerMass(int ai_IonizationType, int ai_DP, int ai_DA, int ai_Charge, int ai_Isotope, int ai_Label, double ad_Loss)
{
	int li_A = ai_DA;
	int li_D = ai_DP - ai_DA;
	return (MASS_A * li_A + MASS_D * li_D + MASS_WATER - ad_Loss + (ai_IonizationType == 0 ? MASS_HYDROGEN : MASS_SODIUM) * ai_Charge + MASS_NEUTRON * ai_Isotope + 2.0 * MASS_HYDROGEN * ai_Label) / ai_Charge;
}


r_OligoHit 
k_ChitoScanner::matchOligomer(double ad_Mz, 
							   double ad_MassAccuracy,
							   tk_IntPair ak_IonizationType,
							   tk_IntPair ak_DP,
							   tk_IntPair ak_Charge,
							   tk_IntPair ak_Isotope,
							   tk_IntPair ak_Label,
							   bool ab_CheckWaterLoss)
{
	// if ai_Charge < 1, then charge is not defined
	QString ls_Result = "(no match)";
	tk_IntPair lk_IonizationType(0, 1);
	tk_IntPair lk_DP(1, 10);
	tk_IntPair lk_Charge(mi_MinCharge, mi_MaxCharge);
	tk_IntPair lk_Isotope(0, mi_MaxIsotopeCount - 1);
	tk_IntPair lk_Label(0, 1);
	tk_IntPair lk_Loss(0, 0);
	if (ak_IonizationType != tk_IntPair(-1, -1))
		lk_IonizationType = ak_IonizationType;
	if (ak_DP != gk_Undefined)
		lk_DP = ak_DP;
	if (ak_Charge != gk_Undefined)
		lk_Charge = ak_Charge;
	if (ak_Isotope != gk_Undefined)
		lk_Isotope = ak_Isotope;
	if (ak_Label != gk_Undefined)
		lk_Label = ak_Label;
	if (ab_CheckWaterLoss)
		lk_Loss = tk_IntPair(0, 1);
	
	QString ls_HashKey = QString("it%1-%2/dp%3-%4/ch%5-%6/iso%7-%8/la%9-%10/ls%11-%12")
		.arg(lk_IonizationType.first).arg(lk_IonizationType.second)
		.arg(lk_DP.first).arg(lk_DP.second)
		.arg(lk_Charge.first).arg(lk_Charge.second)
		.arg(lk_Isotope.first).arg(lk_Isotope.second)
		.arg(lk_Label.first).arg(lk_Label.second)
		.arg(lk_Loss.first).arg(lk_Loss.second);

	if (!mk_TargetCache.contains(ls_HashKey))
	{
		for (int li_IonizationType = lk_IonizationType.first; li_IonizationType <= lk_IonizationType.second; ++li_IonizationType)
		{
			for (int li_DP = lk_DP.first; li_DP <= lk_DP.second; ++li_DP)
			{
				for (int li_DA = 0; li_DA <= li_DP; ++li_DA)
				{
					for (int li_Charge = lk_Charge.first; li_Charge <= lk_Charge.second; ++li_Charge)
					{
						for (int li_Isotope = lk_Isotope.first; li_Isotope <= lk_Isotope.second; ++li_Isotope)
						{
							for (int li_MsnLabel = lk_Label.first; li_MsnLabel <= lk_Label.second; ++li_MsnLabel)
							{
								for (int li_Loss = lk_Loss.first; li_Loss <= lk_Loss.second; ++li_Loss)
								{
									int li_A = li_DA;
									int li_D = li_DP - li_DA;
									double ld_Loss = li_Loss == 0 ? 0.0 : MASS_WATER;
									double ld_Mz = calculateOligomerMass(li_IonizationType, li_DP, li_DA, li_Charge, li_Isotope, li_MsnLabel, ld_Loss);
									mk_TargetCache[ls_HashKey][ld_Mz] = r_OligoHit(li_A, li_D, li_Isotope, li_Charge, li_IonizationType, li_MsnLabel, ld_Mz, ld_Loss);
								}
							}
						}
					}
				}
			}
		}
	}
	
	r_OligoHit lr_Hit;
	QMap<double, r_OligoHit> lk_Targets = mk_TargetCache[ls_HashKey];
	double ld_MinError = 1e20;
	QMap<double, r_OligoHit>::const_iterator lk_BestIter;
	QMap<double, r_OligoHit>::const_iterator lk_Iter = lk_Targets.constBegin();
	for (; lk_Iter != lk_Targets.constEnd(); ++lk_Iter)
	{
		double ld_Mz = lk_Iter.key();
		double ld_Error = fabs(ld_Mz - ad_Mz) / ld_Mz * 1000000.0;
		if (ld_Error < ld_MinError)
		{
			ld_MinError = ld_Error;
			lk_BestIter = lk_Iter;
		}
	}
	
	if (ld_MinError <= ad_MassAccuracy)
	{
		double ld_NextBestError = 1e20;
		if (lk_BestIter != lk_Targets.constBegin())
		{
			double ld_Mz = (lk_BestIter - 1).key();
			double ld_Error = fabs(ld_Mz - ad_Mz) / ld_Mz * 1000000.0;
			ld_NextBestError = std::min<double>(ld_NextBestError, ld_Error);
		}
		if (lk_BestIter != lk_Targets.constEnd())
		{
			double ld_Mz = (lk_BestIter + 1).key();
			double ld_Error = fabs(ld_Mz - ad_Mz) / ld_Mz * 1000000.0;
			ld_NextBestError = std::min<double>(ld_NextBestError, ld_Error);
		}
		lr_Hit = lk_BestIter.value();
		lr_Hit.md_MassAccuracy = ld_MinError;
		lr_Hit.md_NextMassAccuracy = ld_NextBestError;
		// this is only a good hit if the next best hit is far away
		if (lr_Hit.md_NextMassAccuracy > ad_MassAccuracy)
			lr_Hit.mb_IsGood = true;
	}
	return lr_Hit;
}


void k_ChitoScanner::calculateMeanAndStandardDeviation(QList<double> ak_Values, double* ad_Mean_, double* ad_StandardDeviation_)
{
	double ld_Mean = 0.0;
	foreach (double ld_Value, ak_Values)
		ld_Mean += ld_Value;
		
	ld_Mean /= ak_Values.size();
	
	double ld_StandardDeviation = 0.0;
	foreach (double ld_Value, ak_Values)
		ld_StandardDeviation += pow(ld_Value - ld_Mean, 2.0);
		
	ld_StandardDeviation /= ak_Values.size();
	ld_StandardDeviation = sqrt(ld_StandardDeviation);
	
	*ad_Mean_ = ld_Mean;
	*ad_StandardDeviation_ = ld_StandardDeviation;
}
