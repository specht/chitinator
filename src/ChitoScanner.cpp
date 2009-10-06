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


QHash<r_IonizationType::Enumeration, double> gk_IonizationTypeMass;
QHash<r_IonizationType::Enumeration, QString> gk_IonizationTypeLabel;
QHash<r_LabelType::Enumeration, double> gk_LabelMass;
QHash<r_LabelType::Enumeration, QString> gk_LabelLabel;


k_ChitoScanner::k_ChitoScanner(r_ScanType::Enumeration ae_ScanType,
								QList<tk_IntPair> ak_MsLevels,
								double ad_MinSnr,
								double ad_CropLower,
								QSet<r_IonizationType::Enumeration> ak_IonizationType,
								int ai_MinDP, int ai_MaxDP,
								int ai_MinCharge, int ai_MaxCharge, 
								int ai_IsotopeCount, 
								QSet<r_LabelType::Enumeration> ak_VariableLabel,
								QSet<r_LabelType::Enumeration> ak_FixedLabel,
								double ad_PrecursorMassAccuracy,
								double ad_ProductMassAccuracy,
								QTextStream* ak_CompositionFingerprintStream_)
	: k_ScanIterator(ae_ScanType, ak_MsLevels)
	, md_MinSnr(ad_MinSnr)
	, md_CropLower(ad_CropLower)
	, mk_IonizationType(ak_IonizationType)
	, mi_MinDP(ai_MinDP)
	, mi_MaxDP(ai_MaxDP)
	, mi_MinCharge(ai_MinCharge)
	, mi_MaxCharge(ai_MaxCharge)
	, mi_IsotopeCount(ai_IsotopeCount)
	, mk_VariableLabel(ak_VariableLabel)
	, mk_FixedLabel(ak_FixedLabel)
	, md_PrecursorMassAccuracy(ad_PrecursorMassAccuracy)
	, md_ProductMassAccuracy(ad_ProductMassAccuracy)
	, mk_CompositionFingerprintStream_(ak_CompositionFingerprintStream_)
{
	// subtract fixed labels from var labels,
	// because they're not really variable, huh?
	mk_VariableLabel -= mk_FixedLabel;
	
	gk_IonizationTypeMass[r_IonizationType::Proton] = MASS_HYDROGEN;
	gk_IonizationTypeMass[r_IonizationType::Sodium] = MASS_SODIUM;
	gk_IonizationTypeLabel[r_IonizationType::Proton] = "H+";
	gk_IonizationTypeLabel[r_IonizationType::Sodium] = "Na+";
	gk_LabelMass[r_LabelType::ReducingEndLabel2Da] = MASS_HYDROGEN * 2.0;
	gk_LabelMass[r_LabelType::WaterLoss] = -MASS_WATER;
	gk_LabelLabel[r_LabelType::ReducingEndLabel2Da] = "RE label (+2 Da)";
	gk_LabelLabel[r_LabelType::WaterLoss] = "water loss";
}


k_ChitoScanner::~k_ChitoScanner()
{
	if (mk_CompositionFingerprintStream_)
		mk_CompositionFingerprintStream_->flush();
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
		mk_DPs.clear();
		mi_ExtractedMs1PeakCount = 0;
		mi_MatchedMs1PeakCount = 0;
	
		// parse spot
		this->parseFile(ls_Path);
		
		printf(" done.\n");
		
		qSort(mk_DPs.begin(), mk_DPs.end());
		qSort(mk_MS1MassAccuracies.begin(), mk_MS1MassAccuracies.end());

		if (mk_CompositionFingerprintStream_)
		{
			*mk_CompositionFingerprintStream_ << "Amount,A,D\n";
			foreach (tk_IntPair lk_Pair, mk_DAAmounts.keys())
			{
				int li_A = lk_Pair.first;
				int li_D = lk_Pair.second;
				double ld_Amount = mk_DAAmounts[lk_Pair];
				*mk_CompositionFingerprintStream_ << ld_Amount << "," << li_A << "," << li_D << "\n";
			}
		}
		
		printf("From %d MS1 peaks, %d matched to a product (%1.2f%%)\n",
				mi_ExtractedMs1PeakCount, mi_MatchedMs1PeakCount,
				(double)mi_MatchedMs1PeakCount / mi_ExtractedMs1PeakCount * 100.0);
				
/*		printf("MS1 mass accuracies (%d entries):\n", mk_MS1MassAccuracies.size());
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
		
		printf("Product DP median : %d\n", (int)(mk_DPs[mk_DPs.size() / 2]));
		double ld_DPMean = 0.0;
		double ld_DPSD = 0.0;
		calculateMeanAndStandardDeviation(mk_DPs, &ld_DPMean, &ld_DPSD);
		printf("Product DP mean/SD: %1.2f / %1.2f\n\n", ld_DPMean, ld_DPSD);
		
		printf("MS2 peak hit with isotope: %1.2f%%.\n\n", (double)mi_MS2PeakHitCountWithIsotope / mi_MS2PeakHitCount * 100.0);*/
	}
}


void k_ChitoScanner::handleScan(r_Scan& ar_Scan)
{
	// find all peaks in this spectrum
	QList<r_Peak> lk_AllPeaks = k_ScanIterator::findAllPeaks(ar_Scan.mr_Spectrum);
	
	QList<r_Peak> lk_TempPeaks;
	
	// filter out lower (md_CropLower) %
	if (md_CropLower > 0.0)
	{
		double ld_MaxIntensity = 0.0;
		foreach (r_Peak lr_Peak, lk_AllPeaks)
			ld_MaxIntensity = std::max<double>(ld_MaxIntensity, lr_Peak.md_PeakIntensity);
		QList<r_Peak> lk_TempPeaks;
		lk_TempPeaks.clear();
		foreach (r_Peak lr_Peak, lk_AllPeaks)
			if (lr_Peak.md_PeakIntensity >= ld_MaxIntensity * md_CropLower)
				lk_TempPeaks << lr_Peak;
		lk_AllPeaks = lk_TempPeaks;
	}
	
	// filter by SNR
	if (md_MinSnr > 0.0)
	{
		lk_TempPeaks.clear();
		foreach (r_Peak lr_Peak, lk_AllPeaks)
			if (lr_Peak.md_Snr >= md_MinSnr)
				lk_TempPeaks << lr_Peak;
		lk_AllPeaks = lk_TempPeaks;
	}

	if (ar_Scan.mi_MsLevel == 1)
	{
		mi_ExtractedMs1PeakCount += lk_AllPeaks.size();
		foreach (r_Peak lr_Peak, lk_AllPeaks)
		{
			r_OligoHit lr_Hit = 
				matchOligomer(lr_Peak.md_PeakMz, 
							   md_PrecursorMassAccuracy, 
							   mk_IonizationType, mi_MinDP, mi_MaxDP,
							   mi_MinCharge, mi_MaxCharge, mi_IsotopeCount,
							   mk_VariableLabel, mk_FixedLabel);
			if (lr_Hit.mb_IsGood)
			{
				mi_MatchedMs1PeakCount += 1;
				mk_MS1MassAccuracies << lr_Hit.md_MassAccuracy;
				tk_IntPair lk_Pair(lr_Hit.mi_A, lr_Hit.mi_D);
				
				if (!mk_DAAmounts.contains(lk_Pair))
					mk_DAAmounts[lk_Pair] = 0.0;
				mk_DAAmounts[lk_Pair] += lr_Peak.md_PeakIntensity;
				
				mk_DPs << (lr_Hit.mi_A + lr_Hit.mi_D);
			}
//  			printf("%1.4f %s\n", lr_Peak.md_PeakMz, lr_Hit.toString().toStdString().c_str());
		}
	}
	
// 	if (ar_Scan.mi_MsLevel == 2)
// 	{
// 		r_Precursor lr_Precursor = ar_Scan.mk_Precursors.first();
// 		r_OligoHit lr_PrecursorHit = 
// 			matchOligomer(lr_Precursor.md_Mz, md_PrecursorMassAccuracy,
// 						   gk_Undefined, gk_Undefined,
// 						   tk_IntPair(lr_Precursor.mi_ChargeState,
// 									   lr_Precursor.mi_ChargeState));
// 		if (lr_PrecursorHit.mb_IsGood)
// 		{
// 			foreach (r_Peak lr_Peak, lk_AllPeaks)
// 			{
// 				r_OligoHit lr_Hit = matchOligomer(lr_Peak.md_PeakMz, 
// 					md_ProductMassAccuracy, tk_IntPair(0, 0), 
// 					gk_Undefined, tk_IntPair(1, 1), tk_IntPair(0, 1), 
// 					gk_Undefined, true);
// 				QString ls_Fragment;
// 				if (lr_Hit.mb_IsGood)
// 				{
// 					++mi_MS2PeakHitCount;
// 					if (lr_Hit.mi_Isotope > 0)
// 						++mi_MS2PeakHitCountWithIsotope;
// 					ls_Fragment = "[";
// 					for (int i = 0; i < lr_Hit.mi_A; ++i)
// 						ls_Fragment += "A";
// 					for (int i = 0; i < lr_Hit.mi_D; ++i)
// 						ls_Fragment += "D";
// 					ls_Fragment += "]";
// 					if (lr_Hit.mi_Label == 1)
// 						ls_Fragment += "*";
// 				}
// 			}
// 		}
// 	}
}


void k_ChitoScanner::progressFunction(QString as_ScanId, bool)
{
 	printf("\r%s: scan #%s...", ms_CurrentSpot.toStdString().c_str(), as_ScanId.toStdString().c_str());
}

								   
double k_ChitoScanner::calculateOligomerMass(r_IonizationType::Enumeration ae_Ionization, 
											  int ai_DP, int ai_DA, 
											  int ai_Charge, int ai_Isotope, 
											  QSet<r_LabelType::Enumeration> ak_ActiveLabel)
{
	int li_A = ai_DA;
	int li_D = ai_DP - ai_DA;
	double ld_Mass = MASS_A * li_A + MASS_D * li_D + MASS_WATER + 
		gk_IonizationTypeMass[ae_Ionization] * ai_Charge + 
		MASS_NEUTRON * ai_Isotope;
	foreach (r_LabelType::Enumeration le_Label, ak_ActiveLabel)
		ld_Mass += gk_LabelMass[le_Label];
	ld_Mass /= ai_Charge;
	return ld_Mass;
}


r_OligoHit 
k_ChitoScanner::matchOligomer(double ad_Mz,
							   double ad_MassAccuracy,
							   QSet<r_IonizationType::Enumeration> ak_IonizationType,
							   int ai_MinDP, int ai_MaxDP,
							   int ai_MinCharge, int ai_MaxCharge,
							   int ai_IsotopeCount,
							   QSet<r_LabelType::Enumeration> ak_VariableLabel,
							   QSet<r_LabelType::Enumeration> ak_FixedLabel)
{
	// if ai_Charge < 1, then charge is not defined
	QString ls_Result = "(no match)";
	
	QList<r_IonizationType::Enumeration> lk_IonizationTypeList = ak_IonizationType.toList();
	qSort(lk_IonizationTypeList.begin(), lk_IonizationTypeList.end());
	QStringList lk_IonizationTypeCsvList;
	foreach (r_IonizationType::Enumeration le_Ionization, lk_IonizationTypeList)
		lk_IonizationTypeCsvList << QVariant(le_Ionization).toString();
	
	QList<r_LabelType::Enumeration> lk_VariableLabelList = ak_VariableLabel.toList();
	qSort(lk_VariableLabelList.begin(), lk_VariableLabelList.end());
	QStringList lk_VariableLabelCsvList;
	foreach (r_LabelType::Enumeration le_Label, lk_VariableLabelList)
		lk_VariableLabelCsvList << QVariant(le_Label).toString();
	
	QList<r_LabelType::Enumeration> lk_FixedLabelList = ak_FixedLabel.toList();
	qSort(lk_FixedLabelList.begin(), lk_FixedLabelList.end());
	QStringList lk_FixedLabelCsvList;
	foreach (r_LabelType::Enumeration le_Label, lk_FixedLabelList)
		lk_FixedLabelCsvList << QVariant(le_Label).toString();
	
	QString ls_HashKey = QString("it%1/dp%2,%3/ch%4,%5/iso%6/varl%7/fixl%8")
		.arg(lk_IonizationTypeCsvList.join(","))
		.arg(ai_MinDP).arg(ai_MaxDP)
		.arg(ai_MinCharge).arg(ai_MaxCharge)
		.arg(ai_IsotopeCount)
		.arg(lk_VariableLabelCsvList.join(","))
		.arg(lk_FixedLabelCsvList.join(","));
		
	if (!mk_TargetCache.contains(ls_HashKey))
	{
		foreach (r_IonizationType::Enumeration le_Ionization, ak_IonizationType)
		{
			for (int li_DP = ai_MinDP; li_DP <= ai_MaxDP; ++li_DP)
			{
				for (int li_DA = 0; li_DA <= li_DP; ++li_DA)
				{
					int li_A = li_DA;
					int li_D = li_DP - li_DA;
					for (int li_Charge = ai_MinCharge; li_Charge <= ai_MaxCharge; ++li_Charge)
					{
						for (int li_Isotope = 0; li_Isotope < ai_IsotopeCount; ++li_Isotope)
						{
							QList<r_LabelType::Enumeration> lk_VariableLabelList = mk_VariableLabel.toList();
							int li_LabelCombinationCount = 1 << (mk_VariableLabel.size());
							for (int li_LabelCombination = 0; li_LabelCombination < li_LabelCombinationCount; ++li_LabelCombination)
							{
								QSet<r_LabelType::Enumeration> lk_EffectiveLabel;
								for (int i = 0; i < lk_VariableLabelList.size(); ++i)
								{
									if (((li_LabelCombination >> i) & 1) == 1)
										lk_EffectiveLabel << lk_VariableLabelList[i];
								}
								lk_EffectiveLabel |= mk_FixedLabel;
								double ld_Mz = calculateOligomerMass(
									le_Ionization, li_DP, li_DA, 
									li_Charge, li_Isotope, lk_EffectiveLabel);
								mk_TargetCache[ls_HashKey][ld_Mz] = 
									r_OligoHit(le_Ionization, li_A, li_D, 
												li_Charge, li_Isotope, 
												lk_EffectiveLabel, ld_Mz);
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
