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
                                r_AcetylationType::Enumeration ae_AcetylationType,
								int ai_MinDP, int ai_MaxDP,
								int ai_MinCharge, int ai_MaxCharge, 
								int ai_IsotopeCount, 
								double ad_PrecursorMassAccuracy,
								QSet<r_LabelType::Enumeration> ak_Ms1VariableLabel,
								QSet<r_LabelType::Enumeration> ak_Ms1FixedLabel,
								double ad_ProductMassAccuracy,
								QSet<r_LabelType::Enumeration> ak_Ms2VariableLabel,
								QSet<r_LabelType::Enumeration> ak_Ms2FixedLabel,
								QTextStream* ak_CompositionFingerprintStream_,
								QTextStream* ak_MassRangeFingerprintStream_,
								QTextStream* ak_MassCollisionFingerprintStream_,
								QTextStream* ak_Ms2ProductsStream_)
	: k_ScanIterator(ae_ScanType, ak_MsLevels)
	, md_MinSnr(ad_MinSnr)
	, md_CropLower(ad_CropLower)
	, mk_IonizationType(ak_IonizationType)
	, me_AcetylationType(ae_AcetylationType)
	, mi_MinDP(ai_MinDP)
	, mi_MaxDP(ai_MaxDP)
	, mi_MinCharge(ai_MinCharge)
	, mi_MaxCharge(ai_MaxCharge)
	, mi_IsotopeCount(ai_IsotopeCount)
	, md_PrecursorMassAccuracy(ad_PrecursorMassAccuracy)
	, mk_Ms1VariableLabel(ak_Ms1VariableLabel)
	, mk_Ms1FixedLabel(ak_Ms1FixedLabel)
	, md_ProductMassAccuracy(ad_ProductMassAccuracy)
	, mk_Ms2VariableLabel(ak_Ms2VariableLabel)
	, mk_Ms2FixedLabel(ak_Ms2FixedLabel)
	, mk_CompositionFingerprintStream_(ak_CompositionFingerprintStream_)
	, mk_MassRangeFingerprintStream_(ak_MassRangeFingerprintStream_)
	, mk_MassCollisionFingerprintStream_(ak_MassCollisionFingerprintStream_)
	, mk_Ms2ProductsStream_(ak_Ms2ProductsStream_)
	, ms_CurrentSpot(QString())
	, mi_ExtractedMs1PeakCount(0)
	, mi_MatchedMs1PeakCount(0)
	, mi_PrecursorCount(0)
	, mi_MatchedPrecursorCount(0)
	, mi_ExtractedMs2PeakCount(0)
	, mi_MatchedMs2PeakCount(0)
	
{
	// subtract fixed labels from var labels,
	// because they're not really variable, huh?
	mk_Ms1VariableLabel -= mk_Ms1FixedLabel;
	mk_Ms2VariableLabel -= mk_Ms2FixedLabel;
	
	gk_IonizationTypeMass[r_IonizationType::Proton] = MASS_HYDROGEN;
	gk_IonizationTypeMass[r_IonizationType::Sodium] = MASS_SODIUM;
	gk_IonizationTypeLabel[r_IonizationType::Proton] = "H+";
	gk_IonizationTypeLabel[r_IonizationType::Sodium] = "Na+";
	gk_LabelMass[r_LabelType::ReducingEndLabel2Da] = MASS_HYDROGEN * 2.0;
	gk_LabelMass[r_LabelType::WaterLoss] = -MASS_WATER;
	gk_LabelLabel[r_LabelType::ReducingEndLabel2Da] = "RE";
	gk_LabelLabel[r_LabelType::WaterLoss] = "water loss";
}


k_ChitoScanner::~k_ChitoScanner()
{
	if (mk_CompositionFingerprintStream_)
		mk_CompositionFingerprintStream_->flush();
	
	if (mk_MassRangeFingerprintStream_)
		mk_MassRangeFingerprintStream_->flush();

	if (mk_MassCollisionFingerprintStream_)
		mk_MassCollisionFingerprintStream_->flush();
	
	if (mk_Ms2ProductsStream_)
		mk_Ms2ProductsStream_->flush();
}


void k_ChitoScanner::scan(QStringList ak_SpectraFiles)
{
	if (mk_MassCollisionFingerprintStream_ || mk_MassRangeFingerprintStream_)
	{
        // this is for the fingerprints, don't delete it
		r_OligoHit lr_Hit = 
			matchOligomer(200.0, 
							md_PrecursorMassAccuracy, 
							mk_IonizationType, me_AcetylationType, 
							mi_MinDP, mi_MaxDP,
							mi_MinCharge, mi_MaxCharge, mi_IsotopeCount,
							mk_Ms1VariableLabel, mk_Ms1FixedLabel, 1);
		if (mk_MassRangeFingerprintStream_)
		{
			mk_MassRangeFingerprintStream_->setRealNumberNotation(QTextStream::FixedNotation);
			*mk_MassRangeFingerprintStream_ << QString("Amount,A,%1\n").arg(gk_AcetylationChar[(int)me_AcetylationType]);
			foreach (tk_IntPair lk_Pair, mk_MS1MassInRangeFingerprint.keys())
			{
				int li_A = lk_Pair.first;
				int li_D = lk_Pair.second;
				double ld_Amount = mk_MS1MassInRangeFingerprint[lk_Pair];
				*mk_MassRangeFingerprintStream_<< ld_Amount << "," << li_A << "," << li_D << "\n";
			}
		}
		
		if (mk_MassCollisionFingerprintStream_)
		{
			mk_MassCollisionFingerprintStream_->setRealNumberNotation(QTextStream::FixedNotation);
			*mk_MassCollisionFingerprintStream_ << QString("Amount,A,%1\n").arg(gk_AcetylationChar[(int)me_AcetylationType]);
			foreach (tk_IntPair lk_Pair, mk_MS1MassCollisionFingerprint.keys())
			{
				int li_A = lk_Pair.first;
				int li_D = lk_Pair.second;
				double ld_Amount = mk_MS1MassCollisionFingerprint[lk_Pair];
				*mk_MassCollisionFingerprintStream_ << ld_Amount << "," << li_A << "," << li_D << "\n";
			}
		}
		
	}
	
	mk_MS1Fingerprint.clear();
	mk_Ms2Products.clear();
	
	// parse all bands
	foreach (QString ls_Path, ak_SpectraFiles)
	{
		ms_CurrentSpot = QFileInfo(ls_Path).baseName();
		
		mk_MS1MassAccuracies.clear();
		mi_MS2PeakHitCount = 0;
		mi_MS2PeakHitCountWithIsotope = 0;
		mi_ExtractedMs1PeakCount = 0;
		mi_MatchedMs1PeakCount = 0;
		mi_PrecursorCount = 0;
		mi_MatchedPrecursorCount = 0;
		mi_ExtractedMs2PeakCount = 0;
		mi_MatchedMs2PeakCount = 0;
	
		// parse spot
		this->parseFile(ls_Path);
		
		printf(" done.\n");
		
		qSort(mk_MS1MassAccuracies.begin(), mk_MS1MassAccuracies.end());

		if (mk_CompositionFingerprintStream_)
		{
			mk_CompositionFingerprintStream_->setRealNumberNotation(QTextStream::FixedNotation);
			*mk_CompositionFingerprintStream_ << QString("Amount,A,%1\n").arg(gk_AcetylationChar[(int)me_AcetylationType]);
			foreach (tk_IntPair lk_Pair, mk_MS1Fingerprint.keys())
			{
				int li_A = lk_Pair.first;
				int li_D = lk_Pair.second;
				double ld_Amount = mk_MS1Fingerprint[lk_Pair];
				*mk_CompositionFingerprintStream_ << ld_Amount << "," << li_A << "," << li_D << "\n";
			}
		}
		
		if (mk_Ms2ProductsStream_)
		{
			mk_Ms2ProductsStream_->setRealNumberNotation(QTextStream::FixedNotation);
			*mk_Ms2ProductsStream_ << QString("Amount,Precursor A,Precursor %1,Product A,Product %1\n").arg(gk_AcetylationChar[(int)me_AcetylationType]);
			foreach (tk_IntPair lk_PrecursorPair, mk_Ms2Products.keys())
			{
				foreach (tk_IntPair lk_Pair, mk_Ms2Products[lk_PrecursorPair].keys())
					*mk_Ms2ProductsStream_
						<< mk_Ms2Products[lk_PrecursorPair][lk_Pair] << ","
						<< lk_PrecursorPair.first << "," 
						<< lk_PrecursorPair.second << ","
						<< lk_Pair.first << ","
						<< lk_Pair.second << endl;
			}
		}
		
		printf("From %d MS1 peaks, %d matched to a product (%1.2f%%)\n",
				mi_ExtractedMs1PeakCount, mi_MatchedMs1PeakCount,
				(double)mi_MatchedMs1PeakCount / mi_ExtractedMs1PeakCount * 100.0);
		printf("From %d MS2 scan precursors, %d matched to a product (%1.2f%%)\n",
				mi_PrecursorCount, mi_MatchedPrecursorCount,
				(double)mi_MatchedPrecursorCount / mi_PrecursorCount * 100.0);
		printf("From %d MS2 peaks with matched precursors, %d matched to a product (%1.2f%%)\n",
				mi_ExtractedMs2PeakCount, mi_MatchedMs2PeakCount,
				(double)mi_MatchedMs2PeakCount / mi_ExtractedMs2PeakCount * 100.0);
				
/*		printf("MS1 mass accuracies (%d entries):\n", mk_MS1MassAccuracies.size());
			for (int i = 8; i <= 10; ++i)
			printf("%3i%% %9.4f / ", i * 10, mk_MS1MassAccuracies[(mk_MS1MassAccuracies.size() - 1) * i / 10]);
		printf("\n\n");
		
		printf("MS2 peak hit with isotope: %1.2f%%.\n\n", (double)mi_MS2PeakHitCountWithIsotope / mi_MS2PeakHitCount * 100.0);*/
	}
}


void k_ChitoScanner::handleScan(r_Scan& ar_Scan, bool& ab_Continue)
{
    ab_Continue = true;
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
							   mk_IonizationType, me_AcetylationType,
                               mi_MinDP, mi_MaxDP,
							   mi_MinCharge, mi_MaxCharge, mi_IsotopeCount,
							   mk_Ms1VariableLabel, mk_Ms1FixedLabel);
			if (lr_Hit.mb_IsGood)
			{
				mi_MatchedMs1PeakCount += 1;
				mk_MS1MassAccuracies << lr_Hit.md_MassAccuracy;
				tk_IntPair lk_Pair(lr_Hit.mi_A, lr_Hit.mi_D);
				
				if (!mk_MS1Fingerprint.contains(lk_Pair))
					mk_MS1Fingerprint[lk_Pair] = 0.0;
				mk_MS1Fingerprint[lk_Pair] += lr_Peak.md_PeakIntensity;
			}
		}
	}
	
 	if (ar_Scan.mi_MsLevel == 2)
 	{
		if (ar_Scan.mk_Precursors.size() != 1)
		{
			printf("Oops! Scan #%s has %d precursors, where we would expect exactly one --> Go figure out that one!\n", ar_Scan.ms_Id.toStdString().c_str(), ar_Scan.mk_Precursors.size());
			exit(1);
		}
		++mi_PrecursorCount;
 		r_Precursor lr_Precursor = ar_Scan.mk_Precursors.first();
 		r_OligoHit lr_PrecursorHit = 
			matchOligomer(lr_Precursor.md_Mz,
						   md_PrecursorMassAccuracy, 
						   mk_IonizationType, me_AcetylationType,
                           mi_MinDP, mi_MaxDP,
						   mi_MinCharge, mi_MaxCharge, mi_IsotopeCount,
						   mk_Ms1VariableLabel, mk_Ms1FixedLabel);
		if (lr_PrecursorHit.mb_IsGood)
		{
			++mi_MatchedPrecursorCount;
			foreach (r_Peak lr_Peak, lk_AllPeaks)
			{
				++mi_ExtractedMs2PeakCount;
				r_OligoHit lr_Hit = 
					matchOligomer(lr_Peak.md_PeakMz, 
								md_ProductMassAccuracy, 
								QSet<r_IonizationType::Enumeration>() 
									<< lr_PrecursorHit.me_Ionization, 
                                me_AcetylationType,
								1, (lr_PrecursorHit.mi_A + lr_PrecursorHit.mi_D),
								mi_MinCharge, lr_PrecursorHit.mi_Charge, 
								mi_IsotopeCount, mk_Ms2VariableLabel, mk_Ms2FixedLabel);
				if (lr_Hit.mb_IsGood)
				{
					++mi_MatchedMs2PeakCount;
					// oy, we found a MS2 peak hit
					if (lr_Hit.mi_A <= lr_PrecursorHit.mi_A &&
						lr_Hit.mi_D <= lr_PrecursorHit.mi_D)
					{
						tk_IntPair lk_PrecursorKey = tk_IntPair(lr_PrecursorHit.mi_A, lr_PrecursorHit.mi_D);
						tk_IntPair lk_ProductKey = tk_IntPair(lr_Hit.mi_A, lr_Hit.mi_D);
						if (!mk_Ms2Products.contains(lk_PrecursorKey))
							mk_Ms2Products[lk_PrecursorKey] = QHash<tk_IntPair, double>();
						if (!mk_Ms2Products[lk_PrecursorKey].contains(lk_ProductKey))
							mk_Ms2Products[lk_PrecursorKey][lk_ProductKey] = 0.0;
						mk_Ms2Products[lk_PrecursorKey][lk_ProductKey] += lr_Peak.md_PeakIntensity;
					}
				}
			}
		}
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
 	}
}


void k_ChitoScanner::progressFunction(QString as_ScanId, bool)
{
 	printf("\r%s: scan #%s...", ms_CurrentSpot.toStdString().c_str(), as_ScanId.toStdString().c_str());
}

								   
double k_ChitoScanner::calculateOligomerMass(r_IonizationType::Enumeration ae_Ionization, 
                                             r_AcetylationType::Enumeration ae_AcetylationType,
                                             int ai_DP, int ai_DA, 
                                             int ai_Charge, int ai_Isotope, 
                                             QSet<r_LabelType::Enumeration> ak_ActiveLabel)
{
	int li_A = ai_DA;
	int li_D = ai_DP - ai_DA;
	double ld_Mass = MASS_A * li_A + (ae_AcetylationType == r_AcetylationType::Deacetylation ? MASS_D : MASS_R) * li_D + MASS_WATER + 
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
                              r_AcetylationType::Enumeration ae_AcetylationType,
                              int ai_MinDP, int ai_MaxDP,
							  int ai_MinCharge, int ai_MaxCharge,
							  int ai_IsotopeCount,
							  QSet<r_LabelType::Enumeration> ak_VariableLabel,
							  QSet<r_LabelType::Enumeration> ak_FixedLabel,
							  int ai_AdditionalInfoMsLevel)
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
		QHash<tk_IntPair, QList<double> > lk_MzForComposition;
		if (ai_AdditionalInfoMsLevel == 1)
		{
			mk_MS1MassInRangeFingerprint.clear();
			mk_MS1MassCollisionFingerprint.clear();
		}
		for (int li_DP = ai_MinDP; li_DP <= ai_MaxDP; ++li_DP)
		{
			for (int li_DA = 0; li_DA <= li_DP; ++li_DA)
			{
				int li_A = li_DA;
				int li_D = li_DP - li_DA;
				int li_PossibilityCount = 0;
				int li_MassInRangeCount = 0;
				foreach (r_IonizationType::Enumeration le_Ionization, ak_IonizationType)
				{
					for (int li_Charge = ai_MinCharge; li_Charge <= ai_MaxCharge; ++li_Charge)
					{
						for (int li_Isotope = 0; li_Isotope < ai_IsotopeCount; ++li_Isotope)
						{
							QList<r_LabelType::Enumeration> lk_VariableLabelList = ak_VariableLabel.toList();
							int li_LabelCombinationCount = 1 << (ak_VariableLabel.size());
							for (int li_LabelCombination = 0; li_LabelCombination < li_LabelCombinationCount; ++li_LabelCombination)
							{
								QSet<r_LabelType::Enumeration> lk_EffectiveLabel;
								for (int i = 0; i < lk_VariableLabelList.size(); ++i)
								{
									if (((li_LabelCombination >> i) & 1) == 1)
										lk_EffectiveLabel << lk_VariableLabelList[i];
								}
								lk_EffectiveLabel |= ak_FixedLabel;
								double ld_Mz = calculateOligomerMass(
									le_Ionization, ae_AcetylationType, li_DP, li_DA, 
									li_Charge, li_Isotope, lk_EffectiveLabel);
								mk_TargetCache[ls_HashKey].insert(ld_Mz,
									r_OligoHit(le_Ionization, ae_AcetylationType,
                                               li_A, li_D, 
                                               li_Charge, li_Isotope, 
                                               lk_EffectiveLabel, ld_Mz));
								tk_IntPair lk_Composition(li_A, li_D);
								if (!lk_MzForComposition.contains(lk_Composition))
									lk_MzForComposition[lk_Composition] = QList<double>();
								lk_MzForComposition[lk_Composition] << ld_Mz;
								++li_PossibilityCount;
								if (ld_Mz >= 200.0 && ld_Mz <= 3000.0)
								{
									++li_MassInRangeCount;
								}
							}
						}
					}
				}
				double ld_Chance = (double)li_MassInRangeCount / li_PossibilityCount;
				mk_MS1MassInRangeFingerprint[tk_IntPair(li_A, li_D)] = ld_Chance;
			}
		}
/*		printf("%s\n", ls_HashKey.toStdString().c_str());
		QMultiMap<double, r_OligoHit>::const_iterator lk_Iter = mk_TargetCache[ls_HashKey].constBegin();
		for (; lk_Iter != mk_TargetCache[ls_HashKey].constEnd(); ++lk_Iter)
			printf("%10.4f: %s\n", lk_Iter.key(), lk_Iter.value().infoString().toStdString().c_str());*/
		// determine mass collision fingerprint
		if (ai_AdditionalInfoMsLevel == 1)
		{
			foreach (tk_IntPair lk_Composition, lk_MzForComposition.keys())
			{
				double ld_UsableLengthSum = 0.0;
				foreach (double ld_Mz, lk_MzForComposition[lk_Composition])
				{
					QMultiMap<double, r_OligoHit>::const_iterator lk_Iter = mk_TargetCache[ls_HashKey].constFind(ld_Mz);
					QList<double> lk_ProbeMz;
					if (lk_Iter != mk_TargetCache[ls_HashKey].constBegin())
					{
						QMultiMap<double, r_OligoHit>::const_iterator lk_PrevIter = lk_Iter;
						--lk_PrevIter;
						lk_ProbeMz << lk_PrevIter.key();
					}
					if (lk_Iter != mk_TargetCache[ls_HashKey].constEnd())
					{
						QMultiMap<double, r_OligoHit>::const_iterator lk_NextIter = lk_Iter;
						++lk_NextIter;
						lk_ProbeMz << lk_NextIter.key();
					}
					QPair<double, double> lk_TargetSpan = massBounds(ld_Mz, md_PrecursorMassAccuracy);
					double ld_OriginalLength = lk_TargetSpan.second - lk_TargetSpan.first;
					foreach (double ld_ProbeMz, lk_ProbeMz)
					{
						QPair<double, double> lk_ProbeMz = massBounds(ld_ProbeMz, md_PrecursorMassAccuracy);
						if (lk_ProbeMz.first < lk_TargetSpan.second && lk_ProbeMz.second > lk_TargetSpan.second)
							lk_TargetSpan.second = lk_ProbeMz.first;
						if (lk_ProbeMz.second > lk_TargetSpan.first && lk_ProbeMz.first < lk_TargetSpan.first)
							lk_TargetSpan.first = lk_ProbeMz.second;
					}
					double ld_NewLength = lk_TargetSpan.second - lk_TargetSpan.first;
					if (ld_NewLength < 0.0)
						ld_NewLength = 0.0;
					double ld_UsableLength = ld_NewLength / ld_OriginalLength;
					ld_UsableLengthSum += ld_UsableLength;
				}
				mk_MS1MassCollisionFingerprint[lk_Composition] = ld_UsableLengthSum / lk_MzForComposition[lk_Composition].size();
			}
		}
	}
	
	r_OligoHit lr_Hit;
	QMultiMap<double, r_OligoHit> lk_Targets = mk_TargetCache[ls_HashKey];
	double ld_MinError = 1e20;
	QMultiMap<double, r_OligoHit>::const_iterator lk_BestIter;
	QMultiMap<double, r_OligoHit>::const_iterator lk_Iter = lk_Targets.constBegin();
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


QPair<double, double> k_ChitoScanner::massBounds(double ad_Mass, double ad_MassAccuracy)
{
	double ld_Error = ad_Mass * ad_MassAccuracy / 1000000.0;
	return QPair<double, double>(ad_Mass - ld_Error, ad_Mass + ld_Error);
}
