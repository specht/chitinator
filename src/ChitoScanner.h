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


#define MASS_HYDROGEN 1.0078250321
#define MASS_SODIUM 22.98976967
#define MASS_WATER 18.01057
#define MASS_NEUTRON 1.002
#define MASS_D 161.0689
#define MASS_A 203.0794


typedef QPair<int, int> tk_IntPair;
const tk_IntPair gk_Undefined(-1, -1);


struct r_IonizationType
{
	enum Enumeration
	{
		Proton = 0,
		Sodium
	};
};


struct r_LabelType
{
	enum Enumeration
	{
		ReducingEndLabel2Da = 0,
		WaterLoss
	};
};


extern QHash<r_IonizationType::Enumeration, double> gk_IonizationTypeMass;
extern QHash<r_IonizationType::Enumeration, QString> gk_IonizationTypeLabel;
extern QHash<r_LabelType::Enumeration, double> gk_LabelMass;
extern QHash<r_LabelType::Enumeration, QString> gk_LabelLabel;


struct r_OligoHit
{
	r_OligoHit(r_IonizationType::Enumeration ae_Ionization = r_IonizationType::Proton,
				int ai_A = 0, int ai_D = 0, 
				int ai_Charge = 0, int ai_Isotope = 0,
				QSet<r_LabelType::Enumeration> ak_Label = QSet<r_LabelType::Enumeration>(),
				double ad_Mz = 0.0,
				double ad_MassAccuracy = 0.0, 
				double ad_NextMassAccuracy = 0.0,
				bool ab_IsGood = false)
		: me_Ionization(ae_Ionization)
		, mi_A(ai_A)
		, mi_D(ai_D)
		, mi_Charge(ai_Charge)
		, mi_Isotope(ai_Isotope)
		, mk_Label(ak_Label)
		, md_Mz(ad_Mz)
		, md_MassAccuracy(ad_MassAccuracy)
		, md_NextMassAccuracy(ad_NextMassAccuracy)
		, mb_IsGood(ab_IsGood)
	{
	}
	
	r_OligoHit(const r_OligoHit& ar_Other)
		: me_Ionization(ar_Other.me_Ionization)
		, mi_A(ar_Other.mi_A)
		, mi_D(ar_Other.mi_D)
		, mi_Charge(ar_Other.mi_Charge)
		, mi_Isotope(ar_Other.mi_Isotope)
		, mk_Label(ar_Other.mk_Label)
		, md_Mz(ar_Other.md_Mz)
		, md_MassAccuracy(ar_Other.md_MassAccuracy)
		, md_NextMassAccuracy(ar_Other.md_NextMassAccuracy)
		, mb_IsGood(ar_Other.mb_IsGood)
	{
	}
	
	QString toString() const
	{
		QList<r_LabelType::Enumeration> lk_LabelList = mk_Label.toList();
		qSort(lk_LabelList.begin(), lk_LabelList.end());
		QStringList lk_LabelString;
		foreach (r_LabelType::Enumeration le_Label, lk_LabelList)
			lk_LabelString << gk_LabelLabel[le_Label];
		if (lk_LabelString.empty())
			lk_LabelString << "(unlabeled)";
		
		if (mb_IsGood)
			return QString("A%1D%2+%3 (%4%5), %7 (%8 ppm, next: %9 ppm)")
				.arg(mi_A).arg(mi_D).arg(mi_Isotope).arg(mi_Charge)
				.arg(gk_IonizationTypeLabel[me_Ionization])
				.arg(lk_LabelString.join(", "))
				.arg(md_MassAccuracy).arg(md_NextMassAccuracy);
		else
			return QString("[no hit]");
	}
	
	QString compositionString() const
	{
		if (mb_IsGood)
			return QString("A%1D%2%3")
				.arg(mi_A).arg(mi_D).arg(mk_Label.contains(r_LabelType::ReducingEndLabel2Da) ? "*" : "");
		else
			return QString("[no hit]");
	}
	
	QString infoString() const
	{
		QList<r_LabelType::Enumeration> lk_LabelList = mk_Label.toList();
		qSort(lk_LabelList.begin(), lk_LabelList.end());
		QStringList lk_LabelString;
		foreach (r_LabelType::Enumeration le_Label, lk_LabelList)
			lk_LabelString << gk_LabelLabel[le_Label];
		if (lk_LabelString.empty())
			lk_LabelString << "(unlabeled)";
		
		return QString("A%1D%2+%3 (%4%5), %7")
			.arg(mi_A).arg(mi_D).arg(mi_Isotope).arg(mi_Charge)
			.arg(gk_IonizationTypeLabel[me_Ionization])
			.arg(lk_LabelString.join(", "));
	}
	
	r_IonizationType::Enumeration me_Ionization;
	int mi_A;
	int mi_D;
	int mi_Charge;
	int mi_Isotope;
	QSet<r_LabelType::Enumeration> mk_Label;
	double md_Mz;
	double md_MassAccuracy;
	double md_NextMassAccuracy;
	bool mb_IsGood;
};


class k_ChitoScanner: public k_ScanIterator
{
public:
	k_ChitoScanner(r_ScanType::Enumeration ae_ScanType,
					QList<tk_IntPair> ak_MsLevels,
					double ad_MinSnr,
					double ad_CropLower,
					QSet<r_IonizationType::Enumeration> ak_IonizationType,
					int ai_MinDP, int ai_MaxDP,
					int ai_MinCharge, int ai_MaxCharge, int ai_IsotopeCount,
					double ad_PrecursorMassAccuracy,
					QSet<r_LabelType::Enumeration> ak_Ms1VariableLabel,
					QSet<r_LabelType::Enumeration> ak_Ms1FixedLabel,
					double ad_ProductMassAccuracy,
					QSet<r_LabelType::Enumeration> ak_Ms2VariableLabel,
					QSet<r_LabelType::Enumeration> ak_Ms2FixedLabel,
					QTextStream* ak_CompositionFingerprintStream_ = NULL,
					QTextStream* ak_MassRangeFingerprintStream_ = NULL,
					QTextStream* ak_MassCollisionFingerprintStream_ = NULL,
					QTextStream* ak_Ms2ProductsStream_ = NULL);
	virtual ~k_ChitoScanner();
	
	// quantify takes a list of spectra files and a hash of (peptide => protein) entries
	virtual void scan(QStringList ak_SpectraFiles);
	virtual void handleScan(r_Scan& ar_Scan);
	virtual void progressFunction(QString as_ScanId, bool ab_InterestingScan);
	virtual double calculateOligomerMass(r_IonizationType::Enumeration ae_Ionization, 
										  int ai_DP, int ai_DA, 
										  int ai_Charge, int ai_Isotope, 
										  QSet<r_LabelType::Enumeration> ak_ActiveLabel);
	virtual r_OligoHit matchOligomer(double ad_Mz,
									  double ad_MassAccuracy,
									  QSet<r_IonizationType::Enumeration> ak_IonizationType,
									  int ai_MinDP, int ai_MaxDP,
									  int ai_MinCharge, int ai_MaxCharge,
									  int ai_IsotopeCount,
									  QSet<r_LabelType::Enumeration> ak_VariableLabel,
									  QSet<r_LabelType::Enumeration> ak_FixedLabel,
									  int ai_AdditionalInfoMsLevel = 0
								   );
	virtual void calculateMeanAndStandardDeviation(
		QList<double> ak_Values, double* ad_Mean_, double* ad_StandardDeviation_);
	virtual QPair<double, double> massBounds(double ad_Mass, double ad_MassAccuracy);
	
protected:
	double md_MinSnr;
	double md_CropLower;
	QSet<r_IonizationType::Enumeration> mk_IonizationType;
	int mi_MinDP;
	int mi_MaxDP;
	int mi_MinCharge;
	int mi_MaxCharge;
	int mi_IsotopeCount;
	double md_PrecursorMassAccuracy;
	QSet<r_LabelType::Enumeration> mk_Ms1VariableLabel;
	QSet<r_LabelType::Enumeration> mk_Ms1FixedLabel;
	double md_ProductMassAccuracy;
	QSet<r_LabelType::Enumeration> mk_Ms2VariableLabel;
	QSet<r_LabelType::Enumeration> mk_Ms2FixedLabel;
	QTextStream* mk_CompositionFingerprintStream_;
	QTextStream* mk_MassRangeFingerprintStream_;
	QTextStream* mk_MassCollisionFingerprintStream_;
	QTextStream* mk_Ms2ProductsStream_;
	
	QString ms_CurrentSpot;

	// how many peaks were there, how many matched to a product?
	int mi_ExtractedMs1PeakCount;
	int mi_MatchedMs1PeakCount;
	
	int mi_PrecursorCount;
	int mi_MatchedPrecursorCount;
	
	int mi_ExtractedMs2PeakCount;
	int mi_MatchedMs2PeakCount;
	
	QList<double> mk_MS1MassAccuracies;
	// this hash records the abundance of every (#A,#D) combination.
	QHash<tk_IntPair, double> mk_MS1Fingerprint;
	QHash<tk_IntPair, double> mk_MS1MassInRangeFingerprint;
	QHash<tk_IntPair, double> mk_MS1MassCollisionFingerprint;
	QHash<tk_IntPair, QHash<tk_IntPair, double> > mk_Ms2Products;
	
	int mi_MS2PeakHitCount;
	int mi_MS2PeakHitCountWithIsotope;
	
	QHash<QString, QMultiMap<double, r_OligoHit> > mk_TargetCache;
};