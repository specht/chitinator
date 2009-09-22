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


struct r_OligoHit
{
	r_OligoHit(int ai_A = 0, int ai_D = 0, int ai_Isotope = 0,
				int ai_Charge = 0, int ai_IonizationType = 0,
				int ai_Label = 0, double ad_Mz = 0.0,
				double ad_Loss = 0.0,
				double ad_MassAccuracy = 0.0, 
				double ad_NextMassAccuracy = 0.0,
				bool ab_IsGood = false)
		: mi_A(ai_A)
		, mi_D(ai_D)
		, mi_Isotope(ai_Isotope)
		, mi_Charge(ai_Charge)
		, mi_IonizationType(ai_IonizationType)
		, mi_Label(ai_Label)
		, md_Mz(ad_Mz)
		, md_Loss(ad_Loss)
		, md_MassAccuracy(ad_MassAccuracy)
		, md_NextMassAccuracy(ad_NextMassAccuracy)
		, mb_IsGood(ab_IsGood)
	{
	}
	
	r_OligoHit(const r_OligoHit& ar_Other)
		: mi_A(ar_Other.mi_A)
		, mi_D(ar_Other.mi_D)
		, mi_Isotope(ar_Other.mi_Isotope)
		, mi_Charge(ar_Other.mi_Charge)
		, mi_IonizationType(ar_Other.mi_IonizationType)
		, mi_Label(ar_Other.mi_Label)
		, md_Mz(ar_Other.md_Mz)
		, md_Loss(ar_Other.md_Loss)
		, md_MassAccuracy(ar_Other.md_MassAccuracy)
		, md_NextMassAccuracy(ar_Other.md_NextMassAccuracy)
		, mb_IsGood(ar_Other.mb_IsGood)
	{
	}
	
	QString toString()
	{
		if (mb_IsGood)
			return QString("A%1D%2+%3 (%4%5+%6), %7 (%8 ppm, next: %9 ppm)")
				.arg(mi_A).arg(mi_D).arg(mi_Isotope).arg(mi_Charge)
				.arg(mi_IonizationType == 0 ? "H" : "Na")
				.arg(md_Loss == 0.0 ? "" : QString(", -%1").arg(md_Loss))
				.arg(mi_Label == 0 ? "unlabeled" : "labeled")
				.arg(md_MassAccuracy).arg(md_NextMassAccuracy);
		else
			return QString("[no hit]");
	}
	
	int mi_A;
	int mi_D;
	int mi_Isotope;
	int mi_Charge;
	int mi_IonizationType;
	int mi_Label;
	double md_Mz;
	double md_Loss;
	double md_MassAccuracy;
	double md_NextMassAccuracy;
	bool mb_IsGood;
};


class k_ChitoScanner: public k_ScanIterator
{
public:
	k_ChitoScanner(r_ScanType::Enumeration ae_ScanType,
					QList<tk_IntPair> ak_MsLevels,
					int ai_MaxIsotopeCount, 
					int ai_MinCharge, int ai_MaxCharge, 
					double ad_MinSnr, 
					double ad_PrecursorMassAccuracy,
					double ad_ProductMassAccuracy);
	virtual ~k_ChitoScanner();
	
	// quantify takes a list of spectra files and a hash of (peptide => protein) entries
	virtual void scan(QStringList ak_SpectraFiles);
	virtual void handleScan(r_Scan& ar_Scan);
	virtual void progressFunction(QString as_ScanId, bool ab_InterestingScan);
	virtual double calculateOligomerMass(int ai_IonizationType, 
										  int ai_DP, int ai_DA, 
										  int ai_Charge, int ai_Isotope, 
										  int ai_Label = 1, double ad_Loss = 0.0);
	virtual r_OligoHit matchOligomer(double ad_Mz,
									  double ad_MassAccuracy,
									  tk_IntPair ak_IonizationType = gk_Undefined,
									  tk_IntPair ak_DP = gk_Undefined,
									  tk_IntPair ak_Charge = gk_Undefined,
									  tk_IntPair ak_Isotope = gk_Undefined,
									  tk_IntPair ak_Label = gk_Undefined,
									  bool ab_CheckWaterLoss = false
								   );
	virtual void calculateMeanAndStandardDeviation(
		QList<double> ak_Values, double* ad_Mean_, double* ad_StandardDeviation_);
	
protected:
	int mi_MaxIsotopeCount;
	int mi_MinCharge;
	int mi_MaxCharge;
	double md_MinSnr;
	double md_PrecursorMassAccuracy;
	double md_ProductMassAccuracy;
	QString ms_CurrentSpot;
	QList<double> mk_MS1MassAccuracies;
	// this hash records the abundance of every (#A,#D) combination.
	QHash<tk_IntPair, double> mk_DAAmounts;
	// this hash records the abundance of every (#A,#D) combination,
	// only by occurences, not by peak height
	QHash<tk_IntPair, double> mk_DAAmountsUnweighted;
	// this list contains all MS1 product DPs, it's double so we can
	// do mean and SD afterwards with calculateMeanAndStandardDeviation
	QList<double> mk_DPs; 
	
	int mi_MS2PeakHitCount;
	int mi_MS2PeakHitCountWithIsotope;
	
	QHash<QString, QMap<double, r_OligoHit> > mk_TargetCache;
};
