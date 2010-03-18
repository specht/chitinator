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
	printf("  --minSnr [float] (default: 2.0)\n");
	printf("  --cropLower [float] (default: 0.05)\n");
	printf("  --ionizationType [comma-separated ids] (default: 0)\n");
	printf("      (0: H+, 1: Na+)\n");
    printf("  --acetylationType <D|R> (default: D)\n");
	printf("  --minDP [int] (default: 1)\n");
	printf("  --maxDP [int] (default: 10)\n");
	printf("  --minCharge [int] (default: 1)\n");
	printf("  --maxCharge [int] (default: 3)\n");
	printf("  --isotopeCount [int] (default: 2)\n");
	printf("  --precursorMassAccuracy (ppm) [float] (default: 5.0)\n");
	printf("  --ms1VariableLabel [comma-separated ids] (default: empty)\n");
	printf("      (0: MS2 label [+2 Da], 1: water loss)\n");
	printf("  --ms1FixedLabel [comma-separated ids] (default: empty)\n");
	printf("      (0: MS2 label [+2 Da], 1: water loss)\n");
	printf("  --productMassAccuracy (ppm) [float] (default: 700.0)\n");
	printf("  --ms2VariableLabel [comma-separated ids] (default: empty)\n");
	printf("      (0: MS2 label [+2 Da], 1: water loss)\n");
	printf("  --ms2FixedLabel [comma-separated ids] (default: empty)\n");
	printf("      (0: MS2 label [+2 Da], 1: water loss)\n");
	printf("  --writeCompositionFingerprint [path]\n");
	printf("      write CSV composition fingerprint to [path]\n");
	printf("  --writeMassRangeFingerprint [path]\n");
	printf("      write mass range fingerprint to [path]\n");
	printf("  --writeMassCollisionFingerprint [path]\n");
	printf("      write mass collision fingerprint to [path]\n");
	printf("  --writeMs2Products [path]\n");
	printf("      write CSV MS2 products to [path]\n");
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


QSet<r_LabelType::Enumeration> parseLabel(QString as_Label)
{
	QSet<r_LabelType::Enumeration> lk_Label;
	foreach (QString ls_Item, as_Label.split(","))
	{
		bool lb_Ok = false;
		int li_Id = ls_Item.trimmed().toInt(&lb_Ok);
		if (!lb_Ok)
		{
			printf("Error: Unknown label %s.\n", ls_Item.trimmed().toStdString().c_str());
			exit(1);
		}
		switch(li_Id)
		{
			case 0:
				lk_Label << r_LabelType::ReducingEndLabel2Da;
				break;
			case 1:
				lk_Label << r_LabelType::WaterLoss;
				break;
			default:
				printf("Error: Unknown label %s.\n", ls_Item.trimmed().toStdString().c_str());
				exit(1);
		}
	}
	return lk_Label;
}


int main(int ai_ArgumentCount, char** ac_Arguments__)
{
	QStringList lk_Arguments;
	for (int i = 1; i < ai_ArgumentCount; ++i)
		lk_Arguments << ac_Arguments__[i];
	
	if (lk_Arguments.empty())
		printUsageAndExit();

	QSet<r_IonizationType::Enumeration> lk_IonizationType = 
		QSet<r_IonizationType::Enumeration>() << r_IonizationType::Proton;
    r_AcetylationType::Enumeration le_AcetylationType = r_AcetylationType::Deacetylation;
	double ld_MinSnr = 2.0;
	double ld_CropLower = 0.05;
	int li_MinDP = 1;
	int li_MaxDP = 10;
	int li_MinCharge = 1;
	int li_MaxCharge = 3;
	int li_IsotopeCount = 2;
	double ld_PrecursorMassAccuracy = 5.0;
	QSet<r_LabelType::Enumeration> lk_Ms1VariableLabel = 
		QSet<r_LabelType::Enumeration>();
	QSet<r_LabelType::Enumeration> lk_Ms1FixedLabel = 
		QSet<r_LabelType::Enumeration>();
	double ld_ProductMassAccuracy = 700.0;
	QSet<r_LabelType::Enumeration> lk_Ms2VariableLabel = 
		QSet<r_LabelType::Enumeration>();
	QSet<r_LabelType::Enumeration> lk_Ms2FixedLabel = 
		QSet<r_LabelType::Enumeration>();
	RefPtr<QIODevice> lk_pCompositionFingerprintDevice;
	RefPtr<QTextStream> lk_pCompositionFingerprintStream;
	RefPtr<QIODevice> lk_pMassRangeFingerprintDevice;
	RefPtr<QTextStream> lk_pMassRangeFingerprintStream;
	RefPtr<QIODevice> lk_pMassCollisionFingerprintDevice;
	RefPtr<QTextStream> lk_pMassCollisionFingerprintStream;
	RefPtr<QIODevice> lk_pMs2ProductsDevice;
	RefPtr<QTextStream> lk_pMs2ProductsStream;
	
	// consume options
	int li_Index = 0;
	
	li_Index = lk_Arguments.indexOf("--minSnr");
	if (li_Index > -1)
	{
		ld_MinSnr = QVariant(lk_Arguments[li_Index + 1]).toDouble();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--cropLower");
	if (li_Index > -1)
	{
		ld_CropLower = QVariant(lk_Arguments[li_Index + 1]).toDouble();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--ionizationType");
	if (li_Index > -1)
	{
		QString ls_IonizationType = QVariant(lk_Arguments[li_Index + 1]).toString();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
		lk_IonizationType = QSet<r_IonizationType::Enumeration>();
		foreach (QString ls_Item, ls_IonizationType.split(","))
		{
			bool lb_Ok = false;
			int li_Id = ls_Item.trimmed().toInt(&lb_Ok);
			if (!lb_Ok)
			{
				printf("Error: Unknown ionization type %s.\n", ls_Item.trimmed().toStdString().c_str());
				exit(1);
			}
			switch(li_Id)
			{
				case 0:
					lk_IonizationType << r_IonizationType::Proton;
					break;
				case 1:
					lk_IonizationType << r_IonizationType::Sodium;
					break;
				default:
					printf("Error: Unknown ionization type %s.\n", ls_Item.trimmed().toStdString().c_str());
					exit(1);
			}
		}
	}
	
    li_Index = lk_Arguments.indexOf("--acetylationType");
    if (li_Index > -1)
    {
        QString ls_AcetylationType = QVariant(lk_Arguments[li_Index + 1]).toString();
        lk_Arguments.removeAt(li_Index);
        lk_Arguments.removeAt(li_Index);
        if (ls_AcetylationType == "D")
            le_AcetylationType = r_AcetylationType::Deacetylation;
        else if (ls_AcetylationType == "R")
            le_AcetylationType = r_AcetylationType::Reacetylation;
        else
        {
            printf("Error: Unknown acetylation type %s.\n", ls_AcetylationType.trimmed().toStdString().c_str());
            exit(1);
        }
    }
    
	li_Index = lk_Arguments.indexOf("--minDP");
	if (li_Index > -1)
	{
		li_MinDP = QVariant(lk_Arguments[li_Index + 1]).toInt();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--maxDP");
	if (li_Index > -1)
	{
		li_MaxDP = QVariant(lk_Arguments[li_Index + 1]).toInt();
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
	
	li_Index = lk_Arguments.indexOf("--isotopeCount");
	if (li_Index > -1)
	{
		li_IsotopeCount = QVariant(lk_Arguments[li_Index + 1]).toInt();
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
	
	li_Index = lk_Arguments.indexOf("--ms1VariableLabel");
	if (li_Index > -1)
	{
		QString ls_Label = QVariant(lk_Arguments[li_Index + 1]).toString();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
		lk_Ms1VariableLabel = parseLabel(ls_Label);
	}
	
	li_Index = lk_Arguments.indexOf("--ms1FixedLabel");
	if (li_Index > -1)
	{
		QString ls_Label = QVariant(lk_Arguments[li_Index + 1]).toString();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
		lk_Ms1FixedLabel = parseLabel(ls_Label);
	}
	
	li_Index = lk_Arguments.indexOf("--productMassAccuracy");
	if (li_Index > -1)
	{
		ld_ProductMassAccuracy = QVariant(lk_Arguments[li_Index + 1]).toDouble();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--ms2VariableLabel");
	if (li_Index > -1)
	{
		QString ls_Label = QVariant(lk_Arguments[li_Index + 1]).toString();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
		lk_Ms2VariableLabel = parseLabel(ls_Label);
	}
	
	li_Index = lk_Arguments.indexOf("--ms2FixedLabel");
	if (li_Index > -1)
	{
		QString ls_Label = QVariant(lk_Arguments[li_Index + 1]).toString();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
		lk_Ms2FixedLabel = parseLabel(ls_Label);
	}
	
	li_Index = lk_Arguments.indexOf("--writeCompositionFingerprint");
	if (li_Index > -1)
	{
		QString ls_Path = QVariant(lk_Arguments[li_Index + 1]).toString();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
		lk_pCompositionFingerprintDevice = RefPtr<QIODevice>(new QFile(ls_Path));
		lk_pCompositionFingerprintDevice->open(QIODevice::WriteOnly);
		lk_pCompositionFingerprintStream = RefPtr<QTextStream>(new QTextStream(lk_pCompositionFingerprintDevice.get_Pointer()));
	}
	
	li_Index = lk_Arguments.indexOf("--writeMassRangeFingerprint");
	if (li_Index > -1)
	{
		QString ls_Path = QVariant(lk_Arguments[li_Index + 1]).toString();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
		lk_pMassRangeFingerprintDevice = RefPtr<QIODevice>(new QFile(ls_Path));
		lk_pMassRangeFingerprintDevice->open(QIODevice::WriteOnly);
		lk_pMassRangeFingerprintStream = RefPtr<QTextStream>(new QTextStream(lk_pMassRangeFingerprintDevice.get_Pointer()));
	}
	
	li_Index = lk_Arguments.indexOf("--writeMassCollisionFingerprint");
	if (li_Index > -1)
	{
		QString ls_Path = QVariant(lk_Arguments[li_Index + 1]).toString();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
		lk_pMassCollisionFingerprintDevice = RefPtr<QIODevice>(new QFile(ls_Path));
		lk_pMassCollisionFingerprintDevice->open(QIODevice::WriteOnly);
		lk_pMassCollisionFingerprintStream = RefPtr<QTextStream>(new QTextStream(lk_pMassCollisionFingerprintDevice.get_Pointer()));
	}
	
	li_Index = lk_Arguments.indexOf("--writeMs2Products");
	if (li_Index > -1)
	{
		QString ls_Path = QVariant(lk_Arguments[li_Index + 1]).toString();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
		lk_pMs2ProductsDevice = RefPtr<QIODevice>(new QFile(ls_Path));
		lk_pMs2ProductsDevice->open(QIODevice::WriteOnly);
		lk_pMs2ProductsStream = RefPtr<QTextStream>(new QTextStream(lk_pMs2ProductsDevice.get_Pointer()));
	}
	
	QStringList lk_SpectraFiles;
	foreach (QString ls_Path, lk_Arguments)
		lk_SpectraFiles << ls_Path;
	
	k_ChitoScanner lk_ChitoScanner(r_ScanType::All, 
									QList<tk_IntPair>() << tk_IntPair(1, 10000),
									ld_MinSnr, ld_CropLower,
									lk_IonizationType, le_AcetylationType,
									li_MinDP, li_MaxDP,
									li_MinCharge, li_MaxCharge, li_IsotopeCount,
									ld_PrecursorMassAccuracy,
									lk_Ms1VariableLabel, lk_Ms1FixedLabel,
									ld_ProductMassAccuracy,
									lk_Ms2VariableLabel, lk_Ms2FixedLabel,
									lk_pCompositionFingerprintStream.get_Pointer(),
									lk_pMassRangeFingerprintStream.get_Pointer(),
									lk_pMassCollisionFingerprintStream.get_Pointer(),
									lk_pMs2ProductsStream.get_Pointer());
		
	lk_ChitoScanner.scan(lk_SpectraFiles);
}
