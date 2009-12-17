/*
Copyright (c) 2007-2008 Michael Specht

This file is part of Chitinator.

Chitinator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SimQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Chitinator.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <QtCore>
#include "RefPtr.h"
#include "Digestion.h"
#include "Enzyme.h"
#include "Polymer.h"

void printUsageAndExit()
{
	printf("Usage: chitinator [options]\n");
	printf("Options:\n");
	printf("  --runs [int] (default: 1000)\n");
	printf("  --DP [int] (default: 800)\n");
	printf("  --DA [float] (default: 0.66)\n");
	printf("  --enzyme [string] (default: A|D)\n");
	printf("      A, D and x can be used to specify residues (x: don't care)\n");
	printf("  --maxIterations [int] (default: disabled)\n");
	printf("      Stop after n iterations.\n");
	printf("  --writeCompositionFingerprint [path]\n");
	printf("      write CSV composition fingerprint to [path]\n");
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

	int li_Runs = 1000;
	int li_DP = 800;
	float lf_DA = 0.66;
	QString ls_Enzyme = "A|D";
	int li_MaxIterations = -1;
	
	RefPtr<QIODevice> lk_pCompositionFingerprintDevice;
	RefPtr<QTextStream> lk_pCompositionFingerprintStream;
	
	// consume options
	int li_Index = 0;
	
	li_Index = lk_Arguments.indexOf("--runs");
	if (li_Index > -1)
	{
		li_Runs = QVariant(lk_Arguments[li_Index + 1]).toInt();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--DP");
	if (li_Index > -1)
	{
		li_DP = QVariant(lk_Arguments[li_Index + 1]).toInt();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--DA");
	if (li_Index > -1)
	{
		lf_DA = QVariant(lk_Arguments[li_Index + 1]).toDouble();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--enzyme");
	if (li_Index > -1)
	{
		ls_Enzyme = QVariant(lk_Arguments[li_Index + 1]).toString();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--maxIterations");
	if (li_Index > -1)
	{
		li_MaxIterations = QVariant(lk_Arguments[li_Index + 1]).toInt();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
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
	
	srand(time(NULL));

	r_PolymerParameters lr_PolymerParameters((double)li_DP, 10.0, lf_DA);
	
	QString ls_Pattern = ls_Enzyme;
	int li_Length = ls_Pattern.length() - 1;
	int li_CutPosition = ls_Pattern.indexOf("|");
	ls_Pattern.replace("|", "");
	unsigned int lui_Mask = 0;
	unsigned int lui_Pattern = 0;
	for (int i = li_Length - 1; i >= 0; --i)
	{
		char lc_Char = ls_Pattern.at(i).toAscii();
		lui_Mask <<= 1;
		lui_Pattern <<= 1;
		if (lc_Char == 'A')
		{
			lui_Mask |= 1;
		}
		else if (lc_Char == 'D')
		{
			lui_Mask |= 1;
			lui_Pattern |= 1;
		}
		else if (lc_Char == 'x')
		{
		}
		else
		{
			printf("Error!");
			exit(1);
		}
	}
	//printf("determined enzyme for %s: %d %d %d %d.\n", EXP_ENZYME, li_Length, lui_Mask, lui_Pattern, li_CutPosition);
	r_EnzymeParameters lr_EnzymeParameters(li_Length, lui_Mask, lui_Pattern, li_CutPosition);
	k_Enzyme lk_Enzyme(lr_EnzymeParameters);
	k_Digestion lk_Digestion;
	
	typedef QPair<int, int> tk_IntPair;
	QHash<tk_IntPair, int> mk_Fingerprint;
	for (int li_Run = 1; li_Run <= li_Runs; ++li_Run)
	{
		printf("\rDP %d, enzyme %s, DA %1.2f%%: %d...", li_DP, ls_Enzyme.toStdString().c_str(), lf_DA * 100.0, li_Run);
		QList<RefPtr<k_Polymer> > lk_Polymers;
		RefPtr<k_Polymer> lk_pPolymer(new k_Polymer(lr_PolymerParameters));
		lk_Polymers.push_back(lk_pPolymer);
		
		// simple, full digestion
		QList<RefPtr<k_Polymer> > lk_DigestablePolymers = lk_Polymers;
		QList<RefPtr<k_Polymer> > lk_Products;
		
		int li_IterationCount = 0;
		
		while (!lk_DigestablePolymers.empty())
		{
			if (li_MaxIterations >= 0)
			{
				if (li_IterationCount >= li_MaxIterations)
					break;
			}
			++li_IterationCount;
			QList<RefPtr<k_Polymer> >::iterator lk_Iter = lk_DigestablePolymers.begin();
			bool lb_FoundSite = false;
			while (lk_Iter != lk_DigestablePolymers.end())
			{
				RefPtr<k_Polymer> lk_pPolymer = *lk_Iter;
                QList<int> lk_Permutation;
                for (int k = 0; k < lk_pPolymer->mi_Length; ++k)
                    lk_Permutation << k;
                for (int k = 0; k < lk_pPolymer->mi_Length; ++k)
                {
                    int li_Position = rand() % lk_pPolymer->mi_Length;
                    int li_Temp = lk_Permutation[li_Position];
                    lk_Permutation[li_Position] = lk_Permutation[k];
                    lk_Permutation[k] = li_Temp;
                }
                foreach (int k, lk_Permutation)
				{
					if (lk_Enzyme.canCutPolymerAt(lk_pPolymer.get_Pointer(), k))
					{
						QList<RefPtr<k_Polymer> > lk_Parts = lk_pPolymer->cleaveAt(lk_Enzyme.mr_EnzymeParameters.mi_Cut + k);
						// chuck out the cleaved polymer out of lk_DigestablePolymers
						lk_DigestablePolymers.erase(lk_Iter);
						// add parts to lk_DigestablePolymers
						lk_DigestablePolymers += lk_Parts;
						lb_FoundSite = true;
						break;
					}
				}
				if (lb_FoundSite)
					break;
				else
				{
					// we did not find a cleavage site
					lk_Products.push_back(lk_pPolymer);
					lk_Iter = lk_DigestablePolymers.erase(lk_Iter);
				}
			}
		}
		for (int i = 0; i < lk_Products.size(); ++i)
		{
			int li_D = lk_Products[i]->getDCount();
			int li_A = lk_Products[i]->getLength() - li_D;
			tk_IntPair lk_AD = tk_IntPair(li_A, li_D);
			if (!mk_Fingerprint.contains(lk_AD))
				mk_Fingerprint[lk_AD] = 0;
			mk_Fingerprint[lk_AD] += 1;
		}
	}
	if (lk_pCompositionFingerprintStream)
	{
		*(lk_pCompositionFingerprintStream.get_Pointer()) << "Amount,A,D\n";
		foreach (tk_IntPair lk_AD, mk_Fingerprint.keys())
			*(lk_pCompositionFingerprintStream.get_Pointer()) << mk_Fingerprint[lk_AD] << "," << lk_AD.first << "," << lk_AD.second << "\n";
	}
	printf("\n");
}
