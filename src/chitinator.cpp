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

// 
int main(int ai_ArgumentCount, char** ac_Arguments__)
{
	if (ai_ArgumentCount < 5)
	{
		printf("Usage: [runs] [dp] [da] [enzyme (ADx|)]\n");
		exit(1);
	}
	int RUNS = QVariant(ac_Arguments__[1]).toInt();
	QList<int> lk_Dp;
	foreach (QString s, QString(ac_Arguments__[2]).split(","))
		lk_Dp.push_back(QVariant(s).toInt());
	QList<int> lk_Da;
	foreach (QString s, QString(ac_Arguments__[3]).split(","))
		lk_Da.push_back(QVariant(s).toInt());
	char* EXP_ENZYME = ac_Arguments__[4];
	
	srand(time(NULL));
	printf("This is the Chitinator.\n");

	/*
	unsigned int* lui_Count_ = new unsigned int[1001];
	for (int i = 0; i < 1001; ++i)
		lui_Count_[i] = 0;
	for (int i = 0; i < 100000; ++i)
	{
		r_PolymerParameters lr_PolymerParameters((double)1000.0, 10.0, (double)20.0 / 100.0);
		k_Polymer lk_Polymer(lr_PolymerParameters);
		int li_ACount = 0;
		for (int k = 0; k < lk_Polymer.getLength(); ++k)
		{
			if (lk_Polymer.getMonomer(k) == D)
			{
				if (li_ACount > 0)
					++lui_Count_[li_ACount];
				li_ACount = 0;
			}
			else
				++li_ACount;
		}
	}
	for (int i = 0; i < 1001; ++i)
		printf("%d;%d\n", i, lui_Count_[i]);
	delete [] lui_Count_;
	exit(1);

	*/
	
	foreach (int EXP_DP, lk_Dp)
	{
		foreach (int EXP_DA, lk_Da)
		{
			int li_DP = EXP_DP + 1;
			int li_IndexCount = ((li_DP - 1) * li_DP) / 2 + li_DP - 1;
			printf("We have %d bins.\n", li_IndexCount);
			unsigned int* lui_Histogram_ = new unsigned int[li_IndexCount];
			unsigned int* lui_OldHistogram_ = new unsigned int[li_IndexCount];
			memset(lui_Histogram_, 0, li_IndexCount * sizeof(unsigned int));
			memset(lui_OldHistogram_, 0, li_IndexCount * sizeof(unsigned int));
			unsigned int lui_HistogramMax = 0;
			unsigned int lui_OldHistogramMax = 0;
			
			r_PolymerParameters lr_PolymerParameters((double)EXP_DP, 10.0, (double)EXP_DA / 100.0);
			
			QString ls_Pattern(EXP_ENZYME);
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
			
			for (int li_Run = 1; li_Run <= RUNS; ++li_Run)
			{
				printf("\rDP %d, enzyme %s, DA %d: %d...", EXP_DP, EXP_ENZYME, EXP_DA, li_Run);
				QList<RefPtr<k_Polymer> > lk_Polymers;
				RefPtr<k_Polymer> lk_pPolymer(new k_Polymer(lr_PolymerParameters));
				lk_Polymers.push_back(lk_pPolymer);
				
				// simple, full digestion
				QList<RefPtr<k_Polymer> > lk_Products;
				QList<RefPtr<k_Polymer> > lk_DigestablePolymers = lk_Polymers;
				while (!lk_DigestablePolymers.empty())
				{
					/*
					printf("polymers:\n");
					for (int i = 0; i < lk_DigestablePolymers.size(); ++i)
						printf("%s\n", lk_DigestablePolymers[i]->toString().toStdString().c_str());
					printf("\n");
					*/
					QList<RefPtr<k_Polymer> >::iterator lk_Iter = lk_DigestablePolymers.begin();
					bool lb_FoundSite = false;
					while (lk_Iter != lk_DigestablePolymers.end())
					{
						RefPtr<k_Polymer> lk_pPolymer = *lk_Iter;
						for (int k = 0; k < lk_pPolymer->mi_Length; ++k)
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
					int li_DP = lk_Products[i]->getLength();
					int li_DA = lk_Products[i]->getLength() - lk_Products[i]->getDCount();
					int li_Index = ((li_DP - 1) * li_DP) / 2 + li_DP + li_DA - 1;
					++lui_Histogram_[li_Index];
					lui_HistogramMax = std::max<unsigned int>(lui_HistogramMax, lui_Histogram_[li_Index]);
				}
				/*
		<?xml version='1.0' encoding='UTF-8' standalone='no'?>
		<svg xmlns:svg='http://www.w3.org/2000/svg' xmlns='http://www.w3.org/2000/svg' version='1.0' width='1000' height='1000'>
		<rect x='0' y='0' width='1024' height='1024' fill='white' />
		<circle cx='0.5' cy='0.5' r='1.0' fill='red' />
		</svg>
				*/
				double ld_Difference = 0.0;
				for (int i = 0; i < li_IndexCount; ++i)
					ld_Difference += pow(((double)lui_Histogram_[i] / lui_HistogramMax) - ((double)lui_OldHistogram_[i] / lui_OldHistogramMax), 2.0);
				//printf("%le\n", ld_Difference);
				memcpy(lui_OldHistogram_, lui_Histogram_, li_IndexCount * sizeof(unsigned int));
				lui_OldHistogramMax = lui_HistogramMax;

				if (li_Run == RUNS)
				{
					QFile lk_File(QString("out/out-%1-%2-%3-%4.csv").arg(EXP_DP).arg(QString(EXP_ENZYME).replace("|", "_")).arg(EXP_DA).arg(li_Run));
					lk_File.open(QIODevice::WriteOnly);
					QTextStream lk_Stream(&lk_File);
					lk_Stream << "ab rel,rel mass,DP,DA abs,DA rel\n";
					QMultiMap<int, int> lk_ProductsByAbundance;
					for (int i = 0; i < li_IndexCount; ++i)
						lk_ProductsByAbundance.insert(lui_Histogram_[i], i);
					QList<int> lk_Indices = lk_ProductsByAbundance.values();
					QHash<int, int> lk_IndexToDA;
					QHash<int, int> lk_IndexToDP;
					for (int dp = 1; dp <= EXP_DP; ++dp)
						for (int da = 0; da <= dp; ++da)
						{
							int li_Index = ((dp - 1) * dp) / 2 + dp + da - 1;
							lk_IndexToDP[li_Index] = dp;
							lk_IndexToDA[li_Index] = da;
						}

printf("%d\n", lk_Indices.size());
					for (int i = lk_Indices.size() - 1; i >= 0; --i)
					{
						int li_Index = lk_Indices[i];
						if (lui_Histogram_[li_Index] > 0)
						{
							double ld_Value = (double)lui_Histogram_[li_Index] / lui_HistogramMax;
							double ld_RelativeMass = (double)lui_Histogram_[li_Index] * lk_IndexToDP[li_Index] / ((double)EXP_DP * RUNS);
							lk_Stream << ld_Value << "," << ld_RelativeMass << "," << lk_IndexToDP[li_Index] << "," << lk_IndexToDA[li_Index] << "," << (double)lk_IndexToDA[li_Index] / lk_IndexToDP[li_Index] << "\n";
						}
					}
					
					lk_File.close();
					
					lk_File.setFileName(QString("out/out-%1-%2-%3-%4-mass.svg").arg(EXP_DP).arg(QString(EXP_ENZYME).replace("|", "_")).arg(EXP_DA).arg(li_Run));
					lk_File.open(QIODevice::WriteOnly);
					lk_Stream.setDevice(&lk_File);
					lk_Stream << "<?xml version='1.0' encoding='UTF-8' standalone='no'?> <svg xmlns:svg='http://www.w3.org/2000/svg' xmlns='http://www.w3.org/2000/svg' version='1.0' width='1080' height='600'> <rect x='0' y='0' width='1080' height='600' fill='white' />" << endl;
					lk_Stream << "<marker id='arrow' viewBox='0 0 20 10' refX='0' refY='5' markerUnits='strokeWidth' markerWidth='8' markerHeight='6' orient='auto' fill='#000'><path d='M 0 0 L 20 5 L 0 10 z' /></marker>" << endl;
					lk_Stream << "<g transform='translate(40,40)'>" << endl;
					//lk_Stream << "<g transform='scale(0.5,0.5)'>" << endl;
					for (int li_DP = 1; li_DP <= EXP_DP; ++li_DP)
					{
						double ld_DP0 = (double)li_DP / EXP_DP;
						double ld_DP1 = (double)(li_DP + 1.1) / EXP_DP;
						ld_DP0 = pow(ld_DP0, 0.1) * 2.0 - 1.0;
						ld_DP1 = pow(ld_DP1, 0.1) * 2.0 - 1.0;
						// insert log scale here!
						double ld_R0 = ld_DP0 * 500.0;
						if (li_DP == 1)
							ld_R0 = 0.0;
						double ld_R1 = ld_DP1 * 500.0;
						lk_Stream << QString("<mask id='m%1' maskUnits='userSpaceOnUse' x='0' y='0' width='1000' height='500'><circle cx='500' cy='0' r='%2' fill='white'  /><circle cx='500' cy='0' r='%3' fill='black'  /></mask>").arg(li_DP).arg(ld_R1).arg(ld_R0) << endl;
					}
					for (int li_DP = EXP_DP; li_DP >= 1; --li_DP)
						for (int li_DA = 0; li_DA <= li_DP; ++li_DA)
						{
							int li_Index = ((li_DP - 1) * li_DP) / 2 + li_DP + li_DA - 1;
							if (lui_Histogram_[li_Index] > 0)
							{
								double ld_Radius = (double)lui_Histogram_[li_Index] / lui_HistogramMax;
								ld_Radius = pow(ld_Radius, 0.1);
								double da = (double)li_DA;
								double dp = (double)li_DP;
								double ada0 = da / (dp + 1.0);
								double ada1 = (da + 1.0) / (dp + 1.0);
								ada0 *= M_PI;
								ada1 *= M_PI;
								ada1 += 0.01;
								double x0 = 500.0 - cos(ada0) * 1000.0;
								double y0 = sin(ada0) * 1000.0;
								double x1 = 500.0 - cos(ada1) * 1000.0;
								double y1 = sin(ada1) * 1000.0;
								char lc_Color_[16];
								ld_Radius = 1.0 - ld_Radius;
								
								double ld_RelativeMass = (double)lui_Histogram_[li_Index] * lk_IndexToDP[li_Index] / ((double)EXP_DP * RUNS);
								ld_Radius = pow(ld_Radius, 0.3);
								ld_Radius = 1.0 - ld_Radius;
								
								sprintf(lc_Color_, "%02x%02x%02x", (int)(ld_Radius * 255.0),(int)(ld_Radius * 255.0), (int)(ld_Radius * 255.0));
								lk_Stream << QString("<path d='M500,0 L%1,%2 %3,%4 z' fill='#%5' mask='url(#m%6)' />").arg(x0).arg(y0).arg(x1).arg(y1).arg(QString(lc_Color_)).arg(li_DP) << endl;
								//lk_Stream << QString("<rect x='%1' y='%2' width='5' height='5' fill='black' opacity='%3' />").arg(x).arg(y).arg(ld_Radius) << endl;
							}
						}
					lk_Stream << "<line x1='500' y1='0' x2='500' y2='470' fill='none' stroke='#000' stroke-width='2.5' marker-end='url(#arrow)' />" << endl;
					lk_Stream << "<path transform='translate(500, 0) scale(1, -1) translate(-500, 0)' d='M 0,0 a 500,500 -180 0,1 1000,0' fill='none' stroke='#000' stroke-width='2.5' marker-end='url(#arrow)' />" << endl;
					lk_Stream << "<text x='510' y='476' fill='#000' style='font-family: Bitstream Charter' font-size='20px'>DP</text>" << endl;
					lk_Stream << "<text x='965' y='10' fill='#000' style='font-family: Bitstream Charter' font-size='20px'>DA</text>" << endl;
					for (int i = 0; i <= 100; i += 10)
					{
						lk_Stream << QString("<text fill='#000' transform='translate(500, 0) rotate(%1) translate(-500, 0)' text-anchor='middle' x='500' y='522' style='font-family: Bitstream Charter' font-size='20px'>%2%</text>").arg((double)-i / 100.0 * 180.0 + 90.0).arg(i) << endl;
						if (i != 100)
							lk_Stream << QString("<line fill='none' stroke='#000' stroke-width='2.5' transform='translate(500, 0) rotate(%1) translate(-500, 0)' x1='500' y1='495' x2='500' y2='505' />").arg((double)-i / 100.0 * 180.0 + 90.0) << endl;
					}
					lk_Stream << QString("<text x='0' y='548' fill='#000' style='font-family: Bitstream Charter' font-size='20px'>enzyme: %1, DP: 1000, DA: %2, 1000 iterations, error: %3</text>").arg(QString(EXP_ENZYME)).arg(EXP_DA).arg(ld_Difference) << endl;
					lk_Stream << "</g>" << endl;
					lk_Stream << "</svg>" << endl;
					lk_File.close();

					lk_File.setFileName(QString("out/out-%1-%2-%3-%4-counts.svg").arg(EXP_DP).arg(QString(EXP_ENZYME).replace("|", "_")).arg(EXP_DA).arg(li_Run));
					lk_File.open(QIODevice::WriteOnly);
					lk_Stream.setDevice(&lk_File);
					lk_Stream << "<?xml version='1.0' encoding='UTF-8' standalone='no'?> <svg xmlns:svg='http://www.w3.org/2000/svg' xmlns='http://www.w3.org/2000/svg' version='1.0' width='1080' height='600'> <rect x='0' y='0' width='1080' height='600' fill='white' />" << endl;
					lk_Stream << "<marker id='arrow' viewBox='0 0 20 10' refX='0' refY='5' markerUnits='strokeWidth' markerWidth='8' markerHeight='6' orient='auto' fill='#000'><path d='M 0 0 L 20 5 L 0 10 z' /></marker>" << endl;
					lk_Stream << "<g transform='translate(40,40)'>" << endl;
					//lk_Stream << "<g transform='scale(0.5,0.5)'>" << endl;
					for (int li_DP = 1; li_DP <= EXP_DP; ++li_DP)
					{
						double ld_DP0 = (double)li_DP / EXP_DP;
						double ld_DP1 = (double)(li_DP + 1.1) / EXP_DP;
						ld_DP0 = pow(ld_DP0, 0.1) * 2.0 - 1.0;
						ld_DP1 = pow(ld_DP1, 0.1) * 2.0 - 1.0;
						// insert log scale here!
						double ld_R0 = ld_DP0 * 500.0;
						if (li_DP == 1)
							ld_R0 = 0.0;
						double ld_R1 = ld_DP1 * 500.0;
						lk_Stream << QString("<mask id='m%1' maskUnits='userSpaceOnUse' x='0' y='0' width='1000' height='500'><circle cx='500' cy='0' r='%2' fill='white'  /><circle cx='500' cy='0' r='%3' fill='black'  /></mask>").arg(li_DP).arg(ld_R1).arg(ld_R0) << endl;
					}
					for (int li_DP = EXP_DP; li_DP >= 1; --li_DP)
						for (int li_DA = 0; li_DA <= li_DP; ++li_DA)
						{
							int li_Index = ((li_DP - 1) * li_DP) / 2 + li_DP + li_DA - 1;
							if (lui_Histogram_[li_Index] > 0)
							{
								double ld_Radius = (double)lui_Histogram_[li_Index] / lui_HistogramMax;
								ld_Radius = pow(ld_Radius, 0.1);
								double da = (double)li_DA;
								double dp = (double)li_DP;
								double ada0 = da / (dp + 1.0);
								double ada1 = (da + 1.0) / (dp + 1.0);
								ada0 *= M_PI;
								ada1 *= M_PI;
								ada1 += 0.01;
								double x0 = 500.0 - cos(ada0) * 1000.0;
								double y0 = sin(ada0) * 1000.0;
								double x1 = 500.0 - cos(ada1) * 1000.0;
								double y1 = sin(ada1) * 1000.0;
								char lc_Color_[16];
								ld_Radius = 1.0 - ld_Radius;
								
								sprintf(lc_Color_, "%02x%02x%02x", (int)(ld_Radius * 255.0),(int)(ld_Radius * 255.0), (int)(ld_Radius * 255.0));
								lk_Stream << QString("<path d='M500,0 L%1,%2 %3,%4 z' fill='#%5' mask='url(#m%6)' />").arg(x0).arg(y0).arg(x1).arg(y1).arg(QString(lc_Color_)).arg(li_DP) << endl;
								//lk_Stream << QString("<rect x='%1' y='%2' width='5' height='5' fill='black' opacity='%3' />").arg(x).arg(y).arg(ld_Radius) << endl;
							}
						}
					lk_Stream << "<line x1='500' y1='0' x2='500' y2='470' fill='none' stroke='#000' stroke-width='2.5' marker-end='url(#arrow)' />" << endl;
					lk_Stream << "<path transform='translate(500, 0) scale(1, -1) translate(-500, 0)' d='M 0,0 a 500,500 -180 0,1 1000,0' fill='none' stroke='#000' stroke-width='2.5' marker-end='url(#arrow)' />" << endl;
					lk_Stream << "<text x='510' y='476' fill='#000' style='font-family: Bitstream Charter' font-size='20px'>DP</text>" << endl;
					lk_Stream << "<text x='965' y='10' fill='#000' style='font-family: Bitstream Charter' font-size='20px'>DA</text>" << endl;
					for (int i = 0; i <= 100; i += 10)
					{
						lk_Stream << QString("<text fill='#000' transform='translate(500, 0) rotate(%1) translate(-500, 0)' text-anchor='middle' x='500' y='522' style='font-family: Bitstream Charter' font-size='20px'>%2%</text>").arg((double)-i / 100.0 * 180.0 + 90.0).arg(i) << endl;
						if (i != 100)
							lk_Stream << QString("<line fill='none' stroke='#000' stroke-width='2.5' transform='translate(500, 0) rotate(%1) translate(-500, 0)' x1='500' y1='495' x2='500' y2='505' />").arg((double)-i / 100.0 * 180.0 + 90.0) << endl;
					}
					lk_Stream << QString("<text x='0' y='548' fill='#000' style='font-family: Bitstream Charter' font-size='20px'>enzyme: %1, DP: 1000, DA: %2, 1000 iterations, error: %3</text>").arg(QString(EXP_ENZYME)).arg(EXP_DA).arg(ld_Difference) << endl;
					lk_Stream << "</g>" << endl;
					lk_Stream << "</svg>" << endl;
					lk_File.close();
				}
				
			}
			printf("\n");
			delete [] lui_Histogram_;
			delete [] lui_OldHistogram_;
		}
	}
}
