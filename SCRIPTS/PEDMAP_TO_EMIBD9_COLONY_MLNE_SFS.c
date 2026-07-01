
// PEDMAP_TO_EMIBD9_COLONY_MLNE_SFS.c (only ped and map files)

#include "libhdr"

#define II 250 // Maximum N=250 individuals
#define SS 300001 // Maximum 300000 SNPs segregating

int snp_mono, snp_m2, maxCEROS, s1, s2, type_allele;
int i, s, k, x, NIND, NSNP, snp[II][SS][2], all[SS][5];
int SUM_ALLELES[SS], numSNP, calc_SNP[SS], mono_SNP[SS], allele[SS];
int ind[II], sex[II], numHET[SS], FreqClass[II]; 

double MAF, w, freq[SS][3];
double pheL[II], pheW[II]; 

double numLOCI, EHom, EHet, Hom[II], sum_LR[II], sum_LRnum[II], sum_LRden[II], sum_VR1[II], sum_VR2[II], sum_YANG1[II], sum_YANG2[II], sum_LH2[II];
double NL[II], fLR[II], fvr1[II], fvr2[II], fyang1[II], fyang2[II], fLH1[II], fLH2[II], fhom[II];

char ch;

FILE *fped, *fout, *finb, *fgeno, *fEMIBD9, *fcolony_dat, *fmlne, *fcolony_mark, *fsfs, *ffped;

main()
{
	s1 = 1;
	s2 = 2;
	maxCEROS = 1000;

	fout = fopen ("outfileF","w");
	finb = fopen ("inbreeding_file","w");
	fgeno = fopen ("genotype_file","w");
	fEMIBD9  = fopen ("EMIBD9_file","w");
	fcolony_dat  = fopen ("COLONY_data.txt","w");
	fmlne  = fopen ("MLNE_data.txt","w");
	fcolony_mark = fopen ("COLONY_markers.txt","w");
	fsfs = fopen ("SFS_file","w");

	tracelevel = 1;
	getintandskip("NIND :",&NIND,2,250);
	getintandskip("NSNP :",&NSNP,2,infinity);
	getrealandskip("MAF :",&MAF,0.0,1.0);
	getintandskip("Type_allele :",&type_allele,0,1);

	readfiles();
	alleles();
	inbreeding();

	return(0);
}

/* **************************************************************************** */

readfiles()
{
	// ********** Read data.ped file to get genotypes and phenotypes ********** 

	for (i=1; i<=NIND; i++)
	{
		ind[i] = i;
		sex[i] = 1;
		pheL[i] = 100;
		pheW[i] = 100;
  	}

	fped = fopen ("data.ped","r");

	for (i=1; i<=NIND; i++)
	{
		lookfortext("-9");

		if (type_allele ==0)
		{
			for (s=1; s<=NSNP; s++)
			{
				fscanf(fped,"%c", &ch);
				fscanf(fped,"%c", &ch);
//				if (s == 1) if (i==1)	fprintf(fout, "%c", ch);
//				if (s == 1) if (i==1)	fprintf(fout, "%c", ch);

				if (ch == 'A')	snp[i][s][1] = 1;
				else if (ch == 'T')	snp[i][s][1] = 2;
				else if (ch == 'C')	snp[i][s][1] = 3;
				else if (ch == 'G')	snp[i][s][1] = 4;
				else if (ch == '0')	snp[i][s][1] = 0;

				fscanf(fped,"%c", &ch);
				fscanf(fped,"%c", &ch);
//				if (s == 1) if (i==1)	fprintf(fout, "%c", ch);
//				if (s == 1) if (i==1)	fprintf(fout, "%c", ch);

				if (ch == 'A')	snp[i][s][2] = 1;
				else if (ch == 'T')	snp[i][s][2] = 2;
				else if (ch == 'C')	snp[i][s][2] = 3;
				else if (ch == 'G')	snp[i][s][2] = 4;
				else if (ch == '0')	snp[i][s][2] = 0;

				if(i==2)  fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
			}
		}

		if (type_allele ==1)
		{
			for (s=1; s<=NSNP; s++)
			{
				fscanf(fped,"%c", &ch);
				fscanf(fped,"%c", &ch);
//				if (s == 50) if (i==725)	fprintf(fout, "%c", ch);
//				if (s == 50) if (i==728)	fprintf(fout, "%c", ch);

				if (ch == '1')	snp[i][s][1] = 1;
				else if (ch == '2')	snp[i][s][1] = 2;
				else if (ch == '0')	snp[i][s][1] = 0;

				fscanf(fped,"%c", &ch);
				fscanf(fped,"%c", &ch);
//				if (s == 50) if (i==725)	fprintf(fout, "%c", ch);
//				if (s == 50) if (i==728)	fprintf(fout, "%c", ch);

				if (ch == '1')	snp[i][s][2] = 1;
				else if (ch == '2')	snp[i][s][2] = 2;
				else if (ch == '0')	snp[i][s][2] = 0;
			}
		}
	}

	if (tracelevel != 0)
	{
		fprintf(fout, "\n");
		for (s=1; s<=NSNP; s++)
		{
			if (s == s1)
			{
				fprintf(fout, "s=%d  ", s);
				for (i=1; i<=NIND; i++)	 fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
				fprintf(fout, "\n");
			}
			if (s == s2)
			{
				fprintf(fout, "s=%d  ", s);
				for (i=1; i<=NIND; i++)	 fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
				fprintf(fout, "\n");
			}
		}
	}

	fclose(fped);

	ffped = fopen ("datanew.ped","w");
	for (i=1; i<=NIND; i++)
	{
		fprintf(ffped,"1 IND%d 0 0 1 -9 ", i);
		for (s=1; s<=NSNP; s++)
		{
			if (snp[i][s][1]==1)		fprintf(ffped, "A ");
			else if (snp[i][s][1]==2)		fprintf(ffped, "T ");
			else if (snp[i][s][1]==3)		fprintf(ffped, "C ");
			else if (snp[i][s][1]==4)		fprintf(ffped, "G ");
			else if (snp[i][s][1]==0)		fprintf(ffped, "0 ");

			if (snp[i][s][2]==1)		fprintf(ffped, "A ");
			else if (snp[i][s][2]==2)		fprintf(ffped, "T ");
			else if (snp[i][s][2]==3)		fprintf(ffped, "C ");
			else if (snp[i][s][2]==4)		fprintf(ffped, "G ");
			else if (snp[i][s][2]==0)		fprintf(ffped , "0 ");
		}
		fprintf(ffped,"\n");
	}
	fclose (ffped);

	return(0);
}

/* **************************************************************************** */

alleles()
{
	/******* Alleles segregating *******/

	for (s=1; s<=NSNP; s++)
	for (i=1; i<=NIND; i++)
	for (k=1; k<=2; k++)
	{
		if      (snp[i][s][k] == 1) 		all[s][1] ++;
		else if (snp[i][s][k] == 2) 	all[s][2] ++;
		else if (snp[i][s][k] == 3) 	all[s][3] ++;
		else if (snp[i][s][k] == 4) 	all[s][4] ++;
		else if (snp[i][s][k] == 0) 	all[s][5] ++;
	}
	for (s=1; s<=NSNP; s++)	SUM_ALLELES[s] = (all[s][1] + all[s][2] + all[s][3] + all[s][4] + all[s][5]);

	for (s=1; s<=NSNP; s++)
	{
		if (all[s][1] > 0) 	allele[s]++;
		if (all[s][2] > 0) 	allele[s]++;
		if (all[s][3] > 0) 	allele[s]++;
		if (all[s][4] > 0) 	allele[s]++;
	}

	if (tracelevel != 0)		
	{
		fprintf(fout, "\n\n");
		for (s=1; s<=NSNP; s++)	if ((s == s1)||(s == s2)) fprintf(fout, "s=%d  %d %d %d %d %d\n", s, all[s][1], all[s][2], all[s][3], all[s][4], all[s][5]);
//		for (s=1; s<=NSNP; s++)	if (allele[s]>2) fprintf(fout, "s=%d  %d %d %d %d %d\n", s, all[s][1], all[s][2], all[s][3], all[s][4], all[s][5]);
		fprintf(fout, "\n\n");
	}

	/******* Change allele codes to 1 and 2 *******/

	for (s=1; s<=NSNP; s++)
	if (allele[s] > 2)
	{
		snp_m2 ++;
	}

	for (s=1; s<=NSNP; s++)
	if (allele[s] <= 2)
	if (all[s][5] <= maxCEROS)
	{
		if ( (all[s][1] != 0) && (all[s][2] != 0) )
		{
			freq[s][1] = all[s][1]/((2.0*NIND)-all[s][5]);
			freq[s][2] = all[s][2]/((2.0*NIND)-all[s][5]);
			calc_SNP[s] = 1;
			numSNP++;
		}
		else if ( (all[s][1] != 0) && (all[s][3] != 0) )
		{
			freq[s][1] = all[s][1]/((2.0*NIND)-all[s][5]);
			freq[s][2] = all[s][3]/((2.0*NIND)-all[s][5]);

			for (i=1; i<=NIND; i++)
			for (k=1; k<=2; k++)
			if (snp[i][s][k] == 3) 	snp[i][s][k] = 2;
			calc_SNP[s] = 1;
			numSNP++;
		}
		else if ( (all[s][1] != 0) && (all[s][4] != 0) )
		{
			freq[s][1] = all[s][1]/((2.0*NIND)-all[s][5]);
			freq[s][2] = all[s][4]/((2.0*NIND)-all[s][5]);

			for (i=1; i<=NIND; i++)
			for (k=1; k<=2; k++)
			if (snp[i][s][k] == 4) 	snp[i][s][k] = 2;
			calc_SNP[s] = 1;
			numSNP++;
		}
		else if ( (all[s][2] != 0) && (all[s][3] != 0) )
		{
			freq[s][1] = all[s][2]/((2.0*NIND)-all[s][5]);
			freq[s][2] = all[s][3]/((2.0*NIND)-all[s][5]);

			for (i=1; i<=NIND; i++)
			for (k=1; k<=2; k++)
			{
				if (snp[i][s][k] == 2) 	snp[i][s][k] = 1;
				if (snp[i][s][k] == 3) 	snp[i][s][k] = 2;
			}
			calc_SNP[s] = 1;
			numSNP++;
		}
		else if ( (all[s][2] != 0) && (all[s][4] != 0) )
		{
			freq[s][1] = all[s][2]/((2.0*NIND)-all[s][5]);
			freq[s][2] = all[s][4]/((2.0*NIND)-all[s][5]);

			for (i=1; i<=NIND; i++)
			for (k=1; k<=2; k++)
			{
				if (snp[i][s][k] == 2) 	snp[i][s][k] = 1;
				if (snp[i][s][k] == 4) 	snp[i][s][k] = 2;
			}
			calc_SNP[s] = 1;
			numSNP++;
		}
		else if ( (all[s][3] != 0) && (all[s][4] != 0) )
		{
			freq[s][1] = all[s][3]/((2.0*NIND)-all[s][5]);
			freq[s][2] = all[s][4]/((2.0*NIND)-all[s][5]);

			for (i=1; i<=NIND; i++)
			for (k=1; k<=2; k++)
			{
				if (snp[i][s][k] == 3) 	snp[i][s][k] = 1;
				if (snp[i][s][k] == 4) 	snp[i][s][k] = 2;
			}
			calc_SNP[s] = 1;
			numSNP++;
		}
		else
		{
			snp_mono ++;
			mono_SNP[s] = 1;
//			printf("s=%d  ", s);
		}
	}

	for (s=1; s<=NSNP; s++)
//	if (calc_SNP[s] == 1)
	for (i=1; i<=NIND; i++)
	{
		if (((snp[i][s][1] == 1)&&(snp[i][s][2] == 2))||((snp[i][s][1] == 2)&&(snp[i][s][2] == 1)))
		numHET[s] ++;		

		if(i==2)  fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
	}

	for (s=1; s<=NSNP; s++)
//	if (calc_SNP[s] == 1)
	{
		for (i=1; i<=NIND; i++)	fprintf(fgeno, "%d%d ", snp[i][s][1], snp[i][s][2]);
		fprintf(fgeno, "%f %f\n", freq[s][2], (double)numHET[s]/(double)NIND);
	}

/*	if (tracelevel != 0)	
	{
		fprintf(fout, "\MONO, MORE THAN 2 alleles AND MORE THAN Maxceros\n");
		for (s=1; s<=NSNP; s++)
		{
			if (mono_SNP[s] == 1)
			{
				fprintf(fout, "s=%d  ", s);
				for (i=1; i<=NIND; i++)	 fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
				fprintf(fout, "\nqs=%f ps=%f\n\n", freq[s][1], freq[s][2]);
			}
			if (allele[s] > 2)
			{
				fprintf(fout, "s=%d  ", s);
				for (i=1; i<=NIND; i++)	 fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
				fprintf(fout, "\nqs=%f ps=%f\n\n", freq[s][1], freq[s][2]);
			}
			if (all[s][5] > maxCEROS)
			{
				fprintf(fout, "s=%d  ", s);
				for (i=1; i<=NIND; i++)	 fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
				fprintf(fout, "\nqs=%f ps=%f\n\n", freq[s][1], freq[s][2]);
			}
		}
		fprintf(fout, "\nCHANGED\n");

		for (s=1; s<=NSNP; s++)
		{
			if (s == s1)
			{
				fprintf(fout, "s=%d  ", s);
				for (i=1; i<=NIND; i++)	 fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
				fprintf(fout, "\nqs=%f ps=%f\n\n", freq[s][1], freq[s][2]);
			}
			if (s == s2)
			{
				fprintf(fout, "s=%d  ", s);
				for (i=1; i<=NIND; i++)	 fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
				fprintf(fout, "\nqs=%f ps=%f\n\n", freq[s][1], freq[s][2]);
			}
		}
	}
*/
	/******* Major allele is 1 and minor is 2 *******/

	for (s=1; s<=NSNP; s++)
	if (calc_SNP[s] == 1)
	if ((freq[s][2] != 0.0)&&(freq[s][2] != 1.0))
	if (freq[s][2] >= MAF)
	{
		numLOCI ++;

		if (freq[s][1] < freq[s][2])
		{
			for (i=1; i<=NIND; i++)
			for (k=1; k<=2; k++)
			{
				if (snp[i][s][k] == 1)	snp[i][s][k] = 2;
				else if (snp[i][s][k] == 2)	snp[i][s][k] = 1;
			}
			w=freq[s][1]; freq[s][1]=freq[s][2]; freq[s][2]=w;
		}
	}

	if (tracelevel != 0)	
	{
		fprintf(fout, "\nCHANGED major is 1\n");
		for (s=1; s<=NSNP; s++)
		if (calc_SNP[s] == 1) 
		if ((freq[s][2] != 0.0)&&(freq[s][2] != 1.0))
		{
			if (s == s1)
			{
				fprintf(fout, "s=%d  ", s);
				for (i=1; i<=NIND; i++)	 fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
				fprintf(fout, "\nqs=%f ps=%f\n\n", freq[s][1], freq[s][2]);
			}
			if (s == s2)
			{
				fprintf(fout, "s=%d  ", s);
				for (i=1; i<=NIND; i++)	 fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
				fprintf(fout, "\nqs=%f ps=%f\n\n", freq[s][1], freq[s][2]);
			}
		}
	}

//	EMIBD9, COLONY and MLNE input data

	int count_snp[II];
 
	for (i=1; i<=NIND; i++)	fprintf(fEMIBD9, "IND%d ", i);
	fprintf(fEMIBD9, "\n");

	for (i=1; i<=NIND; i++)
	{
		count_snp[i] = 0;

		fprintf(fcolony_dat, "IND%d ", i);
		fprintf(fmlne, "IND%d ", i);

		for (s=1; s<=NSNP; s++)
		if (calc_SNP[s] == 1) 
		if ((freq[s][2] != 0.0)&&(freq[s][2] != 1.0))
		if (freq[s][2] >= MAF)
		{
			if ((snp[i][s][1] == 1)&&(snp[i][s][2] == 1))	fprintf(fEMIBD9, "0");
			else if ((snp[i][s][1] == 1)&&(snp[i][s][2] == 2))	fprintf(fEMIBD9, "1");
			else if ((snp[i][s][1] == 2)&&(snp[i][s][2] == 1))	fprintf(fEMIBD9, "1");
			else if ((snp[i][s][1] == 2)&&(snp[i][s][2] == 2))	fprintf(fEMIBD9, "2");
			else if ((snp[i][s][1] == 0)&&(snp[i][s][2] == 0))	fprintf(fEMIBD9, "3");

			if ((snp[i][s][1] == 1)&&(snp[i][s][2] == 1))	fprintf(fmlne, "0");
			else if ((snp[i][s][1] == 1)&&(snp[i][s][2] == 2))	fprintf(fmlne, "1");
			else if ((snp[i][s][1] == 2)&&(snp[i][s][2] == 1))	fprintf(fmlne, "1");
			else if ((snp[i][s][1] == 2)&&(snp[i][s][2] == 2))	fprintf(fmlne, "2");
			else if ((snp[i][s][1] == 0)&&(snp[i][s][2] == 0))	fprintf(fmlne, "3");

			if ((snp[i][s][1] == 1)&&(snp[i][s][2] == 1))	fprintf(fcolony_dat, "1 1 ");
			else if ((snp[i][s][1] == 1)&&(snp[i][s][2] == 2))	fprintf(fcolony_dat, "1 2 ");
			else if ((snp[i][s][1] == 2)&&(snp[i][s][2] == 1))	fprintf(fcolony_dat, "2 1 ");
			else if ((snp[i][s][1] == 2)&&(snp[i][s][2] == 2))	fprintf(fcolony_dat, "2 2 ");
			else if ((snp[i][s][1] == 0)&&(snp[i][s][2] == 0))	fprintf(fcolony_dat, "0 0 ");

			count_snp[i] ++;
		}
		fprintf(fEMIBD9, "\n");
		fprintf(fcolony_dat, "\n");
		fprintf(fmlne, "\n");
//		fprintf(fEMIBD9, " ***** SNPs=%d\n", count_snp[i]);
	}

	fprintf(fcolony_mark, "mk@\n");
	fprintf(fcolony_mark, "0@\n");
	fprintf(fcolony_mark, "0.01@\n");
	fprintf(fcolony_mark, "0.01@\n");

//	SFS

	int cc, copies;
 
	for (s=1; s<=NSNP; s++)
	{
//		fprintf(fsfs, "NIND=%d s=%d freq=%f 2Nfreq=%f\n", NIND, s, freq[s][2], freq[s][2]*(2.0*NIND));
		copies = round(freq[s][2]*(2*NIND-all[s][5]));

		for (cc=1; cc<=NIND; cc++)
		{
			if (cc == copies)
			{
				FreqClass[cc] ++;
//				fprintf(fsfs, "s=%d  cc=%d  freq=%f  copies=%d all[s][5]=%d FreqClass[%d]=%d\n", s, cc, freq[s][2], copies, all[s][5], cc, FreqClass[cc]);
			}
		}
	}
	for (cc=1; cc<=NIND; cc++)		fprintf(fsfs, " %d ", FreqClass[cc]);
	fprintf(fsfs, "\n");

	return(0);
}

/* **************************************************************************** */

inbreeding()
{
	/******* Calculation of inbreeding *******/

	for (i=1; i<=NIND; i++)
	{
		NL[i] = 0.0;
		Hom[i] = 0.0;
		EHom = 0.0;
		EHet = 0.0;
		sum_LR[i] = 0.0;
		sum_LRnum[i] = 0.0;
		sum_LRden[i] = 0.0;
		sum_VR1[i] = 0.0;
		sum_VR2[i] = 0.0; 
		sum_YANG1[i] = 0.0; 
		sum_YANG2[i] = 0.0; 
		sum_LH2[i] = 0.0;

		for (s=1; s<=NSNP; s++)
		if (calc_SNP[s] == 1)
		if ((freq[s][2] != 0.0)&&(freq[s][2] != 1.0))
		if (freq[s][2] >= MAF)
		{
			if ((snp[i][s][1] == 2) && (snp[i][s][2] == 2))
			{
				Hom[i] ++;
				NL[i] ++;
//				sum_LRnum[i] += ((1.0/(1.0-freq[s][2])) - 1.0) / ( (1.0 - 1.0)*(0.5/freq[s][2] + 0.5/(1.0-freq[s][2]) - 1.0) + (1.0*(1/(1.0-freq[s][2]) - 1.0)) );
//				sum_LRnum[i] += 1.0;
				sum_LR[i] += 1.0;
				sum_VR1[i] += (2.0-(2.0*freq[s][2]))*(2.0-(2.0*freq[s][2]));
				sum_VR2[i] += ( ( (2.0-2.0*freq[s][2])*(2.0-2.0*freq[s][2]) ) / (2.0*freq[s][2]*(1.0-freq[s][2])) ) - 1.0; 
				sum_YANG1[i] += ( (4.0)-((1.0+2.0*freq[s][2])*2.0)+(2.0*freq[s][2]*freq[s][2]) ); 
				sum_YANG2[i] += ( (4.0)-((1.0+2.0*freq[s][2])*2.0)+(2.0*freq[s][2]*freq[s][2]) ) / (2.0*freq[s][2]*(1.0-freq[s][2])); 
				EHom += (1.0 - 2.0*freq[s][2]*(1.0-freq[s][2]));
				EHet += 2.0*freq[s][2]*(1.0-freq[s][2]);
			}
			else if ((snp[i][s][1] == 1) && (snp[i][s][2] == 1))
			{
				Hom[i] ++;
				NL[i] ++;
//				sum_LRnum[i] += ((1.0/(1.0-freq[s][2])) - 1.0) / ( (1.0 - 1.0)*(0.5/freq[s][2] + 0.5/(1.0-freq[s][2]) - 1.0) + (1.0*(1/(1.0-freq[s][2]) - 1.0)) );
//				sum_LRnum[i] += 1.0;
				sum_LR[i] += 1.0;
				sum_VR1[i] += (0.0-(2.0*freq[s][2]))*(0.0-(2.0*freq[s][2]));
				sum_VR2[i] += ( ( (0.0-(2.0*freq[s][2]))*(0.0-(2.0*freq[s][2])) ) / (2.0*freq[s][2]*(1.0-freq[s][2])) ) - 1.0; 
				sum_YANG1[i] += (2.0*freq[s][2]*freq[s][2]); 
				sum_YANG2[i] += (2.0*freq[s][2]*freq[s][2]) / (2.0*freq[s][2]*(1.0-freq[s][2])); 
				EHom += (1.0 - 2.0*freq[s][2]*(1.0-freq[s][2]));
				EHet += 2.0*freq[s][2]*(1.0-freq[s][2]);
			}
			else if ( ((snp[i][s][1] == 1) && (snp[i][s][2] == 2)) || ((snp[i][s][1] == 2) && (snp[i][s][2] == 1)) )
			{
				NL[i] ++;
//				sum_LRnum[i] += ((0.0/(1.0-freq[s][2])) - 1.0) / ( (1.0 - 0.0)*(0.5/freq[s][2] + 0.5/(1.0-freq[s][2]) - 1.0) + (0.0*(1/(1.0-freq[s][2]) - 1.0)) );
//				sum_LRnum[i] += -1.0;
//				sum_LRden[i] += ( 0.5/freq[s][2] + 0.5/(1.0-freq[s][2]) - 1.0);
				sum_LR[i] += -1.0 / ( 0.5/freq[s][2] + 0.5/(1.0-freq[s][2]) - 1.0);
				sum_VR1[i] += (1.0-(2.0*freq[s][2]))*(1.0-(2.0*freq[s][2]));
				sum_VR2[i] += ( ( (1.0-(2.0*freq[s][2]))*(1.0-(2.0*freq[s][2])) ) / (2.0*freq[s][2]*(1.0-freq[s][2])) ) - 1.0; 
				sum_YANG1[i] += ( (1.0)-((1.0+2.0*freq[s][2])*1.0)+(2.0*freq[s][2]*freq[s][2]) ); 
				sum_YANG2[i] += ( (1.0)-((1.0+2.0*freq[s][2])*1.0)+(2.0*freq[s][2]*freq[s][2]) ) / (2.0*freq[s][2]*(1.0-freq[s][2])); 
				sum_LH2[i] += ( 1.0 / (2.0*freq[s][2]*(1.0-freq[s][2])) ); 
				EHom += (1.0 - 2.0*freq[s][2]*(1.0-freq[s][2]));
				EHet += 2.0*freq[s][2]*(1.0-freq[s][2]);
			}
		}

//		fLR[i] = sum_LRnum[i] / sum_LRden[i];
		fLR[i] = sum_LR[i] / NL[i];
		fvr1[i] = (sum_VR1[i]/EHet)-1.0;
		fvr2[i] = sum_VR2[i] / NL[i];
		fyang1[i] = sum_YANG1[i] / EHet;
		fyang2[i] = sum_YANG2[i] / NL[i];
		fLH1[i] = (Hom[i]-EHom) / (NL[i]-EHom);
		fLH2[i] = 1.0 - ( sum_LH2[i] / NL[i]);
		fhom[i] = Hom[i] / NL[i];

		if (tracelevel != 0) if(i==2)	fprintf(fout, "i = %d   fLR = %f   sum_LRnum = %f   sum_LRden = %f   fhom = %f   fLH1 = %f   fvr1 = %f\n", i, fLR[i],  sum_LRnum[i],  sum_LRden[i],  fhom[i],  fLH1[i], fvr1[i]);		
		if (tracelevel != 0) if(i==2)	fprintf(fout, "i = %d   NL = %f   Hom = %f   EHom = %f   fhom = %f   fLH1 = %f   fvr1 = %f\n", i, NL[i], Hom[i], EHom, fhom[i],  fLH1[i], fvr1[i]);		
	}

	if (tracelevel != 0)
	{
		fprintf(fout, "\nInbreeding calculations\n");
		fprintf(fout, "NSNP = %d \n", NSNP);
		fprintf(fout, "numLOCI = %d \n", (int)numLOCI);
		fprintf(fout, "snp_mono = %d \n", snp_mono);
		fprintf(fout, "snp_m2 = %d \n", snp_m2);
	}

	fprintf(finb, "n sex pheL pheW pheFC Fvr1 Fvr2 Fyang1 Fyang2 FLH1 FLH2 Fhom FLR\n");
//	fprintf(finb, "n sex pheL pheW pheFC Fvr1 Fvr2 Fyang1 Fyang2 FLH1 FLH2 Fhom\n");
	for (i=1; i<=NIND; i++)
//	fprintf(finb, "%d  %d  %6.1f  %6.1f  %f  %f  %f  %f  %f  %f  %f  %f  \n", i, sex[i], pheL[i], pheW[i], (pheW[i] / pow(pheL[i],3.0))*100000, fvr1[i], fvr2[i], fyang1[i], fyang2[i], fLH1[i], fLH2[i], fhom[i]);

	fprintf(finb, "%d  %d  %6.1f  %6.1f  %f  %f  %f  %f  %f  %f  %f  %f  %f  \n", i, sex[i], pheL[i], pheW[i], (pheW[i] / pow(pheL[i],3.0))*100000, fvr1[i], fvr2[i], fyang1[i], fyang2[i], fLH1[i], fLH2[i], fhom[i], fLR[i]);

	return(0);
}

/* ********************************************************************************************* */

lookfortext(s)
char *s;
{
   int len, i, curchar;
   char c;

   curchar = 0;
   len = 0;

   for (i=0; i<=100; i++)
   {
      if (s[i] == '\0') break;
      len++;
   }
   do
   {
      c = getc(fped);

      if (c==s[curchar])
      {
         curchar++;
         if (curchar==len) return(0);
      }
      else curchar = 0;
   }
   while (c != EOF);
}

/* ********************************************************************************************* */

