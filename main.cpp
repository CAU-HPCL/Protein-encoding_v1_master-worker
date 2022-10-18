#include <iostream>
#include <cstdlib>
#include<stdio.h>
#include <math.h>
#include <time.h>

using namespace std;

#define _CRT_SECURE_NO_WARNINGS
#define OBJECTIVE_NUM 3		// number of objective function
#define MAX_CODON 6			

enum AA {
	A = 'A', C = 'C', D = 'D', E = 'E', F = 'F', G = 'G', H = 'H', I = 'I', K = 'K', L = 'L', 
	M = 'M', N = 'N', P = 'P', Q = 'Q', R = 'R', S = 'S', T = 'T', V = 'V', W = 'W', Y = 'Y'
};
typedef struct{
	AA name;			 			// amino acid abbreviation
	int num_codons;					// number of codon types in a amino acid
	char codons[MAX_CODON][4];		// codons in a amino acid
	double adaptation[MAX_CODON];	// codons's adaptation weight
}Aminoacids;

/* 20 kinds of amino acids */
/* adaptation weight is ascending order */
const Aminoacids aa[20] = {			
	{A, 4, {"GCG", "GCA", "GCC", "GCU"}, {1854 / 13563.0, 5296 / 13563.0, 7223 / 13563.0, 1}},
	{C, 2, {"UGC", "UGU"}, {1234 / 3052.0, 1}},
	{D, 2, {"GAC", "GAU"}, {8960 / 12731.0, 1}},
	{E, 2, {"GAG", "GAA"}, {6172 / 19532.0, 1}},
	{F, 2, {"UUU", "UUC"},{7773 / 8251.0, 1}},
	{G, 4, {"GGG", "GGA", "GGC", "GGU"},{1852 / 15694.0, 2781 / 15694.0, 3600 / 15694.0, 1}},
	{H, 2, {"CAC", "CAU"}, {3288 / 4320.0, 1}},
	{I, 3, {"AUA", "AUC", "AUU"},{3172 / 12071.0, 8251 / 12071.0, 1}},
	{K, 2, {"AAA", "AAG"},{12845 / 15169.0, 1}},
	{L, 6, {"CUC", "CUG", "CUU", "CUA", "UUA", "UUG"}, {1242 / 13329.0, 2852 / 13329.0, 3207 / 13329.0, 4134 / 13329.0, 8549 / 13329.0, 1}},
	{M, 1, {"AUG"}, {1}},
	{N, 2, {"AAU", "AAC"}, {8613 / 9875.0, 1}},
	{P, 4, {"CCG", "CCC", "CCU", "CCA"}, {1064 / 8965.0, 1656 / 8965.0, 4575 / 8965.0, 1}},
	{Q, 2, {"CAG", "CAA"}, {3312 / 10987.0, 1}},
	{R, 6, {"CGG", "CGA", "CGC", "AGG", "CGU", "AGA"}, {342 / 9784.0, 489 / 9784.0, 658 / 9784.0, 2175 / 9784.0, 3307 / 9784.0, 1}},
	{S, 6, {"UCG", "AGC", "AGU", "UCA", "UCC", "UCU"}, {2112 / 10025.0, 2623 / 10025.0, 3873 / 10025.0, 4583 / 10025.0, 6403 / 10025.0, 1}},
	{T, 4, {"ACG", "ACA", "ACC", "ACU"}, {1938 / 9812.0, 5037 / 9812.0, 6660 / 9812.0, 1}},
	{V, 4, {"GUA", "GUG", "GUC", "GUU"}, {3249 / 11442.0, 3700 / 11442.0, 6911 / 11442.0, 1}},
	{W, 1, {"UGG"}, {1}},
	{Y, 2, {"UAU", "UAC"}, {5768 / 7114.0, 1}}
};

typedef struct{
	char** CDSs;						// CDS's bases sequences
	double obj_val[OBJECTIVE_NUM];		// checking objective function value (0 ~ 1) to Pareto Comparsion
}Solution;
typedef struct{
	int counter;				// checking counter to obsolete this solution
	int rank;					// indicate Pareto front (rank)
	double crowding_distance;	// indicate diversity of solution in same rank
	double sel_prob;			// selection probability
	Solution solution;			 
}Population;
 

/* ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- */
/* calculate objective function value */
/* this function is calculate mininum CAI value */
double mCAI(const char* amino_sequence, char** CDSs, int num_CDSs, int len_amino_sequence)
{
	char codon[3];
	bool check;
	double min = 1;


	double tmp;
	int idx;
	for (int i = 0; i < num_CDSs; i++){
		idx = 0;
		tmp = 1;
		for (int j = 0; j < len_amino_sequence; j++){
			check = false;
			codon[0] = CDSs[i][idx++];
			codon[1] = CDSs[i][idx++];
			codon[2] = CDSs[i][idx++];
			for (int k = 0; k < 20; k++){
				if (aa[k].name == amino_sequence[j]) {
					for (int l = 0; l < aa[k].num_codons; l++) {
						if ((codon[0] == aa[k].codons[l][0]) && (codon[1] == aa[k].codons[l][1]) && (codon[2] == aa[k].codons[l][2])) {
							tmp *= pow(aa[k].adaptation[l], 1.0 / len_amino_sequence);
							check = true;
							break;
						}
					}
				}
				if (check) break;
			}
		}
		if (tmp < min) min = tmp;
	}

	return min;
}

/* this function is calculate minimum Hamming Distance */
double mHD(char** CDSs, int num_CDSs, int len_CDS)
{
	int cnt = 0;
	double result = 1.;

	for (int i = 0; i < num_CDSs - 1; i++)
	{
		for (int j = i + 1; j < num_CDSs; j++)
		{
			for (int k = 0; k < len_CDS; k++)
			{
				if (CDSs[i][k] != CDSs[j][k]) cnt++;
			}
			if ((double)cnt/len_CDS < result) result = (double)cnt / len_CDS;
			cnt = 0;
		}
	}

	return result;
}

/* this function calculate maximun length of common substring */
double MLRCS(const char** CDSs, int num_CDS, int len)
{
	return 1;
}


/* this function caculate selection probability */
double CalSelectionProb()
{
	return 1;
}
/* this function sorting by rank and crowding distance */
void Sort(Population* population, int size)
{
	return;
}


/* this function make random CDS */										// free memory is needed
char* GenCDS(const char* amino_sequence, int len_amino_sequence)
{
	char* str = (char*)malloc(sizeof(char) * 3 * len_amino_sequence + 1);	
	int idx = 0;
	int tmp;

	for (int i = 0; i < len_amino_sequence; i++) {
		for (int j = 0; j < 20; j++) {
			if (aa[j].name == amino_sequence[i]) {
				tmp = rand() % aa[j].num_codons;						// random seed check needed
				str[idx++] = aa[j].codons[tmp][0];
				str[idx++] = aa[j].codons[tmp][1];
				str[idx++] = aa[j].codons[tmp][2];
				break;
			}
		}
	}
	str[idx] = NULL;

	return str;
}

/* this function generate random solution */
void GenSolution(Population *population, int num_CDSs, const char* amino_sequence, int len_amino_sequence)
{
	/* initial value setting */
	population->counter = 0;			
	population->rank = 0;
	population->crowding_distance = 0;
	population->sel_prob = 0;

	/* make random CDSs */
	for (int i = 0; i < num_CDSs; i++) {
		population->solution.CDSs[i] = GenCDS(amino_sequence, len_amino_sequence);
	}

	return;
}

/* this function generate random mutated solution */
//Population Mutate(Population population)
//{
//
//	return;
//}


/* show population's attribute values */
void PrintPopulation(const Population* population, int pop_size, int num_CDSs, int len_amino_sequence)
{
	cout << "\n ----------------------- Print Population ---------------------------- " << endl;
	
	for (int i = 0; i < pop_size; i++) {			// population outer loop
		cout << "population" << "[" << i << "]" << endl;
		cout << "rank : " << population[i].rank << endl;
		cout << "crowding distance : " << population[i].crowding_distance << endl;
		cout << "mCAI : " << population[i].solution.obj_val[0] << endl;
		cout << "mHD : " << population[i].solution.obj_val[1] << endl;
		cout << "MLRCS : " << population[i].solution.obj_val[2] << endl;
		for (int j = 0; j < num_CDSs; j++) {		// population's solution loop
			cout << "CDS" << "[" << j << "] :" << population[i].solution.CDSs[j] << endl;
		}
	}

	return;
}


int main()
{
	srand(time(NULL));
	
	/* amino sequence recieve from FASTA format */
	char file_name[20] = "Q5VZP5.fasta.txt";
	char buffer[256];
	char *amino_sequence;		// amino sequence comprising a CDS
	int len_amino_sequence;		// length of amino sequence
	
	/* file processing */
	//cout << "input file name : "; cin >> file_name;
	FILE* fp;
	fopen_s(&fp, file_name, "r");
	if (fp == NULL) {
		cout << "opening input file failed" << endl;
		return EXIT_FAILURE;
	}
	fseek(fp, 0, SEEK_END);
	len_amino_sequence = ftell(fp);
	fseek(fp, 0, SEEK_SET);
	fgets(buffer, 256, fp);				// emptying first line
	len_amino_sequence -= ftell(fp);
	amino_sequence = (char*)malloc(sizeof(char) * len_amino_sequence);
	
	int idx = 0;
	char tmp;
	while (!feof(fp)) {
		tmp = fgetc(fp);
		if (tmp != '\n')amino_sequence[idx++] = tmp;
	}
	amino_sequence[idx] = NULL;
	len_amino_sequence = idx - 1;
	
	fclose(fp);
	/* -------------------- end file -------------------- */


	Population* population;

	/* user input parameter */
	int colony_size;			// number of solutions in population
	int max_cycle;				// number of generations
	int limit;					// number of solution is not updated
	double mprob;				// mutation probability
	int num_CDSs;				// number of CDSs


	/* input parameter values */
	cout << "input colony size : ";	cin >> colony_size;
	cout << "input max_cycle : ";	cin >> max_cycle;
	cout << "input limit : ";	cin >> limit;
	cout << "input mutation probability (0 ~ 1 value): ";	cin >> mprob;
	cout << "input number of CDSs : ";	cin >> num_CDSs;


	/* memory allocation */
	population = (Population*)malloc(sizeof(Population) * (2 * colony_size));
	for (int i = 0; i < (2 * colony_size); i++) {
		population[i].solution.CDSs = (char**)malloc(sizeof(char*) * num_CDSs);
	}
	

	/* initialize Population */
	for (int i = 0; i < colony_size - 1; i++) {
		GenSolution(&population[i], num_CDSs, amino_sequence, len_amino_sequence);
		/* calculate objective function value */
		population[i].solution.obj_val[0] = mCAI(amino_sequence, population->solution.CDSs, num_CDSs, len_amino_sequence);
		population[i].solution.obj_val[1] = mHD(population->solution.CDSs, num_CDSs, len_amino_sequence * 3);
		//population->solution.obj_val[2] = MLRCS(population->solution.CDSs, num_CDSs, len_amino_sequence * 3);
	}
	/* To boost the optimization of solutions witth high CAI values 
	   remaing solution is generated by selecting highest adaptation */

	
	// temp print
	PrintPopulation(population, colony_size - 1, num_CDSs, len_amino_sequence);




	/* free memory */
	for (int i = 0; i < (2*colony_size); i++) {
		for (int j = 0; j < num_CDSs; j++) {
			free(population[i].solution.CDSs[j]);
		}
		free(population[i].solution.CDSs);
	}
	free(population);
	free(amino_sequence);

	return EXIT_SUCCESS;
}