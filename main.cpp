#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define _CRT_SECURE_NO_WARNINGS
#define MAX_CODON 6			

/* ------------------------------------------------- Amino acids, codons and adaptation weigth definition ---------------------------------------- */
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
/* --------------------------------------------------------- end definition ---------------------------------------------------------------------- */


#define OBJECTIVE_NUM 3				// number of objective function
#define _mCAI 0
#define _mHD 1
#define _MLRCS 2
/* --------------------------------------------------------- Population definition ---------------------------------------------------------------- */
typedef struct{
	char** cds;						// CDS's bases sequences
	int p, q, l;						// this if for MLRCS starting point and length
	int obj_cdsidx[OBJECTIVE_NUM][2];	// CDS's index correspond to objective function
	double obj_val[OBJECTIVE_NUM];		// checking objective function value (0 ~ 1) to Pareto Comparsion
}Solution;
typedef struct{
	int counter;						// checking counter to obsolete this solution
	int rank;							// indicate Pareto front (rank)
	double crowding_distance;			// indicate diversity of solution in same rank
	double sel_prob;					// selection probability
	Solution sol;			 
}Population;
/* ----------------------------------------------------------- end definition --------------------------------------------------------------------- */


/* ----------------------------------------------------------- Find index definition --------------------------------------------------------------- */
/* this function find aminoacid index using binary search */
int FindAminoIndex(const AA aminoacid)
{
	int low = 0;
	int high = 20 - 1;
	int mid;

	while (low <= high) {
		mid = (low + high) / 2;

		if (aa[mid].name == aminoacid)
			return mid;
		else if (aa[mid].name > aminoacid)
			high = mid - 1;
		else
			low = mid + 1;
	}
	
}
/* this function find aminoacid's codon index */
int FindCodonIndex(int amino_idx, const char* codon)
{
	for (int i = 0; i < aa[amino_idx].num_codons; i++) {
		if (aa[amino_idx].codons[i][0] == codon[0] &&
			aa[amino_idx].codons[i][1] == codon[1] &&
			aa[amino_idx].codons[i][2] == codon[2]) {
			return i;
		}
	}
}
/* --------------------------------------------------------------- end definition ------------------------------------------------------------------- */


/* ---------------------------------------------------------- Population memroy management --------------------------------------------------------- */
/* this function memory allocation population */
// 나중에 포인터 매개변수 메모리 할당 체크 필요!!
Population* AllocPopulation(int pop_size, int num_cds, int len_amino_seq)
{
	Population* pop;

	pop = (Population*)malloc(sizeof(Population) * pop_size);
	if (pop == NULL) {
		printf("Memory allocation failed at line %d", __LINE__);
		exit(EXIT_FAILURE);
	}
	
	for (int i = 0; i < pop_size; i++) {
		pop[i].sol.cds = (char**)malloc(sizeof(char*) * num_cds);
		if (pop[i].sol.cds == NULL) {
			printf("Memory allocation failed at line %d", __LINE__);
			exit(EXIT_FAILURE);
		}
		for (int j = 0; j < num_cds; j++) {
			pop[i].sol.cds[j] = (char*)malloc(sizeof(char) * (len_amino_seq * 3 + 1));
			if (pop[i].sol.cds[j] == NULL) {
				printf("Memory allocation failed at line %d", __LINE__);
				exit(EXIT_FAILURE);
			}
		}
	}

	return pop;
}
/* this function free population memory */
void FreePopulation(Population* pop, int pop_size, int num_cds)
{
	for (int i = 0; i < pop_size; i++) {
		for (int j = 0; j < num_cds; j++) {
			free(pop[i].sol.cds[j]);
		}
		free(pop[i].sol.cds);
	}
	free(pop);

	return;
}
/* ---------------------------------------------------------- end memory management ---------------------------------------------------------------- */


/* ------------------------------------------------------- make Population and mutate Population fucntions ----------------------------------------- */
#define RANDOM_GEN 0
#define UPPER_GEN 1
/* this function make random CDS */
void GenCDS(char* cds, const int* amino_seq_idx, int len_amino_seq, int type = RANDOM_GEN)
{
	int idx = 0;
	int rand_idx;
	int codon_idx;

	switch (type) {
	/* random creation */
	case RANDOM_GEN:	
		for (int i = 0; i < len_amino_seq; i++) {
			rand_idx = rand() % aa[amino_seq_idx[i]].num_codons;
			cds[idx++] = aa[amino_seq_idx[i]].codons[rand_idx][0];
			cds[idx++] = aa[amino_seq_idx[i]].codons[rand_idx][1];
			cds[idx++] = aa[amino_seq_idx[i]].codons[rand_idx][2];
		}
		break;
	/* creation maximum CAI values */
	case UPPER_GEN:
		for (int i = 0; i < len_amino_seq; i++) {
			codon_idx = aa[amino_seq_idx[i]].num_codons - 1;
			cds[idx++] = aa[amino_seq_idx[i]].codons[codon_idx][0];
			cds[idx++] = aa[amino_seq_idx[i]].codons[codon_idx][1];
			cds[idx++] = aa[amino_seq_idx[i]].codons[codon_idx][2];
		}
		break;
	}
	cds[idx] = NULL;

	return;
}
/* this function generate random solution */
void GenSolution(Population* pop, int num_cds, const int* amino_seq_idx, int len_amino_seq, int type = RANDOM_GEN)
{
	/* initial value setting */
	pop->counter = 0;
	pop->rank = 0;
	pop->crowding_distance = 0;
	pop->sel_prob = 0;

	/* make random CDSs */
	for (int i = 0; i < num_cds; i++) {
		GenCDS(pop->sol.cds[i], amino_seq_idx, len_amino_seq, type);
	}

	return;
}

/* this function copy population */
void CopyPopulation(const Population* origin, Population* target, int num_cds, int len_amino_seq)
{
	target->counter = origin->counter;
	target->rank = origin->rank;
	target->crowding_distance = origin->crowding_distance;
	target->sel_prob = origin->sel_prob;
	
	for (int i = 0; i < num_cds; i++) {
		for (int j = 0; j < len_amino_seq * 3 + 1; j++) {
			target->sol.cds[i][j] = origin->sol.cds[i][j];
		}
	}

	for (int i = 0; i < OBJECTIVE_NUM; i++) {
		target->sol.obj_val[i] = origin->sol.obj_val[i];
		target->sol.obj_cdsidx[i][0] = origin->sol.obj_cdsidx[i][0];
		target->sol.obj_cdsidx[i][1] = origin->sol.obj_cdsidx[i][1];
	}
	target->sol.p = origin->sol.p;
	target->sol.q = origin->sol.q;
	target->sol.l = origin->sol.l;

	return;
}
#define RANDOM_ADAPTATION 0
#define UPPER_ADAPTATION 1
/* this function change codon into synonymous codon which is not same input codon excluding  number of synonymous codons is one */
void ChageSynonymousCodon(AA amino, char* cds, int cd_idx, int type)
{
	char codon[3];
	codon[0] = cds[cd_idx];
	codon[1] = cds[cd_idx + 1];
	codon[2] = cds[cd_idx + 2];
	int idx;
	int rand_idx;
	int amino_idx;
	
	amino_idx = FindAminoIndex(amino);
	idx = FindCodonIndex(amino_idx, codon);

	switch (type)
	{	
	case RANDOM_ADAPTATION:			// change random synomynous codon
		if (aa[amino_idx].num_codons > 1) {
			while (true) {
				rand_idx = rand() % aa[amino_idx].num_codons;
				if (idx != rand_idx)
					break;
			}
			cds[cd_idx] = aa[amino_idx].codons[rand_idx][0];
			cds[cd_idx + 1] = aa[amino_idx].codons[rand_idx][1];
			cds[cd_idx + 2] = aa[amino_idx].codons[rand_idx][2];
		}
		break;
	case UPPER_ADAPTATION:			// change random synomynous codon which has high adaptation value
		if (idx < aa[amino_idx].num_codons - 1) {
			while (true) {
				rand_idx = rand() % aa[amino_idx].num_codons;
				if (idx < rand_idx)
					break;
			}
			cds[cd_idx] = aa[amino_idx].codons[rand_idx][0];
			cds[cd_idx + 1] = aa[amino_idx].codons[rand_idx][1];
			cds[cd_idx + 2] = aa[amino_idx].codons[rand_idx][2];
		}
		break;
	}

	return;
}
/* this function generate random mutated solution */
Population *Mutation(const Population* pop, int num_cds, const char* amino_seq, int len_amino_seq, double mprob)
{
	/* new population memory allocation */
	Population* new_pop;
	new_pop = AllocPopulation(1, num_cds, len_amino_seq);

	// copy population to new_population
	CopyPopulation(pop, new_pop, num_cds, len_amino_seq);
	new_pop->counter = 0;
	new_pop->rank = 0;
	new_pop->crowding_distance = 0;
	new_pop->sel_prob = 0;


	/* generate (0 ~ 1) random number corresponding to codon in CDS */
	double** random_num;
	random_num = (double**)malloc(sizeof(double*) * num_cds);
	for (int i = 0; i < num_cds; i++) {
		random_num[i] = (double*)malloc(sizeof(double) * len_amino_seq);
	}
	for (int i = 0; i < num_cds; i++) {
		for (int j = 0; j < len_amino_seq; j++) {
			random_num[i][j] = ((double)(rand() % 100001)) / 100000;
		}
	}
	

	int type;
	type = rand() % 4;
	/* four type of variations */
	switch (type)
	{
		/* change each codons with random codon in all CDSs */
	case 0:
		for (int i = 0; i < num_cds; i++) {
			for (int j = 0; j < len_amino_seq; j++) {
				if (random_num[i][j] <= mprob) {
					ChageSynonymousCodon((AA)amino_seq[j], new_pop->sol.cds[i], j * 3, RANDOM_ADAPTATION);
				}
			}
		}
		break;
		/* change each codons with higher adaptation value in CDS which have minimun CAI value */
	case 1:
		for (int i = 0; i < len_amino_seq; i++) {
			if (random_num[new_pop->sol.obj_cdsidx[_mCAI][0]][i] <= mprob) {
				ChageSynonymousCodon((AA)amino_seq[i], new_pop->sol.cds[new_pop->sol.obj_cdsidx[_mCAI][0]], i * 3, UPPER_ADAPTATION);
			}
		}
		break;
		/* change each codons with random codon in pair of CDSs with minimun Hamming Distance */
	case 2:
		for (int i = 0; i < len_amino_seq; i++) {
			if (random_num[new_pop->sol.obj_cdsidx[_mHD][0]][i] <= mprob) {
				ChageSynonymousCodon((AA)amino_seq[i], new_pop->sol.cds[new_pop->sol.obj_cdsidx[_mHD][0]], i * 3, RANDOM_ADAPTATION);
			}
			if (random_num[new_pop->sol.obj_cdsidx[_mHD][1]][i] <= mprob) {
				ChageSynonymousCodon((AA)amino_seq[i], new_pop->sol.cds[new_pop->sol.obj_cdsidx[_mHD][1]], i * 3, RANDOM_ADAPTATION);
			}
		}
		break;
		/* chage each codons with random codon in CDSs with longest common substring */
	case 3:
		for (int i = new_pop->sol.p / 3; i < (new_pop->sol.p + new_pop->sol.l - 1) / 3; i++) {
			if (random_num[new_pop->sol.obj_cdsidx[_MLRCS][0]][i] <= mprob) {
				ChageSynonymousCodon((AA)amino_seq[i], new_pop->sol.cds[new_pop->sol.obj_cdsidx[_MLRCS][0]], i * 3, RANDOM_ADAPTATION);
			}
		}
		for (int i = new_pop->sol.q / 3; i < (new_pop->sol.q + new_pop->sol.l - 1) / 3; i++) {
			if (random_num[new_pop->sol.obj_cdsidx[_MLRCS][1]][i] <= mprob) {
				ChageSynonymousCodon((AA)amino_seq[i], new_pop->sol.cds[new_pop->sol.obj_cdsidx[_MLRCS][1]], i * 3, RANDOM_ADAPTATION);
			}
		}
		break;
	}

	/* free memory */
	for (int i = 0; i < num_cds; i++) {
		free(random_num[i]);
	}
	free(random_num);

	return new_pop;
}
/* -------------------------------------------------------- end Population creation and mutaion ----------------------------------------------------- */


/* --------------------------------------------------------- calculate objective function value ----------------------------------------------------- */
/* this function is calculate mininum CAI value */
void mCAI(Population *pop, int num_cds, const int* amino_seq_idx, int len_amino_seq)
{
	char codon[3];
	int idx;
	int codon_idx;
	double tmp;
	
	pop->sol.obj_val[_mCAI] = 1;
	for (int i = 0; i < num_cds; i++) {
		idx = 0;
		tmp = 1;
		for (int j = 0; j < len_amino_seq; j++) {
			codon[0] = pop->sol.cds[i][idx++];
			codon[1] = pop->sol.cds[i][idx++];
			codon[2] = pop->sol.cds[i][idx++];
			codon_idx = FindCodonIndex(amino_seq_idx[j], codon);
			tmp *= pow(aa[amino_seq_idx[j]].adaptation[codon_idx], 1.0 / len_amino_seq);
		}
		if (tmp <= pop->sol.obj_val[_mCAI]) {
			pop->sol.obj_val[_mCAI] = tmp;
			pop->sol.obj_cdsidx[_mCAI][0] = i;		// CDS's index having mCAI value
		}
	}

	return;
}
/* this function is calculate minimum Hamming Distance */
void mHD(Population *pop, int num_cds, int len_amino_seq)
{
	int len_CDS = len_amino_seq * 3;
	int cnt;
	double tmp;

	pop->sol.obj_val[_mHD] = 1;
	for (int i = 0; i < num_cds - 1; i++)
	{
		for (int j = i + 1; j < num_cds; j++)
		{
			cnt = 0;
			for (int k = 0; k < len_CDS; k++)
			{
				if (pop->sol.cds[i][k] != pop->sol.cds[j][k]) cnt++;
			}
			tmp = (double)cnt / len_CDS;
			if (tmp <= pop->sol.obj_val[_mHD]) {
				pop->sol.obj_val[_mHD] = tmp;
				pop->sol.obj_cdsidx[_mHD][0] = i;
				pop->sol.obj_cdsidx[_mHD][1] = j;
			}
		}
	}

	return;
}
/* this function calculate maximun length of common substring */
void MLRCS(Population *pop, int num_cds, int len_amino_seq)
{
	int** LCS;						// matrix which size [length of CDS + 1][length of CDS + 1]
	int len_cds = 3 * len_amino_seq;
	int max_len = 0;
	int tmp;


	/* memory allocation for LCS matrix */
	LCS = (int**)malloc(sizeof(int*) * (len_cds + 1));
	for (int i = 0; i < len_cds + 1; i++) {
		LCS[i] = (int*)malloc(sizeof(int) * (len_cds + 1));
	}


	for (int i = 0; i < num_cds; i++) {							// i'th CDS
		for (int j = i; j < num_cds; j++) {						// j'th CDS
			for (int k = 0; k < len_cds + 1; k++){				// matrix row
				for (int l = 0; l < len_cds + 1; l++) {			// matrix column
					if (i != j) {
						if (k == 0 || l == 0)
							LCS[k][l] = 0;
						else if (pop->sol.cds[i][l - 1] == pop->sol.cds[j][k - 1])
						{
							tmp = LCS[k - 1][l - 1] + 1;
							LCS[k][l] = tmp;
							if (tmp >= max_len)
							{
								max_len = tmp;
								pop->sol.p = l - max_len;
								pop->sol.q = k - max_len;
								pop->sol.obj_cdsidx[_MLRCS][0] = i;
								pop->sol.obj_cdsidx[_MLRCS][1] = j;
							}
						}
						else
							LCS[k][l] = 0;
					}
					else
					{
						if (k == 0 || l == 0 || k == l)				
							LCS[k][l] = 0;
						else if (pop->sol.cds[i][l - 1] == pop->sol.cds[j][k - 1])
						{
							tmp = LCS[k - 1][l - 1] + 1;
							LCS[k][l] = tmp;
							if ((tmp >= max_len) && (l - 1 < k - tmp || k - 1 < l - tmp))
							{
								max_len = tmp;
								pop->sol.p = l - max_len;
								pop->sol.q = k - max_len;
								pop->sol.obj_cdsidx[_MLRCS][0] = i;
								pop->sol.obj_cdsidx[_MLRCS][1] = j;
							}
						}
						else
							LCS[k][l] = 0;
					}
				}
			}
		}
	}
	pop->sol.l = max_len;
	pop->sol.obj_val[_MLRCS] = (double)max_len / len_cds;

	/* free memory LCS matrix */
	for (int i = 0; i < len_cds + 1; i++) {
		free(LCS[i]);
	}
	free(LCS);

	return;
}
/* ---------------------------------------------------------- end objective function caculation ------------------------------------------------------ */


/* if new_population is dominate return true */
bool ParetoComparison(const Population* new_population, const Population* population)
{
	if (new_population->sol.obj_val[0] >= population->sol.obj_val[0] &&
		new_population->sol.obj_val[1] >= population->sol.obj_val[1] &&
		new_population->sol.obj_val[2] >= population->sol.obj_val[2]) 
		return true;
	else 
		return false;
}

/* this function sorting by rank and crowding distance */
void SortbyRankCrowding(const Population* population, int pop_size)
{
	/* memory allocation */
	Population** pop_set;			// store population to sorting
	pop_set = (Population**)malloc(sizeof(Population*) * pop_size);
	for (int i = 0; i < pop_size; i++) {
		pop_set[i] = (Population*)malloc(sizeof(Population) * pop_size);
	}

	for (int i = 0; i < pop_size; i++) {

	}


	/* free memory */
	for (int i = 0; i < pop_size; i++) {
		free(pop_set[i]);
	}
	free(pop_set);

	return;
}

/* this function caculate selection probability */
/*  */
double CalSelectionProb(Population* pop)
{
	return 1;
}

/* Roulette-wheel selection based on selection probability */
int SelectSolution(const Population* pop, int pop_size)
{
	double sum_fitness = 0;
	double sum = 0;
	double point;

	for (int i = 0; i < pop_size; i++) {
		sum_fitness += pop[i].sel_prob;
	}

	point = rand() % (int)(sum_fitness * 10000000) / 10000000.;

	for (int i = 0; i < pop_size; i++) {
		sum += pop[i].sel_prob;
		if (point < sum)
			return i;			// return selected pop index
	}
}


/* print population's attribute values */
void PrintPopulation(const Population* population, int num_CDSs)
{
	static int cnt = 0;				// for number of function calls

	printf("\n ------------------------------ Print Population ------------------------------ \n ");
	printf("\tcounter : %d\n", population->counter);
	printf("\trank : %d\n", population->rank);
	printf("\tcrowding distance : %lf\n", population->crowding_distance);
	printf("\tselection probabiliy : %lf\n", population->sel_prob);

	printf("\tmCAI value : %lf\n", population->sol.obj_val[_mCAI]);
	printf("\tmHD value : %lf\n", population->sol.obj_val[_mHD]);
	printf("\tMLRCS value : %lf\n", population->sol.obj_val[_MLRCS]);

	for (int i = 0; i < num_CDSs; i++) {			// population's CDSs loop
		printf("\nPopulatin's CDS [%d] : \n", i);
		printf("%s\n", population->sol.cds[i]);
	}

	printf("count : %d\n", ++cnt);

	return;
}


int main()
{
	srand(time(NULL));

	/* amino sequence recieve from FASTA format */
	char file_name[20] = "Q5VZP5.fasta.txt";
	char buffer[256];
	char* amino_seq;		// amino sequence comprising a CDS
	int* amino_seq_idx;		// amino sequence corresponding index value to struct 'aa'
	int len_amino_seq;		// length of amino sequence

	/* ---------------------------------------- file processing ------------------------------------------------- */
	// printf("input file name : );
	// scanf_s("%s", &file_name);
	FILE* fp;
	fopen_s(&fp, file_name, "r");
	if (fp == NULL) {
		printf("Opening input file failed at line : %d", __LINE__);
		return EXIT_FAILURE;
	}
	fseek(fp, 0, SEEK_END);
	len_amino_seq = ftell(fp);
	fseek(fp, 0, SEEK_SET);
	fgets(buffer, 256, fp);				// jump over first line
	len_amino_seq -= ftell(fp);
	amino_seq = (char*)malloc(sizeof(char) * len_amino_seq);		// over memory allocation
	if (amino_seq == NULL) {
		printf("Memory allocation failed at line : %d", __LINE__);
		return EXIT_FAILURE;
	}

	int idx = 0;
	char tmp;
	while (!feof(fp)) {
		tmp = fgetc(fp);
		if (tmp != '\n')amino_seq[idx++] = tmp;
	}
	amino_seq[idx] = NULL;
	len_amino_seq = idx - 1;

	fclose(fp);
	
	amino_seq_idx = (int*)malloc(sizeof(int) * len_amino_seq);		// memory alloc
	if (amino_seq_idx == NULL) {
		printf("Memory allocation failed at line : %d", __LINE__);
		return EXIT_FAILURE;
	}
	for (int i = 0; i < len_amino_seq; i++) {
		amino_seq_idx[i] = FindAminoIndex((AA)amino_seq[i]);
	}
	/* -------------------------------------------- end file proess ------------------------------------------------- */
	


	/* user input parameter */
	int max_cycle;				// number of generations
	int colony_size;			// number of solutions in population
	int num_cds;				// number of CDSs
	int limit;					// number of solution is not updated
	double mprob;				// mutation probability

	/* input parameter values */
	printf("input max cycle value : "); scanf_s("%d", &max_cycle);
	if (max_cycle <= 0) {
		printf("input max cycle value > 0\n");
		return EXIT_FAILURE;
	}
	printf("input colony size : "); scanf_s("%d", &colony_size);
	if (colony_size <= 0) {
		printf("input colony size > 0\n");
		return EXIT_FAILURE;
	}
	printf("input number of CDSs : "); scanf_s("%d", &num_cds);
	if (num_cds <= 1) {
		printf("input number of CDSs > 1\n");
		return EXIT_FAILURE;
	}
	printf("input limit value : "); scanf_s("%d", &limit);
	if (limit <= 0) {
		printf("input limit value > 0\n");
		return EXIT_FAILURE;
	}
	printf("input mutation probability (0 ~ 1 value) : "); scanf_s("%lf", &mprob);
	if (mprob < 0 || mprob > 1 ) {
		printf("input mutation probability (0 ~ 1 value) : \n");
		return EXIT_FAILURE;
	}

	/* Population memory allocation */
	Population* pop;
	pop = AllocPopulation(colony_size * 2, num_cds, len_amino_seq);

	/* --------------------------------------------------- initialize Population ------------------------------------------------------- */
	for (int i = 0; i < colony_size - 1; i++) {
		GenSolution(&pop[i], num_cds, amino_seq_idx, len_amino_seq, RANDOM_GEN);
		/* calculate objective function value */
		mCAI(&pop[i], num_cds, amino_seq_idx, len_amino_seq);
		mHD(&pop[i], num_cds, len_amino_seq);
		MLRCS(&pop[i], num_cds, len_amino_seq);
	}
	/* To boost the optimization of solutions witth high CAI values
	   remaing solution is generated by selecting highest adaptation */
	GenSolution(&pop[colony_size - 1], num_cds, amino_seq_idx, len_amino_seq, UPPER_GEN);
	mCAI(&pop[colony_size - 1], num_cds, amino_seq_idx, len_amino_seq);
	mHD(&pop[colony_size - 1], num_cds, len_amino_seq);
	MLRCS(&pop[colony_size - 1], num_cds, len_amino_seq);
	/* -------------------------------------------------------- initialize end ----------------------------------------------------------- */


	bool check;
	Population* new_sol, * sel_sol;
	Population* tmp_sol;		// for Scout bee step
	tmp_sol = AllocPopulation(1, num_cds, len_amino_seq);
	/* ------------------------------------------------------- start cycle to max cycle -------------------------------------------------- */
	for (int i = 0; i < max_cycle; i++)
	{
		/* --------------------------------------- start Employed bees step ----------------------------------------- */
		for (int j = 0; j < colony_size; j++)
		{
			new_sol = Mutation(&pop[j], num_cds, amino_seq, len_amino_seq, mprob);		// Employed Bee search
			/* Calculate Objective Functions */
			mCAI(new_sol, num_cds, amino_seq_idx, len_amino_seq);
			mHD(new_sol, num_cds, len_amino_seq);
			MLRCS(new_sol, num_cds, len_amino_seq);
			/* Pareto Comparision */
			check = ParetoComparison(new_sol, &pop[j]);
			if (check)
				CopyPopulation(new_sol, &pop[j], num_cds, len_amino_seq);
			else
				pop[j].counter += 1;
			FreePopulation(new_sol, 1, num_cds);
		}
		/* ----------------------------------------- end Employed bees step ----------------------------------------- */
		
		//SortbyRankCrowding(pop, colony_size);
		//CalSelectionProb(pop);
		
		/* -------------------------------------- start Onlooker bees step ------------------------------------------ */
		for (int j = colony_size; j < 2 * colony_size; j++)
		{
			sel_sol = &pop[SelectSolution(pop, colony_size)];							// select solution
			new_sol = Mutation(sel_sol, num_cds, amino_seq, len_amino_seq, mprob);		// Onlooker Bee search
			/* Calculate Objective Function */
			mCAI(new_sol, num_cds, amino_seq_idx, len_amino_seq);
			mHD(new_sol, num_cds, len_amino_seq);
			MLRCS(new_sol, num_cds, len_amino_seq);
			/* Pareto Comparison */
			check = ParetoComparison(new_sol, sel_sol);
			if (check)
				CopyPopulation(new_sol, &pop[j], num_cds, len_amino_seq);
			else {
				CopyPopulation(sel_sol, &pop[j], num_cds, len_amino_seq);
				pop[j].counter += 1;
			}
			FreePopulation(new_sol, 1, num_cds);
		}
		/* ------------------------------------------ end Onlooker bees step ----------------------------------------- */
		

		/* ------------------------------------- start Scout bees step ----------------------------------------------- */
		for (int j = 0; j < 2 * colony_size; j++)
		{
			if (pop[j].counter > limit) 
			{
				GenSolution(tmp_sol, num_cds, amino_seq_idx, len_amino_seq, RANDOM_GEN);
				/* Scout Bee search */
				for (int k = 0; k < i; k++) 
				{
					new_sol = Mutation(tmp_sol, num_cds, amino_seq, len_amino_seq, mprob);
					mCAI(new_sol, num_cds, amino_seq_idx, len_amino_seq);
					mHD(new_sol, num_cds, len_amino_seq);
					MLRCS(new_sol, num_cds, len_amino_seq);
					CopyPopulation(new_sol, tmp_sol, num_cds, len_amino_seq);
					FreePopulation(new_sol, 1, num_cds);
				}
				CopyPopulation(tmp_sol, &pop[j], num_cds, len_amino_seq);

				/* Calculate Objective Function */
				mCAI(&pop[j], num_cds, amino_seq_idx, len_amino_seq);
				mHD(&pop[j], num_cds, len_amino_seq);
				MLRCS(&pop[j], num_cds, len_amino_seq);
				pop[j].counter = 0;
			}
		}
		/* ------------------------------------------ end Scout bees step -------------------------------------------- */
		//SortbyRankCrowding(pop, colony_size * 2);
	
	}
	/* ----------------------------------------------------- end max cyelce ---------------------------------------------------------------- */



	// Print 
	for (int i = 0; i < colony_size * 2; i++) {
		PrintPopulation(&pop[i], num_cds);
	}

	/* free memory */
	FreePopulation(pop, colony_size * 2, num_cds);
	FreePopulation(tmp_sol, 1, num_cds);
	free(amino_seq);
	free(amino_seq_idx);
	

	return EXIT_SUCCESS;
}