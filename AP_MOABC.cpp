#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <omp.h>


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
	float adaptation[MAX_CODON];	// codons's adaptation weight
}Aminoacids;
/* 20 kinds of amino acids */
/* adaptation weight is ascending order */
const Aminoacids aa[20] = {
	{A, 4, {"GCG", "GCA", "GCC", "GCU"}, {1854 / 13563.0f, 5296 / 13563.0f, 7223 / 13563.0f, 1.0f}},
	{C, 2, {"UGC", "UGU"}, {1234 / 3052.0f, 1.0f}},
	{D, 2, {"GAC", "GAU"}, {8960 / 12731.0f, 1.0f}},
	{E, 2, {"GAG", "GAA"}, {6172 / 19532.0f, 1.0f}},
	{F, 2, {"UUU", "UUC"},{7773 / 8251.0f, 1.0f}},
	{G, 4, {"GGG", "GGA", "GGC", "GGU"},{1852 / 15694.0f, 2781 / 15694.0f, 3600 / 15694.0f, 1.0f}},
	{H, 2, {"CAC", "CAU"}, {3288 / 4320.0f, 1.0f}},
	{I, 3, {"AUA", "AUC", "AUU"},{3172 / 12071.0f, 8251 / 12071.0f, 1.0f}},
	{K, 2, {"AAA", "AAG"},{12845 / 15169.0f, 1.0f}},
	{L, 6, {"CUC", "CUG", "CUU", "CUA", "UUA", "UUG"}, {1242 / 13329.0f, 2852 / 13329.0f, 3207 / 13329.0f, 4134 / 13329.0f, 8549 / 13329.0f, 1.0f}},
	{M, 1, {"AUG"}, {1.0f}},
	{N, 2, {"AAU", "AAC"}, {8613 / 9875.0f, 1.0f}},
	{P, 4, {"CCG", "CCC", "CCU", "CCA"}, {1064 / 8965.0f, 1656 / 8965.0f, 4575 / 8965.0f, 1.0f}},
	{Q, 2, {"CAG", "CAA"}, {3312 / 10987.0f, 1.0f}},
	{R, 6, {"CGG", "CGA", "CGC", "AGG", "CGU", "AGA"}, {342 / 9784.0f, 489 / 9784.0f, 658 / 9784.0f, 2175 / 9784.0f, 3307 / 9784.0f, 1.0f}},
	{S, 6, {"UCG", "AGC", "AGU", "UCA", "UCC", "UCU"}, {2112 / 10025.0f, 2623 / 10025.0f, 3873 / 10025.0f, 4583 / 10025.0f, 6403 / 10025.0f, 1.0f}},
	{T, 4, {"ACG", "ACA", "ACC", "ACU"}, {1938 / 9812.0f, 5037 / 9812.0f, 6660 / 9812.0f, 1.0f}},
	{V, 4, {"GUA", "GUG", "GUC", "GUU"}, {3249 / 11442.0f, 3700 / 11442.0f, 6911 / 11442.0f, 1.0f}},
	{W, 1, {"UGG"}, {1.0f}},
	{Y, 2, {"UAU", "UAC"}, {5768 / 7114.0f, 1.0f}}
};
/* --------------------------------------------------------- end definition ---------------------------------------------------------------------- */


#define OBJECTIVE_NUM 3					// three objective function
#define _mCAI 0
#define _mHD 1
#define _MLRCS 2
/* --------------------------------------------------------- Population definition ---------------------------------------------------------------- */
typedef struct {
	char* cds;							// CDSs's sequences
	int p, q, l;						// this if for MLRCS starting point and length
	int obj_cdsidx[OBJECTIVE_NUM][2];	// CDS's index correspond to objective function
	float obj_val[OBJECTIVE_NUM];		// objective function value (0 ~ 1) for Pareto Comparsion
}Solution;
typedef struct {
	int cnt;
	int counter;						// checking counter to obsolete this solution
	int rank;							// indicate Pareto front (rank)
	float crowding_distance;			// indicate diversity of solution in same rank
	float fitness;						// fitness = 1 / rank
	float sel_prob;						// selection probability
	Solution sol;
}Population;
/* ----------------------------------------------------------- end definition --------------------------------------------------------------------- */


/* ----------------------------------------------------------- Find index definition --------------------------------------------------------------- */
/* this function find aminoacid index using binary search */
int FindAminoIndex(const AA aminoacid)
{
	int low = 0;
	int high = 19;
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
/* --------------------------------------------------------------- end definition ------------------------------------------------------------------ */


/* ---------------------------------------------------------- Population memroy management --------------------------------------------------------- */
/* this function memory allocation population */
Population* AllocPopulation(int pop_size, int num_cds, int len_amino_seq)
{
	Population* pop;

	pop = (Population*)malloc(sizeof(Population) * pop_size);
	if (pop == NULL) {
		printf("Memory allocation failed at line %d\n", __LINE__);
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < pop_size; i++) {
		pop[i].sol.cds = (char*)malloc(sizeof(char) * num_cds * len_amino_seq * 3);
		if (pop[i].sol.cds == NULL) {
			printf("Memory allocation failed at line %d\n", __LINE__);
			exit(EXIT_FAILURE);
		}
	}

	return pop;
}
/* this function free population memory */
void FreePopulation(Population* pop, int pop_size, int num_cds)
{
	for (int i = 0; i < pop_size; i++) {
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
void GenCDS(char* cds, int num_cds, const int* amino_seq_idx, int len_amino_seq, int type = RANDOM_GEN)
{
	int idx;
	int rand_idx;
	int codon_idx;

	idx = 0;
	switch (type) {
		/* random creation */
	case RANDOM_GEN:
		for (int i = 0; i < num_cds; i++) {
			for (int j = 0; j < len_amino_seq; j++) {
				rand_idx = rand() % aa[amino_seq_idx[j]].num_codons;
				cds[idx++] = aa[amino_seq_idx[j]].codons[rand_idx][0];
				cds[idx++] = aa[amino_seq_idx[j]].codons[rand_idx][1];
				cds[idx++] = aa[amino_seq_idx[j]].codons[rand_idx][2];
			}
		}
		break;
		/* creation maximum CAI values */
	case UPPER_GEN:
		for(int i = 0; i < num_cds; i++){
			for (int j = 0; j < len_amino_seq; j++) {
				codon_idx = aa[amino_seq_idx[j]].num_codons - 1;
				cds[idx++] = aa[amino_seq_idx[j]].codons[codon_idx][0];
				cds[idx++] = aa[amino_seq_idx[j]].codons[codon_idx][1];
				cds[idx++] = aa[amino_seq_idx[j]].codons[codon_idx][2];
			}
		}
		break;
	}

	return;
}
/* this function generate random solution */
void GenSolution(Population* pop, int num_cds, const int* amino_seq_idx, int len_amino_seq, int type = RANDOM_GEN)
{
	/* initial value setting */
	pop->counter = 0;				// new solution counter value is zero
	pop->rank = 0;
	pop->crowding_distance = 0;
	pop->fitness = 0;
	pop->sel_prob = 0;

	/* make random CDSs */
	GenCDS(pop->sol.cds, num_cds, amino_seq_idx, len_amino_seq, type);

	return;
}
/* this function copy population */
void CopyPopulation(const Population* origin, Population* target, int num_cds, int len_amino_seq)
{
	target->counter = origin->counter;
	target->rank = origin->rank;
	target->crowding_distance = origin->crowding_distance;
	target->fitness = origin->fitness;
	target->sel_prob = origin->sel_prob;
	target->cnt = origin->cnt;

	for (int i = 0; i < num_cds * len_amino_seq * 3; i++) {
		target->sol.cds[i] = origin->sol.cds[i];
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
void ChageSynonymousCodon(int amino_idx, char* cds, int cd_idx, int type)
{
	char codon[3];
	int idx;
	int rand_idx;
	
	codon[0] = cds[cd_idx];
	codon[1] = cds[cd_idx + 1];
	codon[2] = cds[cd_idx + 2];
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
			rand_idx = idx + rand() % (aa[amino_idx].num_codons - 1 - idx) + 1;
			cds[cd_idx] = aa[amino_idx].codons[rand_idx][0];
			cds[cd_idx + 1] = aa[amino_idx].codons[rand_idx][1];
			cds[cd_idx + 2] = aa[amino_idx].codons[rand_idx][2];
		}
		break;
	}

	return;
}
/* this function generate random mutated solution */
Population* Mutation(const Population* pop, int num_cds, const int* amino_seq_idx, int len_amino_seq, float mprob)
{
	/* new population memory allocation */
	Population* new_pop;
	new_pop = AllocPopulation(1, num_cds, len_amino_seq);

	// copy population to new_population
	CopyPopulation(pop, new_pop, num_cds, len_amino_seq);
	new_pop->counter = 0;			// mutated population counter value is zero


	/* generate (0 ~ 1) random number corresponding to codon in CDS */
	float* random_num;
	random_num = (float*)malloc(sizeof(float) * num_cds * len_amino_seq);
	for (int i = 0; i < num_cds * len_amino_seq; i++) {
		random_num[i] = (float) rand() / RAND_MAX;
	}


	int type;
	int len_cds;
	type = rand() % 4;
	len_cds = 3 * len_amino_seq;
	/* four type of variations */
	switch (type)
	{
		/* change each codons with random codon in all CDSs */
	case 0:
		for (int i = 0; i < num_cds; i++) {
			for (int j = 0; j < len_amino_seq; j++) {
				if (random_num[i * len_amino_seq + j] < mprob) {
					ChageSynonymousCodon(amino_seq_idx[j], new_pop->sol.cds, i * len_cds + j * 3, RANDOM_ADAPTATION);
				}
			}
		}
		break;
		/* change each codons with higher adaptation value in CDS which have minimun CAI value */
	case 1:
		for (int i = 0; i < len_amino_seq; i++) {
			if (random_num[new_pop->sol.obj_cdsidx[_mCAI][0] * len_amino_seq + i] < mprob) {
				ChageSynonymousCodon(amino_seq_idx[i], new_pop->sol.cds, new_pop->sol.obj_cdsidx[_mCAI][0] * len_cds + i * 3, UPPER_ADAPTATION);
			}
		}
		break;
		/* change each codons with random codon in pair of CDSs with minimun Hamming Distance */
	case 2:
		for (int i = 0; i < len_amino_seq; i++) {
			if (random_num[new_pop->sol.obj_cdsidx[_mHD][0] * len_amino_seq + i] < mprob) {
				ChageSynonymousCodon(amino_seq_idx[i], new_pop->sol.cds, new_pop->sol.obj_cdsidx[_mHD][0] * len_cds + i * 3, RANDOM_ADAPTATION);
			}
			if (random_num[new_pop->sol.obj_cdsidx[_mHD][1] * len_amino_seq + i] < mprob) {
				ChageSynonymousCodon(amino_seq_idx[i], new_pop->sol.cds, new_pop->sol.obj_cdsidx[_mHD][1] * len_cds + i * 3, RANDOM_ADAPTATION);
			}
		}
		break;
		/* chage each codons with random codon in CDSs with longest common substring */
	case 3:
		for (int i = (new_pop->sol.p - new_pop->sol.obj_cdsidx[_MLRCS][0] * len_cds) / 3; i < (new_pop->sol.p + new_pop->sol.l - 1 - new_pop->sol.obj_cdsidx[_MLRCS][0] * len_cds) / 3; i++) {
			if (random_num[new_pop->sol.obj_cdsidx[_MLRCS][0] * len_amino_seq + i] < mprob) {
				ChageSynonymousCodon(amino_seq_idx[i], new_pop->sol.cds, new_pop->sol.obj_cdsidx[_MLRCS][0] * len_cds + i * 3, RANDOM_ADAPTATION);
			}
		}
		for (int i = (new_pop->sol.q - new_pop->sol.obj_cdsidx[_MLRCS][1] * len_cds) / 3; i < (new_pop->sol.q + new_pop->sol.l - 1 - new_pop->sol.obj_cdsidx[_MLRCS][1] * len_cds) / 3; i++) {
			if (random_num[new_pop->sol.obj_cdsidx[_MLRCS][1] * len_amino_seq + i] < mprob) {
				ChageSynonymousCodon(amino_seq_idx[i], new_pop->sol.cds, new_pop->sol.obj_cdsidx[_MLRCS][1] * len_cds + i * 3, RANDOM_ADAPTATION);
			}
		}
		break;
	}

	/* free memory */
	free(random_num);


	return new_pop;
}
/* -------------------------------------------------------- end Population creation and mutaion ----------------------------------------------------- */



/* --------------------------------------------------------- calculate objective function value ----------------------------------------------------- */
/* this function is calculate mininum CAI value and store cds index */
void mCAI(Population* pop, int num_cds, const int* amino_seq_idx, int len_amino_seq)
{
	char codon[3];
	int idx;
	int codon_idx;
	float tmp;

	idx = 0;
	pop->sol.obj_val[_mCAI] = 1;
	for (int i = 0; i < num_cds; i++) {
		tmp = 1;
		for (int j = 0; j < len_amino_seq; j++) {
			codon[0] = pop->sol.cds[idx++];
			codon[1] = pop->sol.cds[idx++];
			codon[2] = pop->sol.cds[idx++];
			codon_idx = FindCodonIndex(amino_seq_idx[j], codon);
			tmp *= pow(aa[amino_seq_idx[j]].adaptation[codon_idx], 1.0 / len_amino_seq);
		}
		if (tmp <= pop->sol.obj_val[_mCAI]) {
			pop->sol.obj_val[_mCAI] = tmp;
			pop->sol.obj_cdsidx[_mCAI][0] = i;			// CDS's index having mCAI value
		}
	}

	return;
}
/* this function is calculate minimum Hamming Distance */
void mHD(Population* pop, int num_cds, int len_amino_seq)
{
	int len_cds = len_amino_seq * 3;
	int cnt;
	float tmp;

	pop->sol.obj_val[_mHD] = 1;
	for (int i = 0; i < num_cds - 1; i++)
	{
		for (int j = i + 1; j < num_cds; j++)
		{
			cnt = 0;
			for (int k = 0; k < len_cds; k++)
			{
				if (pop->sol.cds[i * len_cds + k] != pop->sol.cds[j * len_cds + k]) 
					cnt++;
			}
			tmp = (float)cnt / len_cds;
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
void MLRCS(Population* pop, int num_cds, int len_amino_seq)
{
	int** LCS;						// matrix which size [length of CDS + 1][length of CDS + 1]
	int len_cds = 3 * len_amino_seq;
	int max_len;


	/* memory allocation for LCS matrix */
	LCS = (int**)malloc(sizeof(int*) * (len_cds + 1));
	for (int i = 0; i < len_cds + 1; i++) {
		LCS[i] = (int*)malloc(sizeof(int) * (len_cds + 1));
	}

	max_len = 0;
	for (int i = 0; i < num_cds; i++) {							// i'th CDS
		for (int j = i; j < num_cds; j++) {						// j'th CDS
			for (int k = 0; k < len_cds + 1; k++) {				// matrix row
				for (int l = 0; l < len_cds + 1; l++) {			// matrix column
					if (i != j) {
						if (k == 0 || l == 0)
							LCS[k][l] = 0;
						else if (pop->sol.cds[i * len_cds + k - 1] == pop->sol.cds[j * len_cds + l - 1])
						{
							LCS[k][l] = LCS[k - 1][l - 1] + 1;
							if (LCS[k][l] >= max_len)
							{
								max_len = LCS[k][l];
								pop->sol.p = i * len_cds + k - max_len;
								pop->sol.q = j * len_cds + l - max_len;
								pop->sol.l = max_len;
								pop->sol.obj_cdsidx[_MLRCS][0] = i;
								pop->sol.obj_cdsidx[_MLRCS][1] = j;
							}
						}
						else
							LCS[k][l] = 0;
					}
					else
					{
						if (k == 0 || l == 0 || (k == l))
							LCS[k][l] = 0;
						else if (pop->sol.cds[i * len_cds + k - 1] == pop->sol.cds[j * len_cds + l - 1])
						{
							LCS[k][l] = LCS[k - 1][l - 1] + 1;
							if (LCS[k][l] >= max_len)
							{
								max_len = LCS[k][l];
								pop->sol.p = i * len_cds + k - max_len;
								pop->sol.q = j * len_cds + l - max_len;
								pop->sol.l = max_len;
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
	pop->sol.obj_val[_MLRCS] = (float)max_len / len_cds;

	/* free memory LCS matrix */
	for (int i = 0; i < len_cds + 1; i++) {
		free(LCS[i]);
	}
	free(LCS);

	return;
}
/* If new population dominate old population return true */
bool ParetoComparison(const Population* new_pop, const Population* old_pop)
{
	if ((new_pop->sol.obj_val[_mCAI] == old_pop->sol.obj_val[_mCAI]) &&
		(new_pop->sol.obj_val[_mHD] == old_pop->sol.obj_val[_mHD]) &&
		(new_pop->sol.obj_val[_MLRCS] == old_pop->sol.obj_val[_MLRCS]))
		return false;
	else if ((new_pop->sol.obj_val[_mCAI] >= old_pop->sol.obj_val[_mCAI]) &&
		(new_pop->sol.obj_val[_mHD] >= old_pop->sol.obj_val[_mHD]) &&
		(new_pop->sol.obj_val[_MLRCS] <= old_pop->sol.obj_val[_MLRCS]))
		return true;
	else
		return false;
}
/* ---------------------------------------------------------- end objective function caculation ------------------------------------------------------ */




#define EMPTY -1
/* this function sorting by rank and crowding distance */
void SortbyRankCrowding(Population* pop, int pop_size, int num_cds, int len_amino_seq)
{
	/* this point out population index value */
	int** Sp, ** F;
	int* np, * Q;
	int Sp_idx, F_front, F_idx, Q_idx;

	/* memory allocation */
	Sp = (int**)malloc(sizeof(int*) * pop_size);
	F = (int**)malloc(sizeof(int*) * pop_size);
	for (int i = 0; i < pop_size; i++) {
		Sp[i] = (int*)malloc(sizeof(int) * pop_size);
		F[i] = (int*)malloc(sizeof(int) * pop_size);
	}
	np = (int*)malloc(sizeof(int) * pop_size);
	Q = (int*)malloc(sizeof(int) * pop_size);


	// F empty initialization
	for (int i = 0; i < pop_size; i++) {
		memset(F[i], EMPTY, sizeof(int) * pop_size);
	}


	/* ------------------------------------------------- fast non-dominated sort -------------------------------------------- */
	F_idx = 0;
	for (int i = 0; i < pop_size; i++)							
	{
		Sp_idx = 0;
		np[i] = 0;
		memset(Sp[i], EMPTY, sizeof(int) * pop_size);
		for (int j = 0; j < pop_size; j++)
		{
			if (i != j) {
				if (ParetoComparison(&pop[i], &pop[j]))
					Sp[i][Sp_idx++] = j;
				else if (ParetoComparison(&pop[j], &pop[i]))
					np[i] += 1;
			}
		}
		if (np[i] == 0) {
			pop[i].rank = 1;
			F[0][F_idx++] = i;									
		}
	}
	/* -------------- 1st front setting complete ---------------- */

	int p;						// indicate which p solution's Sp		p value means n'th Pareto Front solution
	F_front = 0;				// indicate 1st front
	F_idx = 0;
	while (F[F_front][F_idx] != EMPTY)
	{
		Q_idx = 0;
		memset(Q, EMPTY, sizeof(int) * pop_size);
		for (F_idx = 0; F[F_front][F_idx] != EMPTY && F_idx < pop_size; F_idx++) {					// F[F_front][F_idx] == population index
			p = F[F_front][F_idx];
			for (Sp_idx = 0; Sp[p][Sp_idx] != EMPTY && Sp_idx < pop_size; Sp_idx++) {
				np[Sp[p][Sp_idx]]--;
				if (np[Sp[p][Sp_idx]] == 0)
				{
					pop[Sp[p][Sp_idx]].rank = F_front + 2;
					Q[Q_idx++] = Sp[p][Sp_idx];
				}
			}
		}

		F_front++;
		if (F_front == pop_size)
			break;
		F_idx = 0;
		Q_idx = 0;
		while (Q[Q_idx] != EMPTY) {
			F[F_front][Q_idx] = Q[Q_idx];
			Q_idx++;
			if (Q_idx >= pop_size)
				break;
		}
	}
	/* ------------------------------------------------- end non dominated sort -------------------------------------------------- */


	/* ------------------------------------------------- crowding distance assignment -------------------------------------------- */
	int p_len;
	int tmp;
	F_front = 0;
	while (F[F_front][0] != EMPTY)
	{
		p_len = 0;							// number of solutions in pareto rank
		for (F_idx = 0; F[F_front][F_idx] != EMPTY && F_idx < pop_size; F_idx++) {
			pop[F[F_front][F_idx]].crowding_distance = 0;
			p_len++;
		}
		for (int i = 0; i < OBJECTIVE_NUM; i++)
		{
			// sort each objective function ascending order
			for (int j = 0; j < p_len; j++) {
				for (int k = 0; k < p_len - 1 - j; k++) {
					if (pop[F[F_front][k]].sol.obj_val[i] > pop[F[F_front][k + 1]].sol.obj_val[i]) {
						tmp = F[F_front][k];
						F[F_front][k] = F[F_front][k + 1];
						F[F_front][k + 1] = tmp;
					}
				}
			}

			pop[F[F_front][0]].crowding_distance = 0x7fffffff;
			pop[F[F_front][p_len - 1]].crowding_distance = 0x7fffffff;

			for (int j = 1; j < p_len - 1; j++)
				pop[F[F_front][j]].crowding_distance += (pop[F[F_front][j + 1]].sol.obj_val[i] - pop[F[F_front][j - 1]].sol.obj_val[i]) / (1 - 0);
		}

		for (int j = 0; j < p_len; j++) {
			for (int k = 0; k < p_len - 1 - j; k++) {
				if (pop[F[F_front][k]].crowding_distance < pop[F[F_front][k + 1]].crowding_distance) {
					tmp = F[F_front][k];
					F[F_front][k] = F[F_front][k + 1];
					F[F_front][k + 1] = tmp;
				}
			}
		}

		F_front++;
		if (F_front >= pop_size)
			break;
	}
	/* ------------------------------------------------ end crowding distance assignment -------------------------------------------- */

	Population* tmp_pop;
	tmp_pop = AllocPopulation(pop_size, num_cds, len_amino_seq);
	for (int i = 0; i < pop_size; i++) {
		CopyPopulation(&pop[i], &tmp_pop[i], num_cds, len_amino_seq);
	}
	F_front = 0;
	F_idx = 0;
	int p_idx = 0;
	while (F[F_front][F_idx] != EMPTY && F_front < pop_size)
	{
		CopyPopulation(&tmp_pop[F[F_front][F_idx]], &pop[p_idx], num_cds, len_amino_seq);
		p_idx++;
		F_idx++;
		if (F[F_front][F_idx] == EMPTY || F_idx == pop_size) {
			F_front++;
			F_idx = 0;
		}
	}


	// allocate fitness value 
	for (int i = 0; i < pop_size; i++) {
		pop[i].fitness = 1.f / pop[i].rank;
	}


	/* free memory */
	FreePopulation(tmp_pop, pop_size, num_cds);
	for (int i = 0; i < pop_size; i++) {
		free(Sp[i]);
		free(F[i]);
	}
	free(Sp);
	free(F);
	free(np);
	free(Q);

	return;
}


/* this function caculate selection probability */
void CalSelectionProb(Population* pop, int pop_size)
{
	float sum;

	sum = 0;
	for (int i = 0; i < pop_size; i++) {
		sum += pop[i].fitness;
	}
	for (int i = 0; i < pop_size; i++) {
		pop[i].sel_prob = pop[i].fitness / sum;
	}

	return;
}
/* Roulette-wheel selection based on selection probability */
int SelectSolution(const Population* pop, int pop_size)
{
	float sum;
	float point;

	point = (float) rand() / RAND_MAX;		// 0 ~ 1 floating point number

	sum = 0;
	for (int i = 0; i < pop_size; i++) {
		sum += pop[i].sel_prob;
		if (point < sum)
			return i;						// return selected pop index
	}
}



void PrintPopulation(const Population* population, int num_cds, int len_amino_seq);
//void PrintAminoAcids();
//void CompareCdsToAminoAcids(const char* cds, int num_cds, const int* amino_seq_idx, const char* amino_seq, int len_amino_seq);
//void CheckMLRCS(const char* s, int size);
//void CheckMutation(const Population* pop1, const Population* pop2, int num_cds, const int* amino_seq_idx, const char* amino_seq, int len_amino_seq);

/* ---------------------------- For Population queue function -------------------------------------- */
typedef struct Sol
{
	int pos;
	Population* pop;
}Sol;

typedef struct Node
{
	Sol s;
	struct Node* next;
}Node;

typedef struct Queue
{
	Node* front;
	Node* rear;
	int cnt;
}Queue;

void initQueue(Queue* queue)
{
	queue->front = queue->rear = NULL;
	queue->cnt = 0;
}

int isEmpty(Queue* queue)		// if queue is empty return true
{
	return queue->cnt == 0;
}

void enqueue(Queue* queue, Population* pop, int position)
{
	Node* newNode = (Node*)malloc(sizeof(Node)); 
	newNode->s.pop = pop;
	newNode->s.pos = position;
	newNode->next = NULL;

	if (isEmpty(queue))    
		queue->front = newNode;
	else   
		queue->rear->next = newNode;    
	queue->rear = newNode;       
	queue->cnt++;    
}

Sol dequeue(Queue* queue)
{
	Sol res;
	Node* ptr;
	if (isEmpty(queue)){
		printf("Error : Queue is empty!\n");
		//return;
	}
	ptr = queue->front; 
	
	res.pop = ptr->s.pop;
	res.pos = ptr->s.pos;
	queue->front = ptr->next;    
	free(ptr);     
	queue->cnt--;    

	return res;
}
/* --------------------------------- queue function end --------------------------------------- */


/* global variable initialize */
int num_threads = 16;
bool stop = false;



/* ------------------------------- Master thread function definition ------------------------------------*/
void MasterTask(Population* pop, Population* sw_pop, Queue* queue, int colony_size, int max_eval, int num_cds, int len_amino_seq)
{
	int eval;
	bool update;
	Sol s;
	Population* tmp;
	

	eval = 0;
	while (eval < max_eval)
	{
		update = false;
		for (int i = 0; i < num_threads - 1; i++) {
			while (isEmpty(&queue[i]) == false && eval < max_eval) {
				s = dequeue(&queue[i]);
				CopyPopulation(s.pop, &sw_pop[s.pos], num_cds, len_amino_seq);
				FreePopulation(s.pop, 1, num_cds);
				update = true;
				// update evaluation
				eval++;
				sw_pop[s.pos].cnt++;
			}
		}
		if (update) {
			SortbyRankCrowding(sw_pop, colony_size * 2, num_cds, len_amino_seq);
			CalSelectionProb(sw_pop, colony_size);
			tmp = pop;
			pop = sw_pop;
			sw_pop = tmp;
		}
	}
	// send stop signal to worker threads (num_threads - 1)
	stop = true;

	return;
}
/* ------------------------------------ Master thread end definition ------------------------------------*/

/* ------------------------------- Worker thread function definition ------------------------------------*/
void WorkerTask(Population* pop, Queue* queue, int colony_size, int limit, float mprob, int tid, int num_cds, int* amino_seq_idx, int len_amino_seq)
{
	int start, end, pos;
	int em_thread, on_thread;
	bool check;
	Population* new_sol, * sel_sol, * tmp_sol;
	Population* solution;

	tmp_sol = AllocPopulation(1, num_cds, len_amino_seq);

	em_thread = num_threads / 2;
	on_thread = num_threads - em_thread - 1;

	/* work distribution check compelete */
	if (tid < em_thread) {
		start = (colony_size / em_thread) * tid + ((colony_size % em_thread) <= tid ? (colony_size % em_thread) : tid);
		end = start + (colony_size / em_thread) + ((colony_size % em_thread) <= tid ? 0 : 1) - 1;
	}
	else {
		start = colony_size + (colony_size / on_thread) * (tid - em_thread) + ((colony_size % on_thread) <= (tid - em_thread) ? (colony_size % on_thread) : (tid - em_thread));
		end = start + (colony_size / on_thread) + ((colony_size % on_thread) <= (tid - em_thread) ? 0 : 1) - 1;
	}
	pos = start;

	while (stop == false)
	{
		solution = AllocPopulation(1, num_cds, len_amino_seq);
		if (tid < em_thread)			// Perform Employed Bee Processing
		{
			new_sol = Mutation(&pop[pos], num_cds, amino_seq_idx, len_amino_seq, mprob);		
			/* Calculate Objective Functions */
			mCAI(new_sol, num_cds, amino_seq_idx, len_amino_seq);
			mHD(new_sol, num_cds, len_amino_seq);
			MLRCS(new_sol, num_cds, len_amino_seq);
			/* Pareto Comparision */
			check = ParetoComparison(new_sol, &pop[pos]);
			if (check) {
				new_sol->cnt = pop[pos].cnt;
				CopyPopulation(new_sol, solution, num_cds, len_amino_seq);
			}
			else {
				CopyPopulation(&pop[pos], solution, num_cds, len_amino_seq);
				solution->counter++;
			}
			FreePopulation(new_sol, 1, num_cds);
		}
		else								// Perform Onlooker Bee Processing
		{
			sel_sol = &pop[SelectSolution(pop, colony_size)];								// select solution
			new_sol = Mutation(sel_sol, num_cds, amino_seq_idx, len_amino_seq, mprob);		// Onlooker Bee search
			/* Calculate Objective Function */
			mCAI(new_sol, num_cds, amino_seq_idx, len_amino_seq);
			mHD(new_sol, num_cds, len_amino_seq);
			MLRCS(new_sol, num_cds, len_amino_seq);
			/* Pareto Comparison */
			check = ParetoComparison(new_sol, sel_sol);
			if (check) {
				new_sol->cnt = sel_sol->cnt;
				CopyPopulation(new_sol, solution, num_cds, len_amino_seq);
			}
			else {
				CopyPopulation(sel_sol, solution, num_cds, len_amino_seq);
				solution->counter++;
			}
			FreePopulation(new_sol, 1, num_cds);
		}

		if (solution->counter > limit)
		{
			// Perform Scout Bee Processing
			GenSolution(tmp_sol, num_cds, amino_seq_idx, len_amino_seq, RANDOM_GEN);
			mCAI(tmp_sol, num_cds, amino_seq_idx, len_amino_seq);
			mHD(tmp_sol, num_cds, len_amino_seq);
			MLRCS(tmp_sol, num_cds, len_amino_seq);
			/* Scout Bee search */
			for (int j = 0; j < solution->cnt; j++)
			{
				new_sol = Mutation(tmp_sol, num_cds, amino_seq_idx, len_amino_seq, mprob);
				mCAI(new_sol, num_cds, amino_seq_idx, len_amino_seq);
				mHD(new_sol, num_cds, len_amino_seq);
				MLRCS(new_sol, num_cds, len_amino_seq);
				CopyPopulation(new_sol, tmp_sol, num_cds, len_amino_seq);
				FreePopulation(new_sol, 1, num_cds);
			}
			tmp_sol->cnt = solution->cnt;
			CopyPopulation(tmp_sol, solution, num_cds, len_amino_seq);
		}
		// queue »ðÀÔ
		printf("pos : %d \n", pos);
		enqueue(&queue[tid], solution, pos);
		pos++;
		if (pos == end + 1)
			pos = start;
	}
	FreePopulation(tmp_sol, 1, num_cds);
	// non-dominated file update ÇÊ¿ä


	return;
}
/* ------------------------------------ Worker thread end definition ------------------------------------*/



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
	float mprob;				// mutation probability

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
	printf("input mutation probability (0 ~ 1 value) : "); scanf_s("%f", &mprob);
	if (mprob < 0 || mprob > 1) {
		printf("input mutation probability (0 ~ 1 value) : \n");
		return EXIT_FAILURE;
	}


	/* Population memory allocation */
	Population* pop, *sw_pop;
	pop = AllocPopulation(colony_size * 2, num_cds, len_amino_seq);
	sw_pop = AllocPopulation(colony_size * 2, num_cds, len_amino_seq);
	
	// qeueu initialize
	Queue* queue;
	queue = (Queue*)malloc(sizeof(Queue) * num_threads);
	for (int i = 0; i < num_threads - 1; i++) {
		initQueue(&queue[i]);
	}

	int tid;
	omp_set_num_threads(num_threads);
#pragma omp parallel private(tid)
	{
		/* --------------------------------------------------- initialize Population ------------------------------------------------------- */
#pragma omp for
		for (int i = 0; i < colony_size; i++)
		{
			if (i == colony_size - 1)
			{
				GenSolution(&pop[i], num_cds, amino_seq_idx, len_amino_seq, UPPER_GEN);
				mCAI(&pop[i], num_cds, amino_seq_idx, len_amino_seq);
				mHD(&pop[i], num_cds, len_amino_seq);
				MLRCS(&pop[i], num_cds, len_amino_seq);
			}
			else {
				GenSolution(&pop[i], num_cds, amino_seq_idx, len_amino_seq, RANDOM_GEN);
				/* calculate objective function value */
				mCAI(&pop[i], num_cds, amino_seq_idx, len_amino_seq);
				mHD(&pop[i], num_cds, len_amino_seq);
				MLRCS(&pop[i], num_cds, len_amino_seq);
			}
			pop[i].cnt = 0;
		}
#pragma omp for
		for (int i = colony_size; i < colony_size * 2; i++) {
			pop[i].sol.obj_val[_mCAI] = 0;
			pop[i].sol.obj_val[_mHD] = 0;
			pop[i].sol.obj_val[_MLRCS] = 1;
		}
#pragma omp single
		{
		SortbyRankCrowding(pop, colony_size, num_cds, len_amino_seq);
		CalSelectionProb(pop, colony_size); 
		}
#pragma omp for
		for (int i = 0; i < colony_size; i++) {
			CopyPopulation(&pop[i], &sw_pop[i], num_cds, len_amino_seq);
		}
		/* -------------------------------------------------------- initialize end ----------------------------------------------------------- */

		tid = omp_get_thread_num();
		/* Master thread */
		if (tid == num_threads - 1)
			MasterTask(pop, sw_pop, queue, colony_size, colony_size * max_cycle, num_cds, len_amino_seq);
		/* Worker thread */
		else
			WorkerTask(pop, queue, colony_size, limit, mprob, tid, num_cds, amino_seq_idx, len_amino_seq);
	}

	// Print 
	for (int i = 0; i < colony_size * 2; i++) {
		PrintPopulation(&pop[i], num_cds, len_amino_seq);
	}
	
	/* free memory */
	FreePopulation(pop, colony_size * 2, num_cds);
	FreePopulation(sw_pop, colony_size * 2, num_cds);
	free(amino_seq);
	free(amino_seq_idx);
	free(queue);

	return EXIT_SUCCESS;
}



/* ---------------------------------- For function test -------------------------------------  */
/* print population's attribute values */
void PrintPopulation(const Population* population, int num_cds, int len_amino_seq)
{
	int idx;

	printf("\n ------------------------------ Print Population ------------------------------ \n ");
	printf("\trank : %d\n", population->rank);
	printf("\tcrowding distance : %f\n", population->crowding_distance);
	printf("\tfitness : %f\n", population->fitness);
	printf("\tcounter : %d\n", population->counter);
	printf("\tselection probabiliy : %f\n", population->sel_prob);

	printf("\tmCAI value : %f\n", population->sol.obj_val[_mCAI]);
	printf("\tmHD value : %f\n", population->sol.obj_val[_mHD]);
	printf("\tMLRCS value : %f\n", population->sol.obj_val[_MLRCS]);

	idx = 0;
	for (int i = 0; i < num_cds; i++) {					// population's CDSs loop
		printf("\nPopulatin's CDS [%d] : \n", i);
		for (int j = 0; j < len_amino_seq * 3; j++) {
			printf("%c", population->sol.cds[i * len_amino_seq * 3 + j]);
		}
	}

	return;
}
/* print Aminoacids definition check */
void PrintAminoAcids()
{
	char file_name[20] = "amino.txt";
	FILE* fp;

	fopen_s(&fp, file_name, "w");
	if (fp == NULL) {
		printf("%s open failure", file_name);
		return;
	}
	
	fprintf(fp, "Aminoacids codon frequency used in out study\n\n");
	for (int i = 0; i < 20; i++) {
		fprintf(fp, "-------------- Aminoacid : %c -----------\n", aa[i].name);
		for (int j = 0; j < aa[i].num_codons; j++) {
			fprintf(fp, "codon[%d] : %s\tadaptation : %f\n", j + 1, aa[i].codons[j], aa[i].adaptation[j]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);

	/* console window printing */
	/*printf("Aminoacids codon frequency used in out study\n");
	for (int i = 0; i < 20; i++) {
		printf("-------------- Aminoacid : %c ---------\n", aa[i].name);
		for (int j = 0; j < aa[i].num_codons; j++) {
			printf("codon[%d] : %s\tadaptation : %f\n", j + 1, aa[i].codons[j], aa[i].adaptation[j]);
		}
		printf("\n");
	}*/

	return;
}
/* cds to amino sequence and comparison */
void CompareCdsToAminoAcids(const char* cds, int num_cds, const int * amino_seq_idx, const char* amino_seq, int len_amino_seq)
{
	char codon[3];
	char* ch_amino_seq;
	char* idx_amino_seq;
	int idx;
	int c_idx;

	// memory allocation
	ch_amino_seq = (char*)malloc(sizeof(char) * num_cds * len_amino_seq * 3);
	idx_amino_seq = (char*)malloc(sizeof(char) * len_amino_seq * 3);


	/* -------------------------------------- CDS to Amino ----------------------------------------- */
	c_idx = 0;
	for (int i = 0; i < num_cds; i++) {
		for (int j = 0; j < len_amino_seq; j++) {
			codon[0] = cds[i * len_amino_seq * 3 + j * 3];
			codon[1] = cds[i * len_amino_seq * 3 + j * 3 + 1];
			codon[2] = cds[i * len_amino_seq * 3 + j * 3 + 2];
			for (int k = 0; k < 20; k++) {
				for (int l = 0; l < aa[k].num_codons; l++) {
					if (codon[0] == aa[k].codons[l][0] &&
						codon[1] == aa[k].codons[l][1] &&
						codon[2] == aa[k].codons[l][2]) {
						ch_amino_seq[c_idx++] = (char)aa[k].name;
						break;
					}
				}
			}
		}
	}

	/*for (int i = 0; i < len_amino_seq; i++){
		printf("%c", ch_amino_seq[i]);
	}*/

	printf("\nCDS to aminoacids seqeunces comparison ....\n");
	c_idx = 0;
	for (int i = 0; i < num_cds; i++) {
		idx = 0;
		for (int j = 0; j < len_amino_seq; j++) {
			if (ch_amino_seq[c_idx++] != amino_seq[idx++]) {
				printf("Warnings : amino acids sequences are different\n");
				return;
			}
		}
	}
	printf("\nCDS to amino seqeunce comparison is compelete !!\n");
	/* ------------------------------------------------------------------------------------------- */


	/* ----------------------------------- amino index to amino ---------------------------------- */
	idx = 0;
	for (int i = 0; i < len_amino_seq; i++) {
		idx_amino_seq[i] = (char)aa[amino_seq_idx[idx++]].name;
	}

	printf("\nAmino indicies, amino seqeunce compare comparison ....\n");
	for (int i = 0; i < len_amino_seq; i++) {
		if (idx_amino_seq[i] != amino_seq[i]) {
			printf("amino sequences's indices are different\n");
			return;
		}
	}
	printf("Amino indicies, amino seqeunce compare is compelete !!\n");
	/* ------------------------------------------------------------------------------------------- */



	// free memory
	free(ch_amino_seq);
	free(idx_amino_seq);

	return;
}

void CheckMLRCS(const char* s, int size)
{
	int p, q, l;
	int** LCS;
	int max;

	LCS = (int**)malloc(sizeof(int*) * (size + 1));
	for (int i = 0; i < size + 1; i++) {
		LCS[i] = (int*)malloc(sizeof(int) * (size + 1));
	}

	max = 0;
	for (int i = 0; i < size + 1; i++) {
		for (int j = 0; j < size + 1; j++) {
			if (i == 0 || j == 0 || (i == j)) {
				LCS[i][j] = 0;
			}
			else if (s[i - 1] == s[j - 1]) {
				LCS[i][j] = LCS[i - 1][j - 1] + 1;
				if (LCS[i][j] >= max) {
					max = LCS[i][j];
					l = max;
					p = i - max;
					q = j - max;
				}
			}
			else
				LCS[i][j] = 0;
		}
	}

	for (int i = 0; i < size + 1; i++) {
		for (int j = 0; j < size + 1; j++) {
			printf("%d ", LCS[i][j]);
		}
		printf("\n");
	}

	printf("p : %d\n", p);
	printf("q : %d\n", q);
	printf("l : %d\n", l);

	for (int i = 0; i < size + 1; i++) {
		free(LCS[i]);
	}
	free(LCS);

	return;
}

/* original population and muatated population comparison */
void CheckMutation(const Population* pop1, const Population* pop2, int num_cds, const int * amino_seq_idx, const char* amino_seq, int len_amino_seq)
{
	int len_cds;

	printf("CompareCdsToAminoAcids pop1 and pop2 ...\n");
	CompareCdsToAminoAcids(pop1->sol.cds, num_cds, amino_seq_idx, amino_seq, len_amino_seq);
	CompareCdsToAminoAcids(pop2->sol.cds, num_cds, amino_seq_idx, amino_seq, len_amino_seq);

	len_cds = len_amino_seq * 3;
	printf("Check which amino sequneces index are changed ...");
	for (int i = 0; i < num_cds; i++) {
		for (int j = 0; j < len_amino_seq; j++) {
			if (pop1->sol.cds[i * len_cds + j * 3] != pop2->sol.cds[i * len_cds + j * 3] ||
				pop1->sol.cds[i * len_cds + j * 3 + 1] != pop2->sol.cds[i * len_cds + j * 3 + 1] ||
				pop1->sol.cds[i * len_cds + j * 3 + 2] != pop2->sol.cds[i * len_cds + j * 3 + 2]) {
				printf("index : %d is changed\n", i * len_cds + j * 3);
			}
		}
	}
}