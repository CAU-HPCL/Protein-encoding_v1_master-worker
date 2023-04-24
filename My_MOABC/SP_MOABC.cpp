#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <random>
#include <chrono>

// random generator
std::random_device rd;
std::knuth_b gen(rd());
std::uniform_real_distribution<> dist(0, 1); // double

#define _CRT_SECURE_NO_WARNINGS

#define MAX_CODON 6
/* ------------------------------------------------- Amino acids, codons and adaptation weigth definition ---------------------------------------- */
enum AA
{
	A = 'A',
	C = 'C',
	D = 'D',
	E = 'E',
	F = 'F',
	G = 'G',
	H = 'H',
	I = 'I',
	K = 'K',
	L = 'L',
	M = 'M',
	N = 'N',
	P = 'P',
	Q = 'Q',
	R = 'R',
	S = 'S',
	T = 'T',
	V = 'V',
	W = 'W',
	Y = 'Y'
};

typedef struct
{
	AA name;					 // amino acid abbreviation
	int num_codons;				 // number of codon types in a amino acid
	char codons[MAX_CODON][4];	 // codons in a amino acid
	float adaptation[MAX_CODON]; // codons's adaptation weight
} Aminoacids;
/* 20 kinds of amino acids */
/* adaptation weight is ascending order */
const Aminoacids aa[20] = {
	{A, 4, {"GCG", "GCA", "GCC", "GCU"}, {1854 / 13563.0f, 5296 / 13563.0f, 7223 / 13563.0f, 1.0f}},
	{C, 2, {"UGC", "UGU"}, {1234 / 3052.0f, 1.0f}},
	{D, 2, {"GAC", "GAU"}, {8960 / 12731.0f, 1.0f}},
	{E, 2, {"GAG", "GAA"}, {6172 / 19532.0f, 1.0f}},
	{F, 2, {"UUU", "UUC"}, {7773 / 8251.0f, 1.0f}},
	{G, 4, {"GGG", "GGA", "GGC", "GGU"}, {1852 / 15694.0f, 2781 / 15694.0f, 3600 / 15694.0f, 1.0f}},
	{H, 2, {"CAC", "CAU"}, {3288 / 4320.0f, 1.0f}},
	{I, 3, {"AUA", "AUC", "AUU"}, {3172 / 12071.0f, 8251 / 12071.0f, 1.0f}},
	{K, 2, {"AAA", "AAG"}, {12845 / 15169.0f, 1.0f}},
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
	{Y, 2, {"UAU", "UAC"}, {5768 / 7114.0f, 1.0f}}};
/* --------------------------------------------------------- end definition ---------------------------------------------------------------------- */

#define OBJECTIVE_NUM 3 // three objective function
#define _mCAI 0
#define _mHD 1
#define _MLRCS 2
/* --------------------------------------------------------- Population definition ---------------------------------------------------------------- */
typedef struct
{
	char *cds;						  // CDSs's sequences
	int p, q, l;					  // this if for MLRCS starting point and length
	int obj_cdsidx[OBJECTIVE_NUM][2]; // CDS's index correspond to objective function
	float obj_val[OBJECTIVE_NUM];	  // objective function value (0 ~ 1) for Pareto Comparsion
} Solution;
typedef struct
{
	int counter;			 // checking counter to obsolete this solution
	int rank;				 // indicate Pareto front (rank)
	float crowding_distance; // indicate diversity of solution in same rank
	float fitness;			 // fitness = 1 / rank
	float sel_prob;			 // selection probability
	Solution sol;
} Population;
/* ----------------------------------------------------------- end definition --------------------------------------------------------------------- */

/* ----------------------------------------------------------- Find index definition --------------------------------------------------------------- */
/* this function find aminoacid index using binary search */
int FindAminoIndex(const AA aminoacid)
{
	int low = 0;
	int high = 19;
	int mid;

	while (low <= high)
	{
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
int FindCodonIndex(int amino_idx, const char *codon)
{
	for (int i = 0; i < aa[amino_idx].num_codons; i++)
	{
		if (aa[amino_idx].codons[i][0] == codon[0] &&
			aa[amino_idx].codons[i][1] == codon[1] &&
			aa[amino_idx].codons[i][2] == codon[2])
		{
			return i;
		}
	}
}
/* --------------------------------------------------------------- end definition ------------------------------------------------------------------ */

/* ---------------------------------------------------------- Population memroy management --------------------------------------------------------- */
/* this function memory allocation population */
Population *AllocPopulation(int pop_size, int num_cds, int len_amino_seq)
{
	Population *pop;

	pop = (Population *)malloc(sizeof(Population) * pop_size);
	if (pop == NULL)
	{
		printf("Memory allocation failed at line %d\n", __LINE__);
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < pop_size; i++)
	{
		pop[i].sol.cds = (char *)malloc(sizeof(char) * num_cds * len_amino_seq * 3);
		if (pop[i].sol.cds == NULL)
		{
			printf("Memory allocation failed at line %d\n", __LINE__);
			exit(EXIT_FAILURE);
		}
	}

	return pop;
}
/* this function free population memory */
void FreePopulation(Population *pop, int pop_size, int num_cds)
{
	for (int i = 0; i < pop_size; i++)
	{
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
void GenCDS(char *cds, int num_cds, const int *amino_seq_idx, int len_amino_seq, int type = RANDOM_GEN)
{
	int idx;
	int rand_idx;
	int codon_idx;

	idx = 0;
	switch (type)
	{
		/* random creation */
	case RANDOM_GEN:
		for (int i = 0; i < num_cds; i++)
		{
			for (int j = 0; j < len_amino_seq; j++)
			{
				rand_idx = gen() % aa[amino_seq_idx[j]].num_codons;
				cds[idx++] = aa[amino_seq_idx[j]].codons[rand_idx][0];
				cds[idx++] = aa[amino_seq_idx[j]].codons[rand_idx][1];
				cds[idx++] = aa[amino_seq_idx[j]].codons[rand_idx][2];
			}
		}
		break;
		/* creation maximum CAI values */
	case UPPER_GEN:
		for (int i = 0; i < num_cds; i++)
		{
			for (int j = 0; j < len_amino_seq; j++)
			{
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
void GenSolution(Population *pop, int num_cds, const int *amino_seq_idx, int len_amino_seq, int type = RANDOM_GEN)
{
	/* initial value setting */
	pop->counter = 0; // new solution counter value is zero
	pop->rank = 0;
	pop->crowding_distance = 0;
	pop->fitness = 0;
	pop->sel_prob = 0;

	/* make random CDSs */
	GenCDS(pop->sol.cds, num_cds, amino_seq_idx, len_amino_seq, type);

	return;
}
/* this function copy population */
void CopyPopulation(const Population *origin, Population *target, int num_cds, int len_amino_seq)
{
	target->counter = origin->counter;
	target->rank = origin->rank;
	target->crowding_distance = origin->crowding_distance;
	target->fitness = origin->fitness;
	target->sel_prob = origin->sel_prob;

	for (int i = 0; i < num_cds * len_amino_seq * 3; i++)
	{
		target->sol.cds[i] = origin->sol.cds[i];
	}

	for (int i = 0; i < OBJECTIVE_NUM; i++)
	{
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
void ChageSynonymousCodon(int amino_idx, char *cds, int cd_idx, int type)
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
	case RANDOM_ADAPTATION: // change random synomynous codon
		if (aa[amino_idx].num_codons > 1)
		{
			while (true)
			{
				rand_idx = gen() % aa[amino_idx].num_codons;
				if (idx != rand_idx)
					break;
			}
			cds[cd_idx] = aa[amino_idx].codons[rand_idx][0];
			cds[cd_idx + 1] = aa[amino_idx].codons[rand_idx][1];
			cds[cd_idx + 2] = aa[amino_idx].codons[rand_idx][2];
		}
		break;
	case UPPER_ADAPTATION: // change random synomynous codon which has high adaptation value
		if (idx < aa[amino_idx].num_codons - 1)
		{
			rand_idx = idx + gen() % (aa[amino_idx].num_codons - 1 - idx) + 1;
			cds[cd_idx] = aa[amino_idx].codons[rand_idx][0];
			cds[cd_idx + 1] = aa[amino_idx].codons[rand_idx][1];
			cds[cd_idx + 2] = aa[amino_idx].codons[rand_idx][2];
		}
		break;
	}

	return;
}

/* this function generate random mutated solution */
Population *Mutation(const Population *pop, int num_cds, const int *amino_seq_idx, int len_amino_seq, float mprob)
{
	/* new population memory allocation */
	Population *new_pop;
	new_pop = AllocPopulation(1, num_cds, len_amino_seq);

	// copy population to new_population
	CopyPopulation(pop, new_pop, num_cds, len_amino_seq);
	new_pop->counter = 0;

	/* generate (0 ~ 1) random number corresponding to codon in CDS */
	float *random_num;
	random_num = (float *)malloc(sizeof(float) * num_cds * len_amino_seq);
	for (int i = 0; i < num_cds * len_amino_seq; i++)
	{
		random_num[i] = (float)dist(gen);
	}

	int type;
	int len_cds;
	type = gen() % 4;
	len_cds = 3 * len_amino_seq;
	/* four type of variations */
	switch (type)
	{
		/* change each codons with random codon in all CDSs */
	case 0:
		for (int i = 0; i < num_cds; i++)
		{
			for (int j = 0; j < len_amino_seq; j++)
			{
				if (random_num[i * len_amino_seq + j] < mprob)
				{
					ChageSynonymousCodon(amino_seq_idx[j], new_pop->sol.cds, i * len_cds + j * 3, RANDOM_ADAPTATION);
				}
			}
		}
		break;
		/* change each codons with higher adaptation value in CDS which have minimun CAI value */
	case 1:
		for (int i = 0; i < len_amino_seq; i++)
		{
			if (random_num[new_pop->sol.obj_cdsidx[_mCAI][0] * len_amino_seq + i] < mprob)
			{
				ChageSynonymousCodon(amino_seq_idx[i], new_pop->sol.cds, new_pop->sol.obj_cdsidx[_mCAI][0] * len_cds + i * 3, UPPER_ADAPTATION);
			}
		}
		break;
		/* change each codons with random codon in pair of CDSs with minimun Hamming Distance */
	case 2:
		for (int i = 0; i < len_amino_seq; i++)
		{
			if (random_num[new_pop->sol.obj_cdsidx[_mHD][0] * len_amino_seq + i] < mprob)
			{
				ChageSynonymousCodon(amino_seq_idx[i], new_pop->sol.cds, new_pop->sol.obj_cdsidx[_mHD][0] * len_cds + i * 3, RANDOM_ADAPTATION);
			}
			if (random_num[new_pop->sol.obj_cdsidx[_mHD][1] * len_amino_seq + i] < mprob)
			{
				ChageSynonymousCodon(amino_seq_idx[i], new_pop->sol.cds, new_pop->sol.obj_cdsidx[_mHD][1] * len_cds + i * 3, RANDOM_ADAPTATION);
			}
		}
		break;
		/* chage each codons with random codon in CDSs with longest common substring */
	case 3:
		for (int i = (new_pop->sol.p - new_pop->sol.obj_cdsidx[_MLRCS][0] * len_cds) / 3; i <= (new_pop->sol.p + new_pop->sol.l - 1 - new_pop->sol.obj_cdsidx[_MLRCS][0] * len_cds) / 3; i++)
		{
			if (random_num[new_pop->sol.obj_cdsidx[_MLRCS][0] * len_amino_seq + i] < mprob)
			{
				ChageSynonymousCodon(amino_seq_idx[i], new_pop->sol.cds, new_pop->sol.obj_cdsidx[_MLRCS][0] * len_cds + i * 3, RANDOM_ADAPTATION);
			}
		}
		for (int i = (new_pop->sol.q - new_pop->sol.obj_cdsidx[_MLRCS][1] * len_cds) / 3; i <= (new_pop->sol.q + new_pop->sol.l - 1 - new_pop->sol.obj_cdsidx[_MLRCS][1] * len_cds) / 3; i++)
		{
			if (random_num[new_pop->sol.obj_cdsidx[_MLRCS][1] * len_amino_seq + i] < mprob)
			{
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
void mCAI(Population *pop, int num_cds, const int *amino_seq_idx, int len_amino_seq)
{
	char codon[3];
	int idx;
	int codon_idx;
	double tmp;

	idx = 0;
	pop->sol.obj_val[_mCAI] = 1;
	for (int i = 0; i < num_cds; i++)
	{
		tmp = 1;
		for (int j = 0; j < len_amino_seq; j++)
		{
			codon[0] = pop->sol.cds[idx++];
			codon[1] = pop->sol.cds[idx++];
			codon[2] = pop->sol.cds[idx++];
			codon_idx = FindCodonIndex(amino_seq_idx[j], codon);
			tmp *= pow(aa[amino_seq_idx[j]].adaptation[codon_idx], 1.0 / len_amino_seq);
		}
		// tmp = pow(tmp, 1.0 / len_amino_seq);
		if ((float)tmp <= pop->sol.obj_val[_mCAI])
		{
			pop->sol.obj_val[_mCAI] = (float)tmp;
			pop->sol.obj_cdsidx[_mCAI][0] = i; // CDS's index having mCAI value
		}
	}

	return;
}
/* this function is calculate minimum Hamming Distance */
void mHD(Population *pop, int num_cds, int len_amino_seq)
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
			if (tmp <= pop->sol.obj_val[_mHD])
			{
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
	int **LCS; // matrix which size [length of CDS + 1][length of CDS + 1]
	int len_cds = 3 * len_amino_seq;
	int max_len;

	/* memory allocation for LCS matrix */
	LCS = (int **)malloc(sizeof(int *) * (len_cds + 1));
	for (int i = 0; i < len_cds + 1; i++)
	{
		LCS[i] = (int *)malloc(sizeof(int) * (len_cds + 1));
	}

	max_len = 0;
	for (int i = 0; i < num_cds; i++)
	{ // i'th CDS
		for (int j = i; j < num_cds; j++)
		{ // j'th CDS
			for (int k = 0; k < len_cds + 1; k++)
			{ // matrix row
				for (int l = 0; l < len_cds + 1; l++)
				{ // matrix column
					if (i != j)
					{
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
	for (int i = 0; i < len_cds + 1; i++)
	{
		free(LCS[i]);
	}
	free(LCS);

	return;
}
/* If new population dominate old population return true */
bool ParetoComparison(const Population *new_pop, const Population *old_pop)
{
	// weak pareto dominance
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
void SortbyRankCrowding(Population *pop, int pop_size, int num_cds, int len_amino_seq)
{
	/* this point out population index value */
	int **Sp, **F;
	int *np, *Q;
	int Sp_idx, F_front, F_idx, Q_idx;

	/* memory allocation */
	Sp = (int **)malloc(sizeof(int *) * pop_size);
	F = (int **)malloc(sizeof(int *) * pop_size);
	for (int i = 0; i < pop_size; i++)
	{
		Sp[i] = (int *)malloc(sizeof(int) * pop_size);
		F[i] = (int *)malloc(sizeof(int) * pop_size);
	}
	np = (int *)malloc(sizeof(int) * pop_size);
	Q = (int *)malloc(sizeof(int) * pop_size);

	// F empty initialization
	for (int i = 0; i < pop_size; i++)
	{
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
			if (i != j)
			{
				if (ParetoComparison(&pop[i], &pop[j]))
					Sp[i][Sp_idx++] = j;
				else if (ParetoComparison(&pop[j], &pop[i]))
					np[i] += 1;
			}
		}
		if (np[i] == 0)
		{
			pop[i].rank = 1;
			F[0][F_idx++] = i;
		}
	}
	/* -------------- 1st front setting complete ---------------- */

	int p;		 // indicate which p solution's Sp		p value means n'th Pareto Front solution
	F_front = 0; // indicate 1st front
	F_idx = 0;
	while (F[F_front][F_idx] != EMPTY)
	{
		Q_idx = 0;
		memset(Q, EMPTY, sizeof(int) * pop_size);
		for (F_idx = 0; F_idx < pop_size && F[F_front][F_idx] != EMPTY; F_idx++)
		{ // F[F_front][F_idx] == population index
			p = F[F_front][F_idx];
			for (Sp_idx = 0; Sp_idx < pop_size && Sp[p][Sp_idx] != EMPTY; Sp_idx++)
			{
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
		while (Q[Q_idx] != EMPTY)
		{
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
		p_len = 0; // number of solutions in pareto rank
		for (F_idx = 0; F_idx < pop_size && F[F_front][F_idx] != EMPTY; F_idx++)
		{
			pop[F[F_front][F_idx]].crowding_distance = 0;
			p_len++;
		}
		for (int i = 0; i < OBJECTIVE_NUM; i++)
		{
			// sort each objective function ascending order
			for (int j = 0; j < p_len; j++)
			{
				for (int k = 0; k < p_len - 1 - j; k++)
				{
					if (pop[F[F_front][k]].sol.obj_val[i] > pop[F[F_front][k + 1]].sol.obj_val[i])
					{
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

		for (int j = 0; j < p_len; j++)
		{
			for (int k = 0; k < p_len - 1 - j; k++)
			{
				if (pop[F[F_front][k]].crowding_distance < pop[F[F_front][k + 1]].crowding_distance)
				{
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

	Population *tmp_pop;
	tmp_pop = AllocPopulation(pop_size, num_cds, len_amino_seq);
	for (int i = 0; i < pop_size; i++)
	{
		CopyPopulation(&pop[i], &tmp_pop[i], num_cds, len_amino_seq);
	}
	F_front = 0;
	F_idx = 0;
	int p_idx = 0;
	while (F_front < pop_size && F[F_front][F_idx] != EMPTY)
	{
		CopyPopulation(&tmp_pop[F[F_front][F_idx]], &pop[p_idx], num_cds, len_amino_seq);
		p_idx++;
		F_idx++;
		if (F_idx == pop_size || F[F_front][F_idx] == EMPTY)
		{
			F_front++;
			F_idx = 0;
		}
	}

	// allocate fitness value
	for (int i = 0; i < pop_size; i++)
	{
		pop[i].fitness = 1.f / pop[i].rank;
	}

	/* free memory */
	FreePopulation(tmp_pop, pop_size, num_cds);
	for (int i = 0; i < pop_size; i++)
	{
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
void CalSelectionProb(Population *pop, int pop_size)
{
	float sum;

	sum = 0;
	for (int i = 0; i < pop_size; i++)
	{
		sum += pop[i].fitness;
	}
	for (int i = 0; i < pop_size; i++)
	{
		pop[i].sel_prob = pop[i].fitness / sum;
	}

	return;
}
/* Roulette-wheel selection based on selection probability */
int SelectSolution(const Population *pop, int pop_size)
{
	float sum;
	float point;

	point = (float)dist(gen); // 0 ~ 1 floating point number

	sum = 0;
	for (int i = 0; i < pop_size; i++)
	{
		sum += pop[i].sel_prob;
		if (point < sum)
			return i; // return selected pop index
	}
}

#define NUM_THREADS 16


#define IDEAL_MCAI 1
#define IDEAL_MHD 1
#define IDEAL_MLRCS 0
#define EUCLID(val1, val2, val3) (float)sqrt(pow(IDEAL_MCAI - val1, 2) + pow(IDEAL_MHD - val2, 2) + pow(val3, 2))
/* Minimum distance to optimal objective value(point) */
float MinEuclid(Population *pop, int size)
{
	float res;
	float tmp;

	res = 10000;
	for (int i = 0; i < size; i++)
	{
		tmp = EUCLID(pop[i].sol.obj_val[_mCAI], pop[i].sol.obj_val[_mHD], pop[i].sol.obj_val[_MLRCS]);
		if (tmp < res)
			res = tmp;
	}

	return res;
}


int main()
{
	FILE *global_fp;
	FILE *pipe;
	int global_cycle;
	char command[100] = "./hv -r \"0 0 1\" test.txt";
	char buffer[128];
	const char *ffile_name[] = {"Q5VZP5.fasta.txt",
								"A4Y1B6.fasta.txt",
								"B3LS90.fasta.txt",
								"B4TWR7.fasta.txt",
								"Q91X51.fasta.txt",
								"Q89BP2.fasta.txt"};
	int ff_cds_num[] = {2, 3, 4, 5, 6, 7};
	const char *ff_result[] = {
		"/home/jeus5771/CPU_protein/My_MOABC/My_SP_result/My_SP_result1.csv",
		"/home/jeus5771/CPU_protein/My_MOABC/My_SP_result/My_SP_result2.csv",
		"/home/jeus5771/CPU_protein/My_MOABC/My_SP_result/My_SP_result3.csv",
		"/home/jeus5771/CPU_protein/My_MOABC/My_SP_result/My_SP_result4.csv",
		"/home/jeus5771/CPU_protein/My_MOABC/My_SP_result/My_SP_result5.csv",
		"/home/jeus5771/CPU_protein/My_MOABC/My_SP_result/My_SP_result6.csv"};

	printf("OpenMP max thread : %d\n", omp_get_max_threads());
	omp_set_num_threads(omp_get_max_threads());
	
	//printf("OpenMP current using thread : %d\n", NUM_THREADS);

	/* amino sequence recieve from FASTA format */
	char file_name[128];
	char buf[256];
	char *amino_seq;	// amino sequence comprising a CDS
	int *amino_seq_idx; // amino sequence corresponding index value to struct 'aa'
	int len_amino_seq;	// length of amino sequence

	int max_cycle;	 // number of generations
	int colony_size; // number of solutions in population
	int num_cds;	 // number of CDSs
	int limit;		 // number of solution is not updated
	float mprob;	 // mutation probability

	/* ---------------------------------------- file processing ------------------------------------------------- */
	// printf("input file name : ");
	// scanf("%s", &file_name);
	/* input parameter values */
	printf("input max cycle value : ");
	scanf("%d", &max_cycle);
	if (max_cycle <= 0)
	{
		printf("input max cycle value > 0\n");
		return EXIT_FAILURE;
	}
	printf("input colony size : ");
	scanf("%d", &colony_size);
	if (colony_size <= 0)
	{
		printf("input colony size > 0\n");
		return EXIT_FAILURE;
	}
	// printf("input number of CDSs : "); scanf("%d", &num_cds);
	// if (num_cds <= 1) {
	// 	printf("input number of CDSs > 1\n");
	// 	return EXIT_FAILURE;

	// }
	printf("input limit value : ");
	scanf("%d", &limit);
	if (limit <= 0)
	{
		printf("input limit value > 0\n");
		return EXIT_FAILURE;
	}
	printf("input mutation probability (0 ~ 1 value) : ");
	scanf("%f", &mprob);
	if (mprob < 0 || mprob > 1)
	{
		printf("input mutation probability (0 ~ 1 value) : \n");
		return EXIT_FAILURE;
	}

	global_cycle = 3;
	float total_time = 0;
	float min_dist = 0;
	for (int zz = 0; zz < sizeof(ff_cds_num) / sizeof(int); zz++)
	{
		global_fp = fopen(ff_result[zz], "w");
		fprintf(global_fp, "time,short distance,hyper volume,%s\n\n", ffile_name[zz]);
		memset(file_name, 0, sizeof(file_name));
		memcpy(file_name, ffile_name[zz], sizeof(file_name));
		num_cds = ff_cds_num[zz];

		for (int z = 0; z < global_cycle; z++)
		{
			total_time = 0;
			FILE *fp;
			fp = fopen(file_name, "r");
			if (fp == NULL)
			{
				printf("Opening input file failed at line : %d", __LINE__);
				return EXIT_FAILURE;
			}
			fseek(fp, 0, SEEK_END);
			len_amino_seq = ftell(fp);
			fseek(fp, 0, SEEK_SET);
			fgets(buf, 256, fp); // jump over first line
			len_amino_seq -= ftell(fp);
			amino_seq = (char *)malloc(sizeof(char) * len_amino_seq); // over memory allocation
			if (amino_seq == NULL)
			{
				printf("Memory allocation failed at line : %d", __LINE__);
				return EXIT_FAILURE;
			}
			int idx = 0;
			char tmp;
			while (!feof(fp))
			{
				tmp = fgetc(fp);
				if (tmp != '\n')
					amino_seq[idx++] = tmp;
			}
			amino_seq[idx] = NULL;
			len_amino_seq = idx - 1;

			fclose(fp);

			amino_seq_idx = (int *)malloc(sizeof(int) * len_amino_seq); // memory alloc
			if (amino_seq_idx == NULL)
			{
				printf("Memory allocation failed at line : %d", __LINE__);
				return EXIT_FAILURE;
			}
			for (int i = 0; i < len_amino_seq; i++)
			{
				amino_seq_idx[i] = FindAminoIndex((AA)amino_seq[i]);
			}
			/* -------------------------------------------- end file proess ------------------------------------------------- */

			/* Population memory allocation */
			Population *pop;
			float total_time = 0;
			pop = AllocPopulation(colony_size * 2, num_cds, len_amino_seq);

			fp = fopen("test.txt", "w");
			std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
#pragma omp parallel
			{
				Population *new_sol, *sel_sol;
				Population *tmp_sol; // Generating new solution for Scout Bees step
				tmp_sol = AllocPopulation(1, num_cds, len_amino_seq);
				bool check;

				/* --------------------------------------------------- initialize Population ------------------------------------------------------- */
#pragma omp for
				for (int i = 0; i < colony_size; i++)
				{
					/* To boost the optimization of solutions witth high CAI values
					   remaing solution is generated by selecting highest adaptation */
					if (i == colony_size - 1)
					{
						GenSolution(&pop[i], num_cds, amino_seq_idx, len_amino_seq, UPPER_GEN);
						mCAI(&pop[i], num_cds, amino_seq_idx, len_amino_seq);
						mHD(&pop[i], num_cds, len_amino_seq);
						MLRCS(&pop[i], num_cds, len_amino_seq);
					}
					else
					{
						GenSolution(&pop[i], num_cds, amino_seq_idx, len_amino_seq, RANDOM_GEN);
						/* calculate objective function value */
						mCAI(&pop[i], num_cds, amino_seq_idx, len_amino_seq);
						mHD(&pop[i], num_cds, len_amino_seq);
						MLRCS(&pop[i], num_cds, len_amino_seq);
					}
				}
				/* -------------------------------------------------------- initialize end ----------------------------------------------------------- */

				/* ------------------------------------------------------- start cycle to max cycle -------------------------------------------------- */
				for (int cycle = 0; cycle < max_cycle; cycle++)
				{
					/* --------------------------------------- start Employed bees step ----------------------------------------- */
#pragma omp for
					for (int i = 0; i < colony_size; i++)
					{
						new_sol = Mutation(&pop[i], num_cds, amino_seq_idx, len_amino_seq, mprob); // Employed Bee search
						/* Calculate Objective Functions */
						mCAI(new_sol, num_cds, amino_seq_idx, len_amino_seq);
						mHD(new_sol, num_cds, len_amino_seq);
						MLRCS(new_sol, num_cds, len_amino_seq);
						/* Pareto Comparision */
						check = ParetoComparison(new_sol, &pop[i]);
						if (check)
							CopyPopulation(new_sol, &pop[i], num_cds, len_amino_seq);
						else
							pop[i].counter += 1;
						FreePopulation(new_sol, 1, num_cds);
					}
					/* ----------------------------------------- end Employed bees step ----------------------------------------- */

#pragma omp single
					{
						SortbyRankCrowding(pop, colony_size, num_cds, len_amino_seq);
						CalSelectionProb(pop, colony_size);
					}

					/* -------------------------------------- start Onlooker bees step ------------------------------------------ */
#pragma omp for
					for (int i = colony_size; i < 2 * colony_size; i++)
					{
						sel_sol = &pop[SelectSolution(pop, colony_size)];						   // select solution
						new_sol = Mutation(sel_sol, num_cds, amino_seq_idx, len_amino_seq, mprob); // Onlooker Bee search
						/* Calculate Objective Function */
						mCAI(new_sol, num_cds, amino_seq_idx, len_amino_seq);
						mHD(new_sol, num_cds, len_amino_seq);
						MLRCS(new_sol, num_cds, len_amino_seq);
						/* Pareto Comparison */
						check = ParetoComparison(new_sol, sel_sol);
						if (check)
							CopyPopulation(new_sol, &pop[i], num_cds, len_amino_seq);
						else
						{
							sel_sol->counter++;
							CopyPopulation(sel_sol, &pop[i], num_cds, len_amino_seq);
						}
						FreePopulation(new_sol, 1, num_cds);
					}
					/* ------------------------------------------ end Onlooker bees step ----------------------------------------- */

					/* ------------------------------------- start Scout bees step ----------------------------------------------- */
#pragma omp for
					for (int i = 0; i < 2 * colony_size; i++)
					{
						if (pop[i].counter > limit)
						{
							GenSolution(tmp_sol, num_cds, amino_seq_idx, len_amino_seq, RANDOM_GEN);
							mCAI(tmp_sol, num_cds, amino_seq_idx, len_amino_seq);
							mHD(tmp_sol, num_cds, len_amino_seq);
							MLRCS(tmp_sol, num_cds, len_amino_seq);
							/* Scout Bee search */
							for (int j = 0; j < cycle / limit; j++)
							{
								new_sol = Mutation(tmp_sol, num_cds, amino_seq_idx, len_amino_seq, mprob);
								mCAI(new_sol, num_cds, amino_seq_idx, len_amino_seq);
								mHD(new_sol, num_cds, len_amino_seq);
								MLRCS(new_sol, num_cds, len_amino_seq);
								CopyPopulation(new_sol, tmp_sol, num_cds, len_amino_seq);
								FreePopulation(new_sol, 1, num_cds);
							}
							CopyPopulation(tmp_sol, &pop[i], num_cds, len_amino_seq);
							pop[i].counter = 0;
						}
					}
					/* ------------------------------------------ end Scout bees step -------------------------------------------- */
#pragma omp single
					{
						SortbyRankCrowding(pop, colony_size * 2, num_cds, len_amino_seq);

						// non-dominated file update
						for (int i = 0; i < colony_size; i++)
						{
							if (pop[i].rank == 1)
							{
								fprintf(fp, "%f %f %f\n", -pop[i].sol.obj_val[_mCAI], -pop[i].sol.obj_val[_mHD] / 0.4, pop[i].sol.obj_val[_MLRCS]);
							}
						}
					}
				}
				FreePopulation(tmp_sol, 1, num_cds);
			}
			/* ----------------------------------------------------- end max cycle ---------------------------------------------------------------- */
			fclose(fp);
			std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
			total_time += static_cast<float>(sec.count());
			//min_dist = MinEuclid(pop, colony_size * 2);


			pipe = popen(command, "r");
			if (!pipe)
			{
				printf("Failed to execute command\n");
				return 1;
			}

			// Read the output of the command
			while (fgets(buffer, sizeof(buffer), pipe))
			{
				printf("%s", buffer);
			}

			// Close the pipe
			pclose(pipe);
			fprintf(global_fp, "%f,%f,%f\n", total_time, min_dist, atof(buffer));

			/* free memory */
			FreePopulation(pop, colony_size * 2, num_cds);
			free(amino_seq);
			free(amino_seq_idx);
		}
		fclose(global_fp);
	}

	return EXIT_SUCCESS;
}