#include <iostream>
#include <cstdlib>

using namespace std;

// 아미노산 20개 대응 코돈 조합 정의
// 해당 코돈 weight 정의
// 

#define OBJECTIVE_NUM 3

class Solution
{
private:
	unsigned int counter;			
	int obj_val[OBJECTIVE_NUM];		// value of objective function
public:
	int a;

	Solution()
	{
		this->counter = 0;
		for (int i = 0; i < OBJECTIVE_NUM; i++)
		{
			this->counter = 0;
		}
	}

};

typedef struct
{
	
}Population;

int main()
{
	/* user input parameter */
	unsigned int colony_size;	// number of solutions in the population
	unsigned int max_cycles;	// maximun number of generations
	unsigned int limit;			// stagnation limit
	double m_prob;				// mutation probability



	return EXIT_SUCCESS;
}