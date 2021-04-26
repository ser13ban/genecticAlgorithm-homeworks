#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <math.h>

using namespace std;
using namespace std::chrono;
ofstream sch("schefel.txt");

#define NUMBER_OF_GENERATIONS 1000
#define PI 3.14159265
#define CROSS_PROBABILITY 0.2
#define MUTATION_PROBABILITY 0.005
#define POPULATION_SIZE 100

#define SAMPLE_SIZE	100


struct cromozome
{
	double fitness;
	double prob_selectie;
	double q;
	vector<bool> bits;
};
vector<cromozome> population(POPULATION_SIZE);
vector<cromozome> nextGeneration;

double DeJong(vector<double> vec)
{
	double result = 0;
	for (int i = 0; i < vec.size(); i++)
	{
		result += pow(vec[i], 2);
	}
	return result;
}


double Schwefel(vector<double> vec)
{
	double result = 0;
	for (int i = 0; i < vec.size(); i++)
	{
		result -= vec[i] * sin(sqrt(abs(vec[i])));
	}
	return result;
}


double Rastrigin(vector<double> vec)
{
	double result = 0;
	for (int i = 0; i < vec.size(); i++)
	{
		result += (pow(vec[i], 2) - 10 * cos(2 * PI * vec[i]));
	}
	return 10 * vec.size() + result;
}

double Michalewicz(vector<double> vec)
{
	double result = 0;
	for (int i = 0; i < vec.size(); i++)
	{
		result += sin(vec[i]) * pow((sin(i * pow(vec[i], 2) / PI)), 2 * 10);
	}
	return -result;
}


double RandomDouble(const double lower, const double upper) {
	if (lower > upper) {
		throw "ERROR: the second argument of RandomDouble must be bigger than the first argument!";
	}
	random_device generator;
	uniform_real_distribution<double> distribution(lower, upper);
	return distribution(generator);
}

int computeLength(double a, double b, int precisioin)
{
	return ceil(log((b - a) * pow(10, precisioin)) / log(2));
}

double decode_dimension(vector<bool>::iterator start, vector<bool>::iterator end, int l, double a, double b)
{
	unsigned long bi = 0;
	for (auto ii = start; ii != end; ++ii)
	{
		bi *= 2;
		bi += *ii;
	}
	double s = bi / (pow(2, l) - 1);
	return s * (b - a) + a;
}

vector<double> decode(vector<bool>& bits, unsigned n, int l, double a, double b)
{
	vector<double> ret;
	vector<bool>::iterator start, end;
	for (int ii = 0; ii < n; ++ii)
	{
		start = bits.begin() + ii * l;
		end = start + l;
		double x = decode_dimension(start, end, l, a, b);
		//cout << x << endl;
		ret.push_back(x);
	}
	return ret;
}

void PopulateCromozomes(vector<cromozome>& population, int l, int n)
{
	for (int i = 0; i <= population.size() - 1; i++)
	{
		vector<bool> bits(n * l, 0);
		for (int j = 0; j < n * l; j++)
		{
			int x = rand() % 2;
			bits[j] = x;
		}
		population[i].bits = bits;
	}
}


void computeFitness(double(*function)(vector<double>), vector<cromozome>& population, double a, double b, int l, int n)
{
	for (int i = 0; i < population.size() - 1; i++)
	{
		//population[i].fitness = pow(1 / function(decode(population[i].bits, n, l, a, b)), 4);
		if (function(decode(population[i].bits, n, l, a, b)) < 0)
			population[i].fitness = pow(1 / (function(decode(population[i].bits, n, l, a, b))), -4);
		else
			population[i].fitness = pow(1 / (function(decode(population[i].bits, n, l, a, b))), 4);
	}
}


void mutatie1(vector<cromozome>& population, int n, int l)
{
	for (int i = 0; i <= population.size() - 1; i++)
	{
		for (int j = 0; j < n * l; j++)
		{
			int x = rand() % 100 + 1;
			if (x == 1)
			{
				//int y = rand() % population[i].bits.size();
				population[i].bits[j] = !population[i].bits[j];
			}
		}


	}
}

void mutatie(vector<cromozome>& population, int n, int l)
{
	for (int i = 0; i <= population.size() - 1; i++)
	{
		for (int j = 0; j < n * l; j++)
		{
			if ((1.0 * rand() / RAND_MAX) < MUTATION_PROBABILITY)
			{
				//int y = rand() % population[i].bits.size();
				population[i].bits[j] = !population[i].bits[j];
			}
		}


	}
}
void makeCross(vector<cromozome>& population, int n, int l, int index1, int index2)
{ 
	int x = rand() % n * l;
	for (int i = x; i < population[1].bits.size(); i++)
	{
		swap(population[index1].bits[i], population[index2].bits[i]);
	}
}


void CrossOver(vector<cromozome>& population, int n, int l)
{
	int cr = 0;
	int index1, index2;
	for (int i = 0; i < population.size(); i++)
	{
		double x = (1.0 * rand() / RAND_MAX);
		//cout << endl << x << endl << " ---- " << endl;
		if (x < CROSS_PROBABILITY)
		{
			cr++;
			if (cr == 2)
			{
				cr = 0;
				index2 = i;
				makeCross(population, n, l, index1, index2);
			}
			if (cr == 1)
			{
				index1 = i;
			}
		}
	}
}

double RandomDouble()
{
	return static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
}

void  Evaluate(vector<cromozome> population, double& best, int n, int l,double a,double b, double(*function)(vector<double>)) {

	for (int i = 0; i <= population.size() - 1; i++)
	{
		if (function(decode(population[i].bits, n, l, a, b)) < best)
			best = function(decode(population[i].bits, n, l, a, b));
	}
}
void Select(vector<cromozome>& population)
{
	double total_fitness = 0;

	for (int i = 0; i <= population.size() - 1; i++)
		total_fitness += population[i].fitness;
	//cout << "total fitness: -------- " << total_fitness << endl;

	for (int i = 0; i <= population.size() - 1; i++)
	{
		population[i].prob_selectie = population[i].fitness / total_fitness;
		//cout << population[i].prob_selectie << endl;
	}

	population[0].q = 0;
	for (int i = 0; i < population.size() - 1; i++)
	{
		population[i + 1].q = population[i].q + population[i].prob_selectie;
	}
	int k = 0;
	for (int i = 0; i <= population.size() - 1; i++)
	{
		double r = RandomDouble();
		for (int j = 0; j < population.size() - 1; j++)
		{
			if (population[j].q < r && r <= population[j + 1].q)
			{
				nextGeneration.resize(k + 1);
				nextGeneration[k++].bits = population[j].bits;
			}
		}
	}

	population = nextGeneration;
}

double GeneticAlgorithm(double(*function)(vector<double>), vector<cromozome>& population, double a, double b, int l, int n)
{
	int generation = 0;
	population.clear();
	population.resize(POPULATION_SIZE);
	PopulateCromozomes(population, l, n);
	computeFitness(function, population, a, b, l, n);

	double best = INFINITY, best_fitness = -INFINITY;

	Evaluate(population, best, n, l, a, b, function);

	while (generation < NUMBER_OF_GENERATIONS)
	{
		generation++;
	
		Select(population);

		mutatie(population, n, l);

		CrossOver(population, n, l);

		computeFitness(function, population, a, b, l, n);

		Evaluate(population, best, n, l, a, b, function);

		cout << "\nGENERATION: " << generation << " BEST CROMOZOME: " << best;

	}
	return best;
}


int main()
{
	double a = -500, b = 500;
	int n = 5; //5, 10, 30
	int precisioin = 5; //precisioin=5
	int L, l;
	l = computeLength(a, b, precisioin);
	L = l * n;
	srand(time(NULL));
	sch << "\nn=5";
	for (int i = 1; i <= SAMPLE_SIZE; ++i)
	{
		auto start = high_resolution_clock::now();
		sch << "\n" << GeneticAlgorithm(Schwefel, population, a, b, l, n);
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<milliseconds>(stop - start);
		sch << "\t" << duration.count();
	}
	n = 10;
	sch << "\nn=10";
	for (int i = 1; i <= SAMPLE_SIZE; ++i)
	{
		auto start = high_resolution_clock::now();
		sch << "\n" << GeneticAlgorithm(Schwefel, population, a, b, l, n);
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<milliseconds>(stop - start);
		sch << "\t" << duration.count();
	}
	sch << "\nn=30";
	n = 30;
	for (int i = 1; i <= SAMPLE_SIZE; ++i)
	{
		auto start = high_resolution_clock::now();
		sch << "\n" << GeneticAlgorithm(Schwefel, population, a, b, l, n);
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<milliseconds>(stop - start);
		sch << "\t" << duration.count();
	}
}
