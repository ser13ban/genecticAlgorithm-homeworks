// Tema2_GA.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <math.h>
/*
* SELECT
* MUTATE
* CROSSOVER
* O sa am nevoie de cele 3 functii: SELECT() CROSSOVER() MUTATE()
* Nu stiu cum ar trebui implementata functia fitness
* Ma gandesc sa folosesc metoda ruleta, pare sa fie cea mai usoara
*/

#define PI  3.14159265358979323846  
#define MAX 10000
#define POP_SIZE 100

using namespace std;
vector<vector<char>> Population(POP_SIZE);
struct Cromozome
{
	double fitness;
	double selectionProbability;
	double q;	//Acumulated probability
	vector<char> bits;
	Cromozome() {
		fitness = 0;
		selectionProbability = 0;
		q = 0;
	}
};
vector<Cromozome> Pop(100);

double DeJong(vector<double> values) {	//case 1
	double result = 0;
	for (int i = 0; i < values.size(); ++i) {
		double value = pow(values[i], 2);
		result += value;
	}
	return result;
}
double Schwefel(vector<double> values) {	//case 2
	double result = 0;
	for (int i = 0; i < values.size(); ++i) {
		double value = (-values[i]) + sin(sqrt(abs(values[i])));
		result += value;
	}
	return result;
}
double Rastrigin(vector<double> values) { //case 3
	double result = 10 * values.size();
	double sum = 0;
	for (int i = 0; i < values.size(); ++i) {
		double value = pow(values[i], 2) - (10 * cos(2 * PI * values[i]));
		sum += value;
	}
	return result + sum;
}
double  Michalewicz(vector<double> values) {//case 4
	double result = 0;
	for (int i = 0; i < values.size(); ++i) {
		double value = sin(values[i]) * pow((sin(i * pow(values[i], 2) / PI)), 20);
		result += value;
	}
	return -result;
}

int compute_l(double a, double b, int prec)
{

	double number_of_points = (b - a) * pow(10, prec); // calculam numarul de puncte din interval
	int length = ceil(log(number_of_points) / log(2)); //il aducem in biti
	return length;
}
double RandomDouble(const double lower, const double upper) {
	if (lower > upper) {
		throw "ERROR: the second argument of RandomDouble must be bigger than the first argument!";
	}
	random_device generator;
	uniform_real_distribution<double> distribution(lower, upper);
	return distribution(generator);
}

void GenerateRandomBits(vector<char>& v)
{
	for (auto ii = v.begin(); ii != v.end(); ++ii)
	{
		//alegere random: 
		*ii = (RandomDouble(0, 1) < 0.5) ? 0 : 1;
	}
}

double decode(double a, double b, int l, vector<char>::iterator  it_start, vector<char>::iterator it_end)
{

	int x_int = 0;
	double x_scaled;
	for (vector<char>::iterator it = it_start; it != it_end; it++)
	{
		x_int *= 2;
		x_int += *it;
	}

	x_scaled = x_int / (pow(2, l) - 1);

	return x_scaled * (b - a) + a;


}
vector<double> decode_values(double a, double b, int n, int l, vector<char>& bits)
{
	vector<char>::iterator it_start, it_end;
	vector<double> values;

	for (int ii = 0; ii < n; ++ii)
	{
		int start = ii * l;
		int end = start + l;
		it_start = bits.begin() + start;
		it_end = bits.begin() + end;
		double x = decode(a, b, l, it_start, it_end);
		//that cout was for testing purposes
		//cout << x << endl;
		values.push_back(x);
	}
	return values;
}


void Populate(vector<vector<char>>& Population, int a, int b, int n, int l) {
	int L = n * l;
	for (int i = 0; i < 100; ++i) {
		vector<char> bits(L);
		GenerateRandomBits(bits);
		Population[i] = bits;
		/*
		vector<double> values = decode_values(a, b, n, l, Population[i]);
		cout << "next chromozome\n";
		*/
	}
}
void PopulateStruct(vector<Cromozome>& Pop, int L) {
	for (int i = 0; i < POP_SIZE; ++i) {
		vector<char> bits(L);
		GenerateRandomBits(bits);
		Pop[i].bits = bits;
	}
}
void CalulateFitness(vector<Cromozome>& Pop, double (*function)(vector<double>), int a, int b, int n, int l, int L) {
	for (int i = 0; i < POP_SIZE; i++) {
		double fit;
		Pop[i].fitness = pow(1 / function(decode_values(a, b, n, l, Pop[i].bits)), 4);
		//cout << Pop[i].fitness<<'\n';
	}
}
void Mutate(vector<Cromozome>& Pop) {
	for (int i = 0; i < POP_SIZE; i++) {
		int x = rand() % 100 + 1;
		if (x == 1) 
		{
			int p = rand() % Pop[i].bits.size();
			Pop[i].bits[p] = !Pop[i].bits[p];
		}
	}
}
int Select(vector<Cromozome> Pop) 
{
	//make the sum of all of the fitnesses:
	double sumOfFitnesses = 0,partialSumOfFitnesses = 0;
	for (int i = 0; i < POP_SIZE; ++i)
	{
		sumOfFitnesses += Pop[i].fitness;
	}
	
	//generate a random number between 0 and S
	double random = RandomDouble(0, sumOfFitnesses);

	//starting form the top or form the bottom ,you should ckeck this in the raport
	int ii = 0;
	while (partialSumOfFitnesses < sumOfFitnesses && ii<POP_SIZE) {
		partialSumOfFitnesses += Pop[ii].fitness;
		ii++;
	}
	ii--;
	return ii;

}
vector<Cromozome> Select2(vector<Cromozome> Pop)
{
	double sumOfFitnesses = 0;
	vector<Cromozome> childPop(0);
	//Total fitnesses
	
	for (int i = 0; i < POP_SIZE; ++i)
	{
		sumOfFitnesses += Pop[i].fitness;
	}
	//Individaul seleciton probability
	for (int i = 0; i < POP_SIZE; ++i)
	{
		Pop[i].selectionProbability = Pop[i].fitness / sumOfFitnesses;
	}
	//acumulated selection probabilites
	Pop[0].q = 0;
	for (int i = 0; i < POP_SIZE - 1; i++)//you need to see if it goes where wiht POP_SIZE - 1 or POP_SIZE
	{
		Pop[i + 1].q = Pop[i].q + Pop[i].selectionProbability;
	}

	//now the selection:
	int k=0;
	for (int i = 0; i < POP_SIZE; ++i) 
	{
		//get a random form 0 to 1
		double random = RandomDouble(0, 1);
		for (int j = 0 ; j < POP_SIZE-1; ++j) 
		{
			if (Pop[j].q < random && random <= Pop[j + 1].q)
			{
				//here is where i need to push to the new stucture
				childPop.resize(k + 5);
				childPop[k] = Pop[j];
				k++;
			}
		}
	}
	return childPop;
}
double Evaluate(vector<Cromozome> Pop, double (*function)(vector<double>), int a, int b, int n, int l) {
	double best = INFINITY,function_result;
	for (int ii = 0; ii < POP_SIZE; ++ii) 
	{
		function_result = function(decode_values(a, b, n, l, Pop[ii].bits));
		if (best > function_result)
			best = function_result;
	}
	return best;
}
void GA(int a, int b, int n, int l, double (*function)(vector<double>)) {
	//am populat random cu 100 de cromozomi
	int L = n * l;
	int t = 0;
	PopulateStruct(Pop, L); //this generates the first generation of cromozomes
	CalulateFitness(Pop, function,a, b, n, l, L); //calculate the fitness funciton for each cromozome
	
	double best = Evaluate(Pop,function,a,b,n,l);
	while (t < 1000) {
		Pop = Select2(Pop);
		//select [I think I did this]
			//chestia aia de pe site pare super ok https://www.tutorialspoint.com/genetic_algorithms/genetic_algorithms_parent_selection.htm
		//mutate [I also think that i did this]
		//crossover no idea yet...
		t++;
	}
}
int main()
{
	srand(time(NULL));
	int a = -5, b = 5, n = 10, l;
	l = compute_l(a, b, 5);
	GA(a, b, n, l,Schwefel);
}