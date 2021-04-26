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


#define PI  3.14159265358979323846  
#define MAX 10000
#define POP_SIZE 100

using namespace std;
ofstream deJong("DeJong.txt");
ofstream michal("Micalewich.txt");
ofstream schefel("Rastigrin.txt");
ofstream rastigrin("Schefel.txt");

struct Cromozome
{
	double fitness;
	double selectionProbability;
	double q;	//Acumulated probability
	vector<bool> bits;
	
};
vector<Cromozome> Pop(POP_SIZE);

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

void GenerateRandomBits(vector<bool>& v)
{
	for (auto ii = v.begin(); ii != v.end(); ++ii)
	{
		//alegere random: 
		*ii = (RandomDouble(0, 1) < 0.5) ? 0 : 1;
	}
}
void swapp (vector<bool>& v, vector<bool>& v2) {
	for (auto ii = v.begin(), jj = v2.begin(); ii != v.end(),jj !=v2.end(); ++ii,++jj)
	{
		auto aux = *ii;
		*ii = *jj;
		*jj = aux;
	}
}

double decode(double a, double b, int l, vector<bool>::iterator  it_start, vector<bool>::iterator it_end)
{

	int x_int = 0;
	double x_scaled;
	for (vector<bool>::iterator it = it_start; it != it_end; it++)
	{
		x_int *= 2;
		x_int += *it;
	}

	x_scaled = x_int / (pow(2, l) - 1);

	return x_scaled * (b - a) + a;


}
vector<double> decode_values(double a, double b, int n, int l, vector<bool>& bits)
{
	vector<bool>::iterator it_start, it_end;
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



void PopulateStruct(vector<Cromozome>& Pop, int L) {
	for (int i = 0; i < POP_SIZE; ++i) {
		vector<bool> bits(L);
		GenerateRandomBits(bits);
		Pop[i].bits = bits;
	}
}
void CalulateFitness(vector<Cromozome>& Pop, double (*function)(vector<double>), int a, int b, int n, int l, int L) {
	for (int i = 0; i < POP_SIZE; i++) {
		double fit;
		double functionResult = function(decode_values(a, b, n, l, Pop[i].bits));
		if (functionResult < 0) {
			Pop[i].fitness = pow(1 / functionResult, -4);
		}
		else {
			Pop[i].fitness = pow(1 / functionResult, 4);
		}
		//cout << Pop[i].fitness<<'\n';
	}
}
void Mutate(vector<Cromozome>& Pop) {
	for (int i = 0; i < POP_SIZE; i++) {
		for (int j = 0; j < Pop[i].bits.size(); j++) {
			if (RandomDouble(0, 100) < 0.5) {
				Pop[i].bits[j] = !Pop[i].bits[j];
			}
		}
	}
}
int Select(vector<Cromozome> Pop)
{
	//make the sum of all of the fitnesses:
	double sumOfFitnesses = 0, partialSumOfFitnesses = 0;
	for (int i = 0; i < POP_SIZE; ++i)
	{
		sumOfFitnesses += Pop[i].fitness;
	}

	//generate a random number between 0 and S
	double random = RandomDouble(0, sumOfFitnesses);

	//starting form the top or form the bottom ,you should ckeck this in the raport
	int ii = 0;
	while (partialSumOfFitnesses < sumOfFitnesses && ii < POP_SIZE) {
		partialSumOfFitnesses += Pop[ii].fitness;
		ii++;
	}
	ii--;
	return ii;

}
void Select2(vector<Cromozome>& Pop)
{
	double sumOfFitnesses = 0;
	vector<Cromozome> childPop(POP_SIZE);
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
	int k = 0;
	for (int i = 0; i < POP_SIZE; ++i)
	{
		//get a random form 0 to 1
		double random = RandomDouble(0, 1);
		for (int j = 0; j < POP_SIZE - 1; ++j)
		{
			if (Pop[j].q < random && random <= Pop[j + 1].q)
			{
				//here is where i need to push to the new stucture
				childPop[k] = Pop[j];
				k++;
			}
		}
	}
	Pop = childPop;
}
double Evaluate(vector<Cromozome> Pop, double (*function)(vector<double>), int a, int b, int n, int l) 
{
	double best = INFINITY, function_result;
	int L = n * l;
	PopulateStruct(Pop, L); //this generates the first generation of cromozomes
	CalulateFitness(Pop, function, a, b, n, l, L); //calculate the fitness funciton for each cromozome

	for (int ii = 0; ii < POP_SIZE; ++ii)
	{
		function_result = function(decode_values(a, b, n, l, Pop[ii].bits));
		if (best > function_result)
			best = function_result;
	}
	return best;
}

void MakeCross(vector<Cromozome>& Pop, int i, int j,int L) 
{
	int random = int(RandomDouble(0, L));
	bool aux;
	//cout << random << '\n';
	for (int ii = random; ii < L; ++ii)
	{
		aux = Pop[i].bits[ii];
		Pop[i].bits[ii] = Pop[j].bits[ii];
		Pop[j].bits[ii] = aux;
		//swap(Pop[i].bits[ii], Pop[j].bits[ii]);
	}
}

void CrossOver(vector<Cromozome>& Pop,int L)
{
	int cross1 = -1, cross2 = -1, i = 0;
	while (cross1 == -1 && cross2 == -1) {
		if (i >= POP_SIZE) { break; }
		if (RandomDouble(0, 100) < 2 && cross1 == -1) 
		{
			cross1 = i;
		}
		else if (RandomDouble(0, 100) < 2 && cross2 == -1&& cross1 != -1) 
		{
			cross2 = i;
			//MakeCross(Pop, cross1, cross2, L);
			cout << cross1 <<"\t" <<cross2<< "\n";
			cross1 = -1;
			cross2 = -1;
			
		}
		++i;
	}
}
void CrossOver2(vector<Cromozome>& Pop, int L) 
{
	int cross1 = -1, cross2 = -1;
	for (int i = 0; i < Pop.size(); ++i) 
	{
		if (RandomDouble(0, 100) < 0.2 && cross1 == -1)
		{
			cross1 = i;
		}
		if (RandomDouble(0, 100) < 0.2 && cross2 == -1 && cross1 != -1)
		{
			cross2 = i;
			//MakeCross(Pop, cross1, cross2, L);
			swapp(Pop[cross1].bits, Pop[cross2].bits);
			//cout << cross1 << "\t" << cross2 << "\n";
			cross1 = -1;
			cross2 = -1;

		}
	}
}
double GA(int a, int b, int n, int l, double (*function)(vector<double>)) 
{
	//am populat random cu 100 de cromozomi
	int L = n * l;
	int t = 0;
	PopulateStruct(Pop, L);
	double best = Evaluate(Pop, function, a, b, n, l),bb;
	while (t < 1000) {
		Select2(Pop);
		Mutate(Pop);
		//CrossOver(Pop, L);
		CrossOver2(Pop, L); //this gives an error
		bb = Evaluate(Pop, function, a, b, n, l);
		if (bb < best)
			best = bb;
		cout << bb << '\n';
		t++;
	}
	return best;
}

int main()
{
	
	//srand(time(NULL));
	int n = 30, l;
	float a = -5.12, b=5.12;
	l = compute_l(a, b, 5);
	int L = n * l;
	cout << GA(a, b, n, l, DeJong);
	/*deJong << "DEJONG\nn=5";
	for (int ii = 1; ii <= 30; ++ii) {
		cout << GA(a, b, n, l, DeJong);
	}
	deJong << "n=10";
	for (int ii = 1; ii <= 30; ++ii) {
		cout << GA(a, b, n, l, DeJong);
	}
	deJong << "n=10";
	for (int ii = 1; ii <= 30; ++ii) {
		cout << GA(a, b, n, l, DeJong);
	}
	deJong<< "GA RESULT: "<< GA(a, b, n, l, DeJong);*/
	
}