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


using namespace std;
using namespace std::chrono;
ofstream outDeJong("DeJong.txt");
ofstream outSchwefel("Schwefel.txt");
ofstream outRastrigin("Rastrigin.txt");
ofstream outMichalewicz(" Michalewicz.txt");


int prec = 2;

double randomFloat(const double lower, const double upper) {
	if (lower > upper) {
		throw "ERROR: the second argument of randomFloat must be bigger than the first argument!";
	}
	random_device generator;
	uniform_real_distribution<double> distribution(lower, upper);
	return distribution(generator);
}

int compute_l(float a, float b, int prec)
{

	float number_of_points = (b - a) * pow(10, prec); // calculam numarul de puncte din interval
	int length = ceil(log(number_of_points) / log(2)); //il aducem in biti
	return length;
}


float decode(float a, float b, int l, vector<char>::iterator  it_start, vector<char>::iterator it_end)
{

	int x_int = 0;
	float x_scaled;
	for (vector<char>::iterator it = it_start; it != it_end; it++)
	{
		x_int *= 2;
		x_int += *it;
	}

	x_scaled = x_int / (pow(2, l) - 1);

	return x_scaled * (b - a) + a;


}


vector<float> decode_values(float a, float b, int n, int l, vector<char>& bits)
{
	vector<char>::iterator it_start, it_end;
	vector<float> values;

	for (int ii = 0; ii < n; ++ii)
	{
		int start = ii * l;
		int end = start + l;
		it_start = bits.begin() + start;
		it_end = bits.begin() + end;
		float x = decode(a, b, l, it_start, it_end);
		//that cout was for testing purposes
		//cout << x << endl;
		values.push_back(x);
	}
	return values;
}

//Give me a random set for the domain:
void GenerateRandomBits(vector<char>& v)
{
	for (auto ii = v.begin(); ii != v.end(); ++ii)
	{
		//alegere random: 
		*ii = (randomFloat(0, 1) < 0.5) ? 0 : 1;
	}
}

float DeJong(vector<float> values) {	//case 1
	float result = 0;
	for (int i = 0; i < values.size(); ++i) {
		float value = pow(values[i], 2);
		result += value;
	}
	return result;
}
float Schwefel(vector<float> values) {	//case 2
	float result = 0;
	for (int i = 0; i < values.size(); ++i) {
		float value = (-values[i]) + sin(sqrt(abs(values[i])));
		result += value;
	}
	return result;
}
float Rastrigin(vector<float> values) { //case 3
	float result = 10 * values.size();
	float sum = 0;
	for (int i = 0; i < values.size(); ++i) {
		float value = pow(values[i], 2) - (10 * cos(2 * PI * values[i]));
		sum += value;
	}
	return result + sum;
}
float  Michalewicz(vector<float> values) {//case 4
	float result = 0;
	for (int i = 0; i < values.size(); ++i) {
		float value = sin(values[i]) * pow((sin(i * pow(values[i], 2) / PI)), 20);
		result += value;
	}
	return -result;
}


float HC_firstImprovment(float a, float b, int n, float (*function)(vector<float>)) {
	int l = compute_l(a, b, prec);
	int L = n * l;
	vector<char> bits(L);
	GenerateRandomBits(bits);

	int  t = 0;
	float best = INFINITY;
	do {
		float evC;
		bool local = false;
		GenerateRandomBits(bits);
		vector<float> domain = decode_values(a, b, n, l, bits);
		evC = function(domain);
		while (local == false)
		{
			for (int ii = 0; ii < bits.size(); ++ii)
			{
				float evN;
				bits[ii] = !bits[ii];
				vector<float> localDomain = decode_values(a, b, n, l, bits);
				evN = function(localDomain);
				bits[ii] = !bits[ii];
				if (evN < evC) {
					evC = evN;
					local = true;
					break;
				}
			}
		}
		t++;
		if (evC <= best) {
			best = evC;
		}
	} while (t <= MAX);
	return best;
}


float HC_bestImprovment(float a, float b, int n, float (*function)(vector<float>)) {
	int l = compute_l(a, b, prec);
	int L = n * l;
	vector<char> bits(L);
	GenerateRandomBits(bits);

	int  t = 0;
	float best = INFINITY;
	do {
		float evC;
		bool local = false;
		GenerateRandomBits(bits);
		vector<float> domain = decode_values(a, b, n, l, bits);
		evC = function(domain);
		while (local == false)
		{
			float localBest = INFINITY;
			for (int ii = 0; ii < bits.size(); ++ii)
			{
				float evN;
				bits[ii] = !bits[ii];
				vector<float> localDomain = decode_values(a, b, n, l, bits);
				evN = function(localDomain);
				bits[ii] = !bits[ii];
				if (evN < localBest) {
					localBest = evN;
				}
				//cout << localBest << "\n";
			}
			evC = localBest;
			local = true;
			//cout << "\t\t" << evC;
			break;
		}
		t++;
		if (evC < best) {
			best = evC;
		}
		//cout << best << '\n';
	} while (t < MAX);
	return best;
}

float SA(float a, float b, int n, float (*function)(vector<float>)) {
	float T = 30;
	int l = compute_l(a, b, prec);
	int L = n * l;
	float best = INFINITY;
	float evC;
	float evN;
	do {
		vector<char> bits(L);
		GenerateRandomBits(bits);
		vector<float> domain = decode_values(a, b, n, l, bits);
		evC = function(domain);
		int t = 0;
		do {
			//select a random neghbour:
			int ri = randomFloat(0, bits.size());
			bits[ri] = !bits[ri];
			domain = decode_values(a, b, n, l, bits);
			bits[ri] = !bits[ri];
			evN = function(domain);
			if (evN < evC) {
				evC = evN;
			}
			else if (randomFloat(0, 1) < exp(-abs(evN - evC) / T)) {
				evC = evN;
			}
			t++;
		} while (t < 500);
		T = T * 0.99;
		if (evC < best) best = evC;
		//cout << '\t' << best << '\n';
	} while (T > 0.00000001);
	return best;
}

//mai trebuie sa implementezi si functia care calculeaza minumul functiei;
int main()
{
	//cout << "\nDe Jong: " << DeJong(domain);
	//cout << "\nSchwefel: " << Schwefel(domain);
	//cout << "\nRastrigin: " << Rastrigin(domain);
	//cout << "\nRastrigin: " << Rastrigin(domain);
	//cout << "\nMichalewicz: " << Michalewicz(domain);

	/*
	outDeJong << "DEJONG n=5\n";
	for (int i = 0; i < 30; ++i) {
	auto start = high_resolution_clock::now();
	outDeJong << "The iteration: "<<i;
	outDeJong<<"\n\tWith the result: " << HC_bestImprovment(-5.12, 5.12, 5, &DeJong);
	auto stop = high_resolution_clock::now();

	auto duration = duration_cast<seconds>(stop - start);
	outDeJong<<"\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}
	cout << "\n\n";
	outDeJong << "DEJONG n=10\n";
	for (int i = 0; i < 30; ++i) {
		auto start = high_resolution_clock::now();
		outDeJong << "The iteration: " << i;
		outDeJong << "\n\tWith the result: " << HC_bestImprovment(-5.12, 5.12, 10, &DeJong);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outDeJong << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}

	cout << "\n\n";
	outDeJong << "DEJONG n=30\n";
	for (int i = 0; i < 30; ++i) {
		auto start = high_resolution_clock::now();
		outDeJong << "The iteration: " << i;
		outDeJong << "\n\tWith the result: " << HC_bestImprovment(-5.12, 5.12, 30, &DeJong);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outDeJong << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}


	outSchwefel << "Schwefel n=5\n";
	for (int i = 0; i < 30; ++i) {
		auto start = high_resolution_clock::now();
		outSchwefel << "The iteration: " << i;
		outSchwefel << "\n\tWith the result: " << HC_bestImprovment(-500,500,5,&Schwefel);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outSchwefel << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}
	/*
	cout << "\n\n";
	outSchwefel << "Schwefel n=10\n";
	for (int i = 0; i < 30; ++i) {
		auto start = high_resolution_clock::now();
		outSchwefel << "The iteration: " << i;
		outSchwefel << "\n\tWith the result: " << HC_bestImprovment(-500, 500, 10, &Schwefel);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outSchwefel << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}
	cout << "\n\n";
	outSchwefel << "Schwefel n=30\n";
	for (int i = 0; i < 30; ++i) {
		auto start = high_resolution_clock::now();
		outSchwefel << "The iteration: " << i;
		outSchwefel << "\n\tWith the result: " << HC_bestImprovment(-500, 500, 30, &Schwefel);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outSchwefel << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}


	outRastrigin << "Rastrigin n=5\n";
	for (int i = 0; i < 30; ++i) {
		auto start = high_resolution_clock::now();
		outRastrigin << "The iteration: " << i;
		outRastrigin << "\n\tWith the result: " << HC_bestImprovment(-5.12, 5.12, 5, &Rastrigin);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outRastrigin << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}

	cout << "\n\n";
	outRastrigin << "Rastrigin n=10\n";
	for (int i = 0; i < 30; ++i) {
		auto start = high_resolution_clock::now();
		outRastrigin << "The iteration: " << i;
		outRastrigin << "\n\tWith the result: " << HC_bestImprovment(-5.12, 5.12, 10, &Rastrigin);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outRastrigin << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}
	cout << "\n\n";
	outRastrigin << "Rastrigin n=30\n";
	for (int i = 0; i < 30; ++i) {
		auto start = high_resolution_clock::now();
		outRastrigin << "The iteration: " << i;
		outRastrigin << "\n\tWith the result: " << HC_bestImprovment(-5.12, 5.12, 30, &Rastrigin);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outRastrigin << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}

	outMichalewicz << " Michalewicz n=5\n";
	for (int i = 0; i < 30; ++i) {
		auto start = high_resolution_clock::now();
		outMichalewicz << "The iteration: " << i;
		outMichalewicz << "\n\tWith the result: " << HC_bestImprovment(0,PI,5,&Michalewicz);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outMichalewicz << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}

	cout << "\n\n";
	outMichalewicz << " Michalewicz n=10\n";
	for (int i = 0; i < 30; ++i) {
		auto start = high_resolution_clock::now();
		outMichalewicz << "The iteration: " << i;
		outMichalewicz << "\n\tWith the result: " << HC_bestImprovment(0, PI, 10, &Michalewicz);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outMichalewicz << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}
	cout << "\n\n";
	outMichalewicz << " Michalewicz n=30\n";
	for (int i = 0; i < 30; ++i) {
		auto start = high_resolution_clock::now();
		outMichalewicz << "The iteration: " << i;
		outMichalewicz << "\n\tWith the result: " << HC_bestImprovment(0, PI, 30, &Michalewicz);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outMichalewicz << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}

	cout << "Starting the DEJONG with 5\n";
	outDeJong << "DEJONG n=5\n";
	for (int i = 0; i < 30; ++i) {
		auto start = high_resolution_clock::now();
		outDeJong << "The iteration: " << i;
		outDeJong << "\n\tWith the result: " << SA(-5.12, 5.12, 5, &DeJong);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outDeJong << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}
	cout << "Starting the DEJONG with 10\n";
	outDeJong << "DEJONG n=10\n";
	for (int i = 0; i < 30; ++i) {
		auto start = high_resolution_clock::now();
		outDeJong << "The iteration: " << i;
		outDeJong << "\n\tWith the result: " << SA(-5.12, 5.12, 10, &DeJong);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outDeJong << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}

	cout << "Starting the DEJONG with 30\n";
	outDeJong << "DEJONG n=30\n";
	for (int i = 0; i < 30; ++i) {
		auto start = high_resolution_clock::now();
		outDeJong << "The iteration: " << i;
		outDeJong << "\n\tWith the result: " << SA(-5.12, 5.12, 30, &DeJong);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outDeJong << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}
	cout << "Starting the Schwefel with 5\n";
	outSchwefel << "Schwefel n=5\n";
	for (int i = 0; i < 30; ++i) {
		auto start = high_resolution_clock::now();
		outSchwefel << "The iteration: " << i;
		outSchwefel << "\n\tWith the result: " << SA(-500, 500, 5, &Schwefel);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outSchwefel << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}

	cout << "Starting the Schwefel with 10\n";
	outSchwefel << "Schwefel n=10\n";
	for (int i = 0; i < 30; ++i) {
		auto start = high_resolution_clock::now();
		outSchwefel << "The iteration: " << i;
		outSchwefel << "\n\tWith the result: " << SA(-500, 500, 10, &Schwefel);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outSchwefel << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}
	cout << "Starting the Schwefel with 30\n";
	outSchwefel << "Schwefel n=30\n";
	for (int i = 0; i < 30; ++i) {
		auto start = high_resolution_clock::now();
		outSchwefel << "The iteration: " << i;
		outSchwefel << "\n\tWith the result: " << SA(-500, 500, 30, &Schwefel);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outSchwefel << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}
	*/
	cout << "Rastigrin n = 5\n";
	outRastrigin << "Rastrigin n=5\n";
	for (int i = 0; i < 30; ++i) {
		cout <<"iteratia: "<< i << '\n';
		auto start = high_resolution_clock::now();
		outRastrigin << "The iteration: " << i;
		outRastrigin << "\n\tWith the result: " << SA(-5.12, 5.12, 5, &Rastrigin);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outRastrigin << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}

	cout << "Rastigrin n = 10\n";
	outRastrigin << "Rastrigin n=10\n";
	for (int i = 0; i < 30; ++i) {
		cout << "iteratia: " << i << '\n';
		auto start = high_resolution_clock::now();
		outRastrigin << "The iteration: " << i;
		outRastrigin << "\n\tWith the result: " << SA(-5.12, 5.12, 10, &Rastrigin);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outRastrigin << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}
	cout << "Rastigrin n = 10\n";
	outRastrigin << "Rastrigin n=30\n";
	for (int i = 0; i < 30; ++i) {
		cout << "iteratia: " << i << '\n';
		auto start = high_resolution_clock::now();
		outRastrigin << "The iteration: " << i;
		outRastrigin << "\n\tWith the result: " << SA(-5.12, 5.12, 30, &Rastrigin);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outRastrigin << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}

	cout << " Michalewicz n=5\n";
	outMichalewicz << " Michalewicz n=5\n";
	for (int i = 0; i < 30; ++i) {
		cout << "iteratia: " << i << '\n';
		auto start = high_resolution_clock::now();
		outMichalewicz << "The iteration: " << i;
		outMichalewicz << "\n\tWith the result: " << SA(0, PI, 5, &Michalewicz);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outMichalewicz << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}

	cout << " Michalewicz n=10\n";
	outMichalewicz << " Michalewicz n=10\n";
	for (int i = 0; i < 30; ++i) {
		cout << "iteratia: " << i << '\n';
		auto start = high_resolution_clock::now();
		outMichalewicz << "The iteration: " << i;
		outMichalewicz << "\n\tWith the result: " << SA(0, PI, 10, &Michalewicz);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outMichalewicz << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}
	cout<< " Michalewicz n=30\n";
	outMichalewicz << " Michalewicz n=30\n";
	for (int i = 0; i < 30; ++i) {
		cout << "iteratia: " << i << '\n';
		auto start = high_resolution_clock::now();
		outMichalewicz << "The iteration: " << i;
		outMichalewicz << "\n\tWith the result: " << SA(0, PI, 30, &Michalewicz);
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		outMichalewicz << "\n\t\t\ttooked: " << duration.count() << " seconds\n";
	}
}