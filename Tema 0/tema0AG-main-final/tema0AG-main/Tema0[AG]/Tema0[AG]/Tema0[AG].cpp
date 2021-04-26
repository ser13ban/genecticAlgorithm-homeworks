#include <iostream>
#include <random>
#include <iomanip>
#include <chrono> 
using namespace std;
using namespace std::chrono;

//implemetnig a random funciton that can take a interval containing negative numbers
double randomFloat(const double lower, const double upper) {
    if (lower > upper) {
        throw "ERROR: the second argument of randomFloat must be bigger than the first argument!";
    }
    std::random_device generator;
    std:: uniform_real_distribution<double> distribution(lower, upper);
    return distribution(generator);
}


//implemeting the booths function;
long double pow2(long double x) {
    return x * x;
}
long double booth(float x1, float x2) {
    return pow2(x1+2*x2-7)+pow2(2*x1+x2-5);
}
//utility funtion miminum:
long double minimum(long double x, long double y) {
    if (x < y) {
        return x;
    }
    else return y;
}

long double FindTheMinimun_DETERMINSTIC(const long double epsilon) {
    long double  fm=1000000;
    for (long double x=-10; x <= 10; x += epsilon) {
        for (long double y = -10; y <= 10; y+=epsilon) {
            fm = minimum(fm, booth(x, y));
        }
    }
    return fm;
}
long double FindTheMinimun_EURISTIC(const long double esantion) {
    float x, y;
    long double fm=INFINITY;
    for (int i = 0; i <= esantion; ++i) {
        x = randomFloat(-10, 10);
        y = randomFloat(-10, 10);
        fm = minimum(fm, booth(x,y));
    }
    return fm;
}

int main()
{
    auto start = high_resolution_clock::now();
    cout <<"The minimum of booth function, the deterministic algorithm: "<< FindTheMinimun_DETERMINSTIC(0.0001)<<'\n';;
    //cout << "The minimum of booth function, the euristic algorithm:" << FindTheMinimun_EURISTIC(400000000)<<'\n';
    auto stop = high_resolution_clock::now();
  
    auto duration = duration_cast<seconds>(stop - start);
    cout << "Time taken by function: "
        << duration.count() << " seconds" << endl;
}