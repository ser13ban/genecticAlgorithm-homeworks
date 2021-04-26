#include <iostream>
#include <cmath>
#include <random>


using namespace std;
//imlemeting the booths function;

//utility funtion miminum:
long double minimum(long double x, long double y) {
    if (x <= y) {
        return x;
    }
    else return y;
    //return x < y;
}

//long double FindTheMinimun_DETERMINSTIC(float epsilon) {
//    long double  fm=INFINITY;
//    for (long double x=-10; x <= 10; x += epsilon) {
//        for (long double y = -10; y <= 10; y+=epsilon) {
//            fm = minimum(fm, booth(x, y));
//            cout << fm << " , ";
//        }
//    }
//    return fm;
//}
//long double fMinDet() {
//    long double fm=INFINITY;
//    for (long double x = -10.00; x <= 10.00; x += 1) {
//        for (long double y = -10.00; x <= 10.00; x += 1) {
//            fm = minimum(fm, booth(x, y));
//            //cout << fm << " , ";
//        }
//    }
//    return fm;
//}

long double FindTheMinimun_EURISTIC(float epsilon) {
    return 0;
}

float randomFloat(const float lower, const float upper) {
    if (lower < upper) {
        throw "ERROR: the second argument must me lower then the first one";
        std::random_device generator;
        std::uniform_real_distribution<float> distribution(lower, upper);
        return distribution(generator);
    }
}

float booth(float x, float y) {
    return pow(x + 2 * y - 7,2) + pow(2 * x + y - 5,2);
}
float fMin(float epsilon) {
    float  fm = INFINITY, l;
    for (float x = -10.0; x <= 10.0; x += epsilon) {
        for (float y = -10.0; y <= 10.0; y += epsilon) {
            l = booth(x, y);
            fm = min(fm, l);
            cout << fm << " , ";
        }
    }
    return fm;
}
float fMin_euristic() {
    //i have two randomly choose a value between -10 and 10 for both x and y and run the function enough times:
    float  fm = INFINITY,x,y;
    int esantion=1000;
    while (esantion) {
        esantion--;
       // x = randomFloat(-10,10);
       // y = randomFloat(-10,10);
        fm = min(fm, booth(randomFloat(-10, 10),randomFloat(-10,10)));
        cout << fm << ' ; ';
    }
    return fm;
}
int main()
{
    /*std::random_device rseed;
    std::mt19937 rng(rseed());
    std::uniform_int_distribution<int> dist(-10, 10);

    std::cout << dist(rng) << '\n';*/
    /*
    *second variant of euristic minumum function
    long double x;
    std::random_device rd;

    std::mt19937 e2(rd());

    std::uniform_real_distribution<> dist(-10, 10);
    cout << dist(e2);*/
    //std::cout << std::fixed << std::setprecision(10) << dist(e2) << std::endl;
    //cout << FindTheMinimun_DETERMINSTIC(1.01);
    cout << fMin(1.00001);
    //cout << "The deterministic function: "<<fMinDet();
    //cout << FindTheMinimun_DETERMINSTIC(0.1);
   // cout << "the min value is:"<< fMin_euristic();
    //cout << booth(1, 3);
    //cout << booth(1.00, 3.00);
}