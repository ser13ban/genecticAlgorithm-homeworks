#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

vector<char> myarr = { 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0 };

unsigned int l = 4;
float a = -5.12, b = 5.12;

//
//
//auto vec_to_int(vector<char>::iterator beg, vector<char>::iterator end, unsigned const int base = 2) {
//	//agl care decodifica
//	int b = 0;
//	for (auto ii = beg; ii != end; ii++) {
//		b *= base;
//		b += *ii;
//		cout << b << ' ' << int(*ii) << endl;
//	}//ok
//	return b;
//}
//


auto vec_to_float(
	vector<char>::iterator beg, vector<char>::iterator end, 
	const unsigned int l,
	unsigned const int base = 2
	) 
{
	//agl care decodifica
	int q = 0;
	for (auto ii = beg; ii != end; ii++) {
		q *= base;
		q += *ii;
		cout << q << ' ' << int(*ii) << endl;
	}//ok
	float s = q / (pow(base, l) - 1);
	float r = q * (b - a) + a;
	return r;
}


vector<float> decode_dimensions(vector<char>&arr, unsigned const int dim_len) {
	unsigned int dim_nr = floor(arr.size() / dim_len);
	vector<float> return_val;

	for (int ii = 0; ii < dim_nr; ++ii) {
		unsigned int start = ii * dim_len;
		unsigned int end = ii * dim_len + dim_len;
		return_val.push_back(
			vec_to_float(arr.begin() + start, arr.begin() + end, dim_len)
		);
	}
	return return_val;
}




int main()
{
	//cout << vec_to_int(myarr) << endl;
	//cout << vec_to_int(myarr1) << endl;



	//- "fAcI tU dEbuG"
	//- thx :/
	vector<float> v = decode_dimensions(myarr, l);
	for (const auto ii : v) {
		cout << ii << " ";
	}
	cout << endl;
	system("pause");
}