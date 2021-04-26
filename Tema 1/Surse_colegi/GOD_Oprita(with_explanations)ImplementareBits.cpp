#include <vector>
#include <cmath>
int nbParameters = 4;
double a = -5.12, b = 5.12;
int precision = 5;
using namespace std;
int bitsPerComponent;

vector<char> bits;

//de pe pagina profului
int CalculateBitsPerComponent(double a, double b, int nbParameters, int precision)
{
	int power = 1;
	for (int i = 1; i <= precision; i++)
		power *= 10;

	return ceil(log2((b - a) * power));
}

//spawn intr o pozitie oarecare
void GenerateRandomBits(vector<char> &v)
{
	for (auto ii = v.begin(); ii != v.end(); ++ii)
	{
		//aici implementati cum vreti alegerea random
		//eu fac cu ce am :))
		*ii = (randomFloat(0, 1) < 0.5) ? 0 : 1;
	}
}

vector<double> ConvertToDouble(vector<char> &bits, double a, double b, int BitsPerComponent)
{
	double power2 = 1; //2^ nr bits per component
	for (int i = 1; i <= bitsPerComponent; i++)
		power2 *= 2;

	vector<double> Converted(bits.size() / bitsPerComponent);
	cout << "Hai sa convertim! \n";
	for (int i = 0; i < Converted.size(); i++)
	{
		int Xint = 0;
		for (int j = i * bitsPerComponent; j < (i + 1) * (bitsPerComponent); j++)
		{
			cout << int(bits[j]);
			Xint *= 2;
			Xint += bits[j];
		}
		cout << " ";
		cout << Xint << "\n";
		Converted[i] = a + Xint * (b - a) / (power2 - 1);
	}
	cout << "\n\n";
	return Converted;
}

int main()
{
	bitsPerComponent = CalculateBitsPerComponent(a, b, nbParameters, precision);
	bits = vector<char>(bitsPerComponent * nbParameters);
	GenerateRandomBits(bits);
	int counter = 0;
	for (int i = 0; i < bits.size(); i++)
	{
		if (counter == bitsPerComponent)
		{
			counter = 0;
			cout << "\n";
		}

		cout << int(bits[i]);
		counter++;
	}

	cout << "\n";
	vector<double> Converted = ConvertToDouble(bits, a, b, bitsPerComponent);
	for (auto a : Converted)
	{
		cout << a << "\n";
	}

	system("PAUSE");
}
