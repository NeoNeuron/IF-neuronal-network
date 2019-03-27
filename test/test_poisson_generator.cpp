// test file of poisson generators;
#include "../include/poisson_generator.h"
using namespace std;

mt19937 rand_gen(1);

int main () {
	PoissonGenerator pg;
	cout << ">> Created pg..." << endl;
	pg.SetRate(2.0);
	pg.SetStrength(0.618);
	pg.SetOuput("./tmp/pg_text.csv");
	cout << ">> Setting accomplished..." << endl;
	queue<Spike> container;
	pg.GenerateNewPoisson(rand_gen, 100, container);
	cout << ">> Poisson sequence generated." << endl;
	cout << ">> Done." << endl;
	return 0;
}
