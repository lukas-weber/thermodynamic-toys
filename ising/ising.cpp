// compile 
#include <cstdint>
#include <ctime>
#include <sstream>
#include <fstream>
#include <iostream>
#include <random>

class mc {
private:
	std::mt19937 gen; // random numbers
	std::uniform_real_distribution<double> rand;
	
	int L;
	double T;
	std::vector<int8_t> spins;

	std::vector<double> absmag;

public:
	mc(int L, double T);
	void sweep();
	void measure();

	void get_results(double &amag, double &stdamag, int binsize=3000);
	void print_config(const char *filename);
};


mc::mc(int L, double T) : gen(time(0)), rand(0,1), L(L), T(T) {
	spins.resize(L*L);
	absmag.reserve(10000);

	for(int i = 0; i < L*L; i++) {
		spins[i] = 1-2*(rand(gen)>0.5);
	}
}

void mc::sweep() {
	for(int y = 0; y < L; y++) {
		for(int x = 0; x < L; x++) {
			double dE = 2*spins[y*L+x]*(
					spins[y*L+((x+1)%L)]
					+spins[((y+1)%L)*L+x]
					+spins[y*L+((x-1+L)%L)]
					+spins[((y-1+L)%L)*L+x]);
			if(dE < 0 || rand(gen) < exp(-dE/T)) {
				spins[y*L+x]*=-1;
			}
		}
	}
}

void mc::measure() {
	double mag = 0;
	for(int i = 0; i < L*L; i++) {
	       mag += spins[i];
	}
	
	//absmag.push_back(fabs(mag)/L/L);
	absmag.push_back(mag/L/L*mag/L/L);
}

void mc::get_results(double &amag, double &stdamag, int binsize) {
	std::vector<double> bins(absmag.size()/binsize + 1, 0);
	std::vector<double> counts(absmag.size()/binsize + 1, 0);
	for(size_t i = 0; i < absmag.size(); i++) {
		bins[i/binsize] += absmag[i];
		counts[i/binsize]++;
	}
	amag = 0;
	double amag2 = 0;
	for(size_t i = 0; i < bins.size(); i++) {
		amag += bins[i]/counts[i];
		amag2 += bins[i]*bins[i]/counts[i]/counts[i];
	}
	amag/=bins.size();
	amag2/=bins.size();
	
	stdamag = sqrt(1./(bins.size()-1)*(amag2-amag*amag));
}

void mc::print_config(const char *filename) {
	std::ofstream file(filename);
	for(int y = 0; y < L; y++) {
		for(int x = 0; x < L; x++) {
			file << (spins[y*L+x]>0) << " ";
		}
		file << "\n";
	}
	file.close();
}


int main(int argc, char **argv) {
	const int thermsweeps = 2000;
	const int sweeps = 100000;

	const std::vector<int> Ls = {10,20,50,100,150};

	const double T0 = 2;
	const double T1 = 3;
	const int Tsteps = 30;

	for(auto L : Ls) {
		std::vector<double> T(Tsteps);
		std::vector<double> amag(Tsteps);
		std::vector<double> stdamag(Tsteps);

#pragma omp parallel for
		for(int t = 0; t < Tsteps; t++) {
			T[t] = T0 + (T1-T0)/(Tsteps-1)*t;
			mc ising(L,T[t]);

			for(int i = 0; i < thermsweeps/T[t]; i++) {
				ising.sweep();
			}

			for(int i = 0; i < sweeps; i++) {
				ising.sweep();
				ising.measure();
			}

			// if(L == Ls.back()) {
			// 	std::stringstream fname;
			// 	fname << "config_T=" << T[t] << ".csv";
			// 	ising.print_config(fname.str().c_str());
			// }

			ising.get_results(amag[t],stdamag[t]);
			std::cout << "L = " << L << ", T = " << T[t] << " done.\n";
		}

		std::stringstream s;
		s << "data_L=" << L << ".csv";
		std::ofstream out(s.str(), std::ofstream::out);
		out << "T\tAbsolute Magnetization\tError\n";
		for(int t = 0; t < Tsteps; t++) {
			out << T[t] << "\t" << amag[t] << "\t" << stdamag[t] << "\n";
		}
		out.close();
	}
}
