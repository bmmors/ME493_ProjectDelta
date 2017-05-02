//Mors, Britny
//ME493 - Autonomy
//Project Delta - Nonlinear Evolutionary Control

#include <iostream>
#include <assert.h>
#include <vector>
#include <random>
#include <time.h>
#include <stdio.h>
#include <cmath>
#include <time.h>
#include <fstream>
#include "LY_NN.h"


using namespace std;

#define BMMRAND (double)rand()/RAND_MAX;

////----------Simulation Setup----------////
class ship {
public:
	vector<double> x;
	vector<double> y;
	vector<double> theta;
	vector<double> omega;
	double dt = 0.2;
	double v = 3;
	double T = 5;

	void init();
	void updatepos(double u);
	void reset();
};

void ship::init() {
	assert(x.size() == 0);
	assert(y.size() == 0);
	assert(omega.size() == 0);
	assert(theta.size() == 0);

	x.push_back(150);
	y.push_back (150);
	theta.push_back(90);
	omega.push_back(0);
}

void ship::updatepos(double u) {
	x.push_back(x.at(x.size() - 1) + v*sin(theta.at(theta.size() - 2))*dt); //x(t+1)=x(t) + v*sin(theta(t))*dt
	y.push_back(y.at(y.size() - 1) + v*cos(theta.at(theta.size() - 2))*dt); //y(t+1)=y(t) + v*cos(theta(t))*dt
	omega.push_back(omega.at(omega.size()-1) + (u - omega.at(omega.size()-1))*dt / T); // w(t+1)=w(t) + (u-w(t))*dt/T
	theta.push_back(theta.at(theta.size() - 1) + omega.at(omega.size() - 2) / dt); //theta(t+1)=theta(t)+w(t)*dt
}

void ship::reset() {
	omega.clear();
	theta.clear();
	x.clear();
	y.clear();
	cout << "Reset!" << endl;
	init();
}

////----------End Simulation Setup----------////




////----------Goal Setup----------////
class goal {
public:
	///To be randomized later
	vector<double> x;
	vector<double> y;

	void init();
};

void goal::init() {
	x.push_back(900);
	x.push_back(950);

	y.push_back(900);
	y.push_back(900);

	//makes goal start at (900,900) and go to (950,900)
}
////----------End Goal Setup----------////



////----------Policy Setup----------////
class policy {
public:
	vector<double> weights;
	double fitness;

	void init(int num_weights);
	void mutate(double mut_size);
};

void policy::init(int num_weights) {
	for (int i = 0; i < num_weights; i++) {
		double t1 = BMMRAND;
		double t2 = BMMRAND;
		weights.push_back(t1 - t2);
	}
	assert(weights.size() == num_weights);
}

void policy::mutate(double mut_size) {
	int mutate = rand() % weights.size(); //selects how many weights to change
	for (int i = 0; i < mutate; i++) {
		int select = rand() % weights.size();
		weights.at(select) += mut_size * BMMRAND - mut_size * BMMRAND;
	}

}
////----------End Policy Setup----------////




////----------Evolutionary Algorithm Setup----------////
class EA {
public:
	vector<policy> population;

	void replicate(int pop_size,double mut_size);
	void evaluate();
	void downselect();
};

void EA::replicate(int pop_size, double mut_size) {
	assert(population.size() == pop_size);
	int new_pop = pop_size * 2;
	policy temp;
	while (population.size() < new_pop) {
		int select = rand() % pop_size; //select random population
		temp = population.at(select); 
		temp.mutate(mut_size);
		population.push_back(temp);
		//cout << "pop_size:\t" << population.size() << endl;
	}
	assert(population.size() == new_pop);
}
////----------End Evolutionary Algorithm----------////




int main() {
	srand(time(NULL));

	////----------Simulation----------////
	ship s;
	s.init();
	////----------End Simulation----------////


	////----------Evolution Algorithm----------////
	EA e;
	int num_pop = 1;
	int num_weights = 1; //to be determined using NN later
	double mut_size = 0.2;
	////----------End Evolutionary Algorithm----------////


	////----------NN----------////

	////----------End NN----------////


	////----------Start Full Sim with EA and NN----------////
	int max_time = 10; //max number of time for each simulation
	int gen = 1; //number of generations
	int SR = 1; //stat runs

	vector<double> u;
	u.push_back(0.45);
	u.push_back(-0.736);

	cout << "\tx:\t" << s.x[0] << "\ty:\t" << s.y[0] << "\tomega\t" << s.omega[0] << "\ttheta\t" << s.theta[0] << endl;
	s.updatepos(-0.123);
	cout << "\tx:\t" << s.x[1] << "\ty:\t" << s.y[1] << "\tomega\t" << s.omega[1] << "\ttheta\t" << s.theta[1] << endl;
	s.updatepos(-0.123);
	cout << "\tx:\t" << s.x[2] << "\ty:\t" << s.y[2] << "\tomega\t" << s.omega[2] << "\ttheta\t" << s.theta[2] << endl;

	
	////----------End Full Sim with EA and NN----------////
	return 0;
}