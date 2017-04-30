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
	void updatepos(int u);
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

void ship::updatepos(int u) {
	omega.push_back(omega.at(omega.size()-1) + (u - omega.at(omega.size()-1))*dt / T); // w(t+1)=w(t) + (u-w(t))*dt/T
	theta.push_back(theta.at(theta.size() - 1) + omega.at(omega.size() - 2) / dt); //theta(t+1)=theta(t)+w(t)*dt
	x.push_back(x.at(x.size() - 1) + v*sin(theta.at(theta.size() - 2))*dt); //x(t+1)=x(t) + v*sin(theta(t))*dt
	y.push_back(y.at(y.size() - 1) + v*cos(theta.at(theta.size() - 2))*dt); //y(t+1)=y(t) + v*cos(theta(t))*dt
}

void ship::reset() {
	omega.clear();
	theta.clear();
	x.clear();
	y.clear();

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
	double u;

	void init(int num_weights);
	void mutate();
};

void policy::init(int num_weights) {
	for (int i = 0; i < num_weights; i++) {
		int test = rand() % 2; //generate 0 or 1. Zero means + weight 1 means 0 weight
		assert(test == 0 || test == 1);
		if (test == 0) {
			double w = BMMRAND
			weights.push_back(w);
		}
		else if (test == 1) {
			double w = BMMRAND;
			w = w * -1;
			weights.push_back(w);
		}
	}
	assert(weights.size() == num_weights);
}
////----------End Policy Setup----------////




////----------Evolutionary Algorithm Setup----------////
class EA {
public:
	vector<policy> population;

	void replicate();
	void evaluate();
	void downselect();
};

////----------End Evolutionary Algorithm----------////




int main() {

	///Simulation///
	ship s;
	s.init();

	//cout << "x:\t" << s.x[i] << "\ty:\t" << s.y[i] << "\ttheta\t" << s.theta[i] << "\tomega\t" << s.omega[i] << endl;
	///End Simulation///

	///Evolution Algorithm///
	EA e;
	int num_poly = 50;
	int num_weights = 10; //to be determined using NN later

	for (int i = 0; i < num_poly; i++) {
		policy p;
		p.init(num_weights);
		e.population.push_back(p);
	}
	assert(e.population.size() == num_poly);
	///End Evolutionary Algorithm///

	///Start Full Sim with EA and NN///
	int max_time = 1000; //max number of time for each simulation
	///End Full Sim with EA and NN///
	return 0;
}