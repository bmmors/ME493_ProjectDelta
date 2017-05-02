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
	double x_start = 1;
	double y_start = 1;
	double boatx;
	double boaty;
	double w; //omega
	double start_w = 10;
	double start_o = 1.578;
	double o; //theta
	double dt = 0.2;
	double v = 3;
	double T = 5;

	double x; //y-1/b=x
	double m; // y2-y1/x2-x1
	double b; //boaty-m*boatx
	vector<double> state;

	void init();
	void updatepos(int u);
};

void ship::init() {
	boatx = x_start;
	boaty = y_start;
	w = start_w;
	o = start_o;

	state.push_back(boatx);
	state.push_back(boaty);
	state.push_back(o);
 }

void ship::updatepos(int u) {
	double tx;
	double ty;
	double to;
	double tw;

	tx = boatx + v*sin(o)*dt;
	ty = boaty + v*cos(o)*dt;
	to = o + (w*dt);
	tw = w + ((u - w)*dt) / T;

	boatx = tx;
	boaty = ty;
	o = to;
	w = tw;


	cout << boatx << "," << boaty << endl;
	//x(t+1)=x(t) + v*sin(theta(t))*dt
	//y(t+1)=y(t) + v*cos(theta(t))*dt
	//w(t+1)=w(t) + (u-w(t))*dt/T
	//theta(t+1)=theta(t)+w(t)*dt
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
	void mutate(double mm);
};

void policy::init(int num_weights) {
	for (int i = 0; i < num_weights; i++) {
		double w1 = BMMRAND;
		double w2 = BMMRAND;
		weights.push_back(w1-w2);
	}
	assert(weights.size() == num_weights);
}

void policy::mutate(double mm) {
	int num_mutate = rand() % 5; //mutate up to 5 weights?
	for (int i = 0; i < num_mutate; i++) {
		int index = rand() % weights.size();
		double nweight = mm*BMMRAND - mm*BMMRAND;
		weights.at(index) += nweight;
	}
}

////----------End Policy Setup----------////




////----------Evolutionary Algorithm Setup----------////
class EA {
public:
	vector<policy> population;
	double u;
	bool goal_pass = false;

	void replicate(int num_pop,double mm);
	void evaluate(int num_pop, int max_time,ship s,neural_network NN,goal g);
	void downselect();
};

void EA::replicate(int num_pop,double mm) {
	assert(population.size() == num_pop);
	while (population.size() < num_pop * 2) {
		policy temp;
		int select = rand() % population.size();
		temp = population.at(select);
		temp.mutate(mm);
		population.push_back(temp);
	}
}

void EA::evaluate(int num_pop,int max_time,ship s,neural_network NN, goal g) {
	for (int i = 0; i < population.size();i++) {
		population.at(i).fitness = -1;
	}
	for (int k = 0; k < population.size(); k++) {
		NN.set_weights(population.at(k).weights, true);
		//simulation loop 
		for (int sim = 0; sim < max_time; sim++) {
			NN.set_vector_input(s.state);
			NN.execute();
			u = NN.get_output(0);
			cout << "u:" << u << endl;
			s.updatepos(u);

			if (sim == max_time - 1 && goal_pass == false) {
				population.at(k).fitness += -100;
			}
			else if (s.boatx > 1000 || s.boatx < 0) {
				population.at(k).fitness += -100;
			}
			else if (s.boaty > 1000 || s.boaty < 0) {
				population.at(k).fitness += -100;
			}
			else if (goal_pass = true) {
				population.at(k).fitness += 10000;
			}
			//if passes through goal +10000 fitness
			//if sim == max_time - 1 && goal == false -100
			//if x > 1000 fitness -100 || if x < 0 -100
			//if y > 1000 fitness -100 || if y < 0 -100

		}
	}
}
////----------End Evolutionary Algorithm----------////




int main() {
	srand(time(NULL));

	///Simulation///
	ship s;
	s.init();
	goal g;
	g.init();
	///End Simulation///


	///Evolution Algorithm///
	EA e;
	int num_pop = 50;
	double mm = 0.2;
	///End Evolutionary Algorithm///

	///Neural Network///
	neural_network NN;
	NN.setup(3, 5, 1);
	int num_weights = NN.get_number_of_weights();
	NN.set_in_min_max(0,1000); //x
	NN.set_in_min_max(0,1000); //y
	NN.set_in_min_max(0,6.28); //theta

	NN.set_out_min_max(-15,15);
	///End Neural Network


	///Start Full Sim with EA and NN///
	int max_time = 10; //max number of time for each simulation
	int gen = 100; //number of generations
	int SR = 3; //stat runs

	while (e.population.size() < num_pop) {
		policy p;
		p.init(num_pop);
		e.population.push_back(p);
		//cout << e.population.size() << endl;
	}

	e.replicate(num_pop, mm);
	e.evaluate(num_pop, max_time, s, NN,g);

	///End Full Sim with EA and NN///
	return 0;
}