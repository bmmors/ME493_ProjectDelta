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


////----------Goal Setup----------////
class goal {
public:
	///To be randomized later
	vector<double> x;
	vector<double> y;
	vector<double> mid;

	void init();
};

void goal::init() {
	x.push_back(900);
	x.push_back(950);

	y.push_back(900);
	y.push_back(900);

	mid.push_back((x.at(0) + x.at(1)) / 2);
	mid.push_back((y.at(0) + y.at(1)) / 2);

	//makes goal start at (900,900) and go to (950,900)
}
////----------End Goal Setup----------////

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
	bool goal_pass;
	double fitness;

	double x; //y-1/b=x
	double m; // y2-y1/x2-x1
	double b; //boaty-m*boatx
	vector<double> state;

	void init(goal g);
	void updatepos(int u,goal g);
};

void ship::init(goal g) {
	boatx = x_start;
	boaty = y_start;
	w = start_w;
	o = start_o;
	goal_pass = false;

	state.push_back(boatx);
	state.push_back(boaty);
	state.push_back(o);
	fitness = sqrt(((boatx - g.mid.at(0))*(boatx - g.mid.at(0))) + ((boaty - g.mid.at(1))*(boaty - g.mid.at(1)))); //calc distance from boat to midpoint of goal.
 }

void ship::updatepos(int u,goal g) {
	double tx;
	double ty;
	double to;
	double tw;

	tx = boatx + v*sin(o)*dt;
	ty = boaty + v*cos(o)*dt;
	to = o + (w*dt);
	tw = w + ((u - w)*dt) / T;

	m = (ty - boaty) / (tx - boatx);
	b = (ty - m*tx);
	x = (900 - b) / m;

	//minimize fitness
	fitness = fitness - sqrt(((tx - g.mid.at(0))*(tx - g.mid.at(0))) + ((ty - g.mid.at(1))*(ty - g.mid.at(1))));
	if (x >= 900 && x <= 950) {
		if (boatx < 900 && ty > 900) { //crossing goal from bottom to top
			//goal passed 
			fitness += -100; //since we are minimizing fitness give it a negative reward for finding the goal
		}
		if (boatx > 900 && ty <900) { //crossing goal from top to bottom
			//goal passed 
			fitness += -100; //since we are minimizing fitness give it a negative reward for finding the goal
		}
	}
	if (tx > 1000 || tx < 0) {
		fitness += 100; //out of bounds = not ideal soln
	}
	if (ty > 1000 || ty < 0) {
		fitness += 100; //out of bounds = not ideal soln
	}
	
	boatx = tx;
	boaty = ty;
	o = to;
	w = tw;

	//cout << boatx << "," << boaty << endl;
	//x(t+1)=x(t) + v*sin(theta(t))*dt
	//y(t+1)=y(t) + v*cos(theta(t))*dt
	//w(t+1)=w(t) + (u-w(t))*dt/T
	//theta(t+1)=theta(t)+w(t)*dt
}

////----------End Simulation Setup----------////




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

	void replicate(int num_pop,double mm);
	void evaluate(int num_pop, int max_time,ship s,neural_network NN,goal g);
	vector<policy> downselect(int num_pop,int mm);
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
	int sim_count = 0;
	assert(population.size() == num_pop * 2);
	for (int k = 0; k < population.size(); k++) {
		NN.set_weights(population.at(k).weights, true);
		//simulation loop 
		for (int sim = 0; sim < max_time; sim++) {
			NN.set_vector_input(s.state);
			NN.execute();
			u = NN.get_output(0);
			//cout << "u:" << u << endl;
			s.updatepos(u,g);
			//cout << sim << endl;
			sim_count++;
		}
		//cout << "sim count:" << sim_count <<endl;
	}
}

vector<policy> EA::downselect(int num_pop,int mm) {
	vector<policy> new_pop;
	while (new_pop.size() < num_pop*2){
		int rand1 = rand() % population.size();
		int rand2 = rand() % population.size();
		if (rand1 == rand2) {
			int rand1 = rand() % population.size();
		}

		double fit1 = population.at(rand1).fitness;
		double fit2 = population.at(rand2).fitness;

		//while (fit1 == fit2) {
			//population.at(rand1).mutate(mm);
			//fit1 = population.at(rand1).fitness;
		//}
		cout << fit1 << "\t" << fit2 << endl;

		if (fit1 < fit2) {
			//fit 1 wins
			new_pop.push_back(population.at(rand1));
		}
		if (fit2 < fit1) {
			//fit 2 wins
			new_pop.push_back(population.at(rand1));
		}
		//cout << new_pop.size() << endl;
		
	}
	assert(new_pop.size() == num_pop);

	return new_pop;
}
////----------End Evolutionary Algorithm----------////




int main() {
	srand(time(NULL));

	///Simulation///
	goal g;
	g.init();
	ship s;
	s.init(g);
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
	int max_time = 100; //max number of time for each simulation
	int gen = 100; //number of generations
	int SR = 1; //stat runs

	for (int i = 0; i < SR; i++) {
		assert(e.population.size() == 0);
		while (e.population.size() < num_pop) {
			policy p;
			p.init(num_weights);
			e.population.push_back(p);
			//cout << e.population.size() << endl;
		}
		cout << "exit population loop" << endl;
		for (int k = 0; k < gen; k++) {
			e.replicate(num_pop, mm);
			e.evaluate(num_pop, max_time, s, NN, g);
			e.population = e.downselect(num_pop,mm);
		}
		e.population.clear();
		assert(e.population.size() == 0);
	}

	///End Full Sim with EA and NN///
	return 0;
}