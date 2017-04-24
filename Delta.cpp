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

using namespace std;

#define BMMRAND (double)rand()/RAND_MAX;

////Simulation Setup////
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

////End Simulation Setup////

int main() {

	///Simulation///
	double u = 0; //set to zero for straight line movement
	ship s;
	s.init();
	s.updatepos(u);
	s.updatepos(u);
	s.updatepos(u);
	
	for (int i = 0; i < s.x.size(); i++) {
		cout << "x:\t" << s.x[i] << "\ty:\t" << s.y[i] << "\ttheta\t" << s.theta[i] << "\tomega\t" << s.omega[i] << endl;
	}
	///End Simulation///

	return 0;
}