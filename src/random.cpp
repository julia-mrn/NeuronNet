#include "random.h"

RandomNumbers::RandomNumbers(unsigned long int s) 
: seed(s) {
	
		if (s==0){
			std::random_device rd;
			seed = rd();
		};
	
	rng = std::mt19937(seed);
	
	
}


void RandomNumbers::uniform_double(std::vector<double>& res, double lower, double upper){ //attention a ne pas mettre les valeurs par defaut dans le cpp
	
    for (auto I = res.begin(); I != res.end(); I++)
    { *I = uniform_double(lower, upper);};
	
}

double RandomNumbers::uniform_double(double lower, double upper){
	
	std::uniform_real_distribution<> unif(lower, upper);
	return unif(rng);
	
}

void RandomNumbers::normal(std::vector<double>& res, double mean, double sd) {
	
	 
    for (auto I = res.begin(); I != res.end(); I++) 
    {*I = normal(mean, sd);};
	
}

double RandomNumbers::normal(double mean, double sd) {
	
	std::normal_distribution<> norm(mean, sd);
	return norm(rng);
	
}

void RandomNumbers::poisson(std::vector<int>& res, double mean) {
	
	for (auto I = res.begin(); I != res.end(); I++) 
	{*I = poisson(mean);};
	
}

 int RandomNumbers::poisson(double mean){
	 
	 std::poisson_distribution<int> poisson(mean);
	 return poisson(rng);
	 
 }
 


