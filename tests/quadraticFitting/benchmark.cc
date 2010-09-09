#include <vector>
#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <iomanip>

void getBestFitVelocityAndAcceleration(std::vector<double> positions, 
				       const std::vector<double> & times,
                                       double & velocity, 
				       double &acceleration, 
				       double &position0);

std::vector<double> 
mk_fill_obs_moments(std::vector<double> X, std::vector<double> time, 
		    unsigned int num_coefficients);



bool areEqual(double a, double b, double epsilon)
{
     if (fabs(a - b) < epsilon) {
	  return true;
     }
     return false;
     
}




double projectLoc(double time, double p0, double v, double acc)
{
     return p0 + v*time + .5 * time* time * acc;
}

#define NUM_CALLS 100000
#define CHECK_CORRECTNESS false

int main()
{


     std::vector<double> times;
     times.push_back(0.);
     times.push_back(0.05);
     times.push_back(1.);
     times.push_back(1.05);
     times.push_back(5.);
     times.push_back(5.05);

     std::vector<std::vector< double > > positions;
     std::vector<std::vector<double> > answerKey;


     // build a vector of positions with 
     for (unsigned int i = 0; i < NUM_CALLS; i++ ) {
	  double p0, v, acc;
	  p0 = (i % 100) / 50.;
	  v = (i % 13) / 20. - 7.;
	  acc = (i % 7) / 50. - 3.;
	  std::vector<double> motion;

	  motion.push_back(p0);
	  motion.push_back(v);
	  motion.push_back(acc);
	  answerKey.push_back(motion);

	  std::vector<double> p;

	  for (unsigned int j = 0; j < times.size(); j++) {
	       p.push_back(projectLoc(times[j], p0, v, acc));
	  }
	  positions.push_back(p);
     }

     std::cout << "Starting calls to Kubica-style fitting...\n";

     double startTime;
     double elapsed;
     startTime = std::clock();
     for (unsigned int i = 0; i < positions.size(); i++) {
	  std::vector<double> res = mk_fill_obs_moments(positions[i], times, 3);
	  

	  if (CHECK_CORRECTNESS) {
	       for (unsigned int j = 0; j < res.size(); j++) {
		    if (!areEqual(res[j], answerKey[i][j], 1e-10)) {
			 
			 std::cout << j << " element of " << i << " fitting was " << 
			      answerKey[i][j] << " but we got " << res[j] << std::endl;
			 return -1;
		    }
	       }
	  }
     }

     elapsed = (std::clock() - startTime) / (double)CLOCKS_PER_SEC;
     std::cout << NUM_CALLS << " calls to Kubica-style fitting took " 
	       << std::setprecision(10)
	       << elapsed 
	       << " sec." << std::endl;


     startTime = std::clock();
     for (unsigned int i = 0; i < positions.size(); i++) {
	  double p0, v, acc;
	  getBestFitVelocityAndAcceleration(positions[i], times, v, acc, p0);


	  if (CHECK_CORRECTNESS) {
	       std::vector <double> res;
	       res.push_back(p0);
	       res.push_back(v);
	       res.push_back(acc);
	       
	       for (unsigned int j = 0; j < res.size(); j++) {
		    
		    if (!areEqual(res[j], answerKey[i][j], 1e-10)) {
			 
			 std::cout << j << " element of " << i << " fitting was " << 
			      answerKey[i][j] << " but we got " << res[j] << std::endl;
			 return -1;
		    }
	       }
	  }
     }

     elapsed = (std::clock() - startTime) / (double)CLOCKS_PER_SEC;
     std::cout << NUM_CALLS << " calls to GSL-style fitting took " 
	       << std::setprecision(10)
	       << elapsed 
	       << " sec." << std::endl;

     return 0;

}
