/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
  //   x, y, theta and their uncertainties from GPS) and all weights to 1. 
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  num_particles = 100;
  std::default_random_engine gen;
 
  // values near the mean are the most likely
  // standard deviation affects the dispersion of generated values from the mean
  std::normal_distribution<> dx{x, std[0]};
  std::normal_distribution<> dy{y, std[1]};
  std::normal_distribution<> dth{theta, std[2]};
 
  for(int p = 0; p < num_particles; p++)
  {
    double xp = dx(gen);
    double yp = dy(gen);
    double tp = dth(gen);
    double wtp = 1.0;
//    cout << p << " x " << xp <<  " y " << yp <<  " theta " << tp << endl;
    Particle par;
    par.id = p;
    par.x = xp;
    par.y = yp;
    par.theta = tp;
    par.weight = wtp;
    particles.push_back(par);
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/
//  bool zeroyawrate = fabs(yaw_rate) != 0.0;
  if(fabs(yaw_rate) > 10.0 || velocity > 100.0)
  {
    cout << "********************* BAD Sensor data *************************" << endl;
//    return;
  }
  bool zeroyawrate = fabs(yaw_rate) < 0.00001;
  std::default_random_engine gen;
  
  for(int p = 0; p < num_particles; p++)
  {
    Particle par = particles[p];
    
    double xp = par.x;
    double yp = par.y;
    double tp = par.theta;

    if(zeroyawrate)
    {
      double dist = velocity * delta_t;

      particles[p].x = xp + dist * cos(tp);
      particles[p].y = yp + dist * sin(tp);
    }
    else
    {
      double ang = yaw_rate * delta_t;
      double vbyw = velocity / yaw_rate;

      particles[p].x = xp + vbyw * ( sin(tp + ang) - sin(tp));
      particles[p].y = yp + vbyw * (-cos(tp + ang) + cos(tp));
      
      tp = tp + ang;
      if((tp > M_PI || tp < -M_PI) && p == 0)
	cout << "Theta " << tp << endl;
//      while (tp> M_PI) tp-=2.*M_PI;
//      while (tp<-M_PI) tp+=2.*M_PI;
      particles[p].theta = tp;
    }

    std::normal_distribution<> dx{0.0, std_pos[0]};
    std::normal_distribution<> dy{0.0, std_pos[1]};
    std::normal_distribution<> dt{0.0, std_pos[2]};
    particles[p].x += dx(gen);
    particles[p].y += dy(gen);
    particles[p].theta += dt(gen);
//    while (particles[p].theta> M_PI) particles[p].theta-=2.*M_PI;
//    while (particles[p].theta<-M_PI) particles[p].theta+=2.*M_PI;
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
  //   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
				   const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
  // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
  //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
  // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
  //   according to the MAP'S coordinate system. You will need to transform between the two systems.
  //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
  //   The following is a good resource for the theory:
  //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
  //   and the following is a good resource for the actual equation to implement (look at equation 
  //   3.33
  //   http://planning.cs.uiuc.edu/node99.html
  std::vector<Map::single_landmark_s> landmark_list = map_landmarks.landmark_list;
  double denox = 2.0 * std_landmark[0] * std_landmark[0];
  double denoy = 2.0 * std_landmark[1] * std_landmark[1];
  double term1 = 1.0 / (2.0 * M_PI * std_landmark[0] * std_landmark[1]);
  int on = observations.size();
  int ln = landmark_list.size();
  
  double total_wt = 0.0;
  for(int p = 0; p < num_particles; p++)
  {
    Particle par = particles[p];
    
    double xp = par.x;
    double yp = par.y;
    double tp = par.theta;
    double ctp = cos(tp);
    double stp = sin(tp);
    double mvgd = 1.0;
    std::vector<int> lassociations;
    std::vector<double> lsense_x;
    std::vector<double> lsense_y;
    for(int o = 0; o < on; o++)
    {
      double obs_x = observations[o].x;
      double obs_y = observations[o].y;

      double global_obs_x = obs_x * ctp - obs_y * stp + xp;
      double global_obs_y = obs_x * stp + obs_y * ctp + yp;

      double best_dist = 1e6;
      int best_index = -1;
      for(int l = 0; l < ln; l++)
      {
	double range = hypot(landmark_list[l].x_f - xp, landmark_list[l].y_f - yp);
	if(range <= sensor_range)
	{
	  double dist = hypot(landmark_list[l].x_f - global_obs_x, landmark_list[l].y_f - global_obs_y);
	  if(dist < best_dist)
	  {
	    best_dist = dist;
	    best_index = l;
	  }
	}
      }
      if(best_index >= 0)
      {
	double lx = landmark_list[best_index].x_f;
	double ly = landmark_list[best_index].y_f;
	double xd = global_obs_x - lx;
	double yd = global_obs_y - ly;
	double term2 = ((xd * xd) / denox) + ((yd * yd) / denoy);
	mvgd *= term1 * exp(-term2);
      }
      else
	mvgd *= 1e-6;
      lassociations.push_back(best_index+1);
      lsense_x.push_back(global_obs_x);
      lsense_y.push_back(global_obs_y);
    }
    particles[p].weight = mvgd;
    particles[p].associations = lassociations;
    particles[p].sense_x = lsense_x;
    particles[p].sense_y = lsense_y;
    total_wt += particles[p].weight;
  }
  if(total_wt > 0.0)
  {
    for(int p = 0; p < num_particles; p++)
    {
      particles[p].weight /= total_wt;
    }
  }
  
}

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight. 
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  std::vector<double> weights;
  for(int p = 0; p < num_particles; p++)
  {
    weights.push_back(particles[p].weight);
  }
  vector<Particle> new_particles(num_particles);
  random_device rd;
  default_random_engine gen(rd());
  for(int p = 0; p < num_particles; p++)
  {
    discrete_distribution<int> index(weights.begin(), weights.end());
    new_particles[p] = particles[index(gen)];
    new_particles[p].id = p;
  }
  particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
					 const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
  //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates

  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
  return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
  vector<int> v = best.associations;
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
string ParticleFilter::getSenseX(Particle best)
{
  vector<double> v = best.sense_x;
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
string ParticleFilter::getSenseY(Particle best)
{
  vector<double> v = best.sense_y;
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
