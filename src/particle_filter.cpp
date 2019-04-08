/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"


using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  // Set the number of particles
  num_particles = 100;  
  std::default_random_engine gen;
  // Create normal distributions for x, y, and theta
    normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  // Sample from these normal distributions like this: 
  Particle samp;
  //std::cout <<  "GPS x: "  << x << "GPS y: " << y  << std::endl;
  for (int i = 0; i < num_particles; ++i) {
    samp.id = i;
    samp.x = dist_x(gen);
    samp.y = dist_y(gen);
    samp.theta = dist_theta(gen);
    samp.weight = 1;
    
    particles.push_back(samp);
    //std::cout << "particle: " << particles[i].id << "  x: "  << particles[i].x << "   y: " << particles[i].y  << std::endl;
  }
  is_initialized = true;
  
  

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine gen;
  double x_pred;
  double y_pred;
  for (int i = 0; i < num_particles; ++i) {
    if (yaw_rate != 0) {
    	x_pred = particles[i].x + velocity/yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
    	y_pred = particles[i].y + velocity/yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
    } else {
      	x_pred = particles[i].x + velocity * cos(particles[i].theta)*delta_t;
      	y_pred = particles[i].y + velocity * sin(particles[i].theta)*delta_t;
    }
    double theta_pred = particles[i].theta + yaw_rate*delta_t;
    /*if (i < 10) {
      std::cout << "particle: " << particles[i].id << "   i: " << i << std::endl;
      std::cout << "velocity: " << velocity << "   theta: " << particles[i].theta << "  yaw rate: " << yaw_rate << "   delta_t: " << delta_t << std::endl;
      std::cout << "initial x: " << particles[i].x << "   predicted x: " << x_pred << std::endl;
      std::cout << "initial y: " << particles[i].y << "   predicted y: " << y_pred <<  std::endl;
      std::cout << "initial t: " << particles[i].theta <<  "   predicted t: " << theta_pred << std::endl;
    }*/
    
    normal_distribution<double> dist_x(x_pred, std_pos[0]);
  	normal_distribution<double> dist_y(y_pred, std_pos[1]);
  	normal_distribution<double> dist_theta(theta_pred, std_pos[2]);
    
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);

    /*if (i < 10) {
      std::cout << "x with noise: " << particles[i].x << std::endl;
      std::cout << "y with noise: " << particles[i].y << std::endl;
      std::cout << "t with noise: " << particles[i].theta << std::endl;
      std::cout << std::endl;
    }*/
    
  }  
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
	
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  Particle part;	
  double min_dist;
  double min_x;
  double min_y;
  double x_map;
  double y_map;
  int association;
  
  for (int i = 0; i < num_particles; ++i) {
    part = particles[i];
    part.associations = {};
    part.sense_x = {};
    part.sense_y = {};
    part.weight = 1;
    /*if (debugger < 2) {
      std::cout << "Particle: " << part.id << std::endl;
      std::cout << "Weight: " << part.weight << std::endl;
      std::cout << "Associations: " << part.associations.size() << std::endl;
      std::cout << "X observations: " << part.sense_x.size() << std::endl;
      std::cout << "Y observations: " << part.sense_y.size() << std::endl;
    }*/
    
    for (int j = 0; j < observations.size(); ++j) {
      x_map = part.x + (cos(part.theta) * observations[j].x) - (sin(part.theta) * observations[j].y);
      y_map = part.y + (sin(part.theta) * observations[j].x) + (cos(part.theta) * observations[j].y);
      min_dist = sensor_range;
      min_x = sensor_range;
      min_y = sensor_range;
      /*if (debugger < 2) {
        std::cout << "Observation: " << j << "  x map: " << x_map << "   y map: " << y_map << std::endl;
      }*/
      for (int k = 0; k < map_landmarks.landmark_list.size(); ++k) {
        double distance = dist(map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f, x_map, y_map);
        /*if (debugger < 2) {
        	std::cout << "Landmark: " << map_landmarks.landmark_list[k].id_i << "  distance: " << distance << "  landmark x: " << map_landmarks.landmark_list[k].x_f << "  landmark y: " << map_landmarks.landmark_list[k].y_f  << std::endl;
        }*/
        if (distance < min_dist) {
          association = map_landmarks.landmark_list[k].id_i;
          min_dist = distance;
          min_x = map_landmarks.landmark_list[k].x_f;
          min_y = map_landmarks.landmark_list[k].y_f;
        }
      }
      part.associations.push_back(association);
      part.sense_x.push_back(x_map);
      part.sense_y.push_back(y_map);
      double wt = multiv_prob(std_landmark[0], std_landmark[1], x_map, y_map, min_x, min_y);
      part.weight *=  multiv_prob(std_landmark[0], std_landmark[1], x_map, y_map, min_x, min_y);
      /*if (part.weight > 10 & debugger < 10) {
        //std::cout << "Observation: " << j << "  x map: " << x_map << "   y map: " << y_map << std::endl;
        std::cout << "Landmark   : " << association << "  landmark x: " << min_x << "  landmark y: " << min_y << std::endl;
        std::cout << "Weight     : " << wt << "  Cum weight: " << part.weight << std::endl;
        std::cout << std::endl; 
        debugger++;
      }*/
    }
    
  	//std::cout << "Particle: " << part.id << "  Weight: " << part.weight << std::endl;
  	particles[i] = part;
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::vector<double> particle_weights;
  std::vector<Particle> particles_resamp;
  std::default_random_engine gen;
  Particle part;
  
  for (int j = 0; j < num_particles; ++j) {
    particle_weights.push_back(particles[j].weight);
    }
  std::discrete_distribution<> d(particle_weights.begin(), particle_weights.end());
  for(int n=0; n<num_particles; ++n) {
    part = particles[d(gen)];
    particles_resamp.push_back(part);
    }
  particles = particles_resamp;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}