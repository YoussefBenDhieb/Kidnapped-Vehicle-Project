/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <map>
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
	default_random_engine gen;

	ParticleFilter::num_particles = 100;
	ParticleFilter::weights.assign(ParticleFilter::num_particles,1.0);

	//Create normal distributions for x, y and theta.
	normal_distribution<double> dist_x(x, std[0]);		
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	for (int i = 0; i < ParticleFilter::num_particles; ++i) {				
		Particle particle;
		particle.id = i;
        particle.x = dist_x(gen);
        particle.y = dist_y(gen);
        particle.theta = dist_theta(gen);
        particle.weight = weights.at(i);

        particles.push_back(particle);
    }
    ParticleFilter::is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;
	double xf, yf, thetaf; 
	Particle particle;
	for (int i = 0; i < ParticleFilter::num_particles; i++) {
		particle = ParticleFilter::particles.at(i);
		if(fabs(yaw_rate)>0.001){
			xf = particle.x + (velocity/yaw_rate)*(sin(particle.theta+yaw_rate*delta_t)-sin(particle.theta));
			yf = particle.y + (velocity/yaw_rate)*(cos(particle.theta)-cos(particle.theta+yaw_rate*delta_t));
			thetaf = particle.theta + yaw_rate*delta_t;
		}else{
			xf = particle.x + velocity*cos(particle.theta)*delta_t;
			yf = particle.y + velocity*(sin(particle.theta)*delta_t);
			thetaf = particle.theta;
		}

		normal_distribution<double> dist_x(xf, std_pos[0]);		
		normal_distribution<double> dist_y(yf, std_pos[1]);
		normal_distribution<double> dist_theta(thetaf, std_pos[2]);


		
        particle.x = dist_x(gen);
        particle.y = dist_y(gen);
        particle.theta = dist_theta(gen);
        
        ParticleFilter::particles.at(i) = particle;
    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for (int i=0; i<observations.size(); i++){
		double min= 10000.0;
		for (int j=0; j<predicted.size(); j++){
			double distance = dist(predicted.at(j).x,predicted.at(j).y, observations.at(i).x, observations.at(i).y);
			if (distance < min){
				min = distance;
				observations.at(i).id = predicted.at(j).id;
			}
		}
	}


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
	double sum_weights = 0;
	for(int p = 0; p<ParticleFilter::num_particles; p++){
		Particle particle = ParticleFilter::particles.at(p);
		
		//Create a vector for observations in map coordinates and predicted measurements within sensor range
		std::vector<LandmarkObs> observations_map, predicted;
		
		//convert the observations in map coordinates
		for (int i=0; i<observations.size(); i++){
			LandmarkObs observation;
			observation.x = cos(particle.theta)*observations.at(i).x -sin(particle.theta)*observations.at(i).y + particle.x;
			observation.y = sin(particle.theta)*observations.at(i).x +cos(particle.theta)*observations.at(i).y + particle.y;
			observation.id = observations.at(i).id;
			observations_map.push_back(observation);
		}

		//find predicted measurements
		for (int i=0; i<map_landmarks.landmark_list.size(); i++){
			if (dist(map_landmarks.landmark_list[i].x_f,map_landmarks.landmark_list[i].y_f, particle.x, particle.y) <= sensor_range){
				LandmarkObs pred;
			 	pred.x = map_landmarks.landmark_list[i].x_f;
			 	pred.y = map_landmarks.landmark_list[i].y_f;
			 	pred.id = map_landmarks.landmark_list[i].id_i;
			 	predicted.push_back(pred);
			 }
		}

		ParticleFilter::dataAssociation( predicted, observations_map);
		//Compute weight
		double weight =1.0;
		for (int i=0; i<observations_map.size(); ++i){
		double x_obs= observations_map.at(i).x;
		double y_obs= observations_map.at(i).y;
		double mu_x, mu_y;

		//find mu_x and mu_y
		bool id_found = false;
		int j =0;
		while(!id_found){
			if(predicted.at(j).id == observations_map.at(i).id){
				mu_x = predicted.at(j).x;
				mu_y = predicted.at(j).y;
				id_found = true;
			}
			j++;
		}
		
		//calculate normalization term
		double gauss_norm= 1/(2 * M_PI * std_landmark[0] * std_landmark[1]);

		//calculate exponent
		double exponent= pow((x_obs - mu_x),2)/(2 * pow(std_landmark[0],2)) + pow((y_obs - mu_y),2)/(2 * pow(std_landmark[1],2));

		//calculate weight using normalization terms and exponent
		weight= weight * gauss_norm * exp(-exponent);

		}
		particle.weight = weight;
		sum_weights += weight;

		//update the particle
		ParticleFilter::particles.at(p) = particle;
	}

	//normalize weights
	for(int p = 0; p<ParticleFilter::num_particles; p++){
		ParticleFilter::weights.at(p) = ParticleFilter::particles.at(p).weight/ sum_weights;
		ParticleFilter::particles.at(p).weight = ParticleFilter::weights.at(p);
		
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
		std::vector<Particle> new_particles;
		std::random_device rd;
    	std::mt19937 gen(rd());
    	std::discrete_distribution<> d(weights.begin(), weights.end());
    	std::map<int, int> m;
    	for(int n=0; n<ParticleFilter::num_particles; ++n) {
        	new_particles.push_back(ParticleFilter::particles[d(gen)]);
   		}
   		ParticleFilter::particles = new_particles;
    		

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
