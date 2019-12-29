/*******************************************************
 * 
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *
 * Author: Dai Bo (bdairobot@gmail.com)
 *******************************************************/

#pragma once
#include <vector>
#include <map>
#include <iostream>
#include <mutex>
#include <thread>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <ceres/ceres.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include "LocalCartesian.hpp"
#include "tic_toc.h"

using namespace std;

class GlobalOptimization
{
public:
	GlobalOptimization();
	~GlobalOptimization();
	void inputKeyframe(double t, Eigen::Vector3d OdomP, Eigen::Quaterniond OdomQ, double t_dev, double q_dev);
	void inputGPS(double gps_t, vector<double> &GPS_pose);
	void restart();
	void inputBaro(double baro_t, vector<double> &Baro);
	void inputAtt(double att_t, vector<double> &Att);
	void getGlobalOdom(Eigen::Vector3d OdomP, Eigen::Quaterniond OdomQ, Eigen::Vector3d &global_odomP, Eigen::Quaterniond &global_odomQ);
	nav_msgs::Path global_path;
	bool att_init;
	bool pos_init;
	Eigen::Matrix4d WGPS_T_WVIO;
	double sim_scale[1];

private:
	void GPS2XYZ(double latitude, double longitude, double altitude, double* xyz);
	void optimize();
	void updateGlobalPath();

	// format t, tx,ty,tz,qw,qx,qy,qz + t_dev, q_dev
	map<double, vector<double>> localPoseMap;
	map<double, vector<double>> globalPoseMap;

	map<double, vector<double>> GPSPositionMap;
	map<double, vector<double>> baroMap;
	map<double, vector<double>> attMap;
	
	bool newGPS;
	bool newBaro;
	bool newAtt;
	double mag_decl[1][1];
	GeographicLib::LocalCartesian geoConverter;
	std::mutex mPoseMap;
	Eigen::Vector3d lastP;
	Eigen::Quaterniond lastQ;
	std::thread threadOpt;

};