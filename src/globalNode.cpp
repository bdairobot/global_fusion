/*******************************************************
 * 
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *
 * Author: Dai Bo (bdairobot@gmail.com)
 *******************************************************/

#include "ros/ros.h"
#include <ros/console.h>
#include "globalOpt.h"
#include <sensor_msgs/NavSatFix.h>
#include <sensor_msgs/FluidPressure.h>
#include <sensor_msgs/Imu.h>
#include <std_msgs/Header.h>
#include <std_msgs/Bool.h>
#include <std_msgs/Float64MultiArray.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <geometry_msgs/PointStamped.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <iostream>
#include <stdio.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>
#include <fstream>
#include <queue>
#include <mutex>
#include "LocalCartesian.hpp"
#include "utility/utility.h"
using namespace std;

GlobalOptimization globalEstimator;
ros::Publisher pub_global_pose, pub_gps_pose;
ros::Publisher pub_global_kf_path, pub_global_path, pub_gps_path, pub_gp_odom_path;
ros::Publisher pub_baro_height, pub_tmp_point;
ros::Publisher pub_vins_restart;
ros::Subscriber sub_myeye_imu;
nav_msgs::Path *global_kf_path, global_path, gps_path, gp_odom_path;
mutex m_buf, map_buf;

// gps information: time, x, y, z, xy_var, z_var
queue<pair<double, vector<double>>> tmpGPSQueue;
GeographicLib::LocalCartesian geoConverter;
// baro information: time, z, z_var
queue<pair<double, vector<double>>> tmpBaroQueue;
queue<pair<double, vector<double>>> tmpVIOQueue;
queue<pair<double, vector<double>>> tmpAttQueue;

// flag bit 0: gps, 1: baro, 2: att 3: vio
int init_flag = 0;
int error_flag = 0;

static const float BETA_TABLE[5] = {0,
				    2.71,
				    4.61,
				    6.25,
                    7.78
				   };
static double NOISE_BARO = 1.5;
static double NOISE_ATT = 0.01;

static uint BUF_DURATION = 5;

// aligned GPS, Baro and Att data with VIO using time stamp
static map<double, vector<double>> map_GPS; // x,y,z,xy_var,z_var. position
static map<double, vector<double>> map_Baro; // z, z_var. position
static map<double, vector<double>> map_Att; // w, x, y, z. var
static map<double, vector<double>> map_VIO_2G; // x_w, y_w, z_w, xyz_var, q_var
static double myeye_t = -1.0;
static double gap_t = 0.0;
static bool time_aligned = false;
static bool restart_flag = false;

class DataAnalyses {
public:
    DataAnalyses(uint data_len):_data_len(data_len), _sum(0.0), _sum_var(0.0){}
    ~DataAnalyses(){}
    void put_data(double data){
        if (_buf_data.size() == _data_len){
            _sum -= _buf_data.front();
            _sum_var -= _buf_data.front()*_buf_data.front();
            _buf_data.pop();
        }
        _buf_data.push(data);
        _sum += data;
        _sum_var += data*data;
    }
    pair<double,double> result() {
        if (_buf_data.size() < _data_len)
            return make_pair(0.0, -1.0); // return error value
        double mean = _sum / _data_len;
        double dev = sqrt(_sum_var/_data_len - mean*mean);
        return make_pair(mean, dev);
    }
private:
    queue<double> _buf_data;
    uint _data_len;
    double _sum;
    double _sum_var;
};

void GPS_callback(const sensor_msgs::NavSatFixConstPtr &GPS_msg)
{
    if (!time_aligned) return;
    double gps_t = GPS_msg->header.stamp.toSec();
    if (GPS_msg->status.status == -1){
        // error_flag |= 1<<0;
        return;
    }
    static int count = 0;
    if (!(init_flag & 1<<0)){
        // position accuracy less than 2 meters or after 2 seconds
        if (GPS_msg->position_covariance[0] < 4.0 || count == 5){
            geoConverter.Reset(GPS_msg->latitude, GPS_msg->longitude, GPS_msg->altitude);
            init_flag |= 1<<0;
            ROS_INFO("geo init: lat %12.8f,lon %12.8f, alt %12.8f, cov: %8.4f", GPS_msg->latitude,GPS_msg->longitude, GPS_msg->altitude, GPS_msg->position_covariance[0]);
        }
        count++;
    }
    if (init_flag & 1<<0) {
        
        geometry_msgs::PoseStamped pose_gps;
        pose_gps.header.stamp = ros::Time(gps_t);
        pose_gps.header.frame_id = "world";
        double xyz[3];
        geoConverter.Forward(GPS_msg->latitude, GPS_msg->longitude,GPS_msg->altitude, xyz[0],xyz[1],xyz[2]);
           //printf("gps_callback! ");
        
        static double last_gps_dev[3] = {0.0,0.0,0.0}; 
        static double diff_gps_dev_xy = 0;
        static double diff_gps_dev_z = 0;

        if (last_gps_dev[0] > 0.0){
            diff_gps_dev_xy = 0.3*diff_gps_dev_xy + 0.7 *(sqrt(GPS_msg->position_covariance[0]) - last_gps_dev[1])/(gps_t -  last_gps_dev[0]);

            diff_gps_dev_z = 0.3*diff_gps_dev_z + 0.7 *(sqrt(GPS_msg->position_covariance[8]) - last_gps_dev[2])/(gps_t -  last_gps_dev[0]);
        }
       
        vector<double> gps_info = {xyz[0],xyz[1],xyz[2], GPS_msg->position_covariance[0],GPS_msg->position_covariance[8], diff_gps_dev_xy, diff_gps_dev_z};
        m_buf.lock();
        last_gps_dev[0] = gps_t;
        last_gps_dev[1] = sqrt(GPS_msg->position_covariance[0]);
        last_gps_dev[2] = sqrt(GPS_msg->position_covariance[8]);
        tmpGPSQueue.push(make_pair(gps_t, gps_info));
        geometry_msgs::PointStamped test;
        test.header = pose_gps.header;
        test.point.x = sqrt(GPS_msg->position_covariance[0]);
        test.point.y = diff_gps_dev_xy;
        test.point.z = diff_gps_dev_xy / sqrt(GPS_msg->position_covariance[0]) * 10;
        pub_tmp_point.publish(test);
        // store 3s data
        m_buf.unlock();
        if (tmpGPSQueue.size() > 5 * BUF_DURATION)
            tmpGPSQueue.pop();
        pose_gps.pose.position.x = xyz[0];
        pose_gps.pose.position.y = xyz[1];
        pose_gps.pose.position.z = xyz[2];
        pose_gps.pose.orientation.w = 1.0;
        
        gps_path.header = pose_gps.header;
        gps_path.poses.push_back(pose_gps);
        pub_gps_path.publish(gps_path);
    }
}

void baro_callback(const sensor_msgs::FluidPressureConstPtr &baro_msg)
{   
    double baro_t = baro_msg->header.stamp.toSec();
    if (myeye_t >= 0.0 && !time_aligned){
        gap_t = baro_t - myeye_t;
        gap_t = (gap_t< 0.5 && gap_t > -0.5) ? 0 : gap_t; 
        std::cout << "Init gap time between pix and myeye: "<< gap_t << std::endl;
        sub_myeye_imu.shutdown();
        time_aligned = true;
    }
    static double CONSTANTS_ABSOLUTE_NULL_CELSIUS = -273.15;
    static double CONSTANTS_AIR_GAS_CONST = 287.1f;
    static double BARO_MSL = 101.325; /* current pressure at MSL in kPa */
    static double CONSTANTS_ONE_G = 9.80665f;
    static double T1 = 15.0 - CONSTANTS_ABSOLUTE_NULL_CELSIUS;	/* temperature at base height in Kelvin */
    static double a  = -6.5 / 1000.0;	/* temperature gradient in degrees per metre */

    /* measured pressure in kPa */
    const double p = baro_msg->fluid_pressure * 0.001;
    /*
        * Solve:
        *
        *     /               -(aR / g)    \
        *    | (p / BARO_MSL)          . T1 | - T1
        *     \                            /
        * h = -------------------------------  + h1
        *                   a
        */
    double height = (((pow((p / BARO_MSL), (-(a * CONSTANTS_AIR_GAS_CONST) / CONSTANTS_ONE_G))) * T1) - T1) / a;
    static double init_height = 0.0;
    static double last_height = height;
    double filter_factor = 0.5;
    double filter_height = filter_factor*height + (1-filter_factor) * last_height;
    last_height = filter_height;
    m_buf.lock();
    vector<double> baro_info = {filter_height - init_height, NOISE_BARO*NOISE_BARO};
    tmpBaroQueue.push(make_pair(baro_t,baro_info));

    // store 3s data
    if (tmpBaroQueue.size() > 50*BUF_DURATION)
        tmpBaroQueue.pop();
    static unsigned int init_size = 10;
    if (!(init_flag & 1<<1) && tmpBaroQueue.size() == init_size){        
        init_height = filter_height;
        while (!tmpBaroQueue.empty())
            tmpBaroQueue.pop();
        init_flag |= 1<<1;
        ROS_INFO("init baro: height %8.4f", init_height);
    }
    m_buf.unlock();
    if (init_flag & 1<<1){
        geometry_msgs::PointStamped baro_height;
        baro_height.header = baro_msg->header;
        baro_height.point.x = height - init_height;
        baro_height.point.y = filter_height - init_height;
        pub_baro_height.publish(baro_height);

        // static DataAnalyses data_analy(100);
        // data_analy.put_data(height - init_height);
        // if (data_analy.result().second > 0.0)
        //     ROS_INFO("baro mean: %8.4f, dev: %8.4f", data_analy.result().first, data_analy.result().second);
    }
}

void keyframe_callback(const nav_msgs::OdometryConstPtr &pose_msg)
{
    if (!time_aligned) return;
    // if (!(error_flag & 1<<3 && pose_msg->pose.covariance[0]>=0.0)){
    double cov = pose_msg->pose.covariance[0];
    double q_cov = pose_msg->pose.covariance[28];
    if (cov < 0.0) {
        cov = 1.0;
        q_cov = 0.1;
    }
    {
        double kf_t = pose_msg->header.stamp.toSec() + gap_t;
        Eigen::Vector3d vio_p(pose_msg->pose.pose.position.x, pose_msg->pose.pose.position.y, pose_msg->pose.pose.position.z);
        Eigen::Quaterniond vio_q;
        vio_q.w() = pose_msg->pose.pose.orientation.w;
        vio_q.x() = pose_msg->pose.pose.orientation.x;
        vio_q.y() = pose_msg->pose.pose.orientation.y;
        vio_q.z() = pose_msg->pose.pose.orientation.z;
	bool input_flag = false;
        map_buf.lock();
        if(map_GPS.find(kf_t) != map_GPS.end() && !(error_flag & 1<<0)){
            globalEstimator.inputGPS(kf_t, map_GPS[kf_t]);
            geometry_msgs::PoseStamped pose_gps;
            pose_gps.header.stamp = ros::Time(kf_t);
            pose_gps.header.frame_id = "world";
            pose_gps.pose.position.x = map_GPS[kf_t][0];
            pose_gps.pose.position.y = map_GPS[kf_t][1];
            pose_gps.pose.position.z = map_GPS[kf_t][2];
            pose_gps.pose.orientation.w = 1.0;
            pub_gps_pose.publish(pose_gps);
            input_flag = true;
            // std::cout << "map_GPS[kf_t]: " << map_GPS[kf_t][0] << ", " << map_GPS[kf_t][1] << ", " << map_GPS[kf_t][2] << std::endl;
        }
        if(map_Baro.find(kf_t) != map_Baro.end()){
            globalEstimator.inputBaro(kf_t, map_Baro[kf_t]);
            // std::cout << "map_Baro[kf_t]: " << map_Baro[kf_t][0] << ", " << map_Baro[kf_t][1] << ", " << map_Baro[kf_t][2] << std::endl;
            input_flag = true;
        }
        if(map_Att.find(kf_t) != map_Att.end()){
            globalEstimator.inputAtt(kf_t, map_Att[kf_t]);
            input_flag = true;
	}
        if (input_flag){
            globalEstimator.inputKeyframe(kf_t, vio_p, vio_q, cov, q_cov);
            pub_global_kf_path.publish(*global_kf_path);
        }
        map_buf.unlock();
        

    }
}

void vio_callback(const nav_msgs::OdometryConstPtr &pose_msg)
{
    if (!time_aligned) return;
    static double last_update_t = pose_msg->header.stamp.toSec() + gap_t;
    double t = pose_msg->header.stamp.toSec() + gap_t;
    // ROS_INFO("delta time between vision and imu: %8.4f", tmpAttQueue.back().first - t); // about 100ms

    Eigen::Vector3d vio_p(pose_msg->pose.pose.position.x, pose_msg->pose.pose.position.y, pose_msg->pose.pose.position.z);
    Eigen::Quaterniond vio_q;
    vio_q.w() = pose_msg->pose.pose.orientation.w;
    vio_q.x() = pose_msg->pose.pose.orientation.x;
    vio_q.y() = pose_msg->pose.pose.orientation.y;
    vio_q.z() = pose_msg->pose.pose.orientation.z;

    if (t - last_update_t > 0.25) {
       globalEstimator.restart(); // vins restarted
       std::cout << "No vio info for " << t - last_update_t << std::endl;
       map_VIO_2G.clear();
    }
    last_update_t = t;

    if (!globalEstimator.att_init){
        Eigen::Quaterniond att(1.0, 0.0, 0.0, 0.0);
        if (!tmpAttQueue.empty()){
            vector<double> att_info = tmpAttQueue.back().second;
            att = Eigen::Quaterniond(att_info[0], att_info[1], att_info[2], att_info[3]);
        } else {
            static int i = 0; i++;
            if(i >10) {i=0; std::cout << "Need att info to start..." << std::endl;}
	    return;
        }

        globalEstimator.WGPS_T_WVIO.block<3,3>(0,0) = (att*vio_q.inverse()).toRotationMatrix();
        globalEstimator.att_init = true;
    }
    if (!globalEstimator.pos_init){
        Eigen::Vector3d init_pos(0.0,0.0,0.0);
        if (tmpGPSQueue.size()>0 && (init_flag & 1<<0)){
            init_pos(0) = tmpGPSQueue.back().second[0];
            init_pos(1) = tmpGPSQueue.back().second[1];                
        }
        if (tmpBaroQueue.size()>0 && (init_flag & 1 << 1))
            init_pos(2) = tmpBaroQueue.back().second[0]; 
        
        globalEstimator.WGPS_T_WVIO.block<3,1>(0,3) = init_pos - globalEstimator.WGPS_T_WVIO.block<3,3>(0,0)*vio_p;
        globalEstimator.pos_init = true;
    }

    //printf("vio_callback! ");
    // int vio_error_flag = int(-pose_msg->pose.covariance[0]);
    // if ((vio_error_flag & 1<<3) || (vio_error_flag & 1<<1)){
    //     return;
    if (restart_flag){
         std_msgs::Bool restart_flag;
         restart_flag.data = true;
        
         pub_vins_restart.publish(restart_flag);
         globalEstimator.restart();
         ROS_INFO("restarting optimization! ");
         return;
     }
    error_flag &= ~(1<<3);
    
    Eigen::Vector3d global_t;
    Eigen::Quaterniond global_q;
    globalEstimator.getGlobalOdom(vio_p, vio_q, global_t, global_q);
    geometry_msgs::PoseStamped vio_pose;
    vio_pose.header = pose_msg->header;
    vio_pose.header.frame_id = "world";
    vio_pose.pose.position.x = global_t.x();
    vio_pose.pose.position.y = global_t.y();
    vio_pose.pose.position.z = global_t.z();
    vio_pose.pose.orientation.x = global_q.x();
    vio_pose.pose.orientation.y = global_q.y();
    vio_pose.pose.orientation.z = global_q.z();
    vio_pose.pose.orientation.w = global_q.w();
    pub_global_pose.publish(vio_pose);
    global_path.header = vio_pose.header;
    global_path.poses.push_back(vio_pose);
    pub_global_path.publish(global_path);

    static uint duration = 4;
    { // align with baro
        static uint step = 10;
        if (init_flag & 1<<1){
            while(!tmpBaroQueue.empty()){
                double baro_t = tmpBaroQueue.front().first;
                if (t < baro_t - 0.02 && (tmpBaroQueue.size() == 50*BUF_DURATION)){
                    ROS_ERROR("VIO is away behind BARO information! ");                    
                    assert(0);
                } else if (t <= baro_t + 0.02 && t >= baro_t - 0.02){
                    map_buf.lock();
                    if (map_Baro.size() == duration*step)
                        map_Baro.erase(map_Baro.begin());

                    map_Baro[t] = tmpBaroQueue.front().second;
                    map_buf.unlock();
                    tmpBaroQueue.pop();
                    break;
                } else if (t - baro_t > 0.02){
                    tmpBaroQueue.pop();
                    continue;
                } else break;
            }
        }
    }
    //std::cout << "map_Baro.size()" << map_Baro.size() << " tmpBaroQueue.size(): " <<tmpBaroQueue.size() << " t: " << t << " baro_t: " << tmpBaroQueue.back().first << "t-barot: " << t-tmpBaroQueue.back().first << std::endl;

    { // align with att
        static uint step = 10;
        while(!tmpAttQueue.empty()){
            double att_t = tmpAttQueue.front().first;
            if (t < att_t - 0.02 && (tmpAttQueue.size() == 50*BUF_DURATION)){
                ROS_ERROR("VIO is away behind ATT information! ");
                assert(0);
            } else if (t <= att_t + 0.02 && t >= att_t - 0.02){
                map_buf.lock();
                if (map_Att.size() == duration*step)
                    map_Att.erase(map_Att.begin());

                map_Att[t] = tmpAttQueue.front().second;
                map_buf.unlock();
                tmpAttQueue.pop();
                break;
            } else if (t - att_t > 0.02){
                tmpAttQueue.pop();
                continue;
            } else break;
        }
    }

    { // check with gps
        static uint step = 10;
        if (init_flag & 1<<0){
            while(!tmpGPSQueue.empty()){
                pair<double, vector<double>> GPS_info = tmpGPSQueue.front();
                static bool gps_noisy = false;
                static bool pop_flag = false;

                if((GPS_info.second[3] > 16.0 && GPS_info.second[5] > 0.2) || GPS_info.second[3] > 25.0 || (GPS_info.second[3] > 16 && gps_noisy)){
                    gps_noisy = true;
                    pop_flag = false;
                    tmpGPSQueue.pop(); continue;
                }
                gps_noisy = false;
                double gps_t = GPS_info.first;
                Eigen::Vector3d gps_pose(GPS_info.second[0],GPS_info.second[1],GPS_info.second[2]);
                static pair<double, vector<double>> last_pop;
                if (t < gps_t - 0.1 && (tmpGPSQueue.size()==5*BUF_DURATION)){
                    ROS_ERROR("VIO is away behind GPS information! ");
                    assert(0);
                } else if (t <= gps_t + 0.1 && t >= gps_t - 0.1){
                    map_buf.lock();
                    if (map_GPS.size() == duration*step){
                        map_GPS.erase(map_GPS.begin());
                        map_VIO_2G.erase(map_VIO_2G.begin());
                    }

                    map_GPS[t] = GPS_info.second;
                    //map_VIO_2G[t] = {global_t.x(), global_t.y(), global_t.z(), pose_msg->pose.covariance[0], pose_msg->pose.covariance[28]};
                    Eigen::Vector3d pos = globalEstimator.WGPS_T_WVIO.block<3,3>(0,0)*vio_p;
                    map_VIO_2G[t] = {pos.x(), pos.y(), pos.z(), pose_msg->pose.covariance[0], pose_msg->pose.covariance[28]};
                    last_pop = GPS_info;
                    // tmpGPSQueue.pop();
                    pop_flag = true;
                    map_buf.unlock();
                    break;
                // } else tmpGPSQueue.pop();
                }else if(!pop_flag || t > gps_t + 0.1){
                    last_pop = GPS_info;
                    tmpGPSQueue.pop();
                    pop_flag = true;
                    continue;
                }
                if(pop_flag && t < gps_t+0.2 && t > last_pop.first && gps_t - last_pop.first > 0){
                    Eigen::Vector3d insert_gps_pose = gps_pose - (gps_t - t)/(gps_t - last_pop.first)*(gps_pose - Eigen::Vector3d{last_pop.second[0],last_pop.second[1],last_pop.second[2]});
                    map_buf.lock();
                    if (map_GPS.size() == duration*step){
                        map_GPS.erase(map_GPS.begin());
                        map_VIO_2G.erase(map_VIO_2G.begin());
                    }

                    map_GPS[t] = vector<double> {insert_gps_pose(0), insert_gps_pose(1), insert_gps_pose(2), GPS_info.second[3], GPS_info.second[4]};
                    map_VIO_2G[t] = {global_t.x(), global_t.y(), global_t.z(), pose_msg->pose.covariance[0], pose_msg->pose.covariance[28]};
                    last_pop = GPS_info;
                    // tmpGPSQueue.pop();
                    pop_flag = true;
                    map_buf.unlock();
                    break;
                } else {
                    // ROS_INFO("t-gps_t: %8.4f, t - last_pop.first: %8.4f", t-gps_t, t - last_pop.first);
                    break;
                };
            }
            //GPS check
            {  
                if (map_VIO_2G.size() == duration*step){
                    auto iter_back = map_VIO_2G.rbegin();
                    for (int i = 0; i < 10; i++) iter_back++;
                    auto iter_front = map_VIO_2G.begin();
                    int counter[3] = {0, 0, 0};
                    for (int i = 0; i < 10; i++, iter_front++, iter_back--){
                        Eigen::Vector2d vio(iter_back->second[0] - iter_front->second[0], iter_back->second[1]-iter_front->second[1]);
                        Eigen::Vector2d gps(map_GPS[iter_back->first][0] - map_GPS[iter_front->first][0], map_GPS[iter_back->first][1] - map_GPS[iter_front->first][1]);
                        double chi_value = (gps - vio).norm() /(map_GPS[iter_back->first][3]+iter_back->second[3]);
                        if (chi_value >= 1) counter[0]++;
			if (chi_value >= 0.45) counter[1]++;
			if (chi_value >= 0.35) counter[2]++;
                    }
		    if (counter[0] >1 || counter[1] > 5 || counter[2] > 7)
			std::cout << "Restarting" << std::endl;
                }
            }
        }
    }
}

void myeye_imu_callback(const sensor_msgs::ImuConstPtr &imu_msg){
    myeye_t = imu_msg->header.stamp.toSec();
}

void imu_callback(const sensor_msgs::ImuConstPtr &imu_msg){
    if (!time_aligned) return;
    double att_t = imu_msg->header.stamp.toSec();
    m_buf.lock();
    vector<double> att_info = {imu_msg->orientation.w, imu_msg->orientation.x, imu_msg->orientation.y, imu_msg->orientation.z, NOISE_ATT*NOISE_ATT};
    tmpAttQueue.push(make_pair(att_t, att_info));
    m_buf.unlock();

    if(tmpAttQueue.size() > 50*BUF_DURATION)
        tmpAttQueue.pop();
}

void gp_odom_callback(nav_msgs::Odometry::ConstPtr gp_odom)
{
        gp_odom_path.header = gp_odom->header;
	gp_odom_path.header.frame_id = "world";
	geometry_msgs::PoseStamped pose;
	pose.header = gp_odom->header;
	pose.header.frame_id = "world";
	pose.pose = gp_odom->pose.pose;
        gp_odom_path.poses.push_back(pose)	;
        pub_gp_odom_path.publish(gp_odom_path);
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "globalEstimator");
    ros::NodeHandle n("~");
    global_kf_path = &globalEstimator.global_path;

    ros::Subscriber sub_GPS = n.subscribe("/mavros/global_position/raw/fix", 200, GPS_callback);
    ros::Subscriber sub_baro = n.subscribe("/mavros/imu/static_pressure",200, baro_callback);
    ros::Subscriber sub_kf = n.subscribe("/vins_estimator/keyframe_pose", 200, keyframe_callback);
    ros::Subscriber sub_vio = n.subscribe("/vins_estimator/odometry", 200, vio_callback);
    sub_myeye_imu = n.subscribe("/mynteye/imu/data_raw", 1000, myeye_imu_callback);
    ros::Subscriber sub_imu = n.subscribe("/mavros/imu/data", 1000, imu_callback);
    ros::Subscriber sub_gp_odom = n.subscribe("/mavros/local_position/odom", 100, gp_odom_callback);

    pub_global_path = n.advertise<nav_msgs::Path>("global_path", 1000);
    pub_global_kf_path = n.advertise<nav_msgs::Path>("global_kf_path", 1000);
    pub_gps_path = n.advertise<nav_msgs::Path>("gps_path", 1000);
    pub_gp_odom_path = n.advertise<nav_msgs::Path>("gp_odom_path", 1000);
    pub_global_pose = n.advertise<geometry_msgs::PoseStamped>("global_pose", 100);
    pub_gps_pose = n.advertise<geometry_msgs::PoseStamped>("gps_pose", 100);
    pub_baro_height = n.advertise<geometry_msgs::PointStamped>("baro_height", 100);
    pub_vins_restart = n.advertise<std_msgs::Bool>("/vins_restart", 100);
    pub_tmp_point = n.advertise<geometry_msgs::PointStamped> ("/test_point", 100);

    std::cout << "Starting Global Fusion ..." << std::endl;
    ros::spin();
    std::ofstream foutKF("/home/bdai/output/global_kf_path.txt", std::ios::out);
    std::ofstream fout("/home/bdai/output/global_path.txt", std::ios::out);
    std::ofstream foutGPS("/home/bdai/output/gps_path.txt", std::ios::out);
    std::ofstream foutGP("/home/bdai/output/gp_odom_path.txt", std::ios::out);
    if (global_kf_path->poses.size()){
        for (uint i = 0; i < global_kf_path->poses.size(); i ++){
            foutKF.setf(ios::fixed, ios::floatfield);
            foutKF.precision(6);
            foutKF << global_kf_path->poses[i].header.stamp.toSec() << " ";
            foutKF << global_kf_path->poses[i].pose.position.x << " "
                   << global_kf_path->poses[i].pose.position.y << " "
                   << global_kf_path->poses[i].pose.position.z << " "
                   << global_kf_path->poses[i].pose.orientation.x << " "
                   << global_kf_path->poses[i].pose.orientation.y << " "
                   << global_kf_path->poses[i].pose.orientation.z << " "
                   << global_kf_path->poses[i].pose.orientation.w << endl;
        }
    }
    if (global_path.poses.size()){
        for (uint i = 0; i < global_path.poses.size(); i ++){
            fout.setf(ios::fixed, ios::floatfield);
            fout.precision(6);
            fout << global_path.poses[i].header.stamp.toSec()<< " ";
            fout << global_path.poses[i].pose.position.x << " "
                   << global_path.poses[i].pose.position.y << " "
                   << global_path.poses[i].pose.position.z << " "
                   << global_path.poses[i].pose.orientation.x << " "
                   << global_path.poses[i].pose.orientation.y << " "
                   << global_path.poses[i].pose.orientation.z << " "
                   << global_path.poses[i].pose.orientation.w << endl;
        }
    }
    if (gps_path.poses.size()){
        for (uint i = 0; i < gps_path.poses.size(); i ++){
            foutGPS.setf(ios::fixed, ios::floatfield);
            foutGPS.precision(6);
            foutGPS << gps_path.poses[i].header.stamp.toSec()<< " ";
            foutGPS << gps_path.poses[i].pose.position.x << " "
                    << gps_path.poses[i].pose.position.y << " "
                    << gps_path.poses[i].pose.position.z << " "
                    << gps_path.poses[i].pose.orientation.x << " "
                    << gps_path.poses[i].pose.orientation.y << " "
                    << gps_path.poses[i].pose.orientation.z << " "
                    << gps_path.poses[i].pose.orientation.w << endl;
        }
    }
    if (gp_odom_path.poses.size()){
        for (uint i = 0; i < gp_odom_path.poses.size(); i ++){
            foutGP.setf(ios::fixed, ios::floatfield);
            foutGP.precision(6);
            foutGP << gp_odom_path.poses[i].header.stamp.toSec()<< " ";
            foutGP << gp_odom_path.poses[i].pose.position.x << " "
                    << gp_odom_path.poses[i].pose.position.y << " "
                    << gp_odom_path.poses[i].pose.position.z << " "
                    << gp_odom_path.poses[i].pose.orientation.x << " "
                    << gp_odom_path.poses[i].pose.orientation.y << " "
                    << gp_odom_path.poses[i].pose.orientation.z << " "
                    << gp_odom_path.poses[i].pose.orientation.w << endl;
        }
    }

    foutKF.close();
    fout.close();
    foutGPS.close();
    return 0;
}
