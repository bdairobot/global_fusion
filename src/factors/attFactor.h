/*******************************************************
 * 
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *
 * Author: Dai Bo (bdairobot@gmail.com)
 *******************************************************/

#pragma once
#include <ceres/ceres.h>
#include <eigen3/Eigen/Dense>
#include "../utility/utility.h"

class AttFactor : public ceres::SizedCostFunction<3,4>
{
public:
    AttFactor()= delete;
    AttFactor(double w, double x, double y, double z,  double att_dev): q_inv(Eigen::Quaterniond(w, x, y, z).inverse()){
        sqrt_info << 1.0/att_dev, 0.0, 0.0, 0.0, 1.0/att_dev, 0.0, 0.0, 0.0, 1.0/att_dev;
    };
    virtual bool Evaluate(double const *const *paramters, double *residuals, double **jacobians) const
    {
        Eigen::Quaterniond q_w_i(paramters[0][0], paramters[0][1], paramters[0][2], paramters[0][3]);
        Eigen::Map<Eigen::Matrix<double, 3, 1>> resudual(residuals);
        resudual = 2.0*sqrt_info*(q_inv*q_w_i).vec();
        if(jacobians) {
            if(jacobians[0]) {
                Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor>> jaco_q(jacobians[0]);
                jaco_q.setZero();
                jaco_q.block<3,3>(0,1) = sqrt_info*Utility::Qleft(q_inv*q_w_i).bottomRightCorner<3,3>();
            }
        }
        return true;
    }

    void check(double **parameters)
    {
        double *residual = new double[3];
        double **jaco = new double *[1];
        jaco[0] = new double[3*3];
        Evaluate(parameters, residual, jaco);

        std::cout << "check begins" << std::endl << "my: " << std::endl;
        std::cout << Eigen::Map<Eigen::Matrix<double, 3, 1>> (residual).transpose() << std::endl << std::endl;
        std::cout << Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor>> (jaco[0]) << std::endl << std::endl;

        Eigen::Quaterniond q_w_i(parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3]);
        std::cout << "q_w_i is: " << parameters[0][0] << ", "<< parameters[0][1] << ", "<< parameters[0][2] << ", "<< parameters[0][3] << std::endl;

        puts("num: ");
        Eigen::Matrix<double,3,1> res = 2.0*sqrt_info*(q_inv*q_w_i).vec();
        std::cout << res.transpose() << std::endl;
        Eigen::Matrix3d num_jaco;
        const double eps = 1e-6;
        for (int k = 0; k < 3; k++){
            Eigen::Quaterniond q_new;
            Eigen::Vector3d delta = Eigen::Vector3d(k==0, k==1, k==2)*eps;
            q_new = q_w_i*Utility::deltaQ(delta);
            Eigen::Vector3d tmp_res = 2.0*sqrt_info*(q_inv*q_new).vec();
            num_jaco.col(k) = (tmp_res - res) / eps;
        }
        std::cout << num_jaco << std::endl;
    }
private:
    Eigen::Quaterniond q_inv;
    Eigen::Matrix3d sqrt_info;
};