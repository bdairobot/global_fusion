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
static Eigen::Vector3d mag_w_m(0.0000,0.54159435858580351,-0.84063996);
class MagFactor : public ceres::SizedCostFunction<3,4,1>
{
public:
    MagFactor() = delete;
    MagFactor(double x, double y, double z, double m_dev)
        :xyz(x, y, z) {sqrt_info << 1.0/m_dev, 0.0, 0.0, 0.0, 1.0/m_dev, 0.0, 0.0, 0.0, 1.0/m_dev;}
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const
    {
        Eigen::Quaterniond q(parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3]);
        double decl = parameters[1][0];
        Eigen::Map<Eigen::Matrix<double, 3,1>> residual(residuals);
        Eigen::Vector3d mag_w = Eigen::Quaterniond(cos(decl/2.0), 0.0, 0.0, sin(decl/2.0))*mag_w_m;
        residual = sqrt_info*(xyz - q.inverse()*mag_w);
        if (jacobians){
            if (jacobians[0]){
                Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor>> jacobian_q(jacobians[0]);
                jacobian_q.setZero();
                jacobian_q.block<3,3>(0,1) = -sqrt_info*Utility::skewSymmetric(q.inverse()*mag_w);
            }
            if (jacobians[1]){
                Eigen::Map<Eigen::Matrix<double, 3, 1>> jacobian_decl(jacobians[1]);
                jacobian_decl = sqrt_info*(q.inverse()*Eigen::Vector3d(mag_w_m(0)*cos(decl), mag_w_m(0)*sin(decl), 0.0));
            }
        }
        return true;
    }

    void check(double **parameters)
    {
        double *residual = new double[3];
        double **jaco = new double* [2];
        jaco[0] = new double[3*4];
        jaco[1] = new double[3*1];
        Evaluate(parameters, residual, jaco);

        std::cout << "check begins" << std::endl << "my: " << std::endl;
        std::cout << Eigen::Map<Eigen::Matrix<double, 3, 1>> (residual).transpose() << std::endl << std::endl;
        std::cout << Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor>> (jaco[0]) << std::endl << std::endl;
        std::cout << Eigen::Map<Eigen::Matrix<double, 3, 1>> (jaco[1]) << std::endl << std::endl;

        Eigen::Quaterniond q(parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3]);
        std::cout << "q is: " << parameters[0][0] << ", "<< parameters[0][1] << ", "<< parameters[0][2] << ", "<< parameters[0][3] << std::endl;
        double decl = parameters[1][0];
        Eigen::Vector3d mag_w = Eigen::Quaterniond(cos(decl/2.0), 0.0, 0.0, sin(decl/2.0))*mag_w_m;

        Eigen::Matrix<double, 3, 1> res = sqrt_info*(xyz - q.inverse()*mag_w);
        std::cout << res.transpose() << std::endl << std::endl;
        
        std::cout << "num: " << std::endl;
        Eigen::Matrix<double, 3, 4> num_jaco;
        const double eps = 1e-6;
        Eigen::Quaterniond q_new = q*Utility::deltaQ(Eigen::Vector3d(eps, 0.0, 0.0));
        Eigen::Vector3d new_res = sqrt_info*(xyz - q_new.inverse() * mag_w);
        num_jaco.col(0) = (new_res - res) / eps;
        q_new = q*Utility::deltaQ(Eigen::Vector3d(0.0, eps, 0.0));
        new_res = sqrt_info*(xyz - q_new.inverse() * mag_w);
        num_jaco.col(1) = (new_res - res) / eps;
        q_new = q*Utility::deltaQ(Eigen::Vector3d(0.0, 0.0, eps));
        new_res = sqrt_info*(xyz - q_new.inverse() * mag_w);
        num_jaco.col(2) = (new_res - res) / eps;
        mag_w = Eigen::Quaterniond(cos((decl+eps)/2.0), 0.0, 0.0, sin((decl+eps)/2.0))*mag_w_m;
        new_res = sqrt_info*(xyz - q.inverse() * mag_w);
        num_jaco.col(3) = (new_res - res) / eps;

        std::cout << num_jaco << std::endl;

    }

    Eigen::Vector3d xyz;
    Eigen::Matrix3d sqrt_info;
};