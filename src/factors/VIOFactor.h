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

class VIOFactors : public ceres::SizedCostFunction<6,2,1,4,2,1,4> 
{
public:
    VIOFactors() = delete;
    VIOFactors(const Eigen::Vector3d &p_i, const Eigen::Quaterniond &q_i, const Eigen::Vector3d &p_j, const Eigen::Quaterniond &q_j, double p_dev, double q_dev){
        dpos = q_i.inverse()*(p_j - p_i);
        dq_l_inv = (q_i.inverse()*q_j).inverse();
        sqrt_info_p << 1.0/p_dev,0.0,0.0,0.0,1.0/p_dev,0.0,0.0,0.0,1.0/p_dev;
        sqrt_info_q << 1.0/q_dev,0.0,0.0,0.0,1.0/q_dev,0.0,0.0,0.0,1.0/q_dev;
    }

    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const 
    {
        Eigen::Vector3d p_w_i(parameters[0][0], parameters[0][1], parameters[1][0]);
        // Eigen::Quaterniond q_w_i(parameters[1]); // Cannot use this, because Eigen store quaternion as: x, y, z, w. But construct as w, x, y, z
        Eigen::Quaterniond q_w_i(parameters[2][0], parameters[2][1], parameters[2][2], parameters[2][3]);
        Eigen::Vector3d p_w_j(parameters[3][0], parameters[3][1], parameters[4][0]);
        Eigen::Quaterniond q_w_j(parameters[5][0], parameters[5][1], parameters[5][2], parameters[5][3]);
        Eigen::Map<Eigen::Matrix<double,6,1>> residual(residuals);
        residual.block<3,1>(0,0) = sqrt_info_p*(dpos - q_w_i.inverse()*(p_w_j - p_w_i));
        residual.block<3,1>(3,0) = 2.0*sqrt_info_q*(dq_l_inv*(q_w_i.inverse()*q_w_j)).vec();

        if (jacobians)
    	{
            Eigen::Matrix3d jaco_tmp = sqrt_info_p*q_w_i.inverse().toRotationMatrix();
            if (jacobians[0]){
                Eigen::Map<Eigen::Matrix<double, 6, 2, Eigen::RowMajor>> jacobian_p_i_xy(jacobians[0]);
                jacobian_p_i_xy.setZero();
                jacobian_p_i_xy.block<3,2>(0,0) = jaco_tmp.block<3,2>(0,0);
            } 
            if (jacobians[1]){
                Eigen::Map<Eigen::Matrix<double, 6, 1>> jacobian_p_i_z(jacobians[1]);
                jacobian_p_i_z.setZero();
                jacobian_p_i_z.block<3,1>(0,0) = jaco_tmp.block<3,1>(0,2);
            }
            
            if (jacobians[2]){
                Eigen::Map<Eigen::Matrix<double, 6, 4, Eigen::RowMajor>> jacobian_q_i(jacobians[2]);
                jacobian_q_i.setZero();
                jacobian_q_i.block<3,3>(0,1) = (-sqrt_info_p)*Utility::skewSymmetric(q_w_i.inverse()*(p_w_j-p_w_i));
                jacobian_q_i.block<3,3>(3,1) = (-sqrt_info_q)*(Utility::Qright(q_w_i.inverse()*q_w_j)*Utility::Qleft(dq_l_inv)).bottomRightCorner<3,3>();
            }
            if (jacobians[3]){
                Eigen::Map<Eigen::Matrix<double, 6, 2, Eigen::RowMajor>> jacobian_p_j_xy(jacobians[3]);
                jacobian_p_j_xy.setZero();
                jacobian_p_j_xy.block<3,2>(0,0) = -jaco_tmp.block<3,2>(0,0);
            }

            if (jacobians[4]){
                Eigen::Map<Eigen::Matrix<double, 6, 1>> jacobian_p_j_z(jacobians[4]);
                jacobian_p_j_z.setZero();
                jacobian_p_j_z.block<3,1>(0,0) = -jaco_tmp.block<3,1>(0,2);
            }

            if (jacobians[5]){
                Eigen::Map<Eigen::Matrix<double, 6, 4, Eigen::RowMajor>> jacobian_q_j(jacobians[5]);
                jacobian_q_j.setZero();
                jacobian_q_j.block<3,3>(3,1) = sqrt_info_q*Utility::Qleft(dq_l_inv*(q_w_i.inverse()*q_w_j)).bottomRightCorner<3,3>();
            }
        }
        return true;
    }

    void check(double** parameters)
    {
        double* res = new double[6];
        double** jaco = new double* [6];
        jaco[0] = new double[6*2];
        jaco[1] = new double[6*1];
        jaco[2] = new double[6*4];
        jaco[3] = new double[6*2];
        jaco[4] = new double[6*1];
        jaco[5] = new double[6*4];

        Evaluate(parameters, res, jaco);
        // puts("check begins");

        // puts("my");
        if (sqrt(res[0]*res[0] + res[1]*res[1]) > 0.1 || res[3] > 0.02 || res[4] > 0.02 || res[4] > 0.02)
            std::cout << Eigen::Map<Eigen::Matrix<double, 6, 1>> (res).transpose() << std::endl;
        // std::cout << Eigen::Map<Eigen::Matrix<double, 6, 2, Eigen::RowMajor>> (jaco[0]) << std::endl << std::endl;
        // std::cout << Eigen::Map<Eigen::Matrix<double, 6, 1>> (jaco[1]) << std::endl << std::endl;
        // std::cout << Eigen::Map<Eigen::Matrix<double, 6, 4, Eigen::RowMajor>> (jaco[2]) << std::endl << std::endl;
        // std::cout << Eigen::Map<Eigen::Matrix<double, 6, 2, Eigen::RowMajor>> (jaco[3]) << std::endl << std::endl;
        // std::cout << Eigen::Map<Eigen::Matrix<double, 6, 1>> (jaco[4]) << std::endl << std::endl;
        // std::cout << Eigen::Map<Eigen::Matrix<double, 6, 4, Eigen::RowMajor>> (jaco[5]) << std::endl << std::endl;
        
        Eigen::Vector3d p_w_i(parameters[0][0], parameters[0][1], parameters[1][0]);
        // Eigen::Quaterniond q_w_i(parameters[1]); // Cannot use this, because Eigen store quaternion as: x, y, z, w. But construct as w, x, y, z
        Eigen::Quaterniond q_w_i(parameters[2][0], parameters[2][1], parameters[2][2], parameters[2][3]);
        Eigen::Vector3d p_w_j(parameters[3][0], parameters[3][1], parameters[4][0]);
        Eigen::Quaterniond q_w_j(parameters[5][0], parameters[5][1], parameters[5][2], parameters[5][3]);
        // std::cout << "W: pi " << p_w_i.transpose() << " qi "<< q_w_i.vec().transpose() << " pj " << p_w_j << " q_j " << q_w_j.vec().transpose() << std::endl;

        Eigen::Matrix<double,6,1> residual;
        residual.block<3,1>(0,0) = sqrt_info_p*(dpos - q_w_i.inverse()*(p_w_j - p_w_i));
        residual.block<3,1>(3,0) = 2.0*sqrt_info_q*(dq_l_inv*(q_w_i.inverse()*q_w_j)).vec();

        // puts("num");
        // std::cout << residual.transpose() << std::endl;
        const double eps = 1e-6;
        Eigen::Matrix<double, 6, 12> num_jacobian;
        for (int k=0; k<12; k++){
            Eigen::Vector3d p_w_i(parameters[0][0], parameters[0][1], parameters[1][0]);
            // Eigen::Quaterniond q_w_i(parameters[1]); // Cannot use this, because Eigen store quaternion as: x, y, z, w. But construct as w, x, y, z
            Eigen::Quaterniond q_w_i(parameters[2][0], parameters[2][1], parameters[2][2], parameters[2][3]);
            Eigen::Vector3d p_w_j(parameters[3][0], parameters[3][1], parameters[4][0]);
            Eigen::Quaterniond q_w_j(parameters[5][0], parameters[5][1], parameters[5][2], parameters[5][3]);
            int a = k / 3, b = k % 3;
            Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2)*eps;
            if (a == 0)
                p_w_i += delta;
            else if (a == 1)
                q_w_i = q_w_i * Utility::deltaQ(delta);
            else if (a == 2)
                p_w_j += delta;
            else if (a == 3)
                q_w_j = q_w_j * Utility::deltaQ(delta);
            
            Eigen::Matrix<double, 6, 1> tmp_residual;
            tmp_residual.block<3,1>(0,0) = sqrt_info_p*(dpos - q_w_i.inverse()*(p_w_j - p_w_i));
            tmp_residual.block<3,1>(3,0) = 2.0*sqrt_info_q*(dq_l_inv*(q_w_i.inverse()*q_w_j)).vec();
            num_jacobian.col(k) = (tmp_residual - residual) / eps;
        }
        // std::cout << num_jacobian << std::endl;
    }

    Eigen::Vector3d dpos;
    Eigen::Matrix3d sqrt_info_p, sqrt_info_q;
    Eigen::Quaterniond dq_l_inv;
};



// class VIOFactors : public ceres::SizedCostFunction<6,3,4,3,4> {
// public:
//     VIOFactors() = delete;
//     VIOFactors(const Eigen::Vector3d &p_i, const Eigen::Quaterniond &q_i, const Eigen::Vector3d &p_j, const Eigen::Quaterniond &q_j, double p_dev, double q_dev)
//         :p_dev(p_dev), q_dev(q_dev){
//         dpos = q_i.inverse()*(p_j - p_i);
//         dq_l_inv = (q_i.inverse()*q_j).inverse();
//     }

//     virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const 
//     {
//         Eigen::Vector3d p_w_i(parameters[0]);
//         Eigen::Quaterniond q_w_i(parameters[1]);
//         Eigen::Vector3d p_w_j(parameters[2]);
//         Eigen::Quaterniond q_w_j(parameters[3]);
//         Eigen::Map<Eigen::Matrix<double,6,1>> residual(residuals);
//         residual.block<3,1>(0,0) = (dpos - q_w_i.inverse()*(p_w_j - p_w_i))/p_dev;
//         residual.block<3,1>(3,0) = (dq_l_inv*(q_w_i.inverse()*q_w_j)).vec()*2/q_dev;

//         if (jacobians)
//     	{
//             if (jacobians[0]){
//                 Eigen::Map<Eigen::Matrix<double, 6, 3, Eigen::RowMajor>> jacobian_p_i(jacobians[0]);
//                 jacobian_p_i.setZero();
//                 jacobian_p_i.block<3,3>(0,0) = q_w_i.inverse().toRotationMatrix()/p_dev;
//             }
            
//             if (jacobians[1]){
//                 Eigen::Map<Eigen::Matrix<double, 6, 4, Eigen::RowMajor>> jacobian_q_i(jacobians[1]);
//                 jacobian_q_i.setZero();
//                 jacobian_q_i.block<3,3>(0,0) = Utility::skewSymmetric(q_w_i.inverse()*(p_w_j-p_w_i))/p_dev;
//                 jacobian_q_i.block<3,3>(3,0) = (Utility::Qright(q_w_i.inverse()*q_w_j)*Utility::Qleft(dq_l_inv)).bottomRightCorner<3,3>() * -1.0/q_dev;
//             }
//             if (jacobians[2]){
//                 Eigen::Map<Eigen::Matrix<double, 6, 3, Eigen::RowMajor>> jacobian_p_j(jacobians[2]);
//                 jacobian_p_j.setZero();
//                 jacobian_p_j.block<3,3>(0,0) = -1.0/p_dev * q_w_i.inverse().toRotationMatrix();
//             }
//             if (jacobians[3]){
//                 Eigen::Map<Eigen::Matrix<double, 6, 4, Eigen::RowMajor>> jacobian_q_j(jacobians[3]);
//                 jacobian_q_j.setZero();
//                 jacobian_q_j.block<3,3>(3,0) = Utility::Qleft(dq_l_inv*q_w_i.inverse()*q_w_j).bottomRightCorner<3,3>() * -1.0/q_dev;
//             }
//         }

//         return true;
//     }

//     void check(double** parameters)
//     {
//         double* res = new double[6];
//         double** jaco = new double* [4];
//         jaco[0] = new double[6*3];
//         jaco[1] = new double[6*4];
//         jaco[2] = new double[6*3];
//         jaco[3] = new double[6*4];

//         Evaluate(parameters, res, jaco);
//         puts("check begins");

//         puts("my");
//         std::cout << Eigen::Map<Eigen::Matrix<double, 6, 1>> (res).transpose() << std::endl << std::endl;
//         std::cout << Eigen::Map<Eigen::Matrix<double, 6, 3, Eigen::RowMajor>> (jaco[0]) << std::endl << std::endl;
//         std::cout << Eigen::Map<Eigen::Matrix<double, 6, 4, Eigen::RowMajor>> (jaco[1]) << std::endl << std::endl;
//         std::cout << Eigen::Map<Eigen::Matrix<double, 6, 3, Eigen::RowMajor>> (jaco[2]) << std::endl << std::endl;
//         std::cout << Eigen::Map<Eigen::Matrix<double, 6, 4, Eigen::RowMajor>> (jaco[3]) << std::endl << std::endl;
        
//         Eigen::Vector3d p_w_i(parameters[0]);
//         Eigen::Quaterniond q_w_i(parameters[1]);
//         Eigen::Vector3d p_w_j(parameters[2]);
//         Eigen::Quaterniond q_w_j(parameters[3]);

//         Eigen::Matrix<double,6,1> residual;
//         residual.block<3,1>(0,0) = (dpos - q_w_i.inverse()*(p_w_j - p_w_i))/p_dev;
//         residual.block<3,1>(3,0) = (dq_l_inv*(q_w_i.inverse()*q_w_j)).vec()*2/q_dev;

//         puts("num");
//         std::cout << residual.transpose() << std::endl;
//         const double eps = 1e-6;
//         Eigen::Matrix<double, 6, 12> num_jacobian;
//         for (int k=0; k<12; k++){
//             int a = k / 3, b = k % 3;
//             Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2)*eps;
//             if (a == 0)
//                 p_w_i += delta;
//             else if (a == 1)
//                 q_w_i = q_w_i * Utility::deltaQ(delta);
//             else if (a == 2)
//                 p_w_j += delta;
//             else if (a == 3)
//                 q_w_j = q_w_j * Utility::deltaQ(delta);
            
//             Eigen::Matrix<double, 6, 1> tmp_residual;
//             tmp_residual.block<3,1>(0,0) = (dpos - q_w_i.inverse()*(p_w_j - p_w_i))/p_dev;
//             tmp_residual.block<3,1>(3,0) = (dq_l_inv*(q_w_i.inverse()*q_w_j)).vec()*2/q_dev;
//             num_jacobian.col(k) = (tmp_residual - residual) / eps;
//         }
//         std::cout << num_jacobian << std::endl;
//     }

//     Eigen::Vector3d dpos;
//     Eigen::Quaterniond dq_l_inv;
//     double p_dev, q_dev;
// };