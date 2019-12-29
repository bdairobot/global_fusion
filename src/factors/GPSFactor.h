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
#include <Eigen/Dense>

class GPSFactor : public ceres::SizedCostFunction<3,2,1>
{
public:
    GPSFactor() = delete;
    GPSFactor(double x, double y, double z, double xy_dev, double z_dev)
        :xyz(x,y,z){
            sqrt_info << 1.0/xy_dev, 0.0, 0.0, 0.0, 1.0/xy_dev, 0.0, 0.0, 0.0, 1.0/z_dev;
        }
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const
    {
        Eigen::Vector3d p_w_i(parameters[0]);
        Eigen::Map<Eigen::Matrix<double, 3, 1>> residual(residuals);
        residual.block<3,1>(0,0) = sqrt_info*(xyz - p_w_i);

        if (jacobians){            
            if (jacobians[0]){
                Eigen::Map<Eigen::Matrix<double, 3, 2, Eigen::RowMajor>> jacobian_pose(jacobians[0]);
                jacobian_pose.block<3,2>(0,0) = -sqrt_info.block<3,2>(0,0);
            }
            if (jacobians[1]){
                Eigen::Map<Eigen::Matrix<double, 3, 1>> jacobian_pose(jacobians[1]);
                jacobian_pose = -sqrt_info.block<3,1>(0,2);
            }
        }
        return true;
    }
    void check(double **parameters) {
        double **jaco = new double*[2];
        jaco[0] = new double[3*2];
        jaco[1] = new double[3*1];
        double *res = new double[3];
        Evaluate(parameters, res, jaco);
        puts("check begins");

        puts("my");
        std::cout << Eigen::Map<Eigen::Matrix<double, 3, 1>> (res).transpose() << std::endl << std::endl;
        std::cout << Eigen::Map<Eigen::Matrix<double, 3, 2, Eigen::RowMajor>> (jaco[0]) << std::endl << std::endl;
        std::cout << Eigen::Map<Eigen::Matrix<double, 3, 1>> (jaco[1]) << std::endl << std::endl;

    }

    Eigen::Vector3d xyz;
    Eigen::Matrix3d sqrt_info;
    double xy_dev, z_dev;
};

// class GPSFactor : public ceres::SizedCostFunction<4,2,2,2>
// {
// public:
//     GPSFactor() = delete;
//     GPSFactor(double x, double y, double m_dev, double b_dev)
//         :x(x), y(y), m_dev(m_dev), b_dev(b_dev){}
//     virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const
//     {
//         Eigen::Vector2d xy(x, y);
//         Eigen::Vector2d bias_i(parameters[0][0], parameters[0][1]);
//         Eigen::Vector2d bias_j(parameters[1][0], parameters[1][1]);
//         Eigen::Vector2d p_w_i(parameters[2][0], parameters[2][1]);
//         Eigen::Map<Eigen::Matrix<double, 4, 1>> residual(residuals);
//         residual.block<2,1>(0,0) = (xy - bias_j - p_w_i)/m_dev;
//         residual.block<2,1>(2,0) = (bias_j - bias_i)/b_dev;

//         if (jacobians){            
//             if (jacobians[0]){
//                 Eigen::Map<Eigen::Matrix<double, 4, 2, Eigen::RowMajor>> jacobian_bias_i(jacobians[0]);
//                 jacobian_bias_i.setZero();
//                 jacobian_bias_i.block<2,2>(2,0) = Eigen::Matrix2d::Identity()/m_dev;
//             }
//             if (jacobians[1]){
//                 Eigen::Map<Eigen::Matrix<double, 4, 2, Eigen::RowMajor>> jacobian_bias_j(jacobians[1]);
//                 jacobian_bias_j.setZero();
//                 jacobian_bias_j.block<2,2>(0,0) = -Eigen::Matrix2d::Identity()/m_dev;
//                 jacobian_bias_j.block<2,2>(2,0) = Eigen::Matrix2d::Identity()/b_dev;
//             }
//             if (jacobians[3]){
//                 Eigen::Map<Eigen::Matrix<double, 4, 2, Eigen::RowMajor>> jacobian_pos(jacobians[1]);
//                 jacobian_pos.setZero();
//                 jacobian_pos.block<2,2>(0,0) = -Eigen::Matrix2d::Identity()/b_dev;
//             }
//         }
//         return true;
//     }

//     void check(double** parameters)
//     {
//         double* res = new double[4];
//         double** jaco = new double* [3];
//         jaco[0] = new double[4*2];
//         jaco[1] = new double[4*2];
//         jaco[2] = new double[4*2];

//         Evaluate(parameters, res, jaco);
//         puts("check begins");
//         puts("my");
//         std::cout << Eigen::Map<Eigen::Matrix<double, 2, 1>>(res).transpose() << std::endl << std::endl;
//         std::cout << Eigen::Map<Eigen::Matrix<double, 4, 2, Eigen::RowMajor>>(jaco[0]) << std::endl << std::endl;
//         std::cout << Eigen::Map<Eigen::Matrix<double, 4, 2, Eigen::RowMajor>>(jaco[1]) << std::endl << std::endl;
//         std::cout << Eigen::Map<Eigen::Matrix<double, 4, 2, Eigen::RowMajor>>(jaco[2]) << std::endl << std::endl;

//         Eigen::Vector2d xy(x, y);
//         Eigen::Vector2d bias_i(parameters[0][0], parameters[0][1]);
//         Eigen::Vector2d bias_j(parameters[1][0], parameters[1][1]);
//         Eigen::Vector2d p_w_i(parameters[2][0], parameters[2][1]);
//         Eigen::Matrix<double, 4, 1> residual;
        
//         residual.block<2,1>(0,0) = (xy - bias_j - p_w_i)/m_dev;
//         residual.block<2,1>(2,0) = (bias_j - bias_i)/b_dev;
//         puts("num");
//         std::cout << residual.transpose() << std::endl;
//         const double eps = 1e-6;
//         Eigen::Matrix<double, 4, 6> num_jacobian;
//         for (int k=0; k<6; k++){
//             int a = k / 2, b = k % 2;
//             Eigen::Vector2d delta = Eigen::Vector2d(b==0, b==1)*eps;
//             if (a == 0)
//                 bias_i += delta;
//             else if(a == 1)
//                 bias_j += delta;
//             else if(a == 2)
//                 p_w_i += delta;
//             Eigen::Matrix<double, 4, 1> tmp_residual;
//             tmp_residual.block<2,1>(0,0) = (xy - bias_j - p_w_i)/m_dev;
//             tmp_residual.block<2,1>(2,0) = (bias_j - bias_i)/b_dev;
//             num_jacobian.col(k) = (tmp_residual - residual) / eps;
//         }
//         std::cout << num_jacobian << std::endl;
//     }

//     double x, y, m_dev, b_dev;
// };