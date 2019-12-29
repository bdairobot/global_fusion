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

class BaroFactor : public ceres::SizedCostFunction<1,1>
{
public:
    BaroFactor() = delete;
    BaroFactor(double z, double m_dev)
        :z(z), m_dev(m_dev){}
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const{
        double pos_w_i = parameters[0][0];
        residuals[0] = (z - pos_w_i) / m_dev;
        if (jacobians){
            if (jacobians[0]){
                jacobians[0][0] = -1.0 / m_dev;
            }
        }
        return true;
    }
    double z, m_dev;
};

// class BaroFactor : public ceres::SizedCostFunction<2,1,1,1>{
// public:
//     BaroFactor() = delete;
//     BaroFactor(double z, double m_dev, double b_dev)
//         :z(z), m_dev(m_dev), b_dev(b_dev){}
//     virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const{
//         double bias_i = parameters[0][0];
//         double bias_j = parameters[1][0];
//         double pos_w_j = parameters[2][0];
//         residuals[0] = (z - bias_j - pos_w_j) / m_dev;
//         residuals[1] = (bias_j - bias_i) / b_dev;
//         if (jacobians){
//             if (jacobians[0]){
//                 jacobians[0][0] = 0.0;
//                 jacobians[0][1] = -1.0/b_dev;
//             }
//             if (jacobians[1]){
//                 jacobians[1][0] = -1.0/m_dev;
//                 jacobians[1][0] = 1.0/b_dev;
//             }
//             if (jacobians[2]){
//                 jacobians[2][0] = -1.0/m_dev;
//                 jacobians[2][0] = 0.0;
//             }
//         }
//         return true;
//     }
//     double z, m_dev, b_dev;
// };