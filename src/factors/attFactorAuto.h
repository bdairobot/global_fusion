/*******************************************************
 * 
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *
 * Author: Dai Bo (bdairobot@gmail.com)
 *******************************************************/

#pragma once
#include <ceres/ceres.h>
#include <ceres/rotation.h>

struct attFactorAuto
{
	attFactorAuto(double q_w, double q_x, double q_y, double q_z, double q_dev)
				  :q_w(q_w), q_x(q_x), q_y(q_y), q_z(q_z),q_dev(q_dev){}

	template <typename T>
	bool operator()(const T* const q_w_i, const T* const decl, T* residuals) const
	{
        T q_inv[4] = {T(q_w), T(-q_x), T(-q_y), T(-q_z)};
        T decl_q[4] = {cos(decl[0]/T(2)), T(0.0), T(0.0), -sin(decl[0]/T(2))};
        T real_q_inv[4];
        ceres::QuaternionProduct(q_inv, decl_q, real_q_inv);

		T error_q[4];
		ceres::QuaternionProduct(real_q_inv, q_w_i, error_q);

		residuals[0] = T(2) * error_q[1] / T(q_dev);
		residuals[1] = T(2) * error_q[2] / T(q_dev);
		residuals[2] = T(2) * error_q[3] / T(q_dev);

		return true;
	}

	static ceres::CostFunction* Create(const double q_w, const double q_x, const double q_y, const double q_z, 
                                    const double q_dev) 
	{
	  return (new ceres::AutoDiffCostFunction<
	          attFactorAuto, 3, 4, 1>(
	          	new attFactorAuto(q_w, q_x, q_y, q_z, q_dev)));
	}

	double q_w, q_x, q_y, q_z;
	double q_dev;

};