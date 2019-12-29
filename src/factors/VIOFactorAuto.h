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

struct VIOFactorAuto
{
	VIOFactorAuto(double t_x, double t_y, double t_z, 
					double q_w, double q_x, double q_y, double q_z,
					double p_dev, double q_dev)
				  :t_x(t_x), t_y(t_y), t_z(t_z), 
				   q_w(q_w), q_x(q_x), q_y(q_y), q_z(q_z),
				   p_dev(p_dev), q_dev(q_dev){}

	template <typename T>
	bool operator()(const T* q_w_i, const T* p_i_xy, const T* p_i_z, const T* q_w_j, const T* p_j_xy, const T* p_j_z, T* residuals) const
	{
		T t_w_ij[3];
		t_w_ij[0] = p_j_xy[0] - p_i_xy[0];
		t_w_ij[1] = p_j_xy[1] - p_i_xy[1];
		t_w_ij[2] = p_j_z[0] - p_i_z[0];

		T q_i_w[4] = {q_w_i[0], -q_w_i[1], -q_w_i[2], -q_w_i[3]};
		T t_i_ij[3];
		ceres::QuaternionRotatePoint(q_i_w, t_w_ij, t_i_ij);

		residuals[0] = (t_i_ij[0] - T(t_x)) / T(p_dev);
		residuals[1] = (t_i_ij[1] - T(t_y)) / T(p_dev);
		residuals[2] = (t_i_ij[2] - T(t_z)) / T(p_dev);

		T relative_q[4] = {T(q_w), T(q_x), T(q_y), T(q_z)};
		T q_i_j[4];
		ceres::QuaternionProduct(q_i_w, q_w_j, q_i_j);
		T relative_q_inv[4] = {relative_q[0], -relative_q[1], -relative_q[2], -relative_q[3]};
		T error_q[4];
		ceres::QuaternionProduct(relative_q_inv, q_i_j, error_q); 

		residuals[3] = T(2) * error_q[1] / T(q_dev);
		residuals[4] = T(2) * error_q[2] / T(q_dev);
		residuals[5] = T(2) * error_q[3] / T(q_dev);

		return true;
	}

	static ceres::CostFunction* Create(const double t_x, const double t_y, const double t_z,
									   const double q_w, const double q_x, const double q_y, const double q_z,
									   const double p_dev, const double q_dev) 
	{
	  return (new ceres::AutoDiffCostFunction<
	          VIOFactorAuto, 6, 4, 2, 1, 4, 2, 1>(
	          	new VIOFactorAuto(t_x, t_y, t_z, q_w, q_x, q_y, q_z, p_dev, q_dev)));
	}

	double t_x, t_y, t_z;
	double q_w, q_x, q_y, q_z;
	double p_dev, q_dev;

};

struct VIOSimFactorAuto
{
	VIOSimFactorAuto(double t_x, double t_y, double t_z, 
					double q_w, double q_x, double q_y, double q_z,
					double p_dev, double q_dev)
				  :t_x(t_x), t_y(t_y), t_z(t_z), 
				   q_w(q_w), q_x(q_x), q_y(q_y), q_z(q_z),
				   p_dev(p_dev), q_dev(q_dev){}

	template <typename T>
	bool operator()(const T* q_w_i, const T* p_i_xy, const T* p_i_z, const T* q_w_j, const T* p_j_xy, const T* p_j_z, const T* scale, T* residuals) const
	{
		T t_w_ij[3];
		t_w_ij[0] = p_j_xy[0] - p_i_xy[0];
		t_w_ij[1] = p_j_xy[1] - p_i_xy[1];
		t_w_ij[2] = p_j_z[0] - p_i_z[0];

		T q_i_w[4] = {q_w_i[0], -q_w_i[1], -q_w_i[2], -q_w_i[3]};
		T t_i_ij[3];
		ceres::QuaternionRotatePoint(q_i_w, t_w_ij, t_i_ij);

		residuals[0] = (t_i_ij[0] - scale[0] * T(t_x)) / T(p_dev);
		residuals[1] = (t_i_ij[1] - scale[0] * T(t_y)) / T(p_dev);
		residuals[2] = (t_i_ij[2] - scale[0] * T(t_z)) / T(p_dev);

		T relative_q[4] = {T(q_w), T(q_x), T(q_y), T(q_z)};
		T q_i_j[4];
		ceres::QuaternionProduct(q_i_w, q_w_j, q_i_j);
		T relative_q_inv[4] = {relative_q[0], -relative_q[1], -relative_q[2], -relative_q[3]};
		T error_q[4];
		ceres::QuaternionProduct(relative_q_inv, q_i_j, error_q); 

		residuals[3] = T(2) * error_q[1] / T(q_dev);
		residuals[4] = T(2) * error_q[2] / T(q_dev);
		residuals[5] = T(2) * error_q[3] / T(q_dev);

		return true;
	}

	static ceres::CostFunction* Create(const double t_x, const double t_y, const double t_z,
									   const double q_w, const double q_x, const double q_y, const double q_z,
									   const double p_dev, const double q_dev) 
	{
	  return (new ceres::AutoDiffCostFunction<
	          VIOSimFactorAuto, 6, 4, 2, 1, 4, 2, 1, 1>(
	          	new VIOSimFactorAuto(t_x, t_y, t_z, q_w, q_x, q_y, q_z, p_dev, q_dev)));
	}

	double t_x, t_y, t_z;
	double q_w, q_x, q_y, q_z;
	double p_dev, q_dev;

};