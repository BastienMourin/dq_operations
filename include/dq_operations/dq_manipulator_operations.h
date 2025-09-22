#ifndef DQ_MANIPULATOR_OPERATIONS_H
#define DQ_MANIPULATOR_OPERATIONS_H

#define _USE_MATH_DEFINES
#include <dq_operations/dq.h>

using namespace Eigen;

/**
 * @brief      This class describes dual quaternion manipulator kinematics.
 */
class DQManipulatorOperations
{
	public:

		/**
		 * @brief compute fkm matrix
		 * @details compute fkm matrix from current joint values(with joint types) from initial screw axes
		 * 
		 * @param joint_type [description]
		 * @return [description]
		 */
		static std::vector<Matrix<double,8,1> > fkmSingleWithJointType(std::vector<RowVector3d> u,std::vector<RowVector3d> m, RowVectorXd q /*(can be revolute or translational joint)*/, std::vector<int> joint_type);


		/**
		 * @brief compute fkm matrix
		 * @details compute fkm matrix from current joint values(dual) from initial screw axes
		 * @return fkm matrix
		 */
		static std::vector<Matrix<double,8,1> > fkmDualAngles(RowVectorXd q /*q is defined as [q1 d1 q2 d2 ...]*/, std::vector<RowVector3d> u, std::vector<RowVector3d> m);

		/**
		 * @brief compute fkm matrix for rotational only serial mechanism
		 * @details compute fkm matrix from current joint values(rotational only) from initial screw axes
		 * @return fkm matrix
		 */
		static std::vector<Matrix<double,8,1> > fkmRevoluteOnly(RowVectorXd q /*q is defined as [q1 q2 ...]*/, std::vector<RowVector3d> u, std::vector<RowVector3d> m);

		/**
		 * @brief get Jacobian for revolute joint manipulator 
		 * @details get Jacobian for revolute joint manipulator from joint displacement q,  and initial screw axes params u and m
		 * 
		 * @param joint displacement q, screw axis direction u and screw axis moment m
		 * @return jacobian matrix
		 */
		static MatrixXd jacobianRevoluteOnly(RowVectorXd q /*q is defined as [q1 q2 ...]*/, std::vector<RowVector3d> u, std::vector<RowVector3d> m);

		/**
		 * @brief 8d jacobians to 6d jacobian
		 * @details [long description]
		 * 
		 * @param jacob_8d [description]
		 * @return 6d jacobian
		 */
		static MatrixXd jacob8dTo6d(MatrixXd jacob_8d);

		/**
		 * @brief converts 6d Jacobian matrix to 8d
		 * @details [long description]
		 * 
		 * @param joint_size [description]
		 * @param jacobian [description]
		 * 
		 * @return [description]
		 */
		static MatrixXd jacobian6dTo8d(int joint_size, MatrixXd jacobian);

		static MatrixXd jacobian_dual_vm(RowVectorXd q /*q is defined as [q1 d1 q2 d2 ...]*/, std::vector<int> joint_type, std::vector<RowVector3d> u, std::vector<RowVector3d> p, std::vector<RowVector3d> m, Matrix<double,8,1> pose_ee_init);


		/**
		 * @brief get Jacobian from current fkm matrix and initial screw axes
		 * @details [long description]
		 * 
		 * @param fkm_current [description]
		 * @return [description]
		 */
		static MatrixXd jacobianDual(std::vector<Matrix<double,8,1> > s, std::vector<Matrix<double,8,1> > fkm_current);

		/**
		 * @brief get link velocities from jacobian6d and joint velocities
		 * @details [long description]
		 * 
		 * @param jacobian_6d [description]
		 * @param joint_velocities [description]
		 * 
		 * @return [description]
		 */
		static MatrixXd linkVelocites(MatrixXd jacobian_6d, RowVectorXd joint_velocities);

		/**
		 * @brief get jacobian derivative from link velocity and jacobian6d
		 * @details [long description]
		 * 
		 * @param link_velocity [description]
		 * @param jacobian_6d [description]
		 * 
		 * @return [description]
		 */
		static MatrixXd getJacobianDot(MatrixXd link_velocity, MatrixXd jacobian_6d);

};
#endif