#ifndef DQ_H
#define DQ_H

 #define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <stdexcept> //for range_error
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 
#include <complex>
#include <cmath>
#include <vector>
#include <cstddef>
#include <Eigen/SVD>
#include <time.h>  
using namespace Eigen;

/**
 * @brief      This class describes dual quaternion operations.
 */
class DQ
{


	public:

		/**
		 * @brief      Get translation from unit dual quaternion (UDQ)
		 *
		 * @param  input UDQ
		 *
		 * @return     translation vector 
		 */
		static RowVector3d dqToTranslationAsVector(Matrix<double,8,1> dq);


		/**
		 * @brief      Derive screw unit vector and moment vector from 6-D screw axis in DQ form
		 *
		 * @param  s     6-D screw axis in DQ form 
		 * @param      u     3-D unit screw vector
		 * @param      m     3-D moment vector
		 */
		static void axisMomentFromLineVectorDQ(Matrix<double,8,1> s, RowVector3d& u,RowVector3d& m);

		/**
		 * @brief      Translation vector in quaternion form from UDQ 
		 *
		 * @param  dq     Input pose in UDQ form
		 * @param      trans  Translation vector in quaternion form 
		 */
		static void dqToTranslationAsQuat(Matrix<double,8,1> dq, RowVector4d &trans);


		/**
		 * @brief      3-D Cross product operator 
		 *
		 * @param  vector  The input vector
		 *
		 * @return     3*3 cross product matrix 
		 */
		static Matrix3d crossProductOperator3d(Vector3d vector);

		/**
		 * @brief      6-D Cross product operator related to spatial dynamics 
		 *
		 * @param  vector  The 6-D input vector
		 *
		 * @return     6*6 cross product matrix 
		 */
		static MatrixXd crossProductOperator6D(VectorXd vector);

		/**
		 * @brief      Provides product of two quternions
		 *
		 * @param  p     First quaternion input (order is important!)
		 * @param  q     Second quaternion input  
		 *
		 * @return     Product of quaternions
		 */
		static RowVector4d multiplyQuat(RowVector4d p, RowVector4d q);

		/**
		 * @brief      Provides product of two dual quaternions
		 *
		 * @param  p     First DQ input (order is important!)
		 * @param  q     Second DQ input  
		 *
		 * @return     Product of DQs
		 */
		static Matrix<double,8,1> multiplyDQ(Matrix<double,8,1> p, Matrix<double,8,1> q);

		/**
		 * @brief      Provides conjugate of input quaternion
		 *
		 * @param  p     Input quaternion
		 *
		 * @return     Quaternion Conjugate
		 */
		static RowVector4d conjugateQuat(RowVector4d p);

		/**
		 * @brief      Provides classical conjugate of DQ (conjugate of vector terms)
		 *
		 * @param  dq    Input DQ
		 *
		 * @return     Classical conjugate of DQ
		 */
		static Matrix<double,8,1> classicConjugateDQ(Matrix<double,8,1> dq);

		/**
		 * @brief      Provides dual conjugate of DQ (conjugate of vector terms and then dual term)
		 *
		 * @param  dq    Input DQ
		 *
		 * @return     Dual conjugate of DQ
		 */
		static Matrix<double,8,1> dualConjugateDQ(Matrix<double,8,1> dq);

		/**
		 * @brief      Provides combined conjugate of DQ (conjugate of vector terms and then dual terms)
		 *
		 * @param  dq    Input DQ
		 *
		 * @return     Combined conjugate of DQ
		 */
		static Matrix<double,8,1> combinedConjugateDQ(Matrix<double,8,1> dq);

		/**
		 * @brief      Provides UDQ pose from screw displacement parameters
		 *
		 * @param  theta   Rotation around the screw axis
		 * @param  axis    The axis 3-D vector
		 * @param  d       Translation along the screw axis
		 * @param  moment  The moment 3-D vector
		 *
		 * @return     Pose UDQ related to screw displacement
		 */
		static Matrix<double,8,1> screwDisplacementToDQ(double theta, RowVector3d axis, double d, RowVector3d moment);

		/**
		 * @brief      Returns screw parameters from UDQ. The screwResult DQ contains the screw parameter as [theta d l m]
		 *
		 * @param  dq       Input pose UDQ
		 * @param      theta_e  Rotation around the screw axis
		 * @param      d_e      Translation along the screw axis
		 * @param      l_e      The axis 3-D vector
		 * @param      m_e      The moment 3-D vector
		 */
		static void dqToScrewParameters(Matrix<double,8,1> dq, double &theta_e, double &d_e, RowVector3d &l_e, RowVector3d &m_e); 

		/**
		 * @brief      UDQ to homogeneous transformation matrix (HTM) pose
		 *
		 * @param  Input UDQ pose
		 *
		 * @return     Pose in HTM form
		 */
		static Matrix4d dqToHTM(Matrix<double,8,1> dq);

		/**
		 * @brief      Converts UDQ to rotation and translation (actual distance) in quaternion form
		 *
		 * @param  dq     Input UDQ pose
		 * @param      rot    The rot
		 * @param      trans  The transaction
		 */
		static void dqToRotationAndTranslationAsQuat(Matrix<double,8,1> dq, RowVector4d &rot, RowVector4d &trans);

		/**
		 * @brief      HTM to UDQ pose
		 *
		 * @param  htm   The htm
		 *
		 * @return     UDQ pose
		 */
		static Matrix<double,8,1>  htmToDQ(Matrix4d htm);



		/**
		 * @brief      Converts Rowvector to std::vector
		 *
		 * @param  dq_eigen   The dq eigen
		 * @param      dq_double  The dq double
		 */
		static void dqEigenToDQdoubleVector(RowVectorXd dq_eigen, std::vector<double> &dq_double);

		/**
		 * @brief      Converts std::vector dq to Matrix_8x1 dq 
		 *
		 * @param      dq_eigen   The dq eigen
		 * @param  dq_double  The dq std::vector
		 */
		static void dqDoubleVectorToDQEigen(Matrix<double,8,1> &dq_eigen, std::vector<double> dq_double);

		/**
		 * @brief      Converts dq Rowvector to std::vector
		 *
		 * @param  dq_eigen  The dq eigen
		 *
		 * @return     dq The dq std::vector
		 */
		static std::vector<double> returnDoubleVectorFromRowVector(RowVectorXd dq_eigen);

		/**
		 * @brief      Matrix_8x1 dq to std::vector
		 *
		 * @param  dq_eigen  The dq eigen
		 *
		 * @return     dq std::vector
		 */
		static std::vector<double>  returnDQdoubleVectorFromDQMatrix(Matrix<double,8,1> dq_eigen);

		/**
		 * @brief      Inverse transform dq
		 *
		 * @param  dq    Input pose UDQ
		 *
		 * @return     Output inverse pose UDQ
		 */
		static Matrix<double,8,1> inverseDQtransformation(Matrix<double,8,1> dq); 

		/**
		 * @brief      Transform 3-D line in UDQ form in space with UDQ transformation transform
		 *
		 * @param  line       The line in UDQ form
		 * @param  transform  The transformation UDQ
		 *
		 * @return     Transformed line in UDQ form
		 */
		static Matrix<double,8,1>  transformPluckerLineEigen(Matrix<double,8,1> line, Matrix<double,8,1> transform);

		/**
		 * @brief      Transform 3-D line in Rowvector form in space with UDQ transformation transform
		 *
		 * @param  lineVector  The line in Rowvector form
		 * @param  transform   The transformation UDQ
		 *
		 * @return     Transformed line in 3-D Rowvector form 
		 */
		static RowVectorXd  transformLine3DRowVector(RowVectorXd lineVector, Matrix<double,8,1> transform);

		/**
		 * @brief      Transform 6-D line in Rowvector form in space with UDQ transformation transform
		 *
		 * @param  lineVector  The line vector
		 * @param  transform   The transform
		 *
		 * @return     Transformed line in 6-D Rowvector form 
		 */
		static RowVectorXd  transformPluckerLine6DRowVector(RowVectorXd lineVector, Matrix<double,8,1> transform);

		/**
		 * @brief      Transforms a point in space with UDQ transform
		 *
		 * @param  point      The point
		 * @param  transform  The transform in UDQ form
		 *
		 * @return     Transformed point 
		 */
		static RowVector3d  transformPoint(RowVector3d point, Matrix<double,8,1> transform);

		/**
		 * @brief      normaize angle to value between 0 to To*pi
		 *
		 * @param  theta  The input angle
		 *
		 * @return     Normalized angle
		 */
		static double normalizeAngleZeroTo2Pi(double theta);

		/**
		 * @brief      RowVector to std::vector
		 *
		 * @param  eigenVector  The eigen vector
		 *
		 * @return     RowVector to std::vector
		 */
		static std::vector<double>  returnDQdoubleVectorFromDQRowVector(RowVectorXd eigenVector);
		// static Matrix<double,8,1> sclerp(Matrix<double,8,1> pose_now, Matrix<double,8,1> &pose_intermediate, Matrix<double,8,1> pose_desired, double tau);

		/**
		 * @brief      Matrix_8x1 to RowVector8d
		 *
		 * @param  matrix8d  The matrix 81 UDQ form
		 *
		 * @return     RowVectorX8d
		 */
		static RowVectorXd returnDQRowVectorFromDQMatrix(Matrix<double,8,1> matrix8d);

		/**
		 * @brief      Matrix_8x1 to RowVectorX6D
		 *
		 * @param  matrix8d  The matrix 81 UDQ form
		 *
		 * @return     RowVectorX6D
		 */
		static RowVectorXd returnRowVector6DFromMatrix8D(Matrix<double,8,1> matrix8d);

		/**
		 * @brief       Angular and Linear terms of screw displacement from UDQ
		 *
		 * @param  dq    UDQ transformation
		 *
		 * @return     Angular and Linear terms related to twist [0 w 0 v]
		 */
		static Matrix<double,8,1> dqTotwist(Matrix<double,8,1> dq) ;

		/**
		 * @brief       Linear and Angular terms of screw displacement from UDQ as 6D vector.
		 * 				This 6D screw displacement vectorial form will then be used for lower level Cartesian impedance controller
		 *
		 * @param  dq    UDQ transformation
		 *
		 * @return      Linear and Angular terms related to twist [v w ]
		 */
		static Matrix<double,6,1> dqTotwist_6D(Matrix<double,8,1> dq) ;


		/**
		 * @brief Converts screw displacement dq in twist form to actual dq representation 
		 * 
		 * @param dq_twist screw displacement dq in twist form
		 * @return actual dq representation
		 */
		static Matrix<double,8,1> twistTodq(Matrix<double,8,1> dq_twist);

		/**
		 * @brief      UDQ from Rotation as quaternion and translation vector as quaternion akin to HTM
		 *
		 * @param  rot    The rot
		 * @param  trans  The translation
		 *
		 * @return     UDQ
		 */
		static Matrix<double,8,1> rotationTranslationToDQ(RowVector4d rot, RowVector4d trans);



		/**
		 * @brief      From 6-D Rowvector to 8-D Rowvector
		 *
		 * @param  twist  6-D Rowvector
		 *
		 * @return     8-D Rowvector
		 */
		static RowVectorXd rowVector6DToRowVector8D(RowVectorXd twist);

		/**
		 * @brief      From 8-D Rowvector to 6-D Rowvector
		 *
		 * @param  DQtwist  8-D Rowvector
		 *
		 * @return     6-D Rowvector
		 */
		static RowVectorXd rowVector8DToRowVector6D(RowVectorXd DQtwist);


		/**
		 * @brief      Rowvector 6D to Matrix_8x1 UDQ
		 *
		 * @param  rowVector6D  The row vector 6 d
		 *
		 * @return     Matrix_8x1 UDQ
		 */
		static Matrix<double,8,1> rowVector6DToMatrix8D(RowVectorXd rowVector6D);



		/**
		 * @brief      std::vector UDQ to Matrix_8x1 UDQ
		 *
		 * @param  dq_double  The std::vector UDQ
		 *
		 * @return     Matrix_8x1 UDQ
		 */
		static Matrix<double,8,1> returnMatrixDQFromDoubleVectorDQ(std::vector<double> dq_double);


		/**
		 * @brief      Gets position from UDQ pose.
		 *
		 * @param  pose  The pose
		 *
		 * @return     The position from dq.
		 */
		static RowVector3d getPositionFromDQ(Matrix<double,8,1>  pose);

		/**
		 * @brief gets position and orientaion error from DQ poses for KDL operations
		 * 
		 * @param current_pose 
		 * @param desired_pose 
		 * @return 6d pose error for kdl (omega;v)
		 */
		static RowVectorXd spatial2CartPoseError(Matrix<double,8,1>  desired_pose, Matrix<double,8,1>  current_pose);

		DQ();
		~DQ();
};
#endif