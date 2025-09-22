#include <dq_operations/dq_dual_arm_operations.h>
using namespace Eigen;

DQDualArm::DQDualArm(){}
DQDualArm::~DQDualArm(){}



Matrix<double,8,1> DQDualArm::getAbsolutePose(Matrix<double,8,1> rightPoseNow, Matrix<double,8,1> leftPoseNow)
{
	///
	/// Pseudo code to get absolute pose
	///
	Matrix<double,8,1> p_ee_halfRot, pose_rel, pose_abs;
	p_ee_halfRot << 1, 0, 0, 0, 0, 0, 0, 0;
	pose_rel << 1, 0, 0, 0, 0, 0, 0, 0 ;	
	pose_abs << 1, 0, 0, 0, 0, 0, 0, 0 ;	
	/**
	 * 1. get relative pose in reference arm (right arm) frame
	 */
	pose_rel = DQ::multiplyDQ(DQ::classicConjugateDQ(rightPoseNow), leftPoseNow);

	/**
	 * 2. Now, get the screw parameters for rel_dq
	 * 3. check theta: should be positive
	 */
	double theta_rel, d_rel;
	RowVector3d l_rel, m_rel, position_vec, position_vec_2;
	RowVector4d position_rel, quat_dual;	
	DQ::dqToScrewParameters(pose_rel, theta_rel, d_rel, l_rel, m_rel); 
	// std::cout << "l_rel: " << l_rel << std::endl;
	// std::cout << "theta_rel: " << theta_rel << std::endl;
	/**
	 * 4. get quaternion corresponding to the half rotation around same screw
	 *    axis
	 */
	pose_abs = DQ::screwDisplacementToDQ(theta_rel/2, l_rel, d_rel, m_rel);
	/**
	 * 5. get relative position vector
	 */
	position_vec = DQ::getPositionFromDQ(pose_rel);
	position_rel << 0, position_vec(0), position_vec(1), position_vec(2);
	// std::cout << "position_rel: " << position_rel << std::endl;	
	/**
	 * Matrix4d htm_abs = (DQ::dq2HTM(leftPoseNow) +
	 * DQ::dq2HTM(rightPoseNow)); std::cout << "position_rel: " <<
	 * ", " << htm_abs(0,3) << ", " << htm_abs(1,3) << ", " << htm_abs(2,3) <<
	 * std::endl; / / 6. get the absolute pose from half rotation quaternion and
	 * half relative position vector in reference frame / 
	 */
	quat_dual = 0.5*DQ::multiplyQuat(position_rel/2, pose_abs.head(4));
	// std::cout << "pose_abs: " << pose_abs << std::endl ;
	// std::cout << "quat_dual: " << quat_dual << std::endl; 
	pose_abs << pose_abs.head(4), quat_dual.transpose();
	/**
	 * 7. Get the same in the base frame
	 */
	pose_abs = DQ::multiplyDQ(rightPoseNow, pose_abs);
	// std::cout << "pose_abs_htm: " << DQ::dq2HTM(pose_abs) << std::endl ;
	return pose_abs;
}

Matrix<double,8,1> DQDualArm::getAbsolutePoseInRefFrame(Matrix<double,8,1> pose_rel)
{
	/**
	 * 1. Relative pose already available
	 * 2. Now, get the screw parameters for rel_dq
	 * 3. check theta: should be positive
	 */
	double theta_rel, d_rel;
	RowVector3d l_rel, m_rel, position_vec, position_vec_2;
	RowVector4d position_rel, quat_dual;	
	DQ::dqToScrewParameters(pose_rel, theta_rel, d_rel, l_rel, m_rel); 
	// std::cout << "l_rel: " << l_rel << std::endl;
	// std::cout << "theta_rel: " << theta_rel << std::endl;
	/**
	 * 4. get quaternion corresponding to the half rotation around same screw
	 *    axis
	 */
	Matrix<double,8,1> pose_abs;
	pose_abs = DQ::screwDisplacementToDQ(theta_rel/2, l_rel, d_rel, m_rel);
	/**
	 * 5. get relative position vector
	 */
	position_vec = DQ::getPositionFromDQ(pose_rel);
	position_rel << 0, position_vec(0), position_vec(1), position_vec(2);
	// std::cout << "position_rel: " << position_rel << std::endl;	
	/**
	 * 6. get the absolute pose from half rotation quaternion and
	 * half relative position vector in reference frame / 
	 */
	quat_dual = 0.5*DQ::multiplyQuat(position_rel/2, pose_abs.head(4));
	// std::cout << "pose_abs: " << pose_abs << std::endl ;
	// std::cout << "quat_dual: " << quat_dual << std::endl; 
	pose_abs << pose_abs.head(4), quat_dual.transpose();
	return pose_abs;
}




std::vector<Matrix<double,8,1> > DQDualArm::getRightLeft4mAbsRel(Matrix<double,8,1> absPose, Matrix<double,8,1> relPose)
{
	double theta_rel, d_rel;
	RowVector3d l_rel, m_rel;
	Matrix<double,8,1> p_ee_halfRot, pose_abs_2;
	p_ee_halfRot << 1, 0, 0, 0, 0, 0, 0, 0 ;
	pose_abs_2 << 1, 0, 0, 0, 0, 0, 0, 0 ;
	std::vector<Matrix<double,8,1> > ee_poses;
	ee_poses.resize(2);
	//	First making the scalar part of orientation dq positive. 
	if (relPose(0,0) < 0)
	{
		relPose = -relPose;
		p_ee_halfRot = relPose;
		// p_ee_halfRot = DQ::multiplyDQ(absPose, p_ee_halfRot);
		// std::cout << "p_ee_halfRot: " << p_ee_halfRot << std::endl;		
	}

	if(relPose(0,0) > 0.999995)
	{
		relPose << 1, 0, 0, 0, relPose(4,0), relPose(5,0), relPose(6,0), relPose(7,0);
		p_ee_halfRot = relPose;
		// p_ee_halfRot = DQ::multiplyDQ(absPose, p_ee_halfRot);
		// std::cout << "p_ee_halfRot: " << p_ee_halfRot << std::endl;				
	}
	else 	
	{
		DQ::dqToScrewParameters(relPose, theta_rel, d_rel, l_rel, m_rel); 
		// p_ee_halfRot = DQ::screwDisplacementToDQ(theta_rel/2, l_rel, 0, m_rel);
		// pose_abs_2 = DQ::screwDisplacementToDQ((2*M_PI-theta_rel)/2, -l_rel, 0, m_rel);
		if ((theta_rel/2) > ((2*M_PI-theta_rel)/2))
		{
			theta_rel = 2*M_PI-theta_rel;
			l_rel = -l_rel;
			m_rel = -m_rel;
		}
		p_ee_halfRot = DQ::classicConjugateDQ(DQ::screwDisplacementToDQ(theta_rel/2, l_rel, 0, m_rel));
		// p_ee_halfRot = DQ::multiplyDQ(absPose, p_ee_halfRot);
		// std::cout << "p_ee_halfRot: " << p_ee_halfRot << std::endl;
	}


	RowVector4d distance_rel = 2*DQ::multiplyQuat(relPose.tail(4), DQ::conjugateQuat(relPose.head(4))); 

	RowVector4d position_right = -distance_rel/2;

	position_right = DQ::multiplyQuat(position_right, p_ee_halfRot.head(4))/2;
	ee_poses[0] << p_ee_halfRot(0), p_ee_halfRot(1), p_ee_halfRot(2), p_ee_halfRot(3), position_right(0), position_right(1), position_right(2), position_right(3);
	ee_poses[0] = DQ::multiplyDQ(absPose, ee_poses[0]);
	ee_poses[1] = DQ::multiplyDQ(ee_poses[0], relPose);

	return ee_poses;

}