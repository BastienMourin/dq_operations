#include <dq_operations/dq_manipulator_operations.h>

std::vector<Matrix<double,8,1> > DQManipulatorOperations::fkmDualAngles(RowVectorXd q_dual /*q is defined as [q1 d1 q2 d2 ...]*/, std::vector<RowVector3d> u, std::vector<RowVector3d> m)
{
	int joint_size=q_dual.size()/2;
	std::vector<Matrix<double,8,1> > screwDispalcementArray;
	screwDispalcementArray.clear();
	screwDispalcementArray.resize(joint_size);


	for (int i=0;i<joint_size; i++)
	{
		// double theta_i=q_dual[2*i];
		// RowVector3d u_i=u[i];
		// RowVector3d m_i=m[i];
		// double d_i=q_dual[2*i+1];
		Matrix<double,8,1> screwDispalcementArray_i;
		screwDispalcementArray_i= DQ::screwDisplacementToDQ(q_dual[2*i], u[i], q_dual[2*i+1], m[i]);
		if (i==0)
			screwDispalcementArray[i]=screwDispalcementArray_i;
		else
			screwDispalcementArray[i]=DQ::multiplyDQ(screwDispalcementArray[i-1],screwDispalcementArray_i);
	}

	return screwDispalcementArray;	
}


std::vector<Matrix<double,8,1> > DQManipulatorOperations::fkmSingleWithJointType(std::vector<RowVector3d> u,std::vector<RowVector3d> m, RowVectorXd q, std::vector<int> joint_type)
{
	int joint_size=joint_type.size();
	RowVectorXd q_dual;	
	q_dual=RowVectorXd::Zero(joint_size*2);
	for (int j=0; j< joint_size; j++)
	{
		if(joint_type[j]==0)
			q_dual(2*j)=q(j);
		else
			q_dual(2*j+1)=q(j);
		// ROS_INFO("j=%d", j);
	}
	return DQManipulatorOperations::fkmDualAngles(q_dual, u, m);	

}



std::vector<Matrix<double,8,1> > DQManipulatorOperations::fkmRevoluteOnly(RowVectorXd q /*q is defined as [q1 q2 ...]*/, std::vector<RowVector3d> u, std::vector<RowVector3d> m)
{
	int joint_size=q.size();
	RowVectorXd q_dual;	
	q_dual=RowVectorXd::Zero(joint_size*2);
	for (int j=0; j< joint_size; j++)
	{
		q_dual(2*j)=q(j);
		q_dual(2*j+1)=0;
	}	

	return DQManipulatorOperations::fkmDualAngles(q_dual, u, m);		
}

MatrixXd DQManipulatorOperations::jacobian_dual_vm(RowVectorXd q /*q is defined as [q1 d1 q2 d2 ...]*/, std::vector<int> joint_type, std::vector<RowVector3d> u, std::vector<RowVector3d> p , std::vector<RowVector3d> m, Matrix<double,8,1> pose_ee_init)
{
	int joint_size = u.size();
	RowVectorXd q_dual;	
	q_dual=RowVectorXd::Zero(joint_size*2);
	for (int j=0; j< joint_size; j++)
	{
		if(joint_type[j]==0)
			q_dual(2*j)=q(j);
		else
			q_dual(2*j+1)=q(j);
	}

	MatrixXd jacobian_result(6,joint_size);

	std::vector<Matrix<double,8,1> > screwDispalcementArray;
	screwDispalcementArray.clear();
	screwDispalcementArray.resize(joint_size);

	screwDispalcementArray=	DQManipulatorOperations::fkmDualAngles(q_dual, u, m);
	Matrix<double,8,1> screwDispalcementArray_i, screw_axis_i, pose_i;
	
	screwDispalcementArray_i= screwDispalcementArray[joint_size-1];/*The numbering in C++ starts at 0*/
	
	screwDispalcementArray_i=DQ::multiplyDQ(DQ::multiplyDQ(screwDispalcementArray_i, pose_ee_init), DQ::combinedConjugateDQ(screwDispalcementArray_i));
	
	RowVector3d pose_ee_now;
	Matrix4d htm_ee_now;
	pose_ee_now << screwDispalcementArray_i(5), screwDispalcementArray_i(6), screwDispalcementArray_i(7); 
	if(joint_type[0]==0)
		jacobian_result.col(0)<< (p[0].cross(u[0])- pose_ee_now.cross(u[0])).transpose(), u[0].transpose();/*writing Jacobian seperately for first joint for faster operation*/
	else
		jacobian_result.col(0)<< u[0].transpose(), 0, 0, 0;

	for(int i=1; i<joint_size; i++)
	{
		screwDispalcementArray_i=screwDispalcementArray[i-1];
	
		screw_axis_i<< 0, u[i].transpose(), 0, (p[i].cross(u[i])).transpose();
	
		screw_axis_i=DQ::multiplyDQ(DQ::multiplyDQ(screwDispalcementArray_i, screw_axis_i), DQ::classicConjugateDQ(screwDispalcementArray_i));	
	
		RowVector3d u_i, p_i;
		u_i << screw_axis_i(1), screw_axis_i(2), screw_axis_i(3);
		pose_i << 1, 0, 0, 0, 0, p[i].transpose();
		pose_i= DQ::multiplyDQ(DQ::multiplyDQ(screwDispalcementArray_i, pose_i), DQ::combinedConjugateDQ(screwDispalcementArray_i));
		p_i << pose_i(5), pose_i(6), pose_i(7);
		if(joint_type[i]==0)
			jacobian_result.col(i) << (p_i.cross(u_i)- pose_ee_now.cross(u_i)).transpose(), u_i.transpose(); 	
		else
			jacobian_result.col(i) << u_i.transpose(), 0, 0, 0;
	}

	return jacobian_result;
}

MatrixXd DQManipulatorOperations::jacobianDual(std::vector<Matrix<double,8,1>> s, std::vector<Matrix<double,8,1> > fkm_current)
{
	int joint_size =s.size();
	MatrixXd jacobian, jacobian6d;
	jacobian =MatrixXd::Zero(8,joint_size);
	jacobian6d =MatrixXd::Zero(6,joint_size);
	Matrix<double,8,1> screwDispalcementArray_i, screw_axis_i, pose_i;
	
	jacobian.col(0) = s[0];
	for(int i=1; i<joint_size; i++)
	{
		screwDispalcementArray_i=fkm_current[i-1];	
		screw_axis_i = s[i];
		screw_axis_i=DQ::multiplyDQ(DQ::multiplyDQ(screwDispalcementArray_i, screw_axis_i), DQ::classicConjugateDQ(screwDispalcementArray_i));	
		jacobian.col(i) = screw_axis_i.transpose();
	}
	jacobian6d = DQManipulatorOperations::jacob8dTo6d(jacobian);
	return jacobian6d;
}



MatrixXd DQManipulatorOperations::jacobianRevoluteOnly(RowVectorXd q /*q is defined as [q1 q2 ...]*/, std::vector<RowVector3d> u, std::vector<RowVector3d> m)
{
	int joint_size=q.size();

	RowVectorXd q_dual;	
	q_dual=RowVectorXd::Zero(joint_size*2);
	for (int j=0; j< joint_size; j++)
	{
		q_dual(2*j)=q(j);
		q_dual(2*j+1)=0;
	}
	
	std::vector<Matrix<double,8,1>> s;
	s.resize(joint_size);

	for(int i=1; i<joint_size; i++)
	{
		s[i] << 0, u[i](0), u[i](1), u[i](2), 0, m[i](0), m[i](1), m[i](2); 
	}
	
	std::vector<Matrix<double,8,1> > screwDispalcementArray;
	screwDispalcementArray.clear();
	screwDispalcementArray.resize(joint_size);
	screwDispalcementArray=	DQManipulatorOperations::fkmRevoluteOnly(q, u, m);
	
	return jacobianDual(s, screwDispalcementArray);
}

MatrixXd DQManipulatorOperations::getJacobianDot(MatrixXd link_velocity, MatrixXd jacobian_6d)
{
	int joint_size = jacobian_6d.cols();
	MatrixXd jacobian_6d_dot = MatrixXd::Zero(6, joint_size);
	for (int i=0;i<joint_size; i++)
	{
		jacobian_6d_dot.col(i) = DQ::crossProductOperator6D(link_velocity.col(i))*jacobian_6d.col(i);		
	} 
	return jacobian_6d_dot;
}

MatrixXd DQManipulatorOperations::linkVelocites(MatrixXd jacobian_6d, RowVectorXd joint_velocities)
{
	// link_velocities.col(0) refers to first link velocity and so on, because we assume the joint is fixed with child links
	int joint_size = jacobian_6d.cols();
	MatrixXd link_velocities = MatrixXd::Zero(6, joint_size);
	for (int i=0;i<joint_size; i++)
	{
		if(i == 0)
			link_velocities.col(i) = jacobian_6d.col(i)*joint_velocities(i);
		else
			link_velocities.col(i) = link_velocities.col(i-1) + jacobian_6d.col(i)*joint_velocities(i);
	}
	return link_velocities;
}

MatrixXd DQManipulatorOperations::jacob8dTo6d(MatrixXd jacob_8d)
{
	MatrixXd jacobian_6d;
	jacobian_6d = MatrixXd::Zero(6,jacob_8d.cols());
	jacobian_6d.block(0, 0, 3, jacob_8d.cols()) = jacob_8d.block(1, 0, 3, jacob_8d.cols());
	jacobian_6d.block(3, 0, 3, jacob_8d.cols()) = jacob_8d.block(5, 0, 3, jacob_8d.cols());	
	return jacobian_6d;
}

MatrixXd DQManipulatorOperations::jacobian6dTo8d(int joint_size, MatrixXd jacobian)
{
	MatrixXd jacobian_8d;
	jacobian_8d=MatrixXd::Zero(8,joint_size);
    jacobian_8d.row(1)=jacobian.row(0);
	jacobian_8d.row(2)=jacobian.row(1);
	jacobian_8d.row(3)=jacobian.row(2);
	jacobian_8d.row(5)=jacobian.row(3);
	jacobian_8d.row(6)=jacobian.row(4);
	jacobian_8d.row(7)=jacobian.row(5);
	return jacobian_8d;
}