#include <dq_operations/dq.h>
using namespace Eigen;

RowVector4d DQ::multiplyQuat(RowVector4d p, RowVector4d q)
{
	double s1=	p(0);
	double s2=	q(0);
	RowVector3d v1= RowVector3d(p(1), p(2), p(3));
	RowVector3d v2= RowVector3d(q(1), q(2), q(3));
	double term1= (s1*s2-v1.dot(v2));
	RowVector3d term2=s1*v2+s2*v1 +v1.cross(v2);
	RowVector4d result=RowVector4d(term1, term2(0), term2(1), term2(2));
	return result;
}


Matrix<double,8,1> DQ::multiplyDQ(Matrix<double,8,1> p, Matrix<double,8,1> q)
{
	RowVector4d p1, p2, q1, q2;
	p1<< p(0), p(1), p(2), p(3);	
	p2<< p(4), p(5), p(6), p(7);	
	q1<< q(0), q(1), q(2), q(3);	
	q2<< q(4), q(5), q(6), q(7);	
	
	Matrix<double,8,1> result;
	result << multiplyQuat(p1, q1).transpose(), (multiplyQuat(p1, q2)+multiplyQuat(p2, q1)).transpose();

	return result;
}

RowVector4d DQ::conjugateQuat(RowVector4d p)
{
	RowVector4d result;
	result=RowVector4d(p(0),-p(1),-p(2),-p(3));
	return result;
}

Matrix<double,8,1> DQ::classicConjugateDQ(Matrix<double,8,1> dq)
{
	Matrix<double,8,1> dq_result;
	dq_result << dq(0), -dq(1), -dq(2), -dq(3), dq(4), -dq(5), -dq(6), -dq(7);
	return dq_result;
}

Matrix<double,8,1> DQ::dualConjugateDQ(Matrix<double,8,1> dq)
{
	Matrix<double,8,1> dq_result;	
	dq_result<<dq(0), dq(1), dq(2), dq(3), -dq(4), -dq(5), -dq(6), -dq(7);
	return dq_result;
}


Matrix<double,8,1> DQ::combinedConjugateDQ(Matrix<double,8,1> dq)
{
	Matrix<double,8,1> dq_result;	
	dq_result<<dq(0), -dq(1), -dq(2), -dq(3), -dq(4), dq(5), dq(6), dq(7);	
	return dq_result;
}

RowVector3d DQ::dqToTranslationAsVector(Matrix<double,8,1> dq)
{
	// translation then rotation.
	RowVector3d trans;
	RowVector4d rot, trans4d;
	rot << dq(0), dq(1), dq(2), dq(3);
	trans4d << dq(4), dq(5), dq(6), dq(7);
	trans4d = 2*multiplyQuat(trans4d, conjugateQuat(rot));
	trans << trans4d(1), trans4d(2), trans4d(3); 
	return trans;
}

void DQ::axisMomentFromLineVectorDQ(Matrix<double,8,1> s, RowVector3d& u,RowVector3d& m)
{
	u << s(1), s(2), s(3); 
	m << s(5), s(6), s(7);
}	

void DQ::dqToTranslationAsQuat(Matrix<double,8,1> dq, RowVector4d &trans)
{
	// translation then rotation.
	RowVector4d rot;
	rot << dq(0), dq(1), dq(2), dq(3);
	trans << dq(4), dq(5), dq(6), dq(7);
	trans= 2*multiplyQuat(trans, conjugateQuat(rot));
}


Matrix3d DQ::crossProductOperator3d(Vector3d vector)
{
	Matrix3d crossProductOp = MatrixXd::Zero(3, 3);
	crossProductOp << 0., -vector(2), vector(1),
      				  vector(2), 0., -vector(0),
                      -vector(1), vector(0), 0.;
                      
	return crossProductOp;
}

MatrixXd DQ::crossProductOperator6D(VectorXd vector)
{
	MatrixXd cross_6D = MatrixXd::Zero(6, 6);
	Vector3d w;
	w << vector(0), vector(1), vector(2);
	Vector3d v0;
	v0 << vector(3), vector(4), vector(5);

	cross_6D.block(0, 0, 3, 3) = crossProductOperator3d(w);
	cross_6D.block(3, 0, 3, 3) = crossProductOperator3d(v0);
	cross_6D.block(3, 3, 3, 3) = crossProductOperator3d(w);
	return cross_6D;
}

RowVectorXd DQ::spatial2CartPoseError(Matrix<double,8,1>  desired_pose, Matrix<double,8,1>  current_pose)
{
	Matrix<double,8,1> pose_error_dq;
	pose_error_dq = DQ::multiplyDQ(desired_pose, DQ::classicConjugateDQ(current_pose));
	RowVectorXd kdl_error_twist = RowVectorXd::Zero(6);
	Matrix4d htm_error = DQ::dqToHTM(desired_pose) - DQ::dqToHTM(current_pose);
	RowVector3d v_e, w_e, l_e, m_e;
	double d_e = 0, theta_e = 0;
	DQ::dqToScrewParameters(pose_error_dq, theta_e, d_e, l_e, m_e);
	if (theta_e > M_PI)
	{
		desired_pose=-desired_pose;
		pose_error_dq = DQ::multiplyDQ(desired_pose, DQ::classicConjugateDQ(current_pose));
		DQ::dqToScrewParameters(pose_error_dq, theta_e, d_e, l_e, m_e);			
	}	
	w_e = theta_e*l_e;
	kdl_error_twist << w_e(0), w_e(1), w_e(2), htm_error(0,3), htm_error(1,3), htm_error(2,3);
	return kdl_error_twist;
}


Matrix<double,8,1> DQ::screwDisplacementToDQ(double theta, RowVector3d axis, double d, RowVector3d moment)
{
	RowVector4d q_rot, q_tr;
	q_rot << cos(theta/2),
			sin(theta/2)*axis;

	q_tr << -(d/2)*sin(theta/2),
			(d/2)*cos(theta/2)*axis+sin(theta/2)*moment;

	Matrix<double, 8,1> dq;
	dq << q_rot(0), q_rot(1), q_rot(2), q_rot(3), q_tr(0), q_tr(1), q_tr(2), q_tr(3); 

	return dq;
}

void DQ::dqToScrewParameters(Matrix<double,8,1> dq, double &theta_e, double &d_e, RowVector3d &l_e, RowVector3d &m_e)
{
	double eps_theta=0.0001; /*0.1 degrees*/ 
	double sr, sd, theta;
	RowVector3d vr, vd, vec_translation;
	RowVector4d q_rot, q_trans;
	DQ::dqToRotationAndTranslationAsQuat(dq, q_rot, q_trans);
	//vec_translation << q_trans(0), q_trans(1), q_trans(2);
	vec_translation << q_trans(1), q_trans(2), q_trans(3);

	// if(dq(0,0)<0)
	// {
	// 	dq=-dq;
	// }
	double norm_quat = 0;
	norm_quat = sqrt(dq(0,0)*dq(0,0) + dq(1,0)*dq(1,0) + dq(2,0)*dq(2,0) + dq(3,0)*dq(3,0));	
	if (norm_quat !=1)
	{
		dq = dq/norm_quat;
	}
	// if(dq(0,0) > 1.0)
	// {
	// 	dq << 1, 0, 0, 0, dq(4,0), dq(5,0), dq(6,0), dq(7,0);
	// }	
	sr=dq(0);

	vr << dq(1), dq(2), dq(3);
	sd=dq(4);
	vd << dq(5), dq(6), dq(7);
	theta_e=2*acos(sr); // If no errors occur, the arc cosine of arg (arccos(arg)) in the range [0 , π], is returned.
	// theta_e= DQ::normalizeAngleZeroTo2Pi(theta_e);
	double absTheta=fabs(theta_e);
	if (absTheta > eps_theta && absTheta < (2*M_PI -eps_theta))
	{
		l_e=vr/vr.norm();
		//d_e=-sd*(2/(vr.norm()));
		// d_e = -2*sd/sin(theta_e/2);
		d_e = vec_translation(0)*l_e(0) + vec_translation(1)*l_e(1) + vec_translation(2)*l_e(2) ;
		//m_e=(vd-sr*0.5*d_e*l_e)/vr.norm();
		m_e=(vd - cos(theta_e/2)*(d_e/2)*l_e)/sin(theta_e/2);
	}	
	else
	{
		// theta_e = 0;
		RowVector4d qrr, qdd, tt;
		RowVector3d t;
		qrr << dq(0), dq(1), dq(2), dq(3);	
		qdd << dq(4), dq(5), dq(6), dq(7);
		tt=2*multiplyQuat(qdd, conjugateQuat(qrr));
		t << tt(1), tt(2), tt(3);  
		d_e=t.norm();
		if (d_e == 0)
		{
			l_e << 0,0,0;
		}
		else l_e=t/d_e;
		m_e << 0,0,0;			
	}
}


Matrix4d DQ::dqToHTM(Matrix<double,8,1> dq)
{
    if(dq(0) < 0)
    {
	    dq = -dq;
	}	
    RowVector4d qrr, qtt;
	RowVector3d u;

    qrr=RowVector4d(dq(0),dq(1),dq(2),dq(3));
    qtt=RowVector4d(dq(4),dq(5),dq(6),dq(7));
    qtt=2*multiplyQuat(qtt, conjugateQuat(qrr));
    if(qrr(0) > 1)
    	qrr(0) =1;


    double theta=2*acos(qrr(0));
    if (theta < -0.00001 || theta > 0.00001)
    	u=RowVector3d(qrr(1),qrr(2),qrr(3))/sin(theta/2);
    else
    	u=RowVector3d(0, 0, 1);

    Matrix3d skw, rot;
    skw<< 0, -u(2), u(1),
    		u(2), 0., -u(0),
    		-u(1), u(0), 0;

    rot=MatrixXd::Identity(3,3)+sin(theta)*skw+skw*skw*(1-cos(theta));
    Matrix4d htm;
    htm << rot(0,0), rot(0,1), rot(0,2), qtt(1),
    		rot(1,0), rot(1,1), rot(1,2), qtt(2),
    		rot(2,0), rot(2,1), rot(2,2), qtt(3),
    		0, 0, 0, 1;	
	return htm;
}

void DQ::dqToRotationAndTranslationAsQuat(Matrix<double,8,1> dq, RowVector4d &rot, RowVector4d &trans)
{
	// translation then rotation.
	rot << dq(0), dq(1), dq(2), dq(3);
	trans << dq(4), dq(5), dq(6), dq(7);
	trans= 2*multiplyQuat(conjugateQuat(rot), trans);
}

Matrix<double,8,1>  DQ::htmToDQ(Matrix4d htm_original)
{
	Matrix<double,8,1> dq;
	Matrix3d htm_, htm_eye;
	htm_=htm_original.block<3,3>(0,0);
	htm_eye =htm_-MatrixXd::Identity(3,3);
    EigenSolver<MatrixXd> eigensolver(htm_eye);

    std::vector<double> eigenVal, eigenVec;
    eigenVal.clear();
    eigenVal.resize(htm_.rows());
    eigenVec.clear();

    int index=0;
    for(int i=0;i<htm_.rows(); i++)
    {
      eigenVal[i]=abs(eigensolver.eigenvalues()[i]);
      if(eigenVal[i]<eigenVal[index])
      	index=i;
    }

    if (eigenVal[index]>0.001)
      std::cerr << "Rotation Matrix seems dubious\n";
	RowVector3d vecAxis=eigensolver.eigenvectors().col(index).real();
	if (abs(vecAxis.norm()-1)>0.0001)
		std::cerr << "Non-unit rotation axis"<< std::endl;

	double twoCosTheta=htm_.trace()-1;

	RowVector3d twoSinThetaV= RowVector3d ((htm_(2,1)-htm_(1,2)), (htm_(0,2)-htm_(2,0)), (htm_(1,0)-htm_(0,1)));
	double twoSinTheta=vecAxis*twoSinThetaV.transpose();

	double theta= std::atan2(twoSinTheta,twoCosTheta);

	RowVector4d rot_q=RowVector4d(cos(theta/2), sin(theta/2)*vecAxis(0), sin(theta/2)*vecAxis(1), sin(theta/2)*vecAxis(2));

	RowVector4d trans_q=RowVector4d(0., htm_original(0,3), htm_original(1,3), htm_original(2,3));

	RowVector4d prodRotTrans=0.5*multiplyQuat(trans_q, rot_q);

	dq<< rot_q(0),
		rot_q(1),
		rot_q(2),
		rot_q(3),
		prodRotTrans(0),
		prodRotTrans(1),
		prodRotTrans(2),
		prodRotTrans(3);

	return dq;
}



void DQ::dqEigenToDQdoubleVector(RowVectorXd dq_eigen, std::vector<double> &dq_double)
{
	dq_double.clear();

	for (int i=0; i< dq_eigen.size(); i++)
	{
		dq_double.push_back(dq_eigen(i));
	}
}

void DQ::dqDoubleVectorToDQEigen(Matrix<double,8,1> &dq_eigen, std::vector<double> dq_double)
{

	for (int i=0; i< dq_double.size(); i++)
	{
		dq_eigen(i)=dq_double[i];
	}
}

std::vector<double> DQ::returnDQdoubleVectorFromDQRowVector(RowVectorXd dq_eigen)
{
	std::vector<double> dq_double;
	dq_double.clear();
	for (int i=0; i< dq_eigen.size(); i++)
	{
		dq_double.push_back(dq_eigen(i));
	}
	return dq_double;
}

std::vector<double>  DQ::returnDQdoubleVectorFromDQMatrix(Matrix<double,8,1> dq_eigen)
{
	std::vector<double> dq_double;
	dq_double.clear();
	for (int i=0; i< 8; i++)
	{
		dq_double.push_back(dq_eigen(i,0));
	}
	return dq_double;
}


Matrix<double,8,1> DQ::inverseDQtransformation(Matrix<double,8,1> dq)
{
	RowVector4d rot; 
	RowVector4d trans;
	dqToRotationAndTranslationAsQuat(dq, rot, trans);
	dq << conjugateQuat(rot), 0.5*multiplyQuat(conjugateQuat(rot), -trans);
	return dq;
}

Matrix<double,8,1>  DQ::transformPluckerLineEigen(Matrix<double,8,1> line, Matrix<double,8,1> transform)
{		
	line=multiplyDQ(transform, multiplyDQ(line, classicConjugateDQ(transform)));
	return line;
}

RowVectorXd  DQ::transformLine3DRowVector(RowVectorXd lineVector, Matrix<double,8,1> transform)
{
	Matrix<double,8,1> line_dq;
	line_dq << 0, lineVector(0), lineVector(1), lineVector(2), 0, 0, 0, 0;
	line_dq = multiplyDQ(transform, multiplyDQ(line_dq, classicConjugateDQ(transform)));
	lineVector.resize(3);
	lineVector << line_dq(1), line_dq(2), line_dq(3);
	return lineVector;
}

RowVectorXd  DQ::transformPluckerLine6DRowVector(RowVectorXd lineVector, Matrix<double,8,1> transform)
{
	Matrix<double,8,1> line_dq;
	line_dq << 0, lineVector(0), lineVector(1), lineVector(2), 0, lineVector(3), lineVector(4), lineVector(5);
	line_dq = multiplyDQ(transform, multiplyDQ(line_dq, classicConjugateDQ(transform)));
	lineVector.resize(6);
	lineVector << line_dq(1), line_dq(2), line_dq(3), line_dq(5), line_dq(6), line_dq(7);
	return lineVector;
}

RowVector3d  DQ::transformPoint(RowVector3d point, Matrix<double,8,1> transform)
{
	Matrix<double,8,1> point_dq;
	point_dq << 1, 0, 0, 0, 0, point(0), point(1), point(2);	
	point_dq=multiplyDQ(multiplyDQ(transform, point_dq), combinedConjugateDQ(transform));
	point.resize(3);
	point << point_dq(5), point_dq(6), point_dq(7);
	return point;
}

double DQ::normalizeAngleZeroTo2Pi(double theta)
{
	if (theta < 0)
		while (theta<0)
			theta= theta+ 2*M_PI;
	else if(theta >=2*M_PI)
		while (theta>=2*M_PI)
			theta=theta-2*M_PI;
	return theta;
}

std::vector<double>  DQ::returnDoubleVectorFromRowVector(RowVectorXd eigenVector)
{
	std::vector<double> dq_double;
	dq_double.clear();
	for (int i=0; i< eigenVector.size(); i++)
	{
		dq_double.push_back(eigenVector(i));
	}
	return dq_double;
}


RowVectorXd DQ::returnDQRowVectorFromDQMatrix(Matrix<double,8,1> matrix8D)
{
	RowVectorXd rowVector(8);
	rowVector << matrix8D(0), matrix8D(1), matrix8D(2), matrix8D(3), matrix8D(4), matrix8D(5), matrix8D(6), matrix8D(7);
	return rowVector;
}

RowVectorXd DQ::returnRowVector6DFromMatrix8D(Matrix<double,8,1> matrix8D)
{
	RowVectorXd rowVector(6);
	rowVector << matrix8D(1), matrix8D(2), matrix8D(3), matrix8D(5), matrix8D(6), matrix8D(7);
	return rowVector;
}

Matrix<double,8,1> DQ::dqTotwist(Matrix<double,8,1> dq) 
{
	Matrix<double,8,1> dq_twist;
	double theta_e, d_e; 
	RowVector3d l_e, m_e, v, w;
	dqToScrewParameters(dq, theta_e, d_e, l_e, m_e) ;
	w = l_e*theta_e;
	v =  l_e*d_e + m_e*theta_e;
	dqToScrewParameters(dq, theta_e, d_e, l_e, m_e);
	dq_twist << 0, w.transpose(), 0 , v.transpose(); 
	return dq_twist;
}
// Creating a 6- dimensional twist to be used in modified cartesian impedance controller

Matrix<double,6,1> DQ::dqTotwist_6D(Matrix<double,8,1> dq) 
{
	Matrix<double,6,1> dq_twist_6D;  // linear velocity and then angular velocity
	double theta_e, d_e; 
	RowVector3d l_e, m_e, v, w;
	dqToScrewParameters(dq, theta_e, d_e, l_e, m_e) ;
	w = l_e*theta_e;
	v =  l_e*d_e + m_e*theta_e;
	dqToScrewParameters(dq, theta_e, d_e, l_e, m_e);
	dq_twist_6D <<  v.transpose(), w.transpose(); 
	return dq_twist_6D;
}

Matrix<double,8,1> DQ::twistTodq(Matrix<double,8,1> dq_twist)
{
	Matrix<double,8,1> dq;
	double theta_e, d_e;
	RowVector3d l_e, m_e, w, v;

	w << dq_twist(1,0), dq_twist(2,0), dq_twist(3,0);
	v << dq_twist(5,0), dq_twist(6,0), dq_twist(7,0);
	theta_e = w.norm();
	if (theta_e != 0)
	{
		l_e = w/theta_e;
		d_e = (l_e.transpose()*v)(0,0);
		m_e = (v - d_e*l_e)/theta_e;
		dq = screwDisplacementToDQ(theta_e, l_e, d_e, m_e);
	} 
	else
	{
		dq << 1, 0,  0,  0,  0,  0,  0,  0;
	}
	return dq;
}


Matrix<double,8,1> DQ::rotationTranslationToDQ(RowVector4d rot, RowVector4d trans)
{
	Matrix<double,8,1> dq;
	RowVector4d trans_quat = multiplyQuat(trans, rot)/2;
	dq << rot(0), rot(1), rot(2), rot(3), trans_quat(0), trans_quat(1), trans_quat(2), trans_quat(3);
	return dq;
}



RowVectorXd DQ::rowVector6DToRowVector8D(RowVectorXd twist)
{
	RowVectorXd rowVector(8);
	rowVector << 0, twist(0), twist(1), twist(2), 0, twist(3), twist(4), twist(5);
	return rowVector;
}

RowVectorXd DQ::rowVector8DToRowVector6D(RowVectorXd DQtwist)
{
	RowVectorXd rowVector(6);
	rowVector << DQtwist(1), DQtwist(2), DQtwist(3), DQtwist(5), DQtwist(6), DQtwist(7);
	return rowVector;
}

Matrix<double,8,1> DQ::rowVector6DToMatrix8D(RowVectorXd rowVector6d)
{
	Matrix<double,8,1> matrix8d;
	matrix8d << 0, rowVector6d(0), rowVector6d(1), rowVector6d(2), 0, rowVector6d(3), rowVector6d(4), rowVector6d(5);
	return matrix8d;
}


Matrix<double,8,1> DQ::returnMatrixDQFromDoubleVectorDQ(std::vector<double> dq_double)
{
	Matrix<double,8,1> dq_eigen;
	for (int i=0; i< dq_double.size(); i++)
	{
		dq_eigen(i)=dq_double[i];
	}
	return dq_eigen;
}

RowVector3d DQ::getPositionFromDQ(Matrix<double,8,1>  pose)
{
	Matrix4d htm_desired;
	htm_desired = dqToHTM(pose);
	RowVector3d ee_position;
	ee_position << htm_desired(0,3), htm_desired(1,3), htm_desired(2,3); 
	return ee_position;
}