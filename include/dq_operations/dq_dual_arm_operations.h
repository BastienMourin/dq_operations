#ifndef DQDualArm_H
#define DQDualArm_H

#include <dq_operations/dq.h>

class DQDualArm
{

	public:
		/**
		 * @brief      Gets the absolute pose in reference frame.
		 *
		 * @param[in]  pose_rel  The pose relative
		 *
		 * @return     The absolute pose in reference frame.
		 */
		static Matrix<double,8,1> getAbsolutePoseInRefFrame(Matrix<double,8,1> pose_rel);	

		/**
		 * @brief      Gets the right(ref) left(tool) from absolute relative.
		 *
		 * @param[in]  absPose  The absolute pose
		 * @param[in]  relPose  The relative pose
		 *
		 * @return     The [right left] vector of UDQ pose.
		 */
		static std::vector<Matrix<double,8,1> > getRightLeft4mAbsRel(Matrix<double,8,1> absPose, Matrix<double,8,1> relPose);
			


		/**
		 * @brief      Gets the absolute pose from right(ref) and left(tool) arms pose.
		 *
		 * @param[in]  rightPoseNow  The right pose UDQ
		 * @param[in]  leftPoseNow   The left pose UDQ
		 *
		 * @return     The absolute pose UDQ.
		 */		
		static Matrix<double,8,1> getAbsolutePose(Matrix<double,8,1> rightPoseNow, Matrix<double,8,1> leftPoseNow);

		DQDualArm();
		~DQDualArm();
};
#endif