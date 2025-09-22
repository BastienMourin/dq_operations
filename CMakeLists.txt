cmake_minimum_required(VERSION 3.5)
project(dq_operations)

# Use C++11 and suppress deprecated warnings
add_compile_options(-std=c++11 -Wno-deprecated-declarations)

# Find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(orocos_kdl REQUIRED)

# Include directories
include_directories(
  include
  ${EIGEN3_INCLUDE_DIR}
)

# Libraries
add_library(dq SHARED src/dq.cpp)
ament_target_dependencies(dq rclcpp Eigen3 orocos_kdl)

add_library(dq_manipulator_operations SHARED src/dq_manipulator_operations.cpp)
ament_target_dependencies(dq_manipulator_operations rclcpp Eigen3 orocos_kdl)
target_link_libraries(dq_manipulator_operations dq)

add_library(dq_dual_arm_operations SHARED src/dq_dual_arm_operations.cpp)
ament_target_dependencies(dq_dual_arm_operations rclcpp Eigen3 orocos_kdl)
target_link_libraries(dq_dual_arm_operations dq)

add_library(jacDotSolver SHARED src/kdlJacDotSolver.cpp)
ament_target_dependencies(jacDotSolver rclcpp Eigen3 orocos_kdl)
target_link_libraries(jacDotSolver dq)

# Install targets and headers
install(TARGETS dq dq_manipulator_operations dq_dual_arm_operations jacDotSolver
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

install(DIRECTORY include/
  DESTINATION include/
)

# Export package
ament_export_include_directories(include)
ament_export_libraries(dq dq_manipulator_operations dq_dual_arm_operations jacDotSolver)
ament_export_dependencies(rclcpp Eigen3 orocos_kdl)

ament_package()