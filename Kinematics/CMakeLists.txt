add_library(kinematics SHARED particles.cpp reactions.cpp)
target_link_libraries(kinematics ${ROOT_LIBRARIES})