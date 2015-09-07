ADD_EXECUTABLE(optimize-particle-system optimize-particle-system.cpp)

#added by P.Petrov 2014 (required by c++11 std::thread on LINUX)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
TARGET_LINK_LIBRARIES(optimize-particle-system system pthread)

