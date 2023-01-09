cmake_minimum_required(VERSION 3.10)

project(singRC)

#set(CMAKE_C_FLAGS "-ggdb" )
#set(CMAKE_CXX_FLAGS "-ggdb" )

include_directories(include external/matplotlib)
include_directories(external/mpcc/spline)
include_directories(
  external/acado_generated
  external/qpoases/
  external/qpoases/INCLUDE/
  external/qpoases/SRC/
  external/gurobi
  qp_solver_interface)
include_directories(external/Polyfit)

file(GLOB SOURCES "src/*.cpp"
  "external/mpcc/spline/*.cpp")
file(GLOB ACADO_EXPORTED_SRC
  "external/acado_generated/*.cpp"
  "external/acado_generated/*.c"
  "external/qpoases/SRC/*.cpp"
  "qp_solver_interface/*.cpp")

add_library(acado_exported_rti STATIC ${ACADO_EXPORTED_SRC})

add_executable(main main.cpp ${SOURCES})
target_link_libraries(main acado_exported_rti)
target_link_libraries(main /opt/gurobi910/linux64/lib/libgurobi91.so)

find_package(PythonLibs 2.7)
target_include_directories(main PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(main ${PYTHON_LIBRARIES})
