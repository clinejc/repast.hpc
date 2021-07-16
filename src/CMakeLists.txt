cmake_minimum_required(VERSION 3.0)
project(repast_hpc_project
	VERSION 2.3.1
)

set (rhpc_lib_name repast_hpc)
set (relogo_lib_name relogo)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../cmake/Modules/")

set(rhpc_headers
	repast_hpc/AgentId.h
	repast_hpc/AgentImporterExporter.h
	repast_hpc/AgentRequest.h
	repast_hpc/AgentStatus.h
	repast_hpc/BaseGrid.h
	repast_hpc/CartesianTopology.h
	repast_hpc/Context.h
	repast_hpc/DataSet.h
	repast_hpc/DiffusionLayerND.h
	repast_hpc/DirectedVertex.h
	repast_hpc/Edge.h
	repast_hpc/Graph.h
	repast_hpc/Grid.h
	repast_hpc/Grid2DQuery.h
	repast_hpc/GridComponents.h
	repast_hpc/GridDimensions.h
	repast_hpc/Moore2DGridQuery.h
	repast_hpc/MultipleOccupancy.h
	repast_hpc/NCDataSet.h
	repast_hpc/NCDataSetBuilder.h
	repast_hpc/NCDataSource.h
	repast_hpc/NCReducibleDataSource.h
	repast_hpc/NetworkBuilder.h
	repast_hpc/Point.h
	repast_hpc/Projection.h
	repast_hpc/Properties.h
	repast_hpc/Random.h
	repast_hpc/ReducibleDataSource.h
	repast_hpc/RelativeLocation.h
	repast_hpc/RelativeLocationMap.h
	repast_hpc/RepastErrors.h
	repast_hpc/RepastProcess.h
	repast_hpc/SRManager.h
	repast_hpc/SVDataSet.h
	repast_hpc/SVDataSetBuilder.h
	repast_hpc/SVDataSource.h
	repast_hpc/Schedule.h
	repast_hpc/SharedBaseGrid.h
	repast_hpc/SharedContext.h
	repast_hpc/SharedContinuousSpace.h
	repast_hpc/SharedDiscreteSpace.h
	repast_hpc/SharedNetwork.h
	repast_hpc/SharedSpaces.h
	repast_hpc/SingleOccupancy.h
	repast_hpc/Spaces.h
	repast_hpc/TDataSource.h
	repast_hpc/UndirectedVertex.h
	repast_hpc/Utilities.h
	repast_hpc/VN2DGridQuery.h
	repast_hpc/ValueLayer.h
	repast_hpc/ValueLayerND.h
	repast_hpc/Variable.h
	repast_hpc/Vertex.h
	repast_hpc/WeightedRandomSelector.h
	repast_hpc/initialize_random.h
	repast_hpc/io.h
	repast_hpc/logger.h
	repast_hpc/matrix.h
	repast_hpc/mpi_constants.h
	repast_hpc/spatial_math.h
	)

set(rhpc_src
	repast_hpc/AgentId.cpp
	repast_hpc/AgentId.h
	repast_hpc/AgentImporterExporter.cpp
	repast_hpc/AgentImporterExporter.h
	repast_hpc/AgentRequest.cpp
	repast_hpc/AgentRequest.h
	repast_hpc/AgentStatus.cpp
	repast_hpc/AgentStatus.h
	repast_hpc/BaseGrid.h
	repast_hpc/CartesianTopology.cpp
	repast_hpc/CartesianTopology.h
	repast_hpc/Context.h
	repast_hpc/DataSet.h
	repast_hpc/DirectedVertex.h
	repast_hpc/Edge.h
	repast_hpc/Graph.cpp
	repast_hpc/Graph.h
	repast_hpc/Grid.h
	repast_hpc/Grid2DQuery.h
	repast_hpc/GridComponents.cpp
	repast_hpc/GridComponents.h
	repast_hpc/GridDimensions.cpp
	repast_hpc/GridDimensions.h
	repast_hpc/initialize_random.cpp
	repast_hpc/initialize_random.h
	repast_hpc/io.cpp
	repast_hpc/io.h
	repast_hpc/logger.cpp
	repast_hpc/logger.h
	repast_hpc/matrix.h
	repast_hpc/Moore2DGridQuery.h
	repast_hpc/mpi_constants.h
	repast_hpc/MultipleOccupancy.h
	repast_hpc/NCDataSet.cpp
	repast_hpc/NCDataSet.h
	repast_hpc/NCDataSetBuilder.cpp
	repast_hpc/NCDataSetBuilder.h
	repast_hpc/NCDataSource.cpp
	repast_hpc/NCDataSource.h
	repast_hpc/NCReducibleDataSource.h
	repast_hpc/NetworkBuilder.cpp
	repast_hpc/NetworkBuilder.h
	repast_hpc/Point.h
	repast_hpc/Projection.h
	repast_hpc/Properties.cpp
	repast_hpc/Properties.h
	repast_hpc/Random.cpp
	repast_hpc/Random.h
	repast_hpc/ReducibleDataSource.h
	repast_hpc/RelativeLocation.cpp
	repast_hpc/RelativeLocation.h
	repast_hpc/RelativeLocationMap.h
	repast_hpc/RepastErrors.cpp
	repast_hpc/RepastErrors.h
	repast_hpc/RepastProcess.cpp
	repast_hpc/RepastProcess.h
	repast_hpc/Schedule.cpp
	repast_hpc/Schedule.h
	repast_hpc/SharedBaseGrid.cpp
	repast_hpc/SharedBaseGrid.h
	repast_hpc/SharedContext.cpp
	repast_hpc/SharedContext.h
	repast_hpc/SharedContinuousSpace.h
	repast_hpc/SharedDiscreteSpace.h
	repast_hpc/SharedNetwork.h
	repast_hpc/SharedSpaces.h
	repast_hpc/SingleOccupancy.h
	repast_hpc/Spaces.h
	repast_hpc/spatial_math.cpp
	repast_hpc/spatial_math.h
	repast_hpc/SRManager.cpp
	repast_hpc/SRManager.h
	repast_hpc/SVDataSet.cpp
	repast_hpc/SVDataSet.h
	repast_hpc/SVDataSetBuilder.cpp
	repast_hpc/SVDataSetBuilder.h
	repast_hpc/SVDataSource.h
	repast_hpc/TDataSource.h
	repast_hpc/UndirectedVertex.h
	repast_hpc/Utilities.cpp
	repast_hpc/Utilities.h
	repast_hpc/ValueLayer.cpp
	repast_hpc/ValueLayer.h
	repast_hpc/ValueLayerND.cpp
	repast_hpc/ValueLayerND.h
	repast_hpc/Variable.cpp
	repast_hpc/Variable.h
	repast_hpc/Vertex.h
	repast_hpc/VN2DGridQuery.h
	repast_hpc/WeightedRandomSelector.h
)

set(relogo_headers
	relogo/AbstractRelogoAgent.h
	relogo/AgentSet.h
	relogo/Observer.h
	relogo/Patch.h
	relogo/RandomMove.h
	relogo/RelogoAgent.h
	relogo/RelogoContinuousSpaceAdder.h
	relogo/RelogoDiscreteSpaceAdder.h
	relogo/RelogoErrors.h
	relogo/RelogoLink.h
	relogo/RelogoSharedContinuousSpace.h
	relogo/RelogoSharedDiscreteSpace.h
	relogo/SimulationRunner.h
	relogo/SimulationRunnerPlus.h
	relogo/Turtle.h
	relogo/WorldCreator.h
	relogo/WorldDefinition.h
	relogo/agent_set_functions.h
	relogo/creators.h
	relogo/grid_types.h
	relogo/relogo.h
	relogo/utility.h
	)

set (relogo_src
	relogo/AbstractRelogoAgent.cpp
	relogo/AbstractRelogoAgent.h
	relogo/agent_set_functions.h
	relogo/AgentSet.cpp
	relogo/AgentSet.h
	relogo/creators.cpp
	relogo/creators.h
	relogo/grid_types.h
	relogo/Observer.cpp
	relogo/Observer.h
	relogo/Patch.cpp
	relogo/Patch.h
	relogo/RandomMove.cpp
	relogo/RandomMove.h
	relogo/relogo.h
	relogo/RelogoAgent.cpp
	relogo/RelogoAgent.h
	relogo/RelogoContinuousSpaceAdder.cpp
	relogo/RelogoContinuousSpaceAdder.h
	relogo/RelogoDiscreteSpaceAdder.cpp
	relogo/RelogoDiscreteSpaceAdder.h
	relogo/RelogoErrors.h
	relogo/RelogoLink.cpp
	relogo/RelogoLink.h
	relogo/RelogoSharedContinuousSpace.h
	relogo/RelogoSharedDiscreteSpace.h
	relogo/SimulationRunner.h
	relogo/SimulationRunnerPlus.h
	relogo/Turtle.cpp
	relogo/Turtle.h
	relogo/utility.cpp
	relogo/utility.h
	relogo/WorldCreator.cpp
	relogo/WorldCreator.h
	relogo/WorldDefinition.cpp
	relogo/WorldDefinition.h
)

set (rumor_src
	rumor_model/main.cpp
	rumor_model/RumorModel.cpp
	rumor_model/RumorModel.h
)

set (zombie_src
	zombies/AgentPackage.h
	zombies/Human.cpp
	zombies/Human.h
	zombies/InfectionSum.cpp
	zombies/InfectionSum.h
	zombies/main.cpp
	zombies/Zombie.cpp
	zombies/Zombie.h
	zombies/ZombieObserver.cpp
	zombies/ZombieObserver.h
)

set (core_ut_src
        ../test/core/context_test.cpp
        ../test/core/error_test.cpp
        ../test/core/grid_comp_test.cpp
        ../test/core/grid_test.cpp
        ../test/core/main.cpp
        ../test/core/properties_test.cpp
        ../test/core/random_test.cpp
        ../test/core/schedule_test.cpp
        ../test/core/test.h
        ../test/core/value_layer_tests.cpp
)

set (relogo_ut_src
        ../test/relogo/agent_set_tests.cpp
        ../test/relogo/main.cpp
        ../test/relogo/Objects.h
        ../test/relogo/ObserverTests.cpp
)

message("If Boost cannot be found, point CMake to your Boost installation with:")
message("-DBOOST_ROOT=<path to Boost installation>")
set(Boost_USE_MULTITHREADED ON)
find_package(Boost 
	1.61.0 
	REQUIRED 
	COMPONENTS 
		system 
		filesystem 
		mpi 
		serialization
)

if(NOT TARGET Boost::system)
    add_library(Boost::system IMPORTED INTERFACE)
    set_property(TARGET Boost::system PROPERTY
        INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
    set_property(TARGET Boost::system PROPERTY
        INTERFACE_LINK_LIBRARIES ${Boost_LIBRARIES})
endif()

if(NOT TARGET Boost::filesystem)
    add_library(Boost::filesystem IMPORTED INTERFACE)
    set_property(TARGET Boost::filesystem PROPERTY
        INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
    set_property(TARGET Boost::filesystem PROPERTY
        INTERFACE_LINK_LIBRARIES ${Boost_LIBRARIES})
endif()

if(NOT TARGET Boost::mpi)
    add_library(Boost::mpi IMPORTED INTERFACE)
    set_property(TARGET Boost::mpi PROPERTY
        INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
    set_property(TARGET Boost::mpi PROPERTY
        INTERFACE_LINK_LIBRARIES ${Boost_LIBRARIES})
endif()

if(NOT TARGET Boost::serialization)
    add_library(Boost::serialization IMPORTED INTERFACE)
    set_property(TARGET Boost::serialization PROPERTY
        INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
    set_property(TARGET Boost::serialization PROPERTY
        INTERFACE_LINK_LIBRARIES ${Boost_LIBRARIES})
endif()

# Find netCDF
include (FindNetCDF)

set (NETCDF_CXX "YES")
find_package(NetCDF REQUIRED)
# include_directories(${NETCDF_INCLUDE_DIR})

find_package(CURL REQUIRED)

find_package(MPI REQUIRED)

# For supporting CMake < 3.9:

if(NOT TARGET MPI::MPI_CXX)
    add_library(MPI::MPI_CXX IMPORTED INTERFACE)

    set_property(TARGET MPI::MPI_CXX
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${MPI_CXX_COMPILE_FLAGS})
    set_property(TARGET MPI::MPI_CXX
                 PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${MPI_CXX_INCLUDE_PATH}")
    set_property(TARGET MPI::MPI_CXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${MPI_CXX_LINK_FLAGS} ${MPI_CXX_LIBRARIES})
endif()

message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")


set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS} -std=c++11")
include_directories(${MPI_CXX_INCLUDE_PATH})

add_library(${rhpc_lib_name} SHARED ${rhpc_src})
set_target_properties(${rhpc_lib_name} PROPERTIES 
	VERSION ${PROJECT_VERSION}
	PUBLIC_HEADER "${rhpc_headers}"
)
target_link_libraries(${rhpc_lib_name} 
	Boost::system 
	Boost::filesystem 
	Boost::mpi 
	Boost::serialization 
	${NETCDF_LIBRARIES} 
	${CURL_LIBRARIES} 
	MPI::MPI_CXX
)
target_include_directories(${rhpc_lib_name}
	PRIVATE
		${NETCDF_INCLUDE_DIRS}
)


add_library(${relogo_lib_name} SHARED ${relogo_src})
target_include_directories(${relogo_lib_name} PUBLIC .)
set_target_properties(${relogo_lib_name} PROPERTIES 
	OUTPUT_NAME "${relogo_lib_name}-${PROJECT_VERSION}"
	PUBLIC_HEADER "${relogo_headers}"
)
add_dependencies(${relogo_lib_name} ${rhpc_lib_name})
target_link_libraries(${relogo_lib_name} 
	Boost::system 
	Boost::filesystem 
	Boost::mpi 
	Boost::serialization 
	${NETCDF_LIBRARIES} 
	${CURL_LIBRARIES} 
	MPI::MPI_CXX 
	${rhpc_lib_name}
)

set (rumor_exec rumor_model)
add_executable(${rumor_exec} ${rumor_src})
set_target_properties(${rumor_exec} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./bin/rumor)
file (COPY ../src/rumor_model/config.props ../src/rumor_model/model.props DESTINATION ./bin/rumor) 

add_dependencies(${rumor_exec} ${rhpc_lib_name})
target_link_libraries(${rumor_exec} 
	${rhpc_lib_name}
	Boost::system 
	Boost::filesystem 
	Boost::mpi 
	Boost::serialization 
	${NETCDF_LIBRARIES} 
	${CURL_LIBRARIES} 
	MPI::MPI_CXX 
)
target_include_directories(${rumor_exec}
	PUBLIC 
		.
		${NETCDF_INCLUDE_DIRS}
)

set (zombie_exec zombie_model)
add_executable(${zombie_exec} ${zombie_src})
set_target_properties( ${zombie_exec} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./bin/zombies)
file (COPY ../src/zombies/config.props ../src/zombies/model.props ../src/zombies/model-1cpu.props DESTINATION ./bin/zombies) 
target_include_directories(${zombie_exec} PUBLIC .)
add_dependencies(${zombie_exec} ${rhpc_lib_name} ${relogo_lib_name})
target_link_libraries(${zombie_exec} 
	Boost::system 
	Boost::filesystem 
	Boost::mpi 
	Boost::serialization 
	${NETCDF_LIBRARIES} 
	${CURL_LIBRARIES} 
	MPI::MPI_CXX 
	${rhpc_lib_name} 
	${relogo_lib_name}
)

# environment variable `GTEST_ROOT` should be set
find_package(GTest REQUIRED)

#set (core_test_exec core_unit_tests)
#add_executable(${core_test_exec} ${core_ut_src})
#target_include_directories(${core_test_exec} PUBLIC . ${GTEST_INCLUDE_DIRS})
#add_dependencies(${core_test_exec} ${rhpc_lib_name})
#file (COPY ../test_data/config.props ../test_data/test.properties DESTINATION .)
#target_link_libraries(${core_test_exec} Boost::system Boost::filesystem Boost::mpi Boost::serialization ${NETCDF_LIBRARIES} ${CURL_LIBRARIES} MPI::MPI_CXX ${rhpc_lib_name} ${GTEST_LIBRARIES})

# set (relogo_test_exec relogo_unit_tests)
# add_executable(${relogo_test_exec} ${relogo_ut_src})
# target_include_directories(${relogo_test_exec} PUBLIC . ${GTEST_INCLUDE_DIRS})
# add_dependencies(${relogo_test_exec} ${rhpc_lib_name} ${relogo_lib_name})
# file (COPY ../test_data/config.props ../test_data/test.properties DESTINATION .)
# target_link_libraries(${relogo_test_exec} 
# 	Boost::system 
# 	Boost::filesystem 
# 	Boost::mpi 
# 	Boost::serialization 
# 	${NETCDF_LIBRARIES} 
# 	${CURL_LIBRARIES} 
# 	MPI::MPI_CXX 
# 	${rhpc_lib_name} 
# 	${relogo_lib_name} 
# 	${GTEST_LIBRARIES}
# )


set(DEST_DIR "${CMAKE_INSTALL_PREFIX}")

set(repast_hpc_INCLUDE_DIR "${DEST_DIR}/repast_hpc")
set(relogo_INCLUDE_DIR "${DEST_DIR}/relogo")

INSTALL(TARGETS ${rhpc_lib_name}
	EXPORT RepastHPCTargets
	DESTINATION lib 
	PUBLIC_HEADER
		DESTINATION ${repast_hpc_INCLUDE_DIR}
)

INSTALL(TARGETS ${relogo_lib_name}
	EXPORT RelogoTargets
	DESTINATION lib 
	PUBLIC_HEADER
		DESTINATION ${relogo_INCLUDE_DIR}
)


set(exec_prefix ${CMAKE_INSTALL_PREFIX})
CONFIGURE_FILE("repast_hpc.pc.in" "repast_hpc.pc" @ONLY)

INSTALL(FILES "${CMAKE_BINARY_DIR}/repast_hpc.pc"
	        DESTINATION lib/pkgconfig)


# include(CMakePackageConfigHelpers)
# write_basic_package_version_file(
# 	"RepastHPCConfigVersion.cmake"
# 	COMPATABILITY AnyNewerVersion
# )

install(EXPORT RepastHPCTargets
	FILE RepastHPCTargets.cmake 
	NAMESPACE RepastHPC::
	DESTINATION lib/cmake/RepastHPC
)

configure_file(../cmake/RepastHPCConfig.cmake.in RepastHPCConfig.cmake @ONLY)
install(
	FILES 
		"${CMAKE_CURRENT_BINARY_DIR}/RepastHPCConfig.cmake"
		# "${CMAKE_CURRENT_BINARY_DIR}/RepastHPCConfigVersion.cmake"
	DESTINATION lib/cmake/RepastHPC
)