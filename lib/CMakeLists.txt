cmake_minimum_required(VERSION 3.0)

# Build the matrix library
include_directories(matrix)

# Build the split library
include_directories(split)
add_library(split split/split.cpp)

# Build the config library
include_directories(config)
add_library(config config/Config.cpp)

# Build the threadpool library
include_directories(threadpool)
add_library(threadpool threadpool/ThreadPool.cpp)

# Build the mydefine library
include_directories(mydefine)
add_library(mydefine mydefine/MyDefine.cpp)
target_link_libraries(mydefine config inputparser paramodel treecaller threadpool)

# Build the inputparser library
include_directories(inputparser)
add_library(inputparser inputparser/InputParser.cpp)
target_link_libraries(inputparser mydefine)

# Build the treecaller library
include_directories(treecaller)
add_library(treecaller treecaller/TreeCaller.cpp)
target_link_libraries(treecaller mydefine split)

# Build the bayesopt library
FIND_LIBRARY(NLOPT nlopt)
IF(NLOPT MATCHES NLOPT-NOTFOUND)
  SET(NLOPT_BUILD ON CACHE BOOL "Build included version of NLOPT?")
ELSE(NLOPT MATCHES NLOPT-NOTFOUND)
  SET(NLOPT_BUILD OFF CACHE BOOL "Build included version of NLOPT?")
ENDIF(NLOPT MATCHES NLOPT-NOTFOUND)

IF(NLOPT_BUILD)
  ADD_SUBDIRECTORY(bayesopt/nlopt2)
  include_directories(bayesopt/nlopt2/api)
  SET(EXT_LIBS nlopt)
ELSE(NLOPT_BUILD)
  SET(EXT_LIBS ${NLOPT})
ENDIF(NLOPT_BUILD)

SET(BAYESOPT_SRCS
  ./bayesopt/src/bayesoptcont.cpp
  ./bayesopt/src/bayesoptdisc.cpp
  ./bayesopt/src/bayesoptbase.cpp
  ./bayesopt/src/bopt_state.cpp
  ./bayesopt/src/posteriormodel.cpp
  ./bayesopt/src/posterior_fixed.cpp
  ./bayesopt/src/posterior_empirical.cpp
  ./bayesopt/src/posterior_mcmc.cpp
  ./bayesopt/src/mcmc_sampler.cpp
  ./bayesopt/src/inneroptimization.cpp
  ./bayesopt/src/dataset.cpp
  ./bayesopt/src/nonparametricprocess.cpp
  ./bayesopt/src/kernelregressor.cpp
  ./bayesopt/src/conditionalbayesprocess.cpp
  ./bayesopt/src/gaussian_process.cpp
  ./bayesopt/src/gaussian_process_ml.cpp
  ./bayesopt/src/gaussian_process_normal.cpp
  ./bayesopt/src/gaussian_process_hierarchical.cpp
  ./bayesopt/src/student_t_process_jef.cpp
  ./bayesopt/src/student_t_process_nig.cpp
  ./bayesopt/src/parameters.cpp
  ./bayesopt/src/kernel_functors.cpp
  ./bayesopt/src/criteria_functors.cpp
  ./bayesopt/src/criteria_hedge.cpp
  ./bayesopt/src/mean_functors.cpp
  ./bayesopt/src/gauss_distribution.cpp
  ./bayesopt/src/student_t_distribution.cpp
  )

SET(UTILS_SRC
  ./bayesopt/utils/fileparser.cpp
  ./bayesopt/utils/param_loader.cpp
  ./bayesopt/utils/parser.cpp
  ./bayesopt/utils/ublas_extra.cpp
  )

SET(WRAPPPERS_SRC 
  ./bayesopt/src/wrappers/bayesoptwpr.cpp
  )
#FILE(GLOB bayesopt_src bayesopt/src/*.cpp)
#FILE(GLOB bayesopt_util_src bayesopt/utils/*.cpp)
#FILE(GLOB bayesopt_wrappers_src bayesopt/src/wrappers/*.cpp)

# create main bayesopt API static library
include_directories(bayesopt/include bayesopt/utils)
add_library(bayesopt ${BAYESOPT_SRCS} ${UTILS_SRC} ${WRAPPPERS_SRC})
target_link_libraries(bayesopt ${EXT_LIBS})

IF(NLOPT_BUILD)
  add_dependencies(bayesopt nlopt)
ENDIF(NLOPT_BUILD)

# Build the paramodel library
include_directories(paramodel)
add_library(paramodel paramodel/ParaModel.cpp)
target_link_libraries(paramodel mydefine bayesopt)


