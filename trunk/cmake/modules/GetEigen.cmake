include(ExternalProject)

ExternalProject_Add(
  eigen
  URL http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  LOG_DOWNLOAD ON
  LOG_CONFIGURE OFF
  LOG_BUILD OFF)

ExternalProject_Get_Property(eigen SOURCE_DIR)
set(EIGEN_DIR ${SOURCE_DIR})


