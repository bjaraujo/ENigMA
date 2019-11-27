include(ExternalProject)

ExternalProject_Add(
    eigen
    URL http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz
    URL_HASH SHA1=743c1dc00c6680229d8cc87d44debe5a71d15c01
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    UPDATE_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    LOG_CONFIGURE OFF
    LOG_BUILD OFF
)

ExternalProject_Get_Property(eigen SOURCE_DIR)
set(EIGEN_DIR ${SOURCE_DIR})
