include(ExternalProject)

ExternalProject_Add(
    eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/5.0.0/eigen-5.0.0.tar.gz
    URL_HASH SHA1=67ff2328ab2c771baf14b0dc1254f7272f712c98
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
