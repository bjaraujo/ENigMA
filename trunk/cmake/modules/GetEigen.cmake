include(ExternalProject)

ExternalProject_Add(
    eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz
    URL_HASH SHA1=3e8ab94ff389ae4b3103f5dda77c826d85b9c1d5
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
