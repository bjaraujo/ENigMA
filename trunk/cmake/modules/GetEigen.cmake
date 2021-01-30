include(ExternalProject)

ExternalProject_Add(
    eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.gz
    URL_HASH SHA1=6a5a43a327b3aaeb7e74dc32bf1d7011cf6f149c
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
