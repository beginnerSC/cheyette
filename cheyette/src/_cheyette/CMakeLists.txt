cmake_minimum_required(VERSION 3.11)

project(_cheyette VERSION "0.1")

# Define CMAKE_INSTALL_xxx: LIBDIR, INCLUDEDIR
include(GNUInstallDirs)
include(FetchContent)

################ Fetching pybind11 from github does not work behind a corporate firewall ################
# FetchContent_Declare(
#   pybind11
#   URL https://github.com/pybind/pybind11/archive/refs/tags/v2.10.4.tar.gz
#   URL_HASH SHA256=832e2f309c57da9c1e6d4542dedd34b24e4192ecb4d62f6f4866a737454c9970
# )
# FetchContent_MakeAvailable(pybind11)

set(PYBIND11_FINDPYTHON ON)
set(pybind11_DIR "D:/Users/zhang/AppData/Local/pypoetry/Cache/virtualenvs/cheyette-5_iMRCa0-py3.12/Lib/site-packages/pybind11/share/cmake/pybind11")
find_package(pybind11 REQUIRED)

include_directories("D:/Users/zhang/AppData/Local/pypoetry/Cache/virtualenvs/cheyette-5_iMRCa0-py3.12/Lib/site-packages/pybind11/include")
include_directories("D:/python3.12.4/include")        # to see Python.h
include_directories("D:/3rdparty/eigen-3.4.0")        # to see Eigen

link_directories("D:/python3.12.4/libs")   # for linker to see python312.lib

add_executable(_cheyette_exe   # build executable to debug; one must define a main function as the entry point
  _cheyette_pybind.cpp 
  example1.cpp
  example2.cpp
)

set(python_module_name _cheyette)
pybind11_add_module(${python_module_name} MODULE
  _cheyette_pybind.cpp 
  example1.cpp
  example2.cpp
)

install(TARGETS ${python_module_name} DESTINATION .)
