cmake_minimum_required(VERSION 3.16)

# Project
project(PSToolbox LANGUAGES CXX)

# ---- Options ----
# You can toggle these if needed
option(PSTOOLBOX_WARNINGS "Enable compiler warnings" ON)
option(PSTOOLBOX_BUILD_TESTS "Build tests" OFF)

# ---- C++ standard & compile options ----
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Replicate optimization/debug flags roughly from the Makefile
# Users can still configure with -DCMAKE_BUILD_TYPE=Release/Debug/RelWithDebInfo/MinSizeRel
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

if(PSTOOLBOX_WARNINGS)
  if(MSVC)
    add_compile_options(/W4)
  else()
    add_compile_options(-Wall -Wextra)
    # The Makefile also had -w which disables warnings; we *do not* set that here by default.
  endif()
endif()

# ---- Dependencies: Eigen ----
# Prefer a proper Eigen3 package if available
find_package(Eigen3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
  message(STATUS "Eigen3 not found via find_package. You can set EIGEN3_INCLUDE_DIR to your Eigen headers.")
  # Fallback: let caller pass -DEIGEN3_INCLUDE_DIR=/path/to/eigen3
  if(DEFINED EIGEN3_INCLUDE_DIR)
    include_directories(${EIGEN3_INCLUDE_DIR})
  else()
    # As a last resort, try Homebrew default on macOS (similar to: -I$(brew --prefix eigen)/include/eigen3)
    if(APPLE)
      set(_brew_eigen_dir "/opt/homebrew/include/eigen3")
      if(EXISTS "${_brew_eigen_dir}")
        include_directories("${_brew_eigen_dir}")
      endif()
    endif()
  endif()
else()
  message(STATUS "Found Eigen3: ${Eigen3_DIR}")
  target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_NO_DEBUG=1)
endif()

# ---- Sources ----
# Mimic: SRCS := $(wildcard src/*.cpp)
file(GLOB PSTOOLBOX_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

# Exclude files like in Makefile:
# filter-out: src/CoolPropGas.cpp src/CoolPropHA.cpp src/PSToolboxPlotter.cpp
set(_exclude_list
  "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolPropGas.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolPropHA.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/PSToolboxPlotter.cpp"
)
foreach(_f IN LISTS _exclude_list)
  list(REMOVE_ITEM PSTOOLBOX_SOURCES "${_f}")
endforeach()

# ---- Library target ----
add_library(PSToolbox STATIC ${PSTOOLBOX_SOURCES})

# If Eigen3::Eigen is available, link it; otherwise just rely on include dirs added above
if(TARGET Eigen3::Eigen)
  target_link_libraries(PSToolbox PUBLIC Eigen3::Eigen)
endif()

# Public include directories (if you have headers under include/ add it here)
target_include_directories(PSToolbox
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Compile options roughly matching Makefile intent
if(NOT MSVC)
  target_compile_options(PSToolbox PRIVATE -O3)
endif()

# ---- Output directories ----
# Match the Makefile's lib/ directory for the static archive
set_target_properties(PSToolbox PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib"
  OUTPUT_NAME "PSToolbox"
)

# ---- Example executable (optional) ----
# If you later add a main program, you can enable this block.
# add_executable(pstoolbox_main src/main.cpp)
# target_link_libraries(pstoolbox_main PRIVATE PSToolbox)
# set_target_properties(pstoolbox_main PROPERTIES
#   RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
# )

# ---- Installation (optional) ----
# install(TARGETS PSToolbox
#   ARCHIVE DESTINATION lib
#   LIBRARY DESTINATION lib
#   RUNTIME DESTINATION bin
# )
# install(DIRECTORY include/ DESTINATION include)

message(STATUS "PSToolbox sources: ${PSTOOLBOX_SOURCES}")
