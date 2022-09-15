include(FetchContent)

# Replaces the extension of a given file
#
# Args:
#     oldext [in]: Old extension
#     newext [in]: New extension
#     fname [in]: File name in which extension should be replaced.
#     newfname [out]: File name after extension replacement.
#
function(fnet_replace_extension oldext newext fname newfname)

  string(REGEX REPLACE "\\.${oldext}$" ".${newext}" _newfname ${fname})
  set(${newfname} ${_newfname} PARENT_SCOPE)

endfunction()


# Registers files for preprocessing
#
# Args:
#     preproc [in]: Preprocessor to use
#     preprocopts [in]:  Preprocessor command line arguments (but not in/out file)
#     oldext [in]: Extension of the unpreprocessed files.
#     newext [in]: Extension of the preprocessed files.
#     oldfiles [in]: List of unpreprocessed file names.
#     newfiles [out]: List of preprocessed file names.
#
function(fnet_preprocess preproc preprocopts oldext newext oldfiles newfiles)

  set(_newfiles)
  foreach(oldfile IN LISTS oldfiles)
    fnet_replace_extension(${oldext} ${newext} ${oldfile} newfile)
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${newfile}
      COMMAND ${preproc} ${preprocopts} ${CMAKE_CURRENT_SOURCE_DIR}/${oldfile} ${CMAKE_CURRENT_BINARY_DIR}/${newfile}
      MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${oldfile})
    list(APPEND _newfiles ${CMAKE_CURRENT_BINARY_DIR}/${newfile})
  endforeach()
  set(${newfiles} ${_newfiles} PARENT_SCOPE)

endfunction()


# Build -D command line arguments for Fypp preprocessor based on current configuration
#
# Args:
#     fyppflags [inout]: Current Fypp flags on enter, with -D options extended flags on exit.
#
function (fnet_add_fypp_defines fyppflags)

  set(_fyppflags "${${fyppflags}}")

  if(WITH_MPI)
    list(APPEND _fyppflags -DWITH_MPI)
  endif()

  if(WITH_SOCKETS)
    list(APPEND _fyppflags -DWITH_SOCKETS)
  endif()

  if(BUILD_SHARED_LIBS)
    list(APPEND _fyppflags -DBUILD_SHARED_LIBS)
  endif()

  set(${fyppflags} ${_fyppflags} PARENT_SCOPE)

endfunction()


# Gets Fortnet release information.
#
# Args:
#   release [out]: Release string.
#
function(fnet_get_release_name release)

  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/RELEASE)
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/RELEASE DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  else()
    execute_process(
      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/utils/build/update_release ${CMAKE_CURRENT_BINARY_DIR}/RELEASE
      RESULT_VARIABLE exitcode)
    if(NOT exitcode EQUAL 0)
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/RELEASE "(UNKNOWN RELEASE)")
    endif()
  endif()
  file(READ ${CMAKE_CURRENT_BINARY_DIR}/RELEASE _release)
  separate_arguments(_release)
  set(${release} "${_release}" PARENT_SCOPE)

endfunction()


# Prepends a given prefix to a list of items.
#
# Args:
#     result [out]: Variable containing the results
#     prefix [in]: Prefix to add to each item
#     * [in]: List of items to be prefixed
#
function(fnet_add_prefix prefix itemlist result)
  set(_result)
  foreach(var IN LISTS itemlist)
    list(APPEND _result "${prefix}${var}")
  endforeach()
  set(${result} "${_result}" PARENT_SCOPE)
endfunction()


# Returns the parameters needed to create a pkg-config export file
#
# Args:
#     pkgconfig_requires [out]: Value for the Requires field.
#     pkgconfig_libs [out]: Value for the Libs field.
#     pkgconfig_libs_private [out]: Value for the Libs.private field.
#     pkgconfig_c_flags [out]: Value for the cflags field.
#     pkgconfig_prefix [out]: Value for the installation prefix.
#
function(fnet_get_pkgconfig_params pkgconfig_requires pkgconfig_libs pkgconfig_libs_private
    pkgconfig_c_flags)

  set(_pkgconfig_libs "-L${CMAKE_INSTALL_FULL_LIBDIR}")
  fnet_add_prefix("-l" "${PKG_CONFIG_LIBS}" _libs)
  list(APPEND _pkgconfig_libs "${_libs}")

  set(_pkgconfig_libs_private "${PKG_CONFIG_LIBS_PRIVATE}")

  if(PKGCONFIG_LANGUAGE STREQUAL "C")

    set(implibdirs "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}")
    list(REMOVE_ITEM implibdirs "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
    fnet_add_prefix("-L" "${implibdirs}" implibdirs)
    list(APPEND _pkgconfig_libs_private "${implibdirs}")

    set(implibs "${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
    list(REMOVE_ITEM implibs "${CMAKE_C_IMPLICIT_LINK_LIBRARIES}")
    fnet_library_linking_flags("${implibs}" implibs)
    list(APPEND _pkgconfig_libs_private "${implibs}")

    set(_pkgconfig_c_flags
      "-I${CMAKE_INSTALL_FULL_INCLUDEDIR}/${INSTALL_INCLUDEDIR} ${CMAKE_C_FLAGS}")

  else()

    set(implibdirs "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
    list(REMOVE_ITEM implibdirs "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}")
    fnet_add_prefix("-L" "${implibdirs}" implibdirs)
    list(APPEND _pkgconfig_libs_private "${implibdirs}")

    set(implibs "${CMAKE_C_IMPLICIT_LINK_LIBRARIES}")
    list(REMOVE_ITEM implibs "${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
    fnet_library_linking_flags("${implibs}" implibs)
    list(APPEND _pkgconfig_libs_private "${implibs}")

    set(_pkgconfig_c_flags
      "-I${CMAKE_INSTALL_FULL_INCLUDEDIR}/${INSTALL_MODULEDIR} ${CMAKE_Fortran_FLAGS}")

  endif()

  string(REPLACE ";" " " _pkgconfig_libs "${_pkgconfig_libs}")
  string(REPLACE ";" " " _pkgconfig_libs_private "${_pkgconfig_libs_private}")

  set(_pkgconfig_libs_private "${_pkgconfig_libs_private} ${CMAKE_EXE_LINKER_FLAGS}")
  string(REPLACE ";" " " _pkgconfig_requires "${PKG_CONFIG_REQUIRES}")

  set(${pkgconfig_prefix} "${CMAKE_INSTALL_PREFIX}" PARENT_SCOPE)
  set(${pkgconfig_requires} "${_pkgconfig_requires}" PARENT_SCOPE)
  set(${pkgconfig_libs} "${_pkgconfig_libs}" PARENT_SCOPE)
  set(${pkgconfig_libs_private} "${_pkgconfig_libs_private}" PARENT_SCOPE)
  set(${pkgconfig_c_flags} "${_pkgconfig_c_flags}" PARENT_SCOPE)

endfunction()


# Returns library linking flags for given libraries.
#
# If the library is a full path, it is linked directly, otherwise with the "-l" option
#
# Args:
#     libraries [in]: List of libraries to check.
#     linkflags [out]: List of flags to link the libraries
#
function(fnet_library_linking_flags libraries linkflags)
  set(_linkflags)
  foreach(library IN LISTS libraries)
    string(FIND "${library}" "-" dashpos)
    if(dashpos EQUAL 1)
      list(APPEND _linkflags "${library}")
    elseif(IS_ABSOLUTE "${library}")
      list(APPEND _linkflags "${library}")
    else()
      list(APPEND _linkflags "-l${library}")
    endif()
  endforeach()
  set(${linkflags} "${_linkflags}" PARENT_SCOPE)
endfunction()


# Checks the build configuration on consistency and stops in case of inconsistencies
function (fnet_ensure_config_consistency)

  # Check minimal compiler versions
  set(fortran_minimal_versions "GNU;8.3" "Intel;20.2")
  fnet_check_minimal_compiler_version("Fortran" "${fortran_minimal_versions}")

  set(pkgconfig_languages C Fortran)
  list(FIND pkgconfig_languages "${PKGCONFIG_LANGUAGE}" pos)
  if(pos EQUAL -1)
    string(REPLACE ";" "\", \"" pkgconfig_languages_str "${pkgconfig_languages}")
    set(pkgconfig_languages_str "\"${pkgconfig_languages_str}\"")
    message(FATAL_ERROR
      "Invalid language \"${PKGCONFIG_LANGUAGE}\" for PKGCONFIG_LANGUAGE (possible choices: ${pkgconfig_languages_str})")
  endif()

endfunction()


# Stops the code if the source and the build folders are identical.
#
function(fnet_ensure_out_of_source_build)

  get_filename_component(srcdir "${CMAKE_CURRENT_SOURCE_DIR}" REALPATH)
  get_filename_component(bindir "${CMAKE_CURRENT_BINARY_DIR}" REALPATH)

  if("${srcdir}" STREQUAL "${bindir}")
    message(FATAL_ERROR
      "It is not allowed to configure and build Fortnet from its source folder. Please, create a \
separate build directory and invoke CMake from that directory. See the INSTALL.rst file for \
detailed build instructions.")
  endif()

endfunction()


# Makes sure, that a compiler has been already defined for a given language
#
# Args:
#     language [in]: The language to look at.
#
function(fnet_ensure_compiler_def language)

  if(NOT DEFINED CMAKE_${language}_COMPILER)
    message(FATAL_ERROR "Undefined ${language} compiler. The automatic detection of compilers, \
flags and libraries is disabled. You must provide configuration parameters explicitely (e.g. in a \
toolchain file). See the INSTALL.rst file for detailed instructions.")
  endif()

endfunction()


# Loads global build settings (either from config.cmake or from user defined file)
#
macro (fnet_load_build_settings)

  if(NOT DEFINED BUILD_CONFIG_FILE)
    if(DEFINED ENV{FNET_BUILD_CONFIG_FILE}
        AND NOT "$ENV{FNET_BUILD_CONFIG_FILE}" STREQUAL "")
      set(BUILD_CONFIG_FILE "$ENV{FNET_BUILD_CONFIG_FILE}")
    else()
      set(BUILD_CONFIG_FILE "${CMAKE_CURRENT_SOURCE_DIR}/config.cmake")
    endif()
  endif()
  message(STATUS "Reading global build config file: ${BUILD_CONFIG_FILE}")
  include(${BUILD_CONFIG_FILE})

endmacro()


# Sets up the build type.
function (fnet_setup_build_type)

  set(default_build_type "RelWithDebInfo")
  if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Build type: ${default_build_type} (default single-config)")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
  elseif(CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS
      "Build type: ${CMAKE_BUILD_TYPE} (manually selected single-config)")
  else()
    message(STATUS "Build type: Multi-Config (build type selected at the build step)")
  endif()

endfunction()


# Tries to guess which toolchain to load based on the environment.
#
# Args:
#     toolchain [out]: Name of the selected toolchain or undefined if it could not be selected
#
function(fnet_guess_toolchain toolchain)

  if("${CMAKE_Fortran_COMPILER_ID}|${CMAKE_C_COMPILER_ID}" STREQUAL "GNU|GNU")
    set(_toolchain "gnu")
  elseif("${CMAKE_Fortran_COMPILER_ID}|${CMAKE_C_COMPILER_ID}" STREQUAL "Intel|Intel")
    set(_toolchain "intel")
  elseif("${CMAKE_Fortran_COMPILER_ID}|${CMAKE_C_COMPILER_ID}" STREQUAL "NAG|GNU")
    set(_toolchain "nag")
  else()
    set(_toolchain "generic")
  endif()

  set(${toolchain} "${_toolchain}" PARENT_SCOPE)

endfunction()


# Loads toolchain settings.
#
macro(fnet_load_toolchain_settings)

  if(NOT DEFINED TOOLCHAIN_FILE AND NOT "$ENV{FNET_TOOLCHAIN_FILE}" STREQUAL "")
    set(TOOLCHAIN_FILE "$ENV{FNET_TOOLCHAIN_FILE}")
  endif()
  if(NOT DEFINED TOOLCHAIN AND NOT "$ENV{FNETLUS_TOOLCHAIN}" STREQUAL "")
    set(TOOLCHAIN "$ENV{FNETLUS_TOOLCHAIN}")
  endif()
  if(NOT DEFINED TOOLCHAIN_FILE OR TOOLCHAIN_FILE STREQUAL "")
    if(NOT DEFINED TOOLCHAIN OR TOOLCHAIN STREQUAL "")
      fnet_guess_toolchain(TOOLCHAIN)
    endif()
    set(TOOLCHAIN_FILE ${CMAKE_CURRENT_SOURCE_DIR}/sys/${TOOLCHAIN}.cmake)
  endif()
  message(STATUS "Reading build environment specific toolchain file: ${TOOLCHAIN_FILE}")
  include(${TOOLCHAIN_FILE})
endmacro()


# Sets up the global compiler flags
#
macro(fnet_setup_global_compiler_flags)

  if(CMAKE_BUILD_TYPE)
    set(_buildtypes ${CMAKE_BUILD_TYPE})
  else()
    set(_buildtypes ${CMAKE_CONFIGURATION_TYPES})
  endif()
  foreach(_buildtype IN LISTS _buildtypes)
    foreach (lang IN ITEMS Fortran)
      string(TOUPPER "${_buildtype}" _buildtype_upper)
      set(CMAKE_${lang}_FLAGS " ${${lang}_FLAGS}")
      set(CMAKE_${lang}_FLAGS_${_buildtype_upper} " ${${lang}_FLAGS_${_buildtype_upper}}")
      message(STATUS "Flags for ${lang}-compiler (build type: ${_buildtype}): "
        "${CMAKE_${lang}_FLAGS} ${CMAKE_${lang}_FLAGS_${_buildtype_upper}}")
    endforeach()
  endforeach()
  unset(_buildtypes)
  unset(_buildtype)
  unset(_buildtype_upper)
endmacro()


# Builds a package or tries to find it
#
# Args:
#     package [in]: Name of the package.
#     buildpkg [in]: Whether package should be built (otherwise it will be searched for).
#     sourcedir [in]: Directory with package source.
#     exclude [in]: Exclude options when adding sub directory
#
function(fnet_build_or_find_external package buildpkg sourcedir exclude)
  if(${buildpkg})
    add_subdirectory(${sourcedir} ${exclude})
  else()
    find_package(${package} REQUIRED)
  endif()
endfunction()


# Handles a hybrid dependency.
#
# Depending on the list items in the config_methods variable, it will try to:
#
# - checkout the source as a submodule within the origin sub-folder ("Submodule")
# - find the package as external dependency ("Find")
# - fetch the source from a git repository ("Fetch") into the build folder
#
# The methods are tried in the order of their appearance until success the first eligible one.
#
# The methods "Submodule" and "Fetch" would call the passed sub-directory with add_subdirectory()
# passing two variables with the source and binary directory.
#
# Args:
#     package [in]: Name of the dependency to look for.
#     config_methods [in]: Config methods to try
#     target [in]: Name of the target, which must be exported after the configuration.
#     findpkgopts [in]: Options to pass to find_package()
#     subdir [in]: Subdirectory with CMakeFiles.txt for integrating package source.
#     subdiropts [in]: Options to pass to the add_subdir() command.
#     git_repository [in]: Git repository to fetch the package from.
#     git_tag [in]: Git tag to use when fetching the source.
#
# Variables:
#     <UPPER_PACKAGE_NAME>_SOURCE_DIR, <UPPER_PACKAGE_NAME>_BINARY_DIR:
#         Source and binary directories for the build (to pass to add_subdirectory())
#
macro(fnet_config_hybrid_dependency package target config_methods findpkgopts subdir subdiropts
    git_repository git_tag)

  set(_allowed_methods "submodule;find;fetch")
  string(TOLOWER "${package}" _package_lower)
  string(TOUPPER "${package}" _package_upper)

  foreach(_config_method IN ITEMS ${config_methods})

    string(TOLOWER "${_config_method}" _config_lower)
    if(NOT ${_config_lower} IN_LIST _allowed_methods)
      message(FATAL_ERROR "${package}: Unknown configuration method '${_config_method}'")
    endif()

    if("${_config_lower}" STREQUAL "find")

      message(STATUS "${package}: Trying to find installed package")
      find_package(${package} ${findpkgopts})
      if(${package}_FOUND)
        message(STATUS "${package}: Installed package found")
        break()
      else()
        message(STATUS "${package}: Installed package could not be found")
      endif()

    elseif("${_config_lower}" STREQUAL "submodule")

      if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/origin/CMakeLists.txt
          AND GIT_WORKING_COPY)
        message(STATUS "${package}: Downloading via git submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init ${subdir}/origin
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
      endif()

      if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/origin/CMakeLists.txt)
        message(STATUS "${package}: Using source in ${subdir}/origin")
        set(${_package_upper}_SOURCE_DIR "origin")
        set(${_package_upper}_BINARY_DIR)
        add_subdirectory(${subdir} ${subdiropts})
        break()
      endif()

    elseif("${_config_lower}" STREQUAL "fetch")

      message(STATUS "${package}: Fetching from repository ${git_repository}@${git_tag}")
      FetchContent_Declare(${_package_lower} GIT_REPOSITORY ${git_repository} GIT_TAG ${git_tag})
      FetchContent_GetProperties(${_package_lower})
      if(NOT ${_package_lower}_POPULATED)
        FetchContent_Populate(${_package_lower})
      endif()
      set(${_package_upper}_SOURCE_DIR "${${_package_lower}_SOURCE_DIR}")
      set(${_package_upper}_BINARY_DIR "${${_package_lower}_BINARY_DIR}")
      add_subdirectory(${subdir} ${subdiropts})
      break()

    endif()

  endforeach()

  if(NOT TARGET ${target})
    message(FATAL_ERROR "Could not configure ${package} to export target '${target}'")
  endif()

  unset(_allowed_methods)
  unset(_package_lower)
  unset(_package_upper)
  unset(_config_method)
  unset(_config_lower)

endmacro()


# Checks whether the current compiler fullfills minimal version requirements.
#
#
# Arguments:
#   lang [in]: Language for which the compiler should be checked (e.g. Fortran, C, CXX)
#   compiler_versions [in]: List with alternating compiler ids and minimal version numbers, e.g.
#       "Intel;19.0;GNU;9.0". If the compiler is amoung the listed ones and its version number is
#       less than the specified one, a fatal error message will be issued. Otherwise the function
#       returns silently.
#
function (fnet_check_minimal_compiler_version lang compiler_versions)
  while(compiler_versions)
    list(POP_FRONT compiler_versions compiler version)
    if("${CMAKE_${lang}_COMPILER_ID}" STREQUAL "${compiler}"
        AND CMAKE_${lang}_COMPILER_VERSION VERSION_LESS "${version}")
      message(FATAL_ERROR "${compiler} ${lang} compiler is too old "
          "(found \"${CMAKE_${lang}_COMPILER_VERSION}\", required >= \"${version}\")")
    endif()
  endwhile()
endfunction()
