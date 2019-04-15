#----------------------------------------------------------------
# Generated CMake target import file for configuration "release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "szip-static" for configuration "release"
set_property(TARGET szip-static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(szip-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libszip.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS szip-static )
list(APPEND _IMPORT_CHECK_FILES_FOR_szip-static "${_IMPORT_PREFIX}/lib/libszip.a" )

# Import target "szip-shared" for configuration "release"
set_property(TARGET szip-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(szip-shared PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libszip.so.2.1"
  IMPORTED_SONAME_RELEASE "libszip.so.2.1"
  )

list(APPEND _IMPORT_CHECK_TARGETS szip-shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_szip-shared "${_IMPORT_PREFIX}/lib/libszip.so.2.1" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
