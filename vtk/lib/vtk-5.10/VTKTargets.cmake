# Generated by CMake 2.8.7

IF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" LESS 2.5)
   MESSAGE(FATAL_ERROR "CMake >= 2.6.0 required")
ENDIF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" LESS 2.5)
CMAKE_POLICY(PUSH)
CMAKE_POLICY(VERSION 2.6)
#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
SET(CMAKE_IMPORT_FILE_VERSION 1)

# Create imported target vtksys
ADD_LIBRARY(vtksys STATIC IMPORTED)

# Create imported target vtkzlib
ADD_LIBRARY(vtkzlib STATIC IMPORTED)

# Create imported target vtkhdf5
ADD_LIBRARY(vtkhdf5 STATIC IMPORTED)

# Create imported target vtkhdf5_hl
ADD_LIBRARY(vtkhdf5_hl STATIC IMPORTED)

# Create imported target vtkjpeg
ADD_LIBRARY(vtkjpeg STATIC IMPORTED)

# Create imported target vtkpng
ADD_LIBRARY(vtkpng STATIC IMPORTED)

# Create imported target vtktiff
ADD_LIBRARY(vtktiff STATIC IMPORTED)

# Create imported target vtkexpat
ADD_LIBRARY(vtkexpat STATIC IMPORTED)

# Create imported target vtkfreetype
ADD_LIBRARY(vtkfreetype STATIC IMPORTED)

# Create imported target vtklibxml2
ADD_LIBRARY(vtklibxml2 STATIC IMPORTED)

# Create imported target vtkDICOMParser
ADD_LIBRARY(vtkDICOMParser STATIC IMPORTED)

# Create imported target vtkproj4
ADD_LIBRARY(vtkproj4 STATIC IMPORTED)

# Create imported target mpistubs
ADD_LIBRARY(mpistubs STATIC IMPORTED)

# Create imported target MapReduceMPI
ADD_LIBRARY(MapReduceMPI STATIC IMPORTED)

# Create imported target vtkverdict
ADD_LIBRARY(vtkverdict STATIC IMPORTED)

# Create imported target vtkNetCDF
ADD_LIBRARY(vtkNetCDF STATIC IMPORTED)

# Create imported target vtkNetCDF_cxx
ADD_LIBRARY(vtkNetCDF_cxx STATIC IMPORTED)

# Create imported target vtkmetaio
ADD_LIBRARY(vtkmetaio STATIC IMPORTED)

# Create imported target vtksqlite
ADD_LIBRARY(vtksqlite STATIC IMPORTED)

# Create imported target vtkexoIIc
ADD_LIBRARY(vtkexoIIc STATIC IMPORTED)

# Create imported target LSDyna
ADD_LIBRARY(LSDyna STATIC IMPORTED)

# Create imported target vtkalglib
ADD_LIBRARY(vtkalglib STATIC IMPORTED)

# Create imported target vtkEncodeString
ADD_EXECUTABLE(vtkEncodeString IMPORTED)

# Create imported target vtkftgl
ADD_LIBRARY(vtkftgl STATIC IMPORTED)

# Create imported target vtkCommon
ADD_LIBRARY(vtkCommon STATIC IMPORTED)

# Create imported target vtkFiltering
ADD_LIBRARY(vtkFiltering STATIC IMPORTED)

# Create imported target vtkImaging
ADD_LIBRARY(vtkImaging STATIC IMPORTED)

# Create imported target vtkGraphics
ADD_LIBRARY(vtkGraphics STATIC IMPORTED)

# Create imported target vtkGenericFiltering
ADD_LIBRARY(vtkGenericFiltering STATIC IMPORTED)

# Create imported target vtkIO
ADD_LIBRARY(vtkIO STATIC IMPORTED)

# Create imported target vtkRendering
ADD_LIBRARY(vtkRendering STATIC IMPORTED)

# Create imported target vtkVolumeRendering
ADD_LIBRARY(vtkVolumeRendering STATIC IMPORTED)

# Create imported target vtkHybrid
ADD_LIBRARY(vtkHybrid STATIC IMPORTED)

# Create imported target vtkWidgets
ADD_LIBRARY(vtkWidgets STATIC IMPORTED)

# Create imported target vtkInfovis
ADD_LIBRARY(vtkInfovis STATIC IMPORTED)

# Create imported target vtkGeovis
ADD_LIBRARY(vtkGeovis STATIC IMPORTED)

# Create imported target vtkViews
ADD_LIBRARY(vtkViews STATIC IMPORTED)

# Create imported target vtkCharts
ADD_LIBRARY(vtkCharts STATIC IMPORTED)

# Load information for each installed configuration.
GET_FILENAME_COMPONENT(_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
FILE(GLOB CONFIG_FILES "${_DIR}/VTKTargets-*.cmake")
FOREACH(f ${CONFIG_FILES})
  INCLUDE(${f})
ENDFOREACH(f)

# Commands beyond this point should not need to know the version.
SET(CMAKE_IMPORT_FILE_VERSION)
CMAKE_POLICY(POP)
