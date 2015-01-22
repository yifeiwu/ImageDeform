
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(IMAGE_DEFORMATION)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")



FIND_PACKAGE(VTK)
IF(VTK_FOUND)
        INCLUDE(${VTK_USE_FILE})
ELSE(VTK_FOUND)
        MESSAGE(FATAL_ERROR "Cannot build MUIGLSPreOp without VTK.  Please set VTK_BUILD_DIR.")
ENDIF(VTK_FOUND)

FIND_PACKAGE(VTKMY)
IF(VTKMY_FOUND)
        INCLUDE(${VTKMY_USE_FILE})
ELSE(VTKMY_FOUND)
        MESSAGE(FATAL_ERROR "Cannot build MUIGLSPreOp without without vtkMy. Please set VTKMY_BUILD_DIR.")
ENDIF(VTKMY_FOUND)



LINK_DIRECTORIES(${VTKMY_DIR}/bin)
INCLUDE_DIRECTORIES (${PETSC_DIR}/conf/variables)
LINK_DIRECTORIES(/usr/local/lib)
INCLUDE(/usr/local/lib)
Include_directories(/usr/local/include/eigen3)
INCLUDE_DIRECTORIES(/usr/local/include/eigen3/Eigen)
Include_directories(/home/caleb/src/eigen/test)
INCLUDE_DIRECTORIES (/usr/local/include)


ADD_EXECUTABLE(TestImageDeformation TestImageDeformation.cxx ImageDeformation.cxx)
Target_link_libraries(TestImageDeformation
     vtkCommon
     vtkGraphics
     vtkHybrid
     vtkRendering
     vtkIO   
     vtkmySPMESHUtils     
     vtkmyIO
)
