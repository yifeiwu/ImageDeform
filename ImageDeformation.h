#ifndef __ImageDeformation__h
#define __ImageDeformation__h

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "vtkAnalyzeWriter.h"
#include "vtkShortArray.h"
#include "vtkImageData.h"
#include <vtkPointData.h>
#include <vtkXMLImageDataReader.h>

class ImageDeformation
{
 public:
  // initialize class
  ImageDeformation();
  // destroy class
  ~ImageDeformation();
  // get the x,y,z coordinates
  int GetNodePositions(char * nodefilename);
  // get the 4 nodes that are in a tetrahedron. in - incidence list
  int GetNodeConnectivityList(char * elementfilename);
  // get the nodal displacements/FEM solution
  int GetNodeDisplacements(char * displacementfilename);
  // get the resection list
  int GetResectionList(char * resectionfilename);
  // read in the analyze7.5 pre op image and store the intensity values in OldIntensity
  int Read_Analyze_7_5_Image(char * imagefilename);
  int Read_Image(char * imagefilename); //this reads in a vti
  
  // read in the analyze7.5 header file for the image
  int Read_Analyze_7_5_Header(char * headerfilename);
  // get image dimensions,pixel spacin
  int GetHeaderInfo(char * headerfilename, struct dsr * hdr);
  // save the deformed image
  int SaveDeformedImage(char * outfilename);
  
  
  vtkImageData * DeformedImage;
  void WriteDeformedImage(char * outfilename);
  // member function that deforms the image
  int Execute();
  // print out info about mesh and image
  void PrintInfo();

  // functions that can be called by the user to get info
  int GetNumberOfNodes();
  int GetNumberOfElements();
  float GetPixelWidth();
  float GetPixelHeight();
  float GetPixelDepth();
  int GetImageWidth();
  int GetImageHeight();
  int GetNumberOfSlices();
  const char * GetImageDataType();
  
 protected:
  // get the 3d basis coefficients
  int GetBasisCoefficients(int element, float * a, float * b, float * c, float * d);
  // get the volume of a tet in the deformed image space
  int Volume_DeformedTetrahedron(int element, float * volume);
  // swap bytes
  int Swap_Bytes_in_hdr(struct dsr * pntr);
  int swap_long(unsigned char * ptr);
  int swap_short(unsigned char * ptr);
  // get the minimum and maximum x,y,z coordinates of the 4 vertices of a tet
  int GetMinMax_XYZ_of_element(int element);
  // given a point(x_cube,y_cube,z_cube) check to see if it lies inside the 'element'
  int PointCheck(float x_cube, float y_cube, float z_cube, int element, int * inside_element);
  // get the volume of the 4 sub tets formed with (x_cube,y_cube,z_cube) as one of its vertices
  int Volume_SubTetrahedron1(float x_cube, float y_cube, float z_cube, float * volume_sub_tetrahedron1, int element);
  int Volume_SubTetrahedron2(float x_cube, float y_cube, float z_cube, float * volume_sub_tetrahedron2, int element);
  int Volume_SubTetrahedron3(float x_cube, float y_cube, float z_cube, float * volume_sub_tetrahedron3, int element);
  int Volume_SubTetrahedron4(float x_cube, float y_cube, float z_cube, float * volume_sub_tetrahedron4, int element);
  // Given a point (dx_voxel,dy_voxel,dz_voxel) this function finds the corners of the cube that is formed by the 8 nearest neighbors to the point.
  int FindCornersOfCube(float * dx_voxel, float * dy_voxel, float * dz_voxel,int * x_upper, int * y_upper, int * z_upper, int * x_lower, int * y_lower, int * z_lower);
  // Given a point (dx_voxel,dy_voxel,dz_voxel) this function finds the indices of the 8 nearest neighbors that can be used to find the old intensity values associated with 'em. 
  int Get8NearestNeighbors(float * dx_voxel, float * dy_voxel, float * dz_voxel, int * idxval1, int * idxval2, int * idxval3, int * idxval4, int * idxval5, int * idxval6, int * idxval7, int * idxval8);
  // use trilinear interpolation to find the new intensity values
  int TriLinearInterpolation(float * dx_voxel, float * dy_voxel, float * dz_voxel, int * new_intensity);
  // given the basis coeffs and the volume of the tet get the voxel displacements
  // voxel displacements : dx_voxel = phi(node1)*dx + phi(node2)*dx + phi(node3)*dx + phi(node4)*dx
  int GetVoxelDisplacements(int element, float * a, float * b, float * c, float * d, float x_cube, float y_cube, float z_cube, float volume_tetrahedron, float * dx_voxel, float * dy_voxel, float * dz_voxel);
  // convert the displacements from the mesh space to pixel coordinates
  int ConvertDisplacementsToImageSpace();
  // for the pixel(x_cube,y_cube,z_cube) set the new intensity value in the NewIntensity array
  int SetPixel(float x_cube, float y_cube, float z_cube, int * new_intensity);
  int Initialize_Intensities_and_Displacements();
 private:
  // x,y,z coordinates
  float * x;
  float * y;
  float * z;
  // x,y,z displacements
  float * dx;
  float * dy;
  float * dz;
  // node connectivity list
  int ** in;
  int number_of_nodes;
  int number_of_elements;
  // resection flag
  // flag = 1, if element is a part of the resection, else flag = 0
  int * resection_flag;
  // image dimensions
  int image_dimension[3];
  // voxel spacin
  float pixel_spacing[3];
  // data type
  // refer to Analyze_Header_Structure.h for various types
  int data_type;
  // image size
  int image_size;
  // deformed mesh coordinates
  float * xdx;
  float * ydy;
  float * zdz;
  // (of the 4 nodes in a tet)min and max x,y,z 
  float x_minimum;
  float y_minimum;
  float z_minimum;
  float x_maximum;
  float y_maximum;
  float z_maximum;
  // error criterion
  float TOLERANCE;
   // old intensity - intensity values in the pre op image
  int * OldIntensity;
  // new intensity - intensity values in the deformed image
  int * NewIntensity;

}; // end class

#endif
