#include <iostream>
#include <time.h>
#include "ImageDeformation.h"

 using namespace std;    // two of two, yay!
int main(int argc, char * argv[])
{

  if (argc < 7)
    {
    cout << "Program requires 7 input files: " << endl;
    cout << "1. Node file 2. Element file  3. FEM displacements 4. TumorElements.elm if you want tumor to be resected or NO_RESECTION_FILE 5. Preop Image file 6. File name where deformed image is to be stored(with the extension)" << endl;
    return -1; // failure
    } // end if

  // image dimensions
  int xdim,ydim,zdim;
  // voxel spacin
  float pix_x,pix_y,pix_z;
  // timer funcs
  time_t start, end;
  double diff;

  // instantiate the image deformation class
  ImageDeformation * image_deformation = new ImageDeformation;
  
  image_deformation->GetNodePositions(argv[1]);
  image_deformation->GetNodeConnectivityList(argv[2]);
  image_deformation->GetNodeDisplacements(argv[3]);
  image_deformation->GetResectionList(argv[4]);
  
//  
//  image_deformation->Read_Analyze_7_5_Header(argv[5]);
//  image_deformation->Read_Analyze_7_5_Image(argv[6]);
  image_deformation->Read_Image(argv[5]);
  cout << "Finished Readin the input files" << endl;
  image_deformation->PrintInfo();
  
  time(&start);
  image_deformation->Execute();
  time(&end);
  
  diff = difftime(end,start);
  cout << "Time : " << diff << "secs" << endl;

  cout << "Saving the deformed image" << endl;
  image_deformation->SaveDeformedImage(argv[6]);
//  image_deformation->WriteDeformedImage(argv[7]);
  
  delete image_deformation;
  
  return 0; // success
}

