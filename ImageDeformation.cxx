#include "ImageDeformation.h"
#include "AnalyzeHeaderStructure.h"
#include <string.h>


using namespace std;    
ImageDeformation::ImageDeformation()
{
  this->x = NULL;
  this->y = NULL;
  this->z = NULL;
  this->dx = NULL;
  this->dy = NULL;
  this->dz = NULL;
  this->OldIntensity = NULL;
  
  for (int i=0; i<3; i++)
    {
    this->image_dimension[i] = 0;
    this->pixel_spacing[i] = 0.0;
    } // end for i
  this->in = NULL;
  this->resection_flag = NULL;
  this->image_size = 0;
  this->TOLERANCE = 1e-2;
  this->NewIntensity = NULL;
  this->DeformedImage = NULL;
  this->DeformedImage = vtkImageData::New();
  
} // end constructor

ImageDeformation::~ImageDeformation()
{
  if (x)
    {
    delete [] this->x;
    } // end if
  if (y)
    {
    delete [] this->y;
    } // end if
  if (z)
    {
    delete [] this->z;
    }
  if (in)
    {
    for (int i=0; i<4; i++)
      {
      delete []  this->in[i];
      } // end for i
    delete [] this->in;
    } // end if
  if (resection_flag)
    {
    delete [] resection_flag;
    }
  if (dx)
    {
    delete [] this->dx;
    } // end if
  if (dy)
    {
    delete [] this->dy;
    } // end if
  if (dz)
    {
    delete [] this->dz;
    } // end if
  if (OldIntensity)
    {
    delete [] this->OldIntensity;
    } // end if  
  if (xdx)
    {
    delete [] this->xdx;
    } // end if
  if (ydy)
    { 
    delete [] this->ydy;
    } // end if
  if (zdz)
    {
    delete [] this->zdz;
    } // end if
  if (NewIntensity)
    {
    delete [] this->NewIntensity;
    }
} // end destructor

int ImageDeformation::GetNodePositions(char * nodefilename)
{
  
  int i,I;
  // temp floats to hold the x,y,z positions
  float value1, value2, value3;
  // temp vector to gold the x,y,z positions
  std::vector<float> X;
  std::vector<float> Y;
  std::vector<float> Z;
  // file input stream for node file
  ifstream nodfile;
  nodfile.open(nodefilename,ios::in);

  // error check for input file
  if (!nodfile)
    {
    cout << "Error openin file: " << nodefilename << " check file name" << endl;
    return -1; // failure
    } // end if
  
  while (nodfile >> I >> value1 >> value2 >> value3)
    {
    X.push_back(value1);
    Y.push_back(value2);
    Z.push_back(value3);
    } // end while
  nodfile.close();
  
  // initialize the x,y,z arrays
  this->x = new float[X.size()];
  this->y = new float[Y.size()];
  this->z = new float[Z.size()];
  
  // get the number of nodes
  this->number_of_nodes = X.size();
  
  // store the node coordiantes in the arrays
  for (i=0; i<this->number_of_nodes; i++)
    {
    this->x[i] = X[i];
    this->y[i] = Y[i];
    this->z[i] = Z[i];
    } // end for i
  
  // destroy the vectors
  X.clear();
  Y.clear();
  Z.clear();

  return 0; // success
} // end getnodepositions
  
int ImageDeformation::GetNodeConnectivityList(char * elementfilename)
{
  
  int i,I;
  // temp ints to hold the 4 nodes that form a tet
  int value1, value2, value3, value4, value5;

  // temp vectors that hold the 4 nodes that form a tet
  std::vector<int> node1;
  std::vector<int> node2;
  std::vector<int> node3;
  std::vector<int> node4;
  std::vector<int> material_type;
  // input file stream for handlin the element file
  ifstream elmfile;
  elmfile.open(elementfilename,ios::in);

  // error check for filename
  if (!elmfile)
    {
    cout << "Error openin file: " << elementfilename << "check file name" << endl;
    return -1; // failure
    } // end if

  while (elmfile >> I >> value1 >> value2 >> value3 >> value4 >> value5)
    {
    // NML files have a one based index. therefore subtract one from each node
    node1.push_back(value1-1);
    node2.push_back(value2-1);
    node3.push_back(value3-1);
    node4.push_back(value4-1);
    material_type.push_back(value5);
    } // end while
  elmfile.close();
    
  // init. the node connectivity list
  this->in = new int*[node1.size()];
  for (i=0; i<node1.size(); i++)
    {
    in[i] = new int[4];
    } // end for i
    
  // get the number of elements
  this->number_of_elements = node1.size();

  // store the 4 nodes of the tet in the incidence list
  for (i=0; i<this->number_of_elements; i++)
    {
    this->in[i][0] = node1[i]; 
    this->in[i][1] = node2[i];
    this->in[i][2] = node3[i];
    this->in[i][3] = node4[i];
    } // end for i
  
  // destoy the vectors
  node1.clear();
  node2.clear();
  node3.clear();
  node4.clear();

  return 0; // success
} // end getnodeconnectivitylist

int ImageDeformation::GetNodeDisplacements(char * displacementfilename)
{
  
  int i,I, number_of_dnodes;
  // temp floats to hold the x,y,z displacements
  float value1, value2, value3;
  // temp vector to gold the x,y,z displacements
  std::vector<float> DX;
  std::vector<float> DY;
  std::vector<float> DZ;
  // file input stream for displacement file
  ifstream dispfile;
  dispfile.open(displacementfilename,ios::in);

  // error check for input file
  if (!dispfile)
    {
    cout << "Error openin file: " << displacementfilename << " check file name" << endl;
    return -1; // failure
    } // end if
  
  while (dispfile >> I >> value1 >> value2 >> value3)
    {
    DX.push_back(value1);
    DY.push_back(value2);
    DZ.push_back(value3);
    } // end while
  dispfile.close();
  
  // initialize the dx,dy,dz arrays
  this->dx = new float[DX.size()];
  this->dy = new float[DY.size()];
  this->dz = new float[DZ.size()];
  
  // get the number of nodes
  number_of_dnodes = DX.size();
  
  // sanity check
  if (number_of_dnodes != number_of_nodes)
    {
    cout << " Displacement file is not of the same size as the node file.Check" << endl;
    return -1; // failure
    } // end if

  // store the node displacements in the arrays
  for (i=0; i<this->number_of_nodes; i++)
    {
    this->dx[i] = DX[i];
    this->dy[i] = DY[i];
    this->dz[i] = DZ[i];
    } // end for i
 
  // destroy the vectors
  DX.clear();
  DY.clear();
  DZ.clear();

  return 0; // success
} // end get node displacements

int ImageDeformation::GetResectionList(char * resectionfilename)
{

  int i, I, tmp_resection_elm;
  // temp vector to hold all the resection elements
  std::vector<float> tmp_resection_elm_list;

  // init. resection flag
  this->resection_flag = new int[this->number_of_elements];

  for (i=0; i < this->number_of_elements; i++)
    {
    this->resection_flag[i] = 0;
    } // end for 'i'

  // flag resection elements if user inputs a resction file
  if (strcmp(resectionfilename,"NO_RESECTION_FILE")!=0)
    {
    // input file stream
    ifstream infile;
    infile.open(resectionfilename,ios::in);

    // sanity check
    if (!infile)
      {
      cout << "Resection File " << resectionfilename << " does not exist. Check" << endl;
      return -1; // failure
      } // end if

    // read in the resection elements
    // am assumin that the indexin and element ids are one-based
    while (infile >> I >> tmp_resection_elm)
      {
      tmp_resection_elm_list.push_back(tmp_resection_elm);
      } // end while
    infile.close();

    // now flag the resection elements to 1
    for (i=0; i < tmp_resection_elm_list.size(); i++)
      {
      tmp_resection_elm = int(tmp_resection_elm_list[i]);
      this->resection_flag[tmp_resection_elm] = 1;
      } // end for 'i'
    } // end if

  return 0; // success
} // end getresectionlist

int ImageDeformation::Read_Analyze_7_5_Header(char * headerfilename)
{
  
  // define a structure of the type dsr
  struct dsr hdr;
    
  // file pointer for the header file
  FILE *fptr;

  if ( (fptr = fopen(headerfilename,"r")) == NULL )
    {
    cout << " Error Openin file: " << headerfilename << "Check file name" << endl;
    return -1; // failure
    } // end if

  
  // read in the header
  fread(&hdr,1,sizeof(struct dsr),fptr);
  
  // swap bytes if the number of dimensions in database is bet 0 and 15
  // usually the no. of dimensions is 4
  if (hdr.dime.dim[0] < 0 || hdr.dime.dim[0] > 15)
    {
    this->Swap_Bytes_in_hdr(&hdr);
    } // end if
  
  GetHeaderInfo(headerfilename,&hdr);
  
  return 0; // success  
} // end readanalyze75header

int ImageDeformation::Swap_Bytes_in_hdr(struct dsr * pntr)
{

  // dim[i]:image dimension is of type short int. 
  // therefore swap 2 bytes
  this->swap_short((unsigned char *)&pntr->dime.dim[1]);
  this->swap_short((unsigned char *)&pntr->dime.dim[2]);
  this->swap_short((unsigned char *)&pntr->dime.dim[3]);
  // bitpix: bits per pixel is of type short int. 
  // therefore swap short 2 bytes
  this->swap_short((unsigned char *)&pntr->dime.datatype);
  this->swap_short((unsigned char *)&pntr->dime.bitpix);
  // pixdim: pixel spacin is of type float.
  // therefore swap 4 bytes
  this->swap_long((unsigned char *)&pntr->dime.pixdim[0]);
  this->swap_long((unsigned char *)&pntr->dime.pixdim[1]);
  this->swap_long((unsigned char *)&pntr->dime.pixdim[2]);
  this->swap_long((unsigned char *)&pntr->dime.pixdim[3]);
  
  return 0; // success
} // end swapbytesinhdr

int ImageDeformation::swap_long(unsigned char * ptr)
{

  unsigned char b0, b1, b2, b3;

  b0 = *ptr;
  b1 = *(ptr+1);
  b2 = *(ptr+2);
  b3 = *(ptr+3);

  *ptr = b3;
  *(ptr+1) = b2;
  *(ptr+2) = b1;
  *(ptr+3) = b0;

  return 0; // success
} // end swaplong

int ImageDeformation::swap_short(unsigned char * ptr)
{

  unsigned char b0, b1;

  b0 = *ptr;
  b1 = *(ptr+1);

  *ptr = b1;
  *(ptr+1) = b0;

  return 0; // success
} // end swapshort

int ImageDeformation::GetHeaderInfo(char * headerfilename, struct dsr * hdr)
{

  int i;
  
  // get the image dimensions and the pixel spacing
  // also convert the pixel spacing to 'm' from 'mm' units
  for (i=0; i<3; i++)
    {
    this->image_dimension[i] = hdr->dime.dim[i+1];
    this->pixel_spacing[i] = (hdr->dime.pixdim[i+1])/1000;
    } // end for i
  
  this->data_type = hdr->dime.datatype;

  return 0; // success
} // end getheaderinfo

int ImageDeformation::Initialize_Intensities_and_Displacements()
{

  int i;

  // init. the Old Intensity
  this->OldIntensity = new int[this->image_size];
  // init. the new intensity array
  this->NewIntensity = new int[this->image_size];
  // fill the array with zeros
  for (i=0; i<this->image_size; i++)
    {
    this->NewIntensity[i] = 0;
    } // end for i

  // init the deformed mesh coordinates
  this->xdx = new float[this->number_of_nodes];
  this->ydy = new float[this->number_of_nodes];
  this->zdz = new float[this->number_of_nodes];
  
  return 0; // success
} // end initializeintensitiesanddisplacements
  
int ImageDeformation::Read_Analyze_7_5_Image(char * imagefilename)
{
  
  int i;

  // get the size of the image
  this->image_size = this->image_dimension[0] * this->image_dimension[1] * this->image_dimension[2];
  // might not be the right place to init displacements.this makes it easier for the parallel code
  this->Initialize_Intensities_and_Displacements();

  // for some reason ifstream did not work with binary file.
  FILE * imgfile;
  imgfile = fopen(imagefilename,"rb");

  if (imgfile == NULL) 
    {
    cout << "Error openin file: " << imagefilename << " Check file name" << endl;
    return -1; // failure
    } // end if

  // we do not know what the data type is to begin with
  // therefore read in one byte at a time and then typecast it
  if (this->data_type == 2)
    {
    for (i=0; i<this->image_size; i++)
      {
      unsigned char intensity = 0;
      // read in the intensity values
      fread(&intensity,sizeof(unsigned char),1,imgfile);
      this->OldIntensity[i] = static_cast<int> (intensity);
      } // end for i
    } // end if

  if (this->data_type == 4)
    {
     for (i=0; i<this->image_size; i++)
      {
      signed short intensity = 0;
      // read in the intensity values
      fread(&intensity,sizeof(signed short),1,imgfile);
      this->OldIntensity[i] = static_cast<int> (intensity);
      } // end for i
    } // end if
     
  if (this->data_type == 8)
    {
     for (i=0; i<this->image_size; i++)
      {
      signed int intensity = 0;
      // read in the intensity values
      fread(&intensity,sizeof(signed int),1,imgfile);
      this->OldIntensity[i] = static_cast<int> (intensity);
      } // end for i
    } // end if

  if (this->data_type == 16)
    {
    for (i=0; i<this->image_size; i++)
      {
      float intensity = 0;
      // read in the intensity values
      fread(&intensity,sizeof(float),1,imgfile);
      this->OldIntensity[i] = static_cast<int> (intensity);
      } // end for i
    } // end if

  if (this->data_type == 64)
    {
    for (i=0; i<this->image_size; i++)
      {
      double intensity = 0;
      // read in the intensity values
      fread(&intensity,sizeof(double),1,imgfile);
      this->OldIntensity[i] = static_cast<int> (intensity);
      } // end for i
    } // end if

  fclose(imgfile);

  return 0; // success
} // end readanalyze75image







int ImageDeformation::Read_Image(char * imagefilename)
{
  
  int i, j, k, index;

  // read image file
  vtkXMLImageDataReader * image = vtkXMLImageDataReader::New();
  image->SetFileName( imagefilename );
  image->Update();

  vtkImageData * imagedata = vtkImageData::New();
  imagedata->ShallowCopy(image->GetOutput());
  imagedata->UpdateInformation();
//  this->img = imagedata;
//  this->img->Update();

  // get the size of the image
  
  
//  this->data_type = imagedata->GetScalarType()-1;
//  cout<<"type " <<this->data_type<<endl;

this->data_type=8;

  int *dims;
  dims = imagedata->GetDimensions();
  this->image_size = dims[0] * dims[1] * dims[2];
  this->image_dimension[0] = dims[0];
  this->image_dimension[1] = dims[1];
  this->image_dimension[2] = dims[2];

  cout<<dims[0]<<" "<<dims[1]<<" "<<dims[2]<<endl;
  cout<<this->image_size<<endl;

  double *spacing;
  spacing = imagedata->GetSpacing();
  this->pixel_spacing[0] = spacing[0]/1000;
  this->pixel_spacing[1] = spacing[1]/1000;
  this->pixel_spacing[2] = spacing[2]/1000;

  this->Initialize_Intensities_and_Displacements();

  for (i=0; i<dims[2]; i++)
  {
    for (j=0; j<dims[1]; j++)
    {
       for (k=0; k<dims[0]; k++)
       {
          index = i*this->image_dimension[0]*this->image_dimension[1]
                  + j*this->image_dimension[1]
                  +  k;
          this->OldIntensity[index] = imagedata->GetScalarComponentAsDouble(k,j,i,0);
          // cout<<"index: "<<index<<" "<<k<<" "<<j<<" "<<i<<" "<<this->OldIntensity[index]<<" "<<imagedata->GetScalarComponentAsDouble(k,j,i,0)<<endl;

       }
    }
  }

  int *extent;
  extent = imagedata->GetExtent();
  cout<<extent[0]<<" "<<extent[1]<<" "<<extent[2]<<" "<<extent[3]<<" "<<extent[4]<<" "<<extent[5]<<endl;

  return 0; // success
} // end readimage










































int ImageDeformation::Execute()
{

  int i,j,k;
  // x,y,z coordinates of the cube used to determine if voxels lie within an element
  float x_cube,y_cube,z_cube;
  // basis coefficients
  float a[4], b[4], c[4], d[4];
  // volume of deformed tetrahedron 
  // i.e., volume of tetrahedron after deformation
  float volume_tetrahedron;
  // check to see if point lies inside an element
  // inside_element = 1 : point within element
  // inside_element = 0 : point outside element
  int inside_element;
  // voxel displacements
  float dx_voxel, dy_voxel, dz_voxel;
  // new intensity : intensity associated with the deformed image
  int new_intensity;
  // used to display the no. of elements done
  int statcount = 0;
    
  this->ConvertDisplacementsToImageSpace();
  
  for (i=0; i<this->number_of_elements; i++)
    {
      
    // counter used to display the progress of work
    statcount++ ;
    if (statcount == (int)ceil(this->number_of_elements/100))
      {
      cout << "Elements done: " << i << endl;
      statcount = 0;
      } // end if

    // proceed if element is not a part of the resection list
    if (this->resection_flag[i] == 0)
      {
    
      // for a given element get the min and max x,y,z coordinates of the nodes
      // that form the element
      this->GetMinMax_XYZ_of_element(i);

      // this loop will move through the voxel cube bounded by (x_minimum,y_minimum,z_minimum) & (x_maximum,y_maximum,z_maximum) determinining which voxels are in anelement
      for (x_cube=this->x_minimum; x_cube<this->x_maximum; x_cube++)
	{
	for(y_cube=this->y_minimum; y_cube<this->y_maximum; y_cube++)
	  {
          for (z_cube=this->z_minimum; z_cube<this->z_maximum; z_cube++)	  
	    {
	    this->PointCheck(x_cube,y_cube,z_cube,i,&inside_element);
	    if (inside_element == 1)
	      {
	      // get the basis coefficients
	      for (j=0; j<4; j++)
		{
		a[j] = 0.0;
		b[j] = 0.0;
		c[j] = 0.0;
		d[j] = 0.0;
		} // end for j
	      this->GetBasisCoefficients(i,a,b,c,d);
  
	      // get the volume of the deformed tetrahedron
	      this->Volume_DeformedTetrahedron(i,&volume_tetrahedron);

	      // init the voxel displacements
	      dx_voxel = 0.0;
	      dy_voxel = 0.0;
	      dz_voxel = 0.0;
	      // given the basis coeffs and the volume of the tet get the voxel displacements
	      // voxel displacements : dx_voxel = phi(node1)*dx + phi(node2)*dx + phi(node3)*dx + phi(node4)*dx
	      this->GetVoxelDisplacements(i,a,b,c,d,x_cube,y_cube,z_cube,volume_tetrahedron,&dx_voxel,&dy_voxel,&dz_voxel);

	      // given the voxel displacement find the intensity associated with the pixel (dx_voxel,dy_voxel,dz_voxel) thru trilinear interpolation
	      this->TriLinearInterpolation(&dx_voxel,&dy_voxel,&dz_voxel,&new_intensity);	    
	      // fill in the intensity value in the new intensity array
	      this->SetPixel(x_cube,y_cube,z_cube,&new_intensity);
	      } // end if
	    } // end for z_cube
	  } // end y_cube
	} // end x_cube
      } // end if - check for resection flag
    } // end for 'element loop'
    
  return 0; // success
} // end execute

int ImageDeformation::GetVoxelDisplacements(int element, float * a, float * b, float * c, float * d, float x_cube, float y_cube, float z_cube, float volume_tetrahedron, float * dx_voxel, float * dy_voxel, float * dz_voxel)
{

  int k;
  // basis coeff phi = a + bx + cy + dz
  float phij;
  
  for (k=0; k<4; k++)
    {
    // phi = a + bx + cy + dz
    phij = (a[k] + b[k]*x_cube + c[k]*y_cube + d[k]*z_cube)/6/volume_tetrahedron;
    // get the voxel displacements
    // dx_voxel = phi(node1)*dx + phi(node2) + dx + phi(node3) * dx + phi(node4) * dx
   
    *dx_voxel += this->dx[this->in[element][k]]*phij;
    *dy_voxel += this->dy[this->in[element][k]]*phij;
    *dz_voxel += this->dz[this->in[element][k]]*phij;
    } // end for k
 
  // project backwards
  *dx_voxel = x_cube - *dx_voxel;
  *dy_voxel = y_cube - *dy_voxel;
  *dz_voxel = z_cube - *dz_voxel;

  return 0; // success
} // end getvoxeldisplacements

int ImageDeformation::ConvertDisplacementsToImageSpace()
{

  int i;
  // get the deformed mesh coordinates and convert it to image space
  for (i=0; i<this->number_of_nodes; i++)
    {
    this->xdx[i] = (this->x[i] + this->dx[i])/this->pixel_spacing[0];
    this->ydy[i] = (this->y[i] + this->dy[i])/this->pixel_spacing[1];
    this->zdz[i] = (this->z[i] + this->dz[i])/this->pixel_spacing[2];
    } // end for i

  // convert the mesh displacements to image space
  for (i=0; i<this->number_of_nodes; i++)
    {
    this->dx[i] = this->dx[i]/this->pixel_spacing[0];
    this->dy[i] = this->dy[i]/this->pixel_spacing[1];
    this->dz[i] = this->dz[i]/this->pixel_spacing[2];
    } // end for i

  return 0; // success
} // end convertdisplacementstoimagespace

int ImageDeformation::SetPixel(float x_cube, float y_cube, float z_cube, int * new_intensity)
{

  int i;
  // index value that is used to find the intensity 
  int nx, ny, nz, index;
    
  nx = (int)floor(x_cube);
  ny = (int)floor(y_cube);
  nz = (int)floor(z_cube);
  index = (nz-1)*this->image_dimension[0]*this->image_dimension[1]
        + (ny-1)*this->image_dimension[1]
        +  nx;

  // sanity check. ensures that the index value is not greater than the array size of new intensity
  if (index > this->image_size)
    {
    cout << "Index " << index << "is greater than the image size " <<this->image_size << endl;
    return -1; // failure
    } // end if
  
  this->NewIntensity[index] = *new_intensity;
  
  return 0; // success
} // end setpixel

int ImageDeformation::TriLinearInterpolation(float * dx_voxel, float * dy_voxel, float * dz_voxel, int * new_intensity)
{
  
  // indices of the 8 nearest neighbors to the point (dx_voxel,dy_voxel,dz_voxel)
  int idxval1, idxval2, idxval3, idxval4, idxval5, idxval6, idxval7, idxval8;
  // corners of the cube formed by the 8 nearest neighbors 
  int x_upper, y_upper, z_upper, x_lower, y_lower, z_lower;

  // Given a point (x,y,z) this function finds the indices of the 8 nearest neighbors
  this->Get8NearestNeighbors(dx_voxel,dy_voxel,dz_voxel,&idxval1,&idxval2,&idxval3,&idxval4,&idxval5,&idxval6,&idxval7,&idxval8);
  
  //Given a point (x,y,z) find the corners of the cube that is formed by the 8 nearest neighbors to the point.
  this->FindCornersOfCube(dx_voxel,dy_voxel,dz_voxel,&x_upper,&y_upper,&z_upper,&x_lower,&y_lower,&z_lower);

  *new_intensity = (int)( (this->OldIntensity[idxval1]*(x_upper-*dx_voxel)*(y_upper-*dy_voxel)*(z_upper-*dz_voxel))
			  + (this->OldIntensity[idxval2]*(*dx_voxel-x_lower)*(y_upper-*dy_voxel)*(z_upper-*dz_voxel))
			  + (this->OldIntensity[idxval3]*(x_upper-*dx_voxel)*(*dy_voxel-y_lower)*(z_upper-*dz_voxel))
			  + (this->OldIntensity[idxval4]*(x_upper-*dx_voxel)*(y_upper-*dy_voxel)*(*dz_voxel-z_lower))
			  + (this->OldIntensity[idxval5]*(*dx_voxel-x_lower)*(y_upper-*dy_voxel)*(*dz_voxel-z_lower))
			  + (this->OldIntensity[idxval6]*(x_upper-*dx_voxel)*(*dy_voxel-y_lower)*(*dz_voxel-z_lower))
			  + (this->OldIntensity[idxval7]*(*dx_voxel-x_lower)*(*dy_voxel-y_lower)*(z_upper-*dz_voxel))
			  + (this->OldIntensity[idxval8]*(*dx_voxel-x_lower)*(*dy_voxel-y_lower)*(*dz_voxel-z_lower)) );
 
  return 0; // success
} // end trilinearinterpolation

int ImageDeformation::FindCornersOfCube(float * dx_voxel, float * dy_voxel, float * dz_voxel,int * x_upper, int * y_upper, int * z_upper, int * x_lower, int * y_lower, int * z_lower)
{
  
  *x_lower = (int)floor(*dx_voxel);
  *y_lower = (int)floor(*dy_voxel);
  *z_lower = (int)floor(*dz_voxel);

  *x_upper = (int)ceil(*dx_voxel);
  *y_upper = (int)ceil(*dy_voxel);
  *z_upper = (int)ceil(*dz_voxel);

  return 0;//success
} // end findcornersofcube

int ImageDeformation::Get8NearestNeighbors(float * dx_voxel, float * dy_voxel, float * dz_voxel, int * idxval1, int * idxval2, int * idxval3, int * idxval4, int * idxval5, int * idxval6, int * idxval7, int * idxval8) 
{
  
  // x,y,z coordinates of the first nearest neighbor
  int neighbor1_x = (int)floor(*dx_voxel);
  int neighbor1_y = (int)floor(*dy_voxel);
  int neighbor1_z = (int)floor(*dz_voxel);
  // index value that gives the intensity value associated with the first neighbor
  *idxval1 = (neighbor1_z-1)*this->image_dimension[0]*this->image_dimension[1] + (neighbor1_y-1)*this->image_dimension[1] + neighbor1_x;

  //repeat the same for neighbors 2 thru 8
  int neighbor2_x = (int) ceil(*dx_voxel);
  int neighbor2_y = (int) floor(*dy_voxel);
  int neighbor2_z = (int) floor(*dz_voxel);

  *idxval2 = (neighbor2_z-1)*this->image_dimension[0]*this->image_dimension[1] + (neighbor2_y-1)*this->image_dimension[1] + neighbor2_x;

  int neighbor3_x = (int) floor(*dx_voxel);
  int neighbor3_y = (int) ceil(*dy_voxel);
  int neighbor3_z = (int) floor(*dz_voxel);

  *idxval3 = (neighbor3_z-1)*this->image_dimension[0]*this->image_dimension[1] + (neighbor3_y-1)*this->image_dimension[1] + neighbor3_x;

  int neighbor4_x = (int) floor(*dx_voxel);
  int neighbor4_y = (int) floor(*dy_voxel);
  int neighbor4_z = (int) ceil(*dz_voxel);

  *idxval4 = (neighbor4_z-1)*this->image_dimension[0]*this->image_dimension[1] + (neighbor4_y-1)*this->image_dimension[1] + neighbor4_x;

  int neighbor5_x = (int) ceil(*dx_voxel);
  int neighbor5_y = (int) floor(*dy_voxel);
  int neighbor5_z = (int) ceil(*dz_voxel);

  *idxval5 = (neighbor5_z-1)*this->image_dimension[0]*this->image_dimension[1] + (neighbor5_y-1)*this->image_dimension[1] + neighbor5_x;

  int neighbor6_x = (int) floor(*dx_voxel);
  int neighbor6_y = (int) ceil(*dy_voxel);
  int neighbor6_z = (int) ceil(*dz_voxel);

  *idxval6 = (neighbor6_z-1)*this->image_dimension[0]*this->image_dimension[1] + (neighbor6_y-1)*this->image_dimension[1] + neighbor6_x;

  int neighbor7_x = (int) ceil(*dx_voxel);
  int neighbor7_y = (int) ceil(*dy_voxel);
  int neighbor7_z = (int) floor(*dz_voxel);

  *idxval7 = (neighbor7_z-1)*this->image_dimension[0]*this->image_dimension[1] + (neighbor7_y-1)*this->image_dimension[1] + neighbor7_x;

  int neighbor8_x = (int) ceil(*dx_voxel);
  int neighbor8_y = (int) ceil(*dy_voxel);
  int neighbor8_z = (int) ceil(*dz_voxel);

  *idxval8 = (neighbor8_z-1)*this->image_dimension[0]*this->image_dimension[1] + (neighbor8_y-1)*this->image_dimension[1] + neighbor8_x;

  // sanity check. ensures that none of the index values are greater than the array size of new intensity
  if (*idxval1>this->image_size || *idxval2>this->image_size || *idxval3>this->image_size || *idxval4>this->image_size || *idxval5>this->image_size || *idxval6>this->image_size || *idxval7>this->image_size || *idxval8>this->image_size)
    {
    cout << "Check your nearest neighbors. The index values are greater than the image size " << endl;
    return -1; // failure
    } // end if

  return 0; // success
} // end get8nearestneighbors
  
int ImageDeformation::GetMinMax_XYZ_of_element(int element)
{

  int i;
  // store the x,y,z coordinates of the 4 vertices of a tet
  float x_element[4], y_element[4], z_element[4];
  
  for (i=0; i<4; i++)
    {
    x_element[i] = this->xdx[this->in[element][i]];
    y_element[i] = this->ydy[this->in[element][i]];
    z_element[i] = this->zdz[this->in[element][i]];
    } // end for i

  this->x_minimum = x_element[0];
  this->y_minimum = y_element[0];
  this->z_minimum = z_element[0];

  this->x_maximum = x_element[0];
  this->y_maximum = y_element[0];
  this->z_maximum = z_element[0];

  for (i=0; i<4; i++)
    {
    // find the x minimum
    if (x_element[i] < this->x_minimum)
      {
      this->x_minimum = x_element[i];
      } // end if
    // find the y minimum
    if (y_element[i] < this->y_minimum)
      {
      this->y_minimum = y_element[i];
      } // end if
    // find the z minimum
    if (z_element[i] < this->z_minimum)
      {
      this->z_minimum = z_element[i];
      } // end if

    // find the x maximum
    if (x_element[i] > this->x_maximum)
      {
      this->x_maximum = x_element[i];
      } // end if
    // find the y maximum
    if (y_element[i] > this->y_maximum)
      {
      this->y_maximum = y_element[i];
      } // end if
    // find the z maximum
    if (z_element[i] > this->z_maximum)
      {
      this->z_maximum = z_element[i];
      } // end if

    } // end for i
  
  // round max and min to the closest integer and add half a voxel 
  // for the startin point
  this->x_minimum = floor(this->x_minimum-1)-0.51;
  this->y_minimum = floor(this->y_minimum-1)-0.51;
  this->z_minimum = floor(this->z_minimum-1)-0.51;
  this->x_maximum = (int)(this->x_maximum+1)+0.51;
  this->y_maximum = (int)(this->y_maximum+1)+0.51;
  this->z_maximum = (int)(this->z_maximum+1)+0.51;
  
  return 0; // sucess
} // end getminmaxxyzforelement


int ImageDeformation::PointCheck(float x_cube, float y_cube, float z_cube, int element, int * inside_element)
{
  
  // volume of tetrahedron
  float volume_tetrahedron;
  // volume of the sub tetrahedron that can be formed with x_cube,y_cube,z_cube as one of its vertex
  float volume_sub_tetrahedron1, volume_sub_tetrahedron2, volume_sub_tetrahedron3, volume_sub_tetrahedron4;
  // sum of volumes of all 4 sub-tets
  float volume_sum_sub_tetrahedrons;

  // get the volume of the tetrahedron in the deformed image space
  this->Volume_DeformedTetrahedron(element,&volume_tetrahedron);

  // check for volume
  /*
  if (volume_tetrahedron < 0.0)
    {
    cout << "Backlash at element: " << element << endl;
    } // end if
  */

  // get the volume of the 4 sub tets 
  this->Volume_SubTetrahedron1(x_cube,y_cube,z_cube,&volume_sub_tetrahedron1,element);
  this->Volume_SubTetrahedron2(x_cube,y_cube,z_cube,&volume_sub_tetrahedron2,element);
  this->Volume_SubTetrahedron3(x_cube,y_cube,z_cube,&volume_sub_tetrahedron3,element);
  this->Volume_SubTetrahedron4(x_cube,y_cube,z_cube,&volume_sub_tetrahedron4,element);
 
  // check to see if point is inside the element
  if (volume_sub_tetrahedron1>0 && volume_sub_tetrahedron2>0 && volume_sub_tetrahedron3>0 && volume_sub_tetrahedron4>0) 
    {
    // sum up the volumes of all 4 sub tets
    volume_sum_sub_tetrahedrons = volume_sub_tetrahedron1 + volume_sub_tetrahedron2 + volume_sub_tetrahedron3 + volume_sub_tetrahedron4;
    if ( (volume_sum_sub_tetrahedrons-volume_tetrahedron)<this->TOLERANCE)
      {
      *inside_element = 1; // point is inside element
      } // end if
    } // end if
  
   if (volume_sub_tetrahedron1<0 || volume_sub_tetrahedron1>volume_tetrahedron)
    {
    *inside_element = 0;// point is outside element
    } // end if
   if (volume_sub_tetrahedron2<0 || volume_sub_tetrahedron2>volume_tetrahedron)
    {
    *inside_element = 0;// point is outside element
    } // end if
   if (volume_sub_tetrahedron3<0 || volume_sub_tetrahedron3>volume_tetrahedron)
    {
    *inside_element = 0;// point is outside element
    } // end if
   if (volume_sub_tetrahedron4<0 || volume_sub_tetrahedron4>volume_tetrahedron)
    {
    *inside_element = 0;// point is outside element
    } // end if

  return 0; // success
} // end pointcheck

// 4 sub-tetrahedrons can be formed with (x_cube,y_cube,z-cube) as one of the vertices of the tet
int ImageDeformation::Volume_SubTetrahedron1(float x_cube, float y_cube, float z_cube, float * volume_sub_tetrahedron1, int element)
{

  // Volume of a Tet : det(a)/6 where a:
  // a:      |x2-x_cube y2-y_cube z2-z_cube|
  //         |x3-x_cube y3-y_cube z3-z_cube|
  //         |x4-x_cube y4-y_cube z4-z_cube|
  float a11 = this->xdx[this->in[element][1]] - x_cube;
  float a12 = this->ydy[this->in[element][1]] - y_cube;
  float a13 = this->zdz[this->in[element][1]] - z_cube;

  float a21 = this->xdx[this->in[element][2]] - x_cube;
  float a22 = this->ydy[this->in[element][2]] - y_cube;
  float a23 = this->zdz[this->in[element][2]] - z_cube;

  float a31 = this->xdx[this->in[element][3]] - x_cube;
  float a32 = this->ydy[this->in[element][3]] - y_cube;
  float a33 = this->zdz[this->in[element][3]] - z_cube;

  float det = a11*a22*a33 - a11*a32*a23
            + a12*a31*a23 - a12*a21*a33
            + a13*a21*a32 - a13*a22*a31;

  *volume_sub_tetrahedron1 = det/6;
  
  return 0; // success
} // end volumesubtetrahedron1
  
int ImageDeformation::Volume_SubTetrahedron2(float x_cube, float y_cube, float z_cube, float * volume_sub_tetrahedron2, int element)
{
  
  // Volume of a Tet : det(a)/6 where a:
  // a:      |x2-x1     y2-y1     z2-z1    |
  //         |x_cube-x1 y_cube-y1 z_cube-z1|
  //         |x4-x1     y4-y1     z4-z1    |

  float a11 = x_cube - this->xdx[this->in[element][0]];
  float a12 = y_cube - this->ydy[this->in[element][0]];
  float a13 = z_cube - this->zdz[this->in[element][0]];

  float a21 = this->xdx[this->in[element][2]] - this->xdx[this->in[element][0]];
  float a22 = this->ydy[this->in[element][2]] - this->ydy[this->in[element][0]];
  float a23 = this->zdz[this->in[element][2]] - this->zdz[this->in[element][0]];

  float a31 = this->xdx[this->in[element][3]] - this->xdx[this->in[element][0]];
  float a32 = this->ydy[this->in[element][3]] - this->ydy[this->in[element][0]];
  float a33 = this->zdz[this->in[element][3]] - this->zdz[this->in[element][0]];

  float det = a11*a22*a33 - a11*a32*a23
            + a12*a31*a23 - a12*a21*a33
            + a13*a21*a32 - a13*a22*a31;

  *volume_sub_tetrahedron2 = det/6;

  return 0; // success
} // end volumesubtetrahedron2

int ImageDeformation::Volume_SubTetrahedron3(float x_cube, float y_cube, float z_cube, float * volume_sub_tetrahedron3, int element)
{
  
  // Volume of a Tet : det(a)/6 where a:
  // a:      |x_cube-x1 y_cube-y1 z_cube-z1|
  //         |x3-x1     y3-y1     z3-z1    |
  //         |x4-x1     y4-y1     z4-z1    |

  float a11 = this->xdx[this->in[element][1]] - this->xdx[this->in[element][0]];
  float a12 = this->ydy[this->in[element][1]] - this->ydy[this->in[element][0]];
  float a13 = this->zdz[this->in[element][1]] - this->zdz[this->in[element][0]];

  float a21 = x_cube - this->xdx[this->in[element][0]];
  float a22 = y_cube - this->ydy[this->in[element][0]];
  float a23 = z_cube - this->zdz[this->in[element][0]];

  float a31 = this->xdx[this->in[element][3]] - this->xdx[this->in[element][0]];
  float a32 = this->ydy[this->in[element][3]] - this->ydy[this->in[element][0]];
  float a33 = this->zdz[this->in[element][3]] - this->zdz[this->in[element][0]];

  float det = a11*a22*a33 - a11*a32*a23
            + a12*a31*a23 - a12*a21*a33
            + a13*a21*a32 - a13*a22*a31;
 
  *volume_sub_tetrahedron3 = det/6;
  
  return 0; // success
} // end volumesubtetrahedron3

int ImageDeformation::Volume_SubTetrahedron4(float x_cube, float y_cube, float z_cube, float * volume_sub_tetrahedron4, int element)
{
  
  // Volume of a Tet : det(a)/6 where a:
  // a:      |x2-x1     y2-y1     z2-z1    |
  //         |x3-x1     y3-y1     z3-z1    |
  //         |x_cube-x1 y_cube-y1 z_cube-z1|

  float a11 = this->xdx[this->in[element][1]] - this->xdx[this->in[element][0]];
  float a12 = this->ydy[this->in[element][1]] - this->ydy[this->in[element][0]];
  float a13 = this->zdz[this->in[element][1]] - this->zdz[this->in[element][0]];

  float a21 = this->xdx[this->in[element][2]] - this->xdx[this->in[element][0]];
  float a22 = this->ydy[this->in[element][2]] - this->ydy[this->in[element][0]];
  float a23 = this->zdz[this->in[element][2]] - this->zdz[this->in[element][0]];

  float a31 =  x_cube - this->xdx[this->in[element][0]];
  float a32 =  y_cube - this->ydy[this->in[element][0]];
  float a33 =  z_cube - this->zdz[this->in[element][0]];

  float det = a11*a22*a33 - a11*a32*a23
            + a12*a31*a23 - a12*a21*a33
            + a13*a21*a32 - a13*a22*a31;

  *volume_sub_tetrahedron4 = det/6;
  
  return 0; // success
} // end volumesubtetrahedron4

int ImageDeformation::GetBasisCoefficients(int element, float * a, float * b, float * c, float * d)
{

  int i, j, k, l;
  float A[3][3], B[3][3], C[3][3], D[3][3];
  float det_A, det_B, det_C, det_D;

  for (i=0; i<4; i++) 
    {
    if (i==0)      j=1, k=2, l=3;
    else if (i==1) j=2, k=3, l=0;
    else if (i==2) j=3, k=0, l=1;
    else if (i==3) j=0, k=1, l=2;
    
    //a(i) : det(A)
    //where A :   |x(j) y(j) z(j)|
    //            |x(k) y(k) z(k)|
    //            |x(l) y(l) z(l)|
    // for a(j),a(k),a(l) do a cyclic interchange of the constants and flip the sign
    A[0][0] = this->xdx[this->in[element][j]];
    A[0][1] = this->ydy[this->in[element][j]];
    A[0][2] = this->zdz[this->in[element][j]];

    A[1][0] = this->xdx[this->in[element][k]];
    A[1][1] = this->ydy[this->in[element][k]];
    A[1][2] = this->zdz[this->in[element][k]];

    A[2][0] = this->xdx[this->in[element][l]];
    A[2][1] = this->ydy[this->in[element][l]];
    A[2][2] = this->zdz[this->in[element][l]];

    det_A   = A[0][0]*A[1][1]*A[2][2] - A[0][0]*A[1][2]*A[2][1]
            + A[0][1]*A[1][2]*A[2][0] - A[0][1]*A[1][0]*A[2][2]
            + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[2][0]*A[1][1];

    a[i] = det_A * pow(-1.0,(double)i+2);

    //b(i) : det(B)
    //where B :   |1 y(j) z(j)|
    //            |1 y(k) z(k)|
    //            |1 y(l) z(l)|
    // for b(j),b(k),b(l) do a cyclic interchange of the constants and flip the sign
    B[0][0] = 1;
    B[0][1] = this->ydy[this->in[element][j]];
    B[0][2] = this->zdz[this->in[element][j]];

    B[1][0] = 1;
    B[1][1] = this->ydy[this->in[element][k]];
    B[1][2] = this->zdz[this->in[element][k]];

    B[2][0] = 1;
    B[2][1] = this->ydy[this->in[element][l]];
    B[2][2] = this->zdz[this->in[element][l]];

    det_B   = B[0][0]*B[1][1]*B[2][2] - B[0][0]*B[1][2]*B[2][1]
            + B[0][1]*B[1][2]*B[2][0] - B[0][1]*B[1][0]*B[2][2]
            + B[0][2]*B[1][0]*B[2][1] - B[0][2]*B[2][0]*B[1][1];

    b[i] = det_B * pow(-1.0,(double)i+1);

    //c(i) : det(C)
    //where C :   |x(j) 1 z(j)|
    //            |x(k) 1 z(k)|
    //            |x(l) 1 z(l)|
    // for c(j),c(k),c(l) do a cyclic interchange of the constants and flip the sign
    C[0][0] = this->xdx[this->in[element][j]];
    C[0][1] = 1;
    C[0][2] = this->zdz[this->in[element][j]];

    C[1][0] = this->xdx[this->in[element][k]];
    C[1][1] = 1;
    C[1][2] = this->zdz[this->in[element][k]];

    C[2][0] = this->xdx[this->in[element][l]];
    C[2][1] = 1;
    C[2][2] = this->zdz[this->in[element][l]];

    det_C   = C[0][0]*C[1][1]*C[2][2] - C[0][0]*C[1][2]*C[2][1]
            + C[0][1]*C[1][2]*C[2][0] - C[0][1]*C[1][0]*C[2][2]
            + C[0][2]*C[1][0]*C[2][1] - C[0][2]*C[2][0]*C[1][1];

    c[i] = det_C * pow(-1.0,(double)i+1);

    //d(i) : det(D)
    //where D :   |x(j) 1 z(j)|
    //            |x(k) 1 z(k)|
    //            |x(l) 1 z(l)|
    // for d(j),d(k),d(l) do a cyclic interdhange of the donstants and flip the sign
    D[0][0] = this->xdx[this->in[element][j]];
    D[0][1] = this->ydy[this->in[element][j]];
    D[0][2] = 1;

    D[1][0] = this->xdx[this->in[element][k]];
    D[1][1] = this->ydy[this->in[element][k]];
    D[1][2] = 1;

    D[2][0] = this->xdx[this->in[element][l]];
    D[2][1] = this->ydy[this->in[element][l]];
    D[2][2] = 1;

    det_D   = D[0][0]*D[1][1]*D[2][2] - D[0][0]*D[1][2]*D[2][1]
            + D[0][1]*D[1][2]*D[2][0] - D[0][1]*D[1][0]*D[2][2]
            + D[0][2]*D[1][0]*D[2][1] - D[0][2]*D[2][0]*D[1][1];

    d[i] = det_D * pow(-1.0,(double)i+1);

    } // end for i

  return 0; //success
} // end getbasiscoefficients

int ImageDeformation::Volume_DeformedTetrahedron(int element, float * volume)
{
  
  // Volume of a Tet : det(a)/6 where a:
  // a:      |x2-x1 y2-y1 z2-z1|
  //         |x3-x1 y3-y1 z3-z1|
  //         |x4-x1 y4-y1 z4-z1|

  float a11 = this->xdx[this->in[element][1]] - this->xdx[this->in[element][0]];
  float a12 = this->ydy[this->in[element][1]] - this->ydy[this->in[element][0]];
  float a13 = this->zdz[this->in[element][1]] - this->zdz[this->in[element][0]];

  float a21 = this->xdx[this->in[element][2]] - this->xdx[this->in[element][0]];
  float a22 = this->ydy[this->in[element][2]] - this->ydy[this->in[element][0]];
  float a23 = this->zdz[this->in[element][2]] - this->zdz[this->in[element][0]];

  float a31 = this->xdx[this->in[element][3]] - this->xdx[this->in[element][0]];
  float a32 = this->ydy[this->in[element][3]] - this->ydy[this->in[element][0]];
  float a33 = this->zdz[this->in[element][3]] - this->zdz[this->in[element][0]];

  float det = a11*a22*a33 - a11*a32*a23
            + a12*a31*a23 - a12*a21*a33
            + a13*a21*a32 - a13*a22*a31;

  *volume = det/6;

  /*
  if (*volume <=0.0) 
    {
    cout << "Backlash at element: " << element << endl;
    } // end if
  */

  return 0; // success
} // end volumedeformedtetrahedron

void ImageDeformation::PrintInfo()
{
  
  cout << "Mesh Info " << endl;
  cout << "Number Of Nodes: " << this->number_of_nodes << endl;
  cout << "Number Of Elements: " << this->number_of_elements << endl;
  cout << endl;
  cout << "Image Info " << endl;
  cout << "Image Dimensions: Xdim: " << this->image_dimension[0] << " Ydim: " << image_dimension[1] << " Zdim: " << this->image_dimension[2] << endl;
  
  if (this->data_type == 2)
    {
    cout << "Image Data Type: Unsigned Char" << endl; 
    } // end if
   if (this->data_type == 4)
    {
    cout << "Image Data Type: Signed Short" << endl; 
    } // end if
   if (this->data_type == 8)
    {
    cout << "Image Data Type: Signed Int" << endl; 
    } // end if
   if (this->data_type == 16)
    {
    cout << "Image Data Type: Float" << endl; 
    } // end if
   if (this->data_type == 64)
    {
    cout << "Image Data Type: Double" << endl; 
    } // end if

  cout << "Voxel Spacing in x y z: " << this->pixel_spacing[0] << " " << this->pixel_spacing[1] << " " << this->pixel_spacing[2] << endl;
} // end printinfo


int ImageDeformation::SaveDeformedImage(char * outfilename)
{
 int i;

 ofstream outfile;
 outfile.open(outfilename,ios::out | ios::binary);

 // save the deformed image based on the image data type
   for (i=0; i<this->image_size; i++)
     {
  //   outfile << (unsigned char)this->NewIntensity[i];
       outfile.write((char*) &this->NewIntensity[i], 2);
     } // end for i

 return 0; // success
} // end savedeformedimage













//int ImageDeformation::SaveDeformedImage(char * outfilename)
//{

//  int i;

//  ofstream outfile;
//  outfile.open(outfilename,ios::out | ios::binary);

//  // save the deformed image based on the image data type
//  if (this->data_type==2)
//    {
//    for (i=0; i<this->image_size; i++)
//      {
//      outfile << (unsigned char)this->NewIntensity[i];
//      } // end for i
//    } // end if

//  if (this->data_type==4)
//    {
//    cout<<"yes we are short"<<endl;
//      for (i=0; i<this->image_size; i++)
//	{
//	  outfile << (signed short)this->NewIntensity[i];
//	} // end for i
//    } // end if

//  if (this->data_type==8)
//    {
//      for (i=0; i<this->image_size; i++)
//	{
//	  outfile << (signed int)this->NewIntensity[i];
//	} // end for i
//    } // end if
//      
//  if (this->data_type==16)
//    {
//      for (i=0; i<this->image_size; i++)
//	{
//	  outfile << (float)this->NewIntensity[i];
//	} // end for i
//    } // end if

//  if (this->data_type==64)
//    {
//      for (i=0; i<this->image_size; i++)
//	{
//	  outfile << (double)this->NewIntensity[i];
//	} // end for i
//    } // end if

//  return 0; // success
//} // end savedeformedimage








//void ImageDeformation::WriteDeformedImage(char * outfilename)
//{
//  
//  int i;
//  

//  // data type = unsigned short
//  if (this->data_type==4)
//    {
//    vtkShortArray * tmp_new_int = vtkShortArray::New();
//    tmp_new_int->SetNumberOfTuples(this->GetImageWidth()*this->GetImageHeight()*this->GetNumberOfSlices());
//    for (i=0; i < tmp_new_int->GetNumberOfTuples(); i++)
//      {
//      tmp_new_int->InsertTuple1(i,(short)this->NewIntensity[i]);
//      } // end for 'i'
//    this->DeformedImage->SetScalarTypeToShort();
//    this->DeformedImage->GetPointData()->SetScalars(tmp_new_int);

//    // clear memory
//    tmp_new_int->Delete();
//    } // end if
// 
//  // write out the deformed image using vtkAnalyzeWriter
//  vtkAnalyzeWriter * writer = vtkAnalyzeWriter::New();
//  writer->SetFileName("mytestdeformed.avw");
//  writer->SetInput(this->DeformedImage);
//  //this->DeformedImage->Print(cout);
//  writer->Write();

//  // clear memory
//  writer->Delete();

//} // end WriteDeformedImage







// functions that can be called by the user to get info
int ImageDeformation::GetNumberOfNodes()
{
  return this->number_of_nodes; 
}
  
int ImageDeformation::GetNumberOfElements()
{ 
  return this->number_of_elements;
}

float ImageDeformation::GetPixelWidth()
{ 
  return this->pixel_spacing[0];
}

float ImageDeformation::GetPixelHeight()
{
  return this->pixel_spacing[1];
}
 
float ImageDeformation::GetPixelDepth()
{
  return this->pixel_spacing[2]; 
}

int ImageDeformation::GetImageWidth()
{
  return this->image_dimension[0]; 
}

int ImageDeformation::GetImageHeight()
{ 
  return this->image_dimension[1]; 
}

int ImageDeformation::GetNumberOfSlices()
{
  return this->image_dimension[2]; 
}

const char * ImageDeformation::GetImageDataType()
{ 
  if (this->data_type==2)
    {
      return "UNSIGNED_CHAR";
    } // end if

  if (this->data_type==4)
    { 
      return "SIGNED_SHORT";
    }
  if (this->data_type==8)
    {
      return "SIGNED_INT";
    } // end if

  if (this->data_type==16)
    { 
      return "FLOAT";
    }// end if

  if (this->data_type==64)
    { 
      return "DOUBLE";
    } // end if
}
