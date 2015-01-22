# ImageDeform
Image deformation driven by displacements calculated from a finite element solve


Image Deformation with VTK
This program reads in a MRI/CT image volume and deforms it according to the displacements specified. 
The volume must be first meshed to produce node and element files to create a finite element model. The displacements are calculated from a registration process and are used to deform the images such that they produce an updated image volume. 

Usage:
Compile and run the TestImageDeformation executable. Argument list format will be displayed if inputs are incorrect. The accepted input is a vti file obtained from running the raw to vti converter on the original mhd/raw tomogram. 
