# ICCV_2011_4D_Face_Tracking
Demo code (Matlab) for ICCV 2011 paper:
A FACS Valid 3D Dynamic Action Unit Database With Applications to 3D Dynamic Morphable Facial Modeling. D Cosker, E Krumhuber, A Hilton  
IEEE International Conference on Computer Vision (ICCV), 2296-2303, 2011.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ICCV 2011 Paper Demo. Darren Cosker, University of Bath. 
% This demo (was originally) packaged together in 2012 (though not widely distributed).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This demo goes from raw unaligned OBJ and RGB image files to OBJ files with a regular mesh
structure tracked through a sequence.

I've included 3 steps of the pipeline, streamlined for demo purposes so they are readable, and run
easily. So this highlights the steps of the algorithm clearly.

Each step has a data folder. So you can run each step independently if you want. The final results
are also in a separate folder.

To run all scripts, go to the 'code' folder.
You will be asked to select the appropriate data directory at each stage.


Each script also contains an explanation, similar to below.
Comments should be clear enough to allow you to follow what's going on.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

STEP 1 - paint_cylindrical_UVs (run from command line)

First, some code is included to paint cylindrical UV maps from noisy stereo images, a set of UV coords for the new UV map, original stereo images, and camera calibration matrics for these images (from the cameras).

There is also a Batch file that calls a C++ .exe for smoothing the obj files (but you don't need to run that unless you want to smooth the objs further - this has already been done). The C++ is modified from something akin to a 3D Bilateral Filter.

Run the script 'paint_cylindrical_UVs' from the matlab command line, select the AU_4+7+17+23 directory, and it will create new UV maps AU_4+7+17+24_[IM_NUMBER]_w.bmp


STEP 2 - create_3D_images (run from command line)

This script creates '3D Images'. The 3D image format is pretty convenient. What it allows is for
non-rigid alignment to be calculated based on the RGB values of the UV
maps - say using optical flow, or AAM tracking (see ICCV paper), or other tracking - 
and then the alignment data to be applied to the 3D image. Applying in
this sense means alignment to a reference image - say a neutral
expression. Once this is done, a canonical vertex mapping - or a single
stable known facial vertex configuration - can be applied to the whole
sequence of OBJs. Originally, the OBJS have different vertex topologies
etc, so this fixes the problem. The quality of the topology and it's
stability is directly proportional to the quality of the non-rigid
alignment as applied to the RGB UV maps. So a better method can be
plugged in at a later date. This script however, just creates the 3D
images, which are used in non-rigid alignment in the next part of the
demo.


STEP 3 - fit_canonical_mesh_to_raw_OBJs (run from command line)
 

There are lots of ways to do the mesh fitting. Here, each UV map is
warped to a reference (frame 1) using tracked points (AAM) and a Thin
Plate Spline warp. You can also use optical flow, or any other tracker
you like, but here I am using TPS as it's stable and simple.
The warps are also applied to the 3D Images, and once they have been
aligned they are then sampled using the topology required for the overall
mesh. This can be artist based, but here Im just using a regular grid for
simplicity - so one vertex for every other pixel. This can be made one
vertex per pixel, or vertices wherever the artist likes.


Enjoy!

