file configure_unix.config: 

it is an example of configuration file for taurus in RELEASE mode

--------------------------------------------------------------------------

file ViewReducedMesh.py:

it is a python script that allows to visualize the reduced mesh used 
in the online phase by the hyper reduction algorithm. How to use it?
First, make sure that after launching the **_GenerateROM** executable 
you have the files ReducedMesh.h5 and ReducedMesh.xmf. Open Paraview,
an then open the file ReducedMesh.xmf. Finally, click 
Tools->Python Shell->Run Script and select ViewReducedMesh.py. 
Triangles and Volumes will have different colors. 

