=========================================================================================

To use the BioMesh3D application on the mac there are two ways of
invoking the application:


Basic version:

The easiest way to operate BioMesh3D on the mac is to drag BioMesh3D application to the
Applications folder and to make a short cut in the dock by dragging and dropping the 
application into the dock. 

The next step is to drag and drop the python configuration files on top of the BioMesh3D
application. This will start the process and will launch SCIRun for the interactive parts
of the pipeline.

Examples of BioMesh3D files are included in the DiskImage. Open the examples and drop the
configuration file from the tooth or the mickey onto the BioMesh3D app and it will launch
the process and produce a mesh in the specified output directory.


Advanced version:

The BioMesh3D.app application folder actually contains all the scripts to run each part of
the pipeline. Drag the BioMesh3D application to the location where you want to install the
application. Then open a Terminal and cd into the BioMesh3D.app folder. Inside the application
folder under Contents/Resources/bin/FEMesher all the BioMesh3D scripts are located.

Once you opened the Terminal and changed the current path to the Contents/Resources/bin/FEMesher
folder inside the application, you can use the python command to run the python scripts.
The main python script that runs the pipeline is the BuildMesh.py script.

This way of accessing the application has the advantage that one can specify which parts of the
pipeline to run. Besides that the python scripts can be modified to fulfill specific needs.


For more information on the scripting of the BioMesh3D pipeline, please see the 
BioMesh3D documentation.

=========================================================================================

