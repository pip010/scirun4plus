
MODEL_CONFIG_FILE="model_config.py"

PATH_SEG="/media/ppetrov/WORK/PhD/rats/CBWJ13_P80_indexed_volume/CBWJ13_P80_indexed_volume.nrrd"

PATH_BIO_MESH_CONFIG="/home/ppetrov/Documents/scirun4plus/bin/FEMesher/BuildModelConfig.py"


if [ ! -f MODEL_CONFIG_FILE ]; then 
	
	echo "*** Creating model_configpy file..."
	python $PATH_BIO_MESH_CONFIG $PATH_SEG $MODEL_CONFIG_FILE;
	echo "*** Finished in current directory!"
fi


