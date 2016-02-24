
SCIRUN4_PATH=$(pwd -P )

echo "SCIRun4 path:"
echo '$SCIRUN4_PATH'

if ! [ -d "cleaver2plus" ]; then
	mkdir cleaver2plus
	echo "exec: git clone https://github.com/pip010/Cleaver2.git ."
    git clone https://github.com/pip010/Cleaver2.git cleaver2plus
fi

if [ -d "cleaver2plus" ]; then

	cd cleaver2plus

	if ! [ -d "build" ]; then
		mkdir build
	fi

	cd build

	cmake -DBUILD_CLEAVER_APP=OFF -DSCIRun4_DIR='$SCIRUN4_PATH' ../src

	make

	cd ..
fi