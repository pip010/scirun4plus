
SCIRUN4_PATH=$(pwd -P )

echo "SCIRun4 path : $SCIRUN4_PATH"


QT5_PATH="$HOME/Qt/5.4/gcc_64/lib/cmake" 

if ! [ -d "cleaver2plus" ]; then
	mkdir cleaver2plus
	echo "exec: git clone https://github.com/pip010/Cleaver2.git ."
    git clone https://github.com/pip010/Cleaver2.git cleaver2plus
fi

if [ -d "cleaver2plus" ]; then

	(cd cleaver2plus;

	if ! [ -d "build" ]; then
		mkdir build
	fi
		
		(cd build;
		
		if [ "$#" -gt 0 ]; 
			then cmake -DBUILD_CLEAVER_APP=ON  -DQt5Widgets_DIR="$QT5_PATH/Qt5Widgets" -DQt5OpenGL_DIR="$QT5_PATH/Qt5OpenGL" -DSCIRun4_DIR="$SCIRUN4_PATH" ../src && make
			else cmake -DBUILD_CLEAVER_APP=OFF -DSCIRun4_DIR="$SCIRUN4_PATH" ../src && make
		fi
		)
	)
fi
