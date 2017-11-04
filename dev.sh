echo 'Building, please wait...'
./build.sh -j4 --debug &> out.txt
echo 'Showing any errors : '
cat out.txt | grep error
