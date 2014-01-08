#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

struct vert
{
  float x,y,z;
};

int main(int argc, char **argv)
{
  if ( argc != 3 )
  {
    cout << "Usage: <input file name.pts> " <<
      "<output file name.node>" << endl;
    exit( 0 );
  }
  
  // open the node file
  ifstream in( argv[1] );
  if ( !in )
  {
    cout << "Error reading " << argv[1] << endl;
    exit( 1 );
  }

  // temp variables
  string buffer;

   // read the vertices
  vector<vert> verts;
  vert v;
  while ( in.peek() != EOF )
  {
    in >> buffer;

    v.x = atof( buffer.c_str() ); in >> buffer;
    v.y = atof( buffer.c_str() ); in >> buffer;
    v.z = atof( buffer.c_str() ); getline( in, buffer );

    verts.push_back( v );
  }

  in.close();

  // open the file for writting
  ofstream out( argv[2] );

  // check that it opened okay
  if ( !out )
  { 
    cout << "Error opening output file" << endl;
    exit ( 1 );
  }

  // write out the number of verts
  int num_verts = verts.size();
  out << num_verts << " 3 0 0" << endl;

  // read the vertices and write out
  for ( int i = 0; i < num_verts; i++ )
    out << i << " " << verts[i].x << " " << verts[i].y << " " <<
      verts[i].z << endl;

  out.close();

  return 0;
}



 




