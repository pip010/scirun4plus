#include <cstdlib>
#include <iostream>
#include <SizingField.h>

using namespace std;

int main(int argc, char **argv)
{

  int arg_offset = 0;
  double constant_value = 0.0;
  if ( argc > 2 ) 
  {
	if ( strncmp( argv[1], "-c", strlen("-c") ) == 0 || 
		 strncmp( argv[1], "--constant_value", strlen("--constant_value") ) == 0 )
	{
		constant_value = atof( argv[2] );
		arg_offset = 2;
	}
  }

  if ( (argc != 4+arg_offset) && (argc != 6+arg_offset) )
  {
    cout << "Usage: <opt: -c/--constant_value sizing_field_value> <volume file> <base name> <redo volume name.vol> "
         << "<opt: epsilon delta (default = 0.5 0.4)>" << endl;
    exit( 0 );
  }

  SizingField *_sf;
  _sf = new SizingField( argv[2+arg_offset], argv[1+arg_offset], argv[3+arg_offset] );

  if ( argc == 6+arg_offset )
  {
    _sf->epsilon( atof(argv[4+arg_offset]) );
    _sf->delta( atof(argv[5+arg_offset]) );
  }

  if ( arg_offset ) 
  {
	_sf->generateConstantSizingField( constant_value );
  } else {
  	_sf->redoInitialSizingField();
  	_sf->generateSizingField();
  }


  delete _sf;
  
  return 0;
}



 




