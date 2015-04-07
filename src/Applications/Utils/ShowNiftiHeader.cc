//Author: p.petrov 2015 

#include <iostream>
#include <fstream>

#include "nifti1.h"

#define nifti_header_size 348

union nifti_header
{
	nifti_1_header fields;
	char bytes[nifti_header_size];
};

using namespace std;

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		cout << "Please specify a valid nifty file as first argument!\n";
	}

	string filename(argv[1]);

	if(filename.size() == 0)
	{
		cout << "Please specify a valid nifty file as first argument!\n";
	}

	//! read the nifti header from a file
	nifti_header header;
	{
		ifstream nifti_file;

		cout << "trying to open file: " << filename.c_str() ;

		nifti_file.open( filename.c_str() , ios::binary );

		if( nifti_file.is_open() )
		{
			nifti_file.read( header.bytes ,nifti_header_size );
			cout << "... DONE!\n";
		}
		else
		{
			cout << "Error opening file, make sure it is not open by another application.\n";
		}
	}

	//! print the header info to standard output
	cout << header.fields.srow_x[0] << " " << header.fields.srow_x[1] << " " << header.fields.srow_x[2] << " " << header.fields.srow_x[3] << endl;
	cout << header.fields.srow_y[0] << " " << header.fields.srow_y[1] << " " << header.fields.srow_y[2] << " " << header.fields.srow_y[3] << endl;
	cout << header.fields.srow_z[0] << " " << header.fields.srow_z[1] << " " << header.fields.srow_z[2] << " " << header.fields.srow_z[3] << endl;
	cout << "0 0 0 1" << endl; // last row of the 4x4 always assumed to be this one

}