#include "matrix.h"

HMatrix::HMatrix(float r,float phi,float theta,
		 float x,float y,float z){
    //r=angle of rotation around the x-axis in deg.
    //phi=angle of rotation around the y-axis in deg.
    //theta=angle of rotation around the z-axis in deg.
    //x,y,z is a translation
    initialize(4,4);
    poke(0,0)=fcos(phi)*fcos(theta);
    poke(0,1)=fsin(r)*fsin(phi)*fcos(theta)+fcos(r)*fsin(theta);
    poke(0,2)=-1.0f*fcos(r)*fsin(phi)*fcos(theta)+fsin(r)*fsin(theta);
    poke(0,3)=x;
    poke(1,0)=-1.0f*fcos(phi)*fsin(theta);
    poke(1,1)=-1.0f*fsin(r)*fsin(phi)*fsin(theta)+fcos(r)*fcos(theta);
    poke(1,2)=fcos(r)*fsin(phi)*fsin(theta)+fsin(r)*fcos(theta);
    poke(2,0)=fsin(phi);
    poke(2,1)=-1.0f*fsin(r)*fcos(phi);
    poke(2,2)=fcos(r)*fsin(phi);
}

HMatrix::HMatrix(const Matrix& other){
    // here we could check for zeros on last row and col
    // and a 1 in the bottom corner.
    if((other.c()==4) && (other.r()==4)){
	Matrix(other);
    }
    else{
	cout << "HMatrix(const Matrix& other)" << endl;
	cout << "You have attempted to create an HMatrix from "
	     << "a matrix that is not 4x4" << endl;
	exit(0);
    }
}
	
HMatrix::HMatrix(char* filename){
    // here we could check for zeros on last row and col
    // and a 1 in the bottom corner.
    ifstream file(filename);
    int r,c;
    file >> r >> c;
    if ((r==4) && (c==4)){
	initialize(r,c);
	int i, j;
	for (i = 0; i < _rows; i++) {
	    for (j = 0; j < _cols; j++) {
		if (!file)
		    cerr << "error reading matrix: ran out of data";
		file >> poke(i,j);
	    }
	}
    }
    else{
	cout << "HMatrix(char* filename)" << endl;
	cout << "You have attempted to create an HMatrix from "
	     << "a matrix that is not 4x4" << endl;
	exit(0);
    }
}

HMatrix& HMatrix::operator=(const Matrix& other){
    if((other.c()==4) && (other.r()==4)){
	Matrix::operator=(other);
	return *this;
    }
    else{
	cout << "HMatrix operator=(const Matrix& other)" << endl;
	cout << "You have attempted to create an HMatrix from "
	     << "a matrix that is not 4x4" << endl;
	exit(0);
    }
}
	    
Vector HMatrix::operator%(const Vector& other) const{
    // This is an operator that rotates a 3D point.  It returns
    // a vector of size three.  This operation is equvalent
    // to doing r*pt where pt is a vector of size 4 where
    // the fourth element of the vector is zero.
    if (other.n() == 3){
	Vector ans(3);
	ans.poke(0) = peek(0,0)*other.peek(0) +
	              peek(0,1)*other.peek(1) +
	              peek(0,2)*other.peek(2);
	ans.poke(1) = peek(1,0)*other.peek(0) +
	              peek(1,1)*other.peek(1) +
	              peek(1,2)*other.peek(2);
	ans.poke(2) = peek(2,0)*other.peek(0) +
	              peek(2,1)*other.peek(1) +
	              peek(2,2)*other.peek(2);
	return ans;
    }
    else{
	cout << "Vector HMatrix::operator*(const Vector& other)" << endl;
	cout << "you have attempted to rotate/translate a 3D point that "
	     << "is not 3x1" << endl;
	exit(0);
    }	    
}

/*Vector HMatrix::operator*(const Vector& other) const{
    // This is an operator that rotates and translates a 3D point.
    // It returns a vector of size three.  This operation is equvalent
    // to doing r*pt where pt is a vector of size 4 where
    // the fourth element of the vector is one.
    if (other.n() == 3){
	Vector ans(3);
	ans.poke(0) = peek(0,0)*other.peek(0) +
	              peek(0,1)*other.peek(1) +
	              peek(0,2)*other.peek(2) +
	              peek(0,3);				
	ans.poke(1) = peek(1,0)*other.peek(0) +
	              peek(1,1)*other.peek(1) +
	              peek(1,2)*other.peek(2) +
	              peek(1,3);
	ans.poke(2) = peek(2,0)*other.peek(0) +
	              peek(2,1)*other.peek(1) +
	              peek(2,2)*other.peek(2) +
	              peek(2,3);
	return ans;
    }
    else{
	cout << "Vector HMatrix::operator*(const Vector& other)" << endl;
	cout << "you have attempted to rotate/translate a 3D point that"
	     << "is not 3x1" << endl;
	exit(0);
    }	    
}*/
