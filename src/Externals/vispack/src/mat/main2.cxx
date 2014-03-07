/*****************************************************************************
 *                              EE573 Project2                               *
 * Module : main2.C                                                          *
 *                                                                           *
 *****************************************************************************/

#include <math.h>
#include "matrix.h"

Vector NonlinearCalibration(Matrix&, Matrix&, Vector&, float);
Matrix MakeMatrixG(Matrix&, Vector&, Vector&);

main (int argc, char **argv)
{
    Matrix m5D, m3D, m2D, nm;
    Vector Vbo, Vb;

    m5D.readfile( "gridpts.txt" ); //matrix NX5
    m3D = m5D.subColMatrix(0, 2);  //3D positions Nx3
    m2D = m5D.subColMatrix(3, 4);  //2D coordinates Nx2
    Vbo.readfile( "beta0.txt" );

    float stddev = 1.0f;
    nm = stddev * noiseMatrix(m2D.r(), m2D.c());   
    nm = m2D + nm;  //noisy matrix NX2

    Vb = NonlinearCalibration(m3D, nm, Vbo, 1.0f);    
}

//NonLinear camera calibration, p3D:NX3 p2D:NX2 Vbo:10X1
Vector NonlinearCalibration(Matrix& p3D, Matrix& p2D, Vector& Vbo, float ftiny)
{
    int N = p2D.r();
    Vector Alpha(2*N), Vgt(2*N);
    Vector Beta, delta;
    Matrix mG, tmpG;
    float eps;

    for(int i=0; i < N; i++){
	Alpha.setValue(2*i, p2D[i][0]);   //a = [...ri ci...]'
	Alpha.setValue(2*i+1, p2D[i][1]);
    }
//    cout << "Alpha: " << Alpha << endl;
    Beta = Vbo;
//    cout << "Beta0: " << Beta << endl;
    int count = 1;
    for(; ;){
	mG = MakeMatrixG(p3D, Beta, Vgt); //matrix G
	tmpG = (mG.T() * mG).inv(); //(G'*G)^-1
	delta = (tmpG * mG.T())*(Alpha - Vgt); //deltaB
	delta = delta * 0.5f;
	Beta = Beta + delta;  //B(t+1)=B(t) + dB
	eps = (Alpha-Vgt).norm() / N; //||a - gt||/N
	//  eps = delta / Beta;
	cout << "Iterations: " << count << "  ";
	cout << "eps = " << eps << endl << endl;
	if(eps < ftiny){  
	    cout << "Beta: "<< Beta << endl;
	    break;
	}
	count++;
    }
    return Beta;
}

Matrix MatrixR(Vector& wpk) //matrix R at (w phi k)
//tan(w)=(-R32)/R33; sin(phi)=R31; tan(k)=(-R21)/R11
{
    float w, p, k; //angles w phi k
    w = wpk[0];    p = wpk[1];    k = wpk[2];
    Matrix R(3, 3);
    R[0][0] = cos(p) * cos(k);
    R[0][1] = sin(w) * sin(p) * cos(k) + cos(w) * sin(k);
    R[0][2] = - cos(w) * sin(p) * cos(k) + sin(w) * sin(k);
    R[1][0] = - cos(p) * sin(k);
    R[1][1] = - sin(w) * sin(p) * sin(k) + cos(w) * cos(k);
    R[1][2] = cos(w) * sin(p) * sin(k) + sin(w) * cos(k);
    R[2][0] = sin(p);
    R[2][1] = - sin(w) * cos(p);
    R[2][2] = cos(w) * cos(p);
//    cout << "determinant | R | = " << R.det() << endl;
//    cout << "rotation matrix R : " << R << endl;
    return R;
}

Matrix MdRdw(Vector& wpk) //matrix dR/dw at (w phi k)
{
    float w, p, k; //angles w phi k
    w = wpk[0];    p = wpk[1];    k = wpk[2];
    Matrix R(3, 3);
    R[0][1] = cos(w)*sin(p)*cos(k) - sin(w)*sin(k);
    R[0][2] = sin(w)*sin(p)*cos(k) + cos(w)*sin(k);
    R[1][1] = - cos(w)*sin(p)*sin(k) - sin(w)*cos(k);
    R[1][2] = cos(w)*cos(k) - sin(w) *sin(p)*sin(k);
    R[2][1] = - cos(w)*cos(p);
    R[2][2] = - sin(w)*cos(p);
    return R;
}

Matrix MdRdp(Vector& wpk) //matrix dR/dphi at (w phi k)
{
    float w, p, k; //angles w phi k
    w = wpk[0];    p = wpk[1];    k = wpk[2];
    Matrix R(3, 3);
    R[0][0] = - sin(p)*cos(k);
    R[0][1] = sin(w)*cos(p)*cos(k);
    R[0][2] = - cos(w)*cos(p)*cos(k);
    R[1][0] = sin(p)*sin(k);
    R[1][1] = - sin(w)*cos(p)*sin(k);
    R[1][2] = cos(w)*cos(p)*sin(k);
    R[2][0] = cos(p);
    R[2][1] = sin(w)*sin(p);
    R[2][2] = - cos(w)*sin(p);
    return R;
}

Matrix MdRdk(Vector& wpk) //matrix dR/dk at (w phi k)
{
    float w, p, k; //angles w phi k
    w = wpk[0];    p = wpk[1];    k = wpk[2];
    Matrix R(3, 3);
    R[0][0] = - cos(p) * sin(k);
    R[0][1] = - sin(w)*sin(p)*sin(k) + cos(w)*cos(k);
    R[0][2] = cos(w)*sin(p)*sin(k) + sin(w)*cos(k);
    R[1][0] = - cos(p) * cos(k);
    R[1][1] = - sin(w)*sin(p)*cos(k) - cos(w)*sin(k);
    R[1][2] = cos(w)*sin(p)*cos(k) - sin(w)*sin(k);
    return R;
}

//Gt=Gi(N):2NX10, Gi=[A*B C]:2X10, A:2x3  B:3x6  C:2x4
Matrix MakeMatrixG(Matrix& p3D, Vector& Vbt, Vector& Vgt)
{
    Vector Vooo = Vbt.getVector(0, 2); //xo yo zo
    Vector Vwpk = Vbt.getVector(3, 5); //w phi k
    Vector Vrc = Vbt.getVector(6, 7);  //ro co
    Vector Vf = Vbt.getVector(8, 9);   //fu fv

    Matrix A(2, 3), B(3, 6), dR(3, 3), C(2, 4);
    Matrix Gi(2, 10), Gt(0, 10);
    Vector Vxyz(3), Vpqs(3);
    for(int i=0; i < p3D.r(); i++){
	Vxyz = p3D.getRow(i).T(); //xn yn zn
	Vpqs = MatrixR(Vwpk) * (Vxyz - Vooo); //pqs = R*(xyzn -xyzo)
	//makes vector gt: to be compared with vector a
	Vgt.setValue(2*i, Vf[0]*Vpqs[0]/Vpqs[2]+Vrc[0]);  //fu*pn/sn +ro
	Vgt.setValue(2*i+1, Vf[1]*Vpqs[1]/Vpqs[2]+Vrc[1]);//fv*qn/sn +co
      
	//make matrix A: 2X3
	A[0][0] = Vf[0] / Vpqs[2]; //fu/sn
	A[1][1] = Vf[1] / Vpqs[2]; //fv/sn
	A[0][2] = - Vf[0]*Vpqs[0] / (Vpqs[2]*Vpqs[2]); //-fu*pn/sn*sn
	A[1][2] = - Vf[1]*Vpqs[1] / (Vpqs[2]*Vpqs[2]); //-fv*qn/sn*sn
	//make matrix B: 3X6
	dR.setCol(0, MdRdw(Vwpk)*(Vxyz - Vooo));//dR=MdRdX()*(xyzn-xyzo)
	dR.setCol(1, MdRdp(Vwpk)*(Vxyz - Vooo));
	dR.setCol(2, MdRdk(Vwpk)*(Vxyz - Vooo));
	B = (-1.0f * MatrixR(Vwpk)).colAddMatrix(dR); //B=[-R dR]
	//make matrix C: 2X4
	C[0][0] = 1.0f;  // 1 0
	C[1][1] = 1.0f;  // 0 1
	C[0][2] = Vpqs[0] / Vpqs[2]; // pn/sn 0
	C[1][3] = Vpqs[1] / Vpqs[2]; // 0 qn/sn

	//make matrix Gi: 2X10, Gt: 2NX10
	Gi = (A*B).colAddMatrix(C); //"col-add": Gi = [A*B C] 
	Gt = Gt.rowAddMatrix(Gi);   //"row-add": [G1..Gi..GN]'
    }
    return Gt;
}
