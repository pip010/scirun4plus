#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "mtxlib.h"
#include "Surface.h"
#include "IOScalarField.h"
#include "ScalarField.h"
#include "SurfaceParameters.h"

using namespace std;

struct vert
{
  float x, y, z;
  vector<int> connected_tets;
  vector<int> connected_faces;
};

struct tet
{
  int v1, v2, v3, v4;
  bool inside;
  float ratio;
  vec<3> circumcircle;
};

struct face
{
  int v1, v2, v3;
  bool t1_inside, t2_inside;
  const face& operator = (const face &that)
  {
    v1 = that.v1; v2 = that.v2; v3 = that.v3;
    t1_inside = that.t1_inside; t2_inside = that.t2_inside;
    return *this;
  }
  bool keep;
};

vector<vert> verts;

void radius_ratio( tet &t )
{
  vec<3> points[4];
  points[0].set(verts[t.v1].x, verts[t.v1].y, verts[t.v1].z);
  points[1].set(verts[t.v2].x, verts[t.v2].y, verts[t.v2].z);
  points[2].set(verts[t.v3].x, verts[t.v3].y, verts[t.v3].z);
  points[3].set(verts[t.v4].x, verts[t.v4].y, verts[t.v4].z);

  float
    a = (points[0] - points[3]).length(),
    b = (points[0] - points[2]).length(),
    c = (points[2] - points[3]).length(),
    d = (points[1] - points[2]).length(),
    e = (points[1] - points[3]).length(),
    f = (points[0] - points[1]).length();
  float s = (a*d + b*e + c*f) / 2.0;

  float V = DotProduct(
    CrossProduct((points[2]-points[0]),(points[1]-points[0])),
    (points[3]-points[0])) / 6.0;
  float S = 0.5*(
    (CrossProduct((points[1]-points[0]),(points[2]-points[0]))).length() +
    (CrossProduct((points[1]-points[0]),(points[3]-points[0]))).length() +
    (CrossProduct((points[2]-points[0]),(points[3]-points[0]))).length() +
    (CrossProduct((points[2]-points[1]),(points[3]-points[1]))).length() );

  V = fabs(V);
  S = fabs(S);

  float r_in  = 3.0 * V / (S+EPSILON);
  float r_cir = sqrt(s*(s-a*d)*(s-b*e)*(s-c*f)) / (6.0*V+EPSILON);

  float ideal_ratio = 0.33333;
  t.ratio = (r_in / r_cir) / ideal_ratio;
}

void centroid( tet &t )
{
  vec<3> a, b, c, d;
  a.set(verts[t.v1].x, verts[t.v1].y, verts[t.v1].z);
  b.set(verts[t.v2].x, verts[t.v2].y, verts[t.v2].z);
  c.set(verts[t.v3].x, verts[t.v3].y, verts[t.v3].z);
  d.set(verts[t.v4].x, verts[t.v4].y, verts[t.v4].z);

  t.circumcircle = (a + b + c +d)/4;     
}


void circumcirclecenter( tet &t )
{
  vec<3> a, b, c, d;
  a.set(verts[t.v1].x, verts[t.v1].y, verts[t.v1].z);
  b.set(verts[t.v2].x, verts[t.v2].y, verts[t.v2].z);
  c.set(verts[t.v3].x, verts[t.v3].y, verts[t.v3].z);
  d.set(verts[t.v4].x, verts[t.v4].y, verts[t.v4].z);

  vec<3> ba = b-a;
  vec<3> ca = c-a;
  vec<3> da = d-a;

  t.circumcircle = a +
    ( da.lengthSqr()*CrossProduct( ba, ca ) +
      ca.lengthSqr()*CrossProduct( da, ba ) +
      ba.lengthSqr()*CrossProduct( ca, da ) ) /
    ( 2.0*DotProduct( ba, CrossProduct(ca,da) ) ); 
}

double check_max_dist( tet &t, Surface* field )
{
  vec<3> a, b, c, d;
  a.set(verts[t.v1].x, verts[t.v1].y, verts[t.v1].z);
  b.set(verts[t.v2].x, verts[t.v2].y, verts[t.v2].z);
  c.set(verts[t.v3].x, verts[t.v3].y, verts[t.v3].z);
  d.set(verts[t.v4].x, verts[t.v4].y, verts[t.v4].z);

  SurfacePointParams params;
  double max_dist=0.0;


  field->computeSurfacePointParams( a, params,false, false );
   
  //get most negative point
  if ( params._F < max_dist )
  { 
    max_dist = params._F;
  }
   field->computeSurfacePointParams( b, params,false, false );
   
  //get most negative point
  if ( params._F < max_dist )
  { 
    max_dist = params._F;
  }
   field->computeSurfacePointParams( c, params,false, false );
   
  //get most negative point
  if ( params._F < max_dist )
  { 
    max_dist = params._F;
  }
   field->computeSurfacePointParams( d, params,false, false );
   
  //get most negative point
  if ( params._F < max_dist )
  { 
    max_dist = params._F;
  }
   
   return(max_dist);
}

int main(int argc, char **argv)
{
  if ( (argc != 4) && (argc != 5) )
  {
    cout << "Usage: <tet file basename> " <<
      "<volume filename> <output file name.m> <main IO function>" << endl;
    exit( 0 );
  }
  
  // open the tet file
  char filename[300]; sprintf( filename, "%s.ele", argv[1] );
  ifstream in1( filename );
  if ( !in1 )
  {
    cout << "Error reading " << filename << endl;
    exit( 1 );
  }

  // open the vertex file
  sprintf( filename, "%s.node", argv[1] );
  ifstream in2( filename );
  if ( !in2 )
  {
    cout << "Error reading " << filename << endl;
    exit( 1 );
  }

  // create the reconstructed surface
  Surface *field;
  if ( argc == 4 )
    field = new ScalarField( argv[2], 0.5 );
  else
    field = new IOScalarField ( argv[2], IOScalarField::INTERPOLATING,
                                atoi( argv[4] ) );
  
  // read in the verts
//  cout << "*** RemoteTets: Reading in the vertices..." << endl;
  string buffer;
  getline( in2, buffer );
 
  vert v;
  int i;
  while ( in2.peek() != EOF )
  {
    in2 >> buffer;
    i = atoi( buffer.c_str() ); in2 >> buffer;

    v.x = atof( buffer.c_str() ); in2 >> buffer;
    v.y = atof( buffer.c_str() ); in2 >> buffer;
    v.z = atof( buffer.c_str() ); getline( in2, buffer );

    verts.push_back( v );
  }

  in2.close();

  // read in the tets
//  cout << "*** RemoteTets: Reading in the tets..." << endl;
  getline( in1, buffer );
 
  vector<tet> tets;
  tet t;
  SurfacePointParams params, params2;
  vec<3> pos;  
   double toofar=false;
  while ( in1.peek() != EOF )
  {
    in1 >> buffer;
    i = atoi( buffer.c_str() ); in1 >> buffer;

    t.v1 = atoi( buffer.c_str() ); in1 >> buffer;
    t.v2 = atoi( buffer.c_str() ); in1 >> buffer;
    t.v3 = atoi( buffer.c_str() ); in1 >> buffer;
    t.v4 = atoi( buffer.c_str() ); getline( in1, buffer );
    
    circumcirclecenter( t );
    if ( field->inBounds(t.circumcircle) )
      field->computeSurfacePointParams( t.circumcircle, params,
                                        false, false );
    else
    {
      if ( argc == 4 )
        params._F = 1.0;
      else
        params._F = -1.0;
      //cout << "Out of bounds " << t.circumcircle << endl;
    }

    centroid( t );
    if ( field->inBounds(t.circumcircle) )
      field->computeSurfacePointParams( t.circumcircle, params2,
                                        false, false );
    else
    {
      if ( argc == 4 )
        params2._F = 1.0;
      else
        params2._F = -1.0;
      //cout << "Out of bounds " << t.circumcircle << endl;
    }



    // check if the centroid is positive or negative in the
    //   distance field
    if ( argc == 4 )
    {
      if ( (params._F > 0.0) /*|| (params2._F > 0.0)*/ )
      {
        t.inside = false;
      }
      else
      {
        toofar=check_max_dist(t,field);
        if (toofar > 0.01)
          t.inside = false;
        else        
          t.inside = true;      }
    }
    else
    {
      if ( (params._F < 0.0) /*|| (params2._F < 0.0)*/ )
      {
        t.inside = false;
      }
      else
      {
//        toofar=check_max_dist(t,field);
//        if (toofar< -0.01)
//          t.inside = false;
//        else        
          t.inside = true;
      }
    }
    tets.push_back( t );
  }

  in1.close();
  delete field;

  // create the faces
//  cout << "Creating the faces..." << endl;
  vector<face> faces;
  face f;
  for ( int i = 0; i < tets.size(); i++ )
  {
    // store the tet connections in the vertices
    verts[tets[i].v1].connected_tets.push_back( i );
    verts[tets[i].v2].connected_tets.push_back( i );
    verts[tets[i].v3].connected_tets.push_back( i );
    verts[tets[i].v4].connected_tets.push_back( i );

    // create the faces
    f.v1 = tets[i].v1; f.v2 = tets[i].v2; f.v3 = tets[i].v3;
    f.t1_inside = f.t2_inside = false;
    faces.push_back( f );

    f.v1 = tets[i].v1; f.v2 = tets[i].v2; f.v3 = tets[i].v4;
    f.t1_inside = f.t2_inside = false;
    faces.push_back( f );

    f.v1 = tets[i].v1; f.v2 = tets[i].v3; f.v3 = tets[i].v4;
    f.t1_inside = f.t2_inside = false;
    faces.push_back( f );

    f.v1 = tets[i].v2; f.v2 = tets[i].v3; f.v3 = tets[i].v4;
    f.t1_inside = f.t2_inside = false;
    faces.push_back( f );
  }

//  cout << "*** RemoteTets: Determining the tet connections..." << endl;
//  cout << "num faces = " << faces.size() << endl;
  for ( int i = 0; i < faces.size(); i++ )
  {
    f = faces[i];

    // go through the tets connected to the first vertex and
    //   find the two connected
    bool at_one = false;
    int ti;
    int tcounter = 0;
    for ( int j = 0; j < verts[f.v1].connected_tets.size(); j++ )
    {
      int num_matches = 0;     
      ti = verts[f.v1].connected_tets[j];
      if ( (tets[ti].v1 == f.v1) ||
           (tets[ti].v2 == f.v1) ||
           (tets[ti].v3 == f.v1) ||
           (tets[ti].v4 == f.v1) )
        ++num_matches;
      if ( (tets[ti].v1 == f.v2) ||
           (tets[ti].v2 == f.v2) ||
           (tets[ti].v3 == f.v2) ||
           (tets[ti].v4 == f.v2) )
        ++num_matches;
      if ( (tets[ti].v1 == f.v3) ||
           (tets[ti].v2 == f.v3) ||
           (tets[ti].v3 == f.v3) ||
           (tets[ti].v4 == f.v3) )
        ++num_matches;

      if ( num_matches == 3 )
      {
        // we have a match
        if ( !at_one )
        {
          f.t1_inside = tets[ti].inside;
          at_one = true;
        }
        else
         f.t2_inside = tets[ti].inside;

        ++tcounter;
      }   
    }

    if ( (tcounter != 1) && (tcounter != 2) )
        cout << tcounter << " " << verts[f.v1].connected_tets.size() << endl;

    // now, if one tet is on one side, and the other is on the
    //   other side, add this face to the list
    if ( f.t1_inside == f.t2_inside )
    {
      faces[i] = faces.back();
      faces.pop_back();
      --i;
    }
  }

  tets.clear();
  for ( int i = 0; i < verts.size(); i++ )
    verts[i].connected_tets.clear();

  // go through the faces and get rid of the duplicates
//  cout << "Finding duplicated faces ..." << endl;
//  cout << faces.size() << endl;

  // first, build the face connections with the vertices
  for ( int i = 0; i < faces.size(); i++ )
  {
    verts[faces[i].v1].connected_faces.push_back(i);
    verts[faces[i].v2].connected_faces.push_back(i);
    verts[faces[i].v3].connected_faces.push_back(i);

    faces[i].keep = true;
  }
  
  for ( int i = 0; i < faces.size(); i++ )
  {
    // so now, pick the first vertex for the face, and compare
    //   all the faces connected here
    for ( int j = 0; j < verts[faces[i].v1].connected_faces.size(); j++ )
    {
      int num_same = 0;
      int index = verts[faces[i].v1].connected_faces[j];
      if ( (index == i) || !faces[index].keep )
         continue;
      if ( ( faces[i].v1 == faces[index].v1) ||
           ( faces[i].v1 == faces[index].v2) ||
           ( faces[i].v1 == faces[index].v3) )
        ++num_same;
      if ( ( faces[i].v2 == faces[index].v1) ||
           ( faces[i].v2 == faces[index].v2) ||
           ( faces[i].v2 == faces[index].v3) )
        ++num_same;
      if ( ( faces[i].v3 == faces[index].v1) ||
           ( faces[i].v3 == faces[index].v2) ||
           ( faces[i].v3 == faces[index].v3) )
        ++num_same;
      if ( num_same == 3 )
      {
        faces[i].keep = false;
        j = verts[faces[i].v1].connected_faces.size();
      }
    }
  }
    
//     for ( int j = (i+1); j < faces.size(); j++ )
//     {
//       int num_same = 0;
//       if ( ( faces[i].v1 == faces[j].v1) ||
//            ( faces[i].v1 == faces[j].v2) ||
//            ( faces[i].v1 == faces[j].v3) )
//         ++num_same;
//       if ( ( faces[i].v2 == faces[j].v1) ||
//            ( faces[i].v2 == faces[j].v2) ||
//            ( faces[i].v2 == faces[j].v3) )
//         ++num_same;
//       if ( ( faces[i].v3 == faces[j].v1) ||
//            ( faces[i].v3 == faces[j].v2) ||
//            ( faces[i].v3 == faces[j].v3) )
//         ++num_same;
//       if ( num_same == 3 )
//       {
//         faces[j] = faces.back();
//         faces.pop_back();
//         j = faces.size();
//       }   
//     }

//     if ( !(i%1000) )
//       cout << i << endl;
//   }

  // write out the faces
//  cout << "Writing out the final mesh..." << endl;
  ofstream out( argv[3] );

  // check that it opened okay
  if ( !out )
  { 
    cout << "Error opening output file" << endl;
    exit ( 1 );
  }

    out << "# extracted isosurface" << endl;

//  cout << "   writing vertices...";
  for ( int i = 0; i < verts.size(); i++ )
    out << "Vertex " << (i+1) << "  " << verts[i].x << ' '
      << verts[i].y << ' ' << verts[i].z << endl;

//  cout << " writing faces..." << endl;
  for ( int i = 0; i < faces.size(); i++)
  {
    if ( !faces[i].keep )
      continue;
    out << "Face " << i+1 << "  " << (faces[i].v1+1) << ' ' <<
      (faces[i].v2+1) << ' ' << (faces[i].v3+1) << endl;
  }

//  cout << "Finished writing file" << endl;
  out.close();
  
  return 0;
}



 




