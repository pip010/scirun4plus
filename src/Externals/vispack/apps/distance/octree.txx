//#include <windows.h>
//#include <GL/glut.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

//-------------------------------------------------------------------------
// function   : constructor and destructor
// description: 
//-------------------------------------------------------------------------
template <class T>
OctreeNode<T>::OctreeNode()
{
  _child = NULL;
  _num_points = 0;
}

template <class T> 
OctreeNode<T>::~OctreeNode()
{
  if ( gotChildren() )
    delete [] _child;

  else 
  {
    for ( int i = 0; i < _num_points; i++ )
      delete _point[i];
  }
}


//-------------------------------------------------------------------------
// function   : createChildren
// description: 
//-------------------------------------------------------------------------
template <class T>
void OctreeNode<T>::createChildren()
{
  // create the children 
  _child = new OctreeNode<T>[8];

  // set the dimensions of each child
  //
  // front slab: 3 2
  //             0 1     back slab: 7 6
  //                                4 5
  //
  vec<3> minc, maxc;
  _diff = 0.5*(_max_coord - _min_coord);
  _mid_coord = _min_coord + _diff;

  minc = _min_coord;
  maxc = _min_coord+_diff;
  _child[0].setDimensions( minc, maxc );

  minc(0) += _diff(0);
  maxc(0) += _diff(0);
  _child[1].setDimensions( minc, maxc );

  minc(1) += _diff(1);
  maxc(1) += _diff(1);
  _child[2].setDimensions( minc, maxc );

  minc(0) -= _diff(0);
  maxc(0) -= _diff(0);
  _child[3].setDimensions( minc, maxc );

  minc = _min_coord;
  maxc = _min_coord+_diff;

  minc(2) += _diff(2);
  maxc(2) += _diff(2);
  _child[4].setDimensions( minc, maxc );

  minc(0) += _diff(0);
  maxc(0) += _diff(0);
  _child[5].setDimensions( minc, maxc );

  minc(1) += _diff(1);
  maxc(1) += _diff(1);
  _child[6].setDimensions( minc, maxc );

  minc(0) -= _diff(0);
  maxc(0) -= _diff(0);
  _child[7].setDimensions( minc, maxc );

  // redistribute points to children
  for ( int i = 0; i < _num_points; i++ )
    _child[determineChild(_point[i])].addPoint(_point[i]);

  // clean up this node's points
  _num_points = 0;
}

//-------------------------------------------------------------------------
// function   : determineChild
// description: 
//-------------------------------------------------------------------------
template <class T>
int OctreeNode<T>::determineChild( const T* point )
{
  vec<3> pos( point->position() );

  if ( pos(0) < _mid_coord(0) )
    if ( pos(1) < _mid_coord(1) )
      if ( pos(2) < _mid_coord(2) )
        return 0;
      else
        return 4;
    else
      if ( pos(2) < _mid_coord(2) )
        return 3;
      else 
        return 7;
  else
    if ( pos(1) < _mid_coord(1) )
      if ( pos(2) < _mid_coord(2) )
        return 1;
      else
        return 5;
    else
      if ( pos(2) < _mid_coord(2) )
        return 2;
      else 
        return 6;
}

//-------------------------------------------------------------------------
// function   : addPoint
// description: 
//-------------------------------------------------------------------------
template <class T>
void OctreeNode<T>::addPoint( T* point )
{
  // if this is not a leaf node, send the point to the children
  if ( gotChildren() )
    _child[determineChild(point)].addPoint(point);

  else
  {
    // check if this leaf node's list of points is full
    if ( _num_points == MAX_NUM_POINTS )
    {
      createChildren();
      _child[determineChild(point)].addPoint(point);
    }
    // otherwise, add this point to the list
    else
    {
      _point[_num_points] = point;
      ++_num_points;
    }
  }
}

//-------------------------------------------------------------------------
// function   : distanceToClosestPoint
// description: 
//-------------------------------------------------------------------------
template <class T>
float OctreeNode<T>::distanceToClosestPoint( const T* point )
{
  float distance = FLT_MAX;

  // find the distance to the closest point in the list
  if ( !gotChildren() )
    for ( int i = 0; i < _num_points; i++ )
      distance = 
        min( dist(_point[i]->position(),point->position()), distance );

  // otherwise, we need to find the shortest distance in the children  
  else
  {
    // first, get distances to each child node
    vector<NodeDistance> nd( 8 );
    vector<NodeDistance>::iterator iter;  

    int i = 0;
    for ( iter = nd.begin(); iter != nd.end(); iter++ )
    {
      (*iter)._node_index = i;
      (*iter)._distance = distanceToChildNode( i, point );

      ++i;
    }

    // sort these distances
    sort( nd.begin(), nd.end() );

    // now, go into each child as long as the distance to the node
    //   is less than the shortest point distance
    iter = nd.begin();
    while ( (distance > (*iter)._distance) && (iter != nd.end()) )
    {
      distance = 
        min( _child[(*iter)._node_index].distanceToClosestPoint(point),
             distance );
      ++iter;
    }
  }
  return distance;
}


//-------------------------------------------------------------------------
// function   : distanceToClosestPoint
// description: 
//-------------------------------------------------------------------------
template <class T>
T OctreeNode<T>::closestPoint( const T *point)
{
  float this_distance, min_distance = FLT_MAX;
  T min_point, this_point;
  min_point.invalid();

  // find the distance to the closest point in the list
  if ( !gotChildren() )
    {
    for ( int i = 0; i < _num_points; i++ )
      {
	this_distance = dist(_point[i]->position(),point->position());
	if (this_distance < min_distance)
	  {
	    min_distance = this_distance; 
	    min_point = *(_point[i]);
	  }
      }
    //    cout << "got min point in children : " << min_point << endl;
    //    cout << "for point : " << point->position() << endl;
    }

  // otherwise, we need to find the shortest distance in the children  
  else
  {
    // first, get distances to each child node
    vector<NodeDistance> nd( 8 );
    vector<NodeDistance>::iterator iter;  

    int i = 0;
    for ( iter = nd.begin(); iter != nd.end(); iter++ )
    {
      (*iter)._node_index = i;
      (*iter)._distance = distanceToChildNode( i, point );
      ++i;
    }

    // sort these distances
    sort( nd.begin(), nd.end() );

    // now, go into each child as long as the distance to the node
    //   is less than the shortest point distance
    iter = nd.begin();
    while ( (min_distance > (*iter)._distance) && (iter != nd.end()) )
      {
	this_point = _child[(*iter)._node_index].closestPoint(point);
	this_distance = dist(this_point.position(),point->position());
	if ((this_distance < min_distance)&&this_point.valid())
	  {
	    min_distance = this_distance; 
	    min_point = this_point;
	  }
	++iter;
    }
  }
  return min_point;
}



//-------------------------------------------------------------------------
// function   : distanceToNode
// description: 
//-------------------------------------------------------------------------
template <class T>
float OctreeNode<T>::distanceToChildNode( int node, const T* point )
{
  float distance = 0.0;
  vec<3> pos = point->position();

  switch ( node )
  {
    case 0:
      distance += 
        power( max(pos(0)-_child[0].getMaxDimension(0), 0.0f), 2);
      distance += 
        power( max(pos(1)-_child[0].getMaxDimension(1), 0.0f), 2);
      distance += 
        power( max(pos(2)-_child[0].getMaxDimension(2), 0.0f), 2);
      break;

    case 1:
      distance += 
        power( max(_child[1].getMinDimension(0)-pos(0), 0.0f), 2);
      distance += 
        power( max(pos(1)-_child[1].getMaxDimension(1), 0.0f), 2);
      distance += 
        power( max(pos(2)-_child[1].getMaxDimension(2), 0.0f), 2);
      break;

    case 2:
      distance += 
        power( max(_child[2].getMinDimension(0)-pos(0), 0.0f), 2);
      distance += 
        power( max(_child[2].getMinDimension(1)-pos(1), 0.0f), 2);
      distance += 
        power( max(pos(2)-_child[2].getMaxDimension(2), 0.0f), 2);
      break;

    case 3:
      distance += 
        power( max(pos(0)-_child[3].getMaxDimension(0), 0.0f), 2);
      distance += 
        power( max(_child[3].getMinDimension(1)-pos(1), 0.0f), 2);
      distance += 
        power( max(pos(2)-_child[3].getMaxDimension(2), 0.0f), 2);
      break;

    case 4:
      distance += 
        power( max(pos(0)-_child[4].getMaxDimension(0), 0.0f), 2);
      distance += 
        power( max(pos(1)-_child[4].getMaxDimension(1), 0.0f), 2);
      distance += 
        power( max(_child[4].getMinDimension(2)-pos(2), 0.0f), 2);
      break;

    case 5:
      distance += 
        power( max(_child[5].getMinDimension(0)-pos(0), 0.0f), 2);
      distance += 
        power( max(pos(1)-_child[5].getMaxDimension(1), 0.0f), 2);
      distance += 
        power( max(_child[5].getMinDimension(2)-pos(2), 0.0f), 2);
      break;

    case 6:
      distance += 
        power( max(_child[6].getMinDimension(0)-pos(0), 0.0f), 2);
      distance += 
        power( max(_child[6].getMinDimension(1)-pos(1), 0.0f), 2);
      distance += 
        power( max(_child[6].getMinDimension(2)-pos(2), 0.0f), 2);
      break;

    case 7:
      distance += 
        power( max(pos(0)-_child[7].getMaxDimension(0), 0.0f), 2);
      distance += 
        power( max(_child[7].getMinDimension(1)-pos(1), 0.0f), 2);
      distance += 
        power( max(_child[7].getMinDimension(2)-pos(2), 0.0f), 2);
      break;
  }

  return (sqrt( distance ));
}

//-------------------------------------------------------------------------
// function   : render
// description: 
//-------------------------------------------------------------------------
/*template <class T>
void OctreeNode<T>::render() const
{
  if ( gotChildren() )
    for ( int i = 0; i < 8; i++ )
      _child[i].render();

  else
  {
      glDisable( GL_LIGHTING );
      glColor3f( 0.0, 0.0, 0.0 );
      glLineWidth( 1.0 );
      glBegin( GL_LINE_STRIP );
      {
        glVertex3f( _min_coord(0), _min_coord(1), _min_coord(2) );
        glVertex3f( _max_coord(0), _min_coord(1), _min_coord(2) );
        glVertex3f( _max_coord(0), _max_coord(1), _min_coord(2) );
        glVertex3f( _min_coord(0), _max_coord(1), _min_coord(2) );
        glVertex3f( _min_coord(0), _min_coord(1), _min_coord(2) );
      } glEnd();
      glBegin( GL_LINE_LOOP );
      {
        glVertex3f( _min_coord(0), _min_coord(1), _max_coord(2) );
        glVertex3f( _max_coord(0), _min_coord(1), _max_coord(2) );
        glVertex3f( _max_coord(0), _max_coord(1), _max_coord(2) );
        glVertex3f( _min_coord(0), _max_coord(1), _max_coord(2) );
      } glEnd();
      glBegin( GL_LINES );
      {
        glVertex3f( _min_coord(0), _min_coord(1), _min_coord(2) );
        glVertex3f( _min_coord(0), _min_coord(1), _max_coord(2) );

        glVertex3f( _max_coord(0), _min_coord(1), _min_coord(2) );
        glVertex3f( _max_coord(0), _min_coord(1), _max_coord(2) );

        glVertex3f( _max_coord(0), _max_coord(1), _min_coord(2) );
        glVertex3f( _max_coord(0), _max_coord(1), _max_coord(2) );

        glVertex3f( _min_coord(0), _max_coord(1), _min_coord(2) );
        glVertex3f( _min_coord(0), _max_coord(1), _max_coord(2) );
      }glEnd();
  }
}*/
