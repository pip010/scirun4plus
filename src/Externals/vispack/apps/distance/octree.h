//----------------------------------------------------------------------//
// OctreeNode creates an adaptive octree to store a set of points,      //
//   which can then be queried for a closest point distance. This class //
//   takes pointers to objects of type T, but will take care of         //
//   releasing the memory for objects. T must instantiate the following //
//   functions:
//      vec<dim> position()
//----------------------------------------------------------------------//
#ifndef __OCTREE_H__
#define __OCTREE_H__

#include <mtxlib.h>

#define MAX_NUM_POINTS 20

class Point: public vec<3>
{
 public:
  const vec<3>& position() const { return(*((vec<3>*)this));};
  Point(const vec<3> v):vec<3>(v) {};
    Point():vec<3>() {};
      Point(const Point &p):vec<3>((vec<3>)p) {};
  const Point& operator=(const Point &p) {vec<3>::operator=(p); return(*this);}
  boolean valid  () const {return(operator[](0) != FLT_MAX);}
  void invalid () {operator[](0) = FLT_MAX;}
 }; 

// class Point
// {
//  public:
//   vec<3> _v;
//   const vec<3>& position() const { return(_v);}
//   Point(const vec<3> v){_v = v;}
//   Point() {}
//   Point(const Point &p) {_v = p._v;}
//   const Point& operator=(const Point &p) {_v = p._v;}
//   //  operator const vec<3>& () {return(_v);}
// };

//std::ostream& operator << (std::ostream& outs, 
//                           const Point& p) {outs << p._v;}

class NodeDistance
{
public:
  NodeDistance() {};
  ~NodeDistance() {};

  int _node_index;
  float _distance;
  bool operator < (const NodeDistance &rhs) const
  { return ( _distance < rhs._distance ); };
};

template <class T>
class OctreeNode
{
public:
  OctreeNode();
  ~OctreeNode();

  inline void setDimensions(const vec<3> &min, const vec<3> &max)
  { _min_coord = min; _max_coord = max; };
    inline float getMaxDimension(int i) const
  { return _max_coord(i); };
  inline float getMinDimension(int i) const
  { return _min_coord(i); };

  void addPoint(T* point);

  inline bool gotChildren() const
  { if ( _child ) return true; else return false; };
  inline bool gotPoints() const
  { if ( _num_points > 0 ) return true; else return false; };

  float distanceToClosestPoint(const T* point);
  T closestPoint(const T* point);

  //void render() const;

 private:
  OctreeNode<T> *_child;
  T* _point[MAX_NUM_POINTS];
  int _num_points;
  vec<3> _min_coord, _max_coord, _diff, _mid_coord;

  void createChildren();
  void distributePointsToChildren();
  int  determineChild(const T* point);
  float distanceToChildNode(int node, const T* point);

  inline float dist(const vec<3> &a, const vec<3> &b) const
  { return ( (a-b).length() ); };
};

#endif

#include <octree.txx>

