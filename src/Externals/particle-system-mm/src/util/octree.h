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

#include <vector>
#include <mtxlib.h>

#define FLT_MAX 3.402823466e+38F

#define MAX_NUM_POINTS 20


class Point
{
public:
  Point() {_got_normal=false;};
  ~Point() {};
  vec<3> x;
  vec<3> n;
  float flux;
  bool _got_normal;
  int num_m;
  int m[4];

  inline bool got_normal() const {return _got_normal;}
  inline void got_normal(bool gn) {_got_normal=gn;}
  inline const vec<3>& position() const {return x;};
  inline bool valid() const {return(x[0] != FLT_MAX);}
  inline void invalid() {x[0] = FLT_MAX;}
};

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

  void pointsWithinDistance(const T* point, double distance,
                            std::vector<T> &out_points);

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

#include <octree.txx>

#endif



