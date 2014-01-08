//----------------------------------------------------------------------//
// FILE        : LFS.h                                                
// DESCRIPTION : Reads in a set of medial axis points and build an
//               octree of the points. Has a querying function for
//               computing the lfs at a position.
//----------------------------------------------------------------------//

#ifndef __LFS_H__
#define __LFS_H__

#include <mtxlib.h>
#include <defines.h>
#include <svector.h>
#include <octree.h>

class LocalFeatureSize
{
public:
  LocalFeatureSize(const char *filename, const vec<3> &s,
                   const vec<3> &e,
                   int xdim, int ydim, int zdim=1);
  ~LocalFeatureSize();

  float lfs(const vector_type &pos) const;

private:
  custom_class::svector<Point*> _medial_axis_pts;
  OctreeNode<Point> *_octree;

  int _xdim, _ydim, _zdim;
  vec<3> _start, _end;

  void readInMAPoints(const char *filename);

  void addMedialAxisPoint(const vector_type &pos,
                          const vector_type &n);

};


#endif // __LFS_H__
