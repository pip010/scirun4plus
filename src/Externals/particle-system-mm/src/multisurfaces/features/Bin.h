// FILE        : Bin.h                                                
// DESCRIPTION : The BinningStructure class is a spatial binning data
//               structure to increase the effiency of Points finding
//               neighboring Points.
//----------------------------------------------------------------------//

#ifndef __BIN_H__
#define __BIN_H__

#include <features/mtxlib.h>
#include <features/svector.h>

#ifdef _WIN32
#pragma warning (disable : 4244)
#endif

template <class T>
class Bin   
{ 
public:
  Bin() {};
  ~Bin() {};

  inline bool empty() const { return _elements.empty(); };
  inline void clear() { _elements.clear(); };
  inline int numElements() const { return _elements.size(); };

  void addElement(const T &element);
  void removeElement(const T &element);    
  inline const T& getElementAt(int index) const
  { return _elements.at(index); };
    
  void addNeighbor(int bin) { _neighbors.push_back(bin); };    
  inline int getNeighborAt(int index) const 
  { return _neighbors[index]; };
  inline int numNeighbors() const { return _neighbors.size(); };
    
private:
  custom_class::svector<int> _neighbors;
  custom_class::svector<T> _elements;
};

template <int dim, class T>
class BinningStructure  
{ 
public: 
  //--------------------------------------------------------------------
  // instantiate a BinningStructure by giving the width of a bin and
  //   the spatial domain
  BinningStructure(float bin_w, 
                    const vec<dim> &start, const vec<dim> &end);
  ~BinningStructure() { delete [] _bins; };

  //--------------------------------------------------------------------
  // send in a position, and determine which bin this point lies in
  int whichBin(const vec<dim> &pos) const;

  inline bool isBinEmpty(int bin) const 
  { return _bins[bin].empty(); };

  void clearBins();

  inline void addElementToBin(int bin, const T &element)
  { _bins[bin].addElement( element ); };
  inline void removeElementFromBin(int bin, const T &element)
  { _bins[bin].removeElement( element ); };

  inline int getBinsNumElements(int bin) const 
  { return _bins[bin].numElements(); };
  inline const T& getBinsElementAt(int bin, int index) const 
  { return _bins[bin].getElementAt(index); };

  inline int getBinsNumNeighbors(int bin) const 
  { return _bins[bin].numNeighbors(); };
  inline int getBinsNeighborAt(int bin, int index) 
  { return _bins[bin].getNeighborAt(index); };

  //--------------------------------------------------------------------
  // set all the neighboring bins around each bin -- this must be called
  //   for the bins' neighbor lists to be populated
  void buildBinNeighbors();

  inline void dimensions(int *dimension) const
  { for (int i=0; i<dim; i++) dimension[i] = _dimension[i]; };

  inline void limits(vector_type& start, vector_type& end) const
  { start = _start; end = _end; }

  void printNeighbors() const;
  int numBins() const; 
  
private:
  Bin<T> *_bins;
  int _dimension[dim];
  vector_type _start, _end;
  float _bin_width;
    
}; 

#include "Bin.T"

#endif // __BIN_H__
