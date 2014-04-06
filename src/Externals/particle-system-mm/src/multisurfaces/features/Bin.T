#include <cstdlib>
#include <iostream>

// 4018 -- '<' signed/unsigned mismatch

#ifdef WIN_32
#pragma warning( disable : 4018 )
#endif

//using namespace std;


template <class T>
Bin<T>::Bin() 
{ 
_neighbors.reserve(30); 
}
 
template <class T>
Bin<T>::~Bin() 
{
	//cout << "_neighbors count : " <<  _neighbors.size() << "   _elements count : " <<    _elements.size() << std::endl;
}


//------------------------------------------------------------------------
// Function    : Bin::addElement
// Description : add an element to the list
//------------------------------------------------------------------------
template <class T>
void Bin<T>::addElement( const T &element )
{
  _elements.push_back( element );
}

//------------------------------------------------------------------------
// Function    : Bin::removeElement
// Description : search the bin's element list for the element and
//               remove it
//------------------------------------------------------------------------
template <class T>
void Bin<T>::removeElement( const T &element )
{
  for ( unsigned int i = 0; i < _elements.size(); i++ )
    if ( _elements[i] == element )
    {
	    _elements[i] = _elements.back();
	    _elements.pop_back();

	    return;
    }

  cout << "  Bin<T>::PROBLEM -- never found element!\n";
}

//------------------------------------------------------------------------
// Function    : BinningStructure::Constructor
// Description : 
//------------------------------------------------------------------------
template <int dim, class T>
BinningStructure<dim,T>::BinningStructure( float bin_w, 
                                           const vec<dim> &start,
                                           const vec<dim> &end )
{
  // the start and end of the domain
  _start = start; _end = end;

  // the size of the bins
  _bin_width = bin_w;

  // number of bins in the xyz directions
  for ( int i = 0; i < dim; i++ )
    _dimension[i] = (int)ceil( (_end(i) - _start(i)) / _bin_width );

  // and create the appropriate number of bins
  _bins = new Bin<T>[numBins()];
  
  cout << "BinningStructure num bins: " << numBins() << std::endl;
}

//------------------------------------------------------------------------
// Function    : BinningStructure::whichBin
// Description : determine which bin a position lies within
//------------------------------------------------------------------------
template <int dim, class T>
int BinningStructure<dim,T>::whichBin( const vec<dim> &pos ) const
{
  int index, bin_index=0, offset;
  for ( int i = 0; i < dim; i++ )
  {
    index = min( (int)floor( (pos[i]-_start[i]) / _bin_width ),
                 (_dimension[i]-1) );
    
    offset = 1;
    for ( int j = 0; j < i; j++ )
      offset *= _dimension[j];
      
    bin_index += index * offset;            
  }
  
  return bin_index;
}

//------------------------------------------------------------------------
// Function    : BinningStructure::clearBins
// Description : remove all the objects stored in the Bins
//------------------------------------------------------------------------
template <int dim, class T>
void BinningStructure<dim,T>::clearBins()
{
  for ( int i = 0; i < numBins(); i++ )
    _bins[i].clear();
}

//------------------------------------------------------------------------
// Function    : BinningStructure::printNeighbors
// Description : print out each bins' neighbors
//------------------------------------------------------------------------
template <int dim, class T>
void BinningStructure<dim,T>::printNeighbors() const
{
  for ( int i = 0; i < numBins(); i++ )
  {
    cout << "BIN " << i << " : ";
    for ( int j = 0; j < _bins[i].numNeighbors(); j++ )
	    cout << _bins[i].getNeighborAt( j ) << " ";

    cout << endl;
  }
}

//------------------------------------------------------------------------
// Function    : BinningStructure::numBins
// Description : return the total number of bins in the structure
//------------------------------------------------------------------------
template <int dim, class T>
int BinningStructure<dim,T>::numBins() const
{
  int num = 1;
  for ( int i = 0; i < dim; i++ )
    num *= _dimension[i];
    
  return num;
}

//------------------------------------------------------------------------
// Function    : BinningStructure::buildBinNeighbors
// Description : build the neighbor info of each bin's neighborhood
//------------------------------------------------------------------------
template <int dim, class T>
void BinningStructure<dim,T>::buildBinNeighbors()
{
  /**********************************************************************/
  //                              2D                                    //
  /**********************************************************************/
  if ( dim == 2 )
  {
    // set the neighboring bins
    int ni, nj;
    int index;
    for ( int i = 0; i < _dimension[0]; i++ )
      for ( int j = 0; j < _dimension[1]; j++ )
	    {
	      index = i + j*_dimension[0];
	      
	      //-------------- bottom row ----------------//
	      nj = j-1;
	      if ( nj >= 0 )
	      {
	        // left bin
	        ni = i-1;
	        if ( ni >= 0 )
		        _bins[index].addNeighbor( ni + nj*_dimension[0] );

	        // center bin
	        ni = i;
	        _bins[index].addNeighbor( ni + nj*_dimension[0] );

	        // right bin
	        ni = i+1;
	        if ( ni < _dimension[0] )
		        _bins[index].addNeighbor( ni + nj*_dimension[0] );
			  }
  	    

	      //-------------- middle row ---------------//
	      nj = j;

	      // left bin
	      ni = i-1;
	      if ( ni >= 0 )
	        _bins[index].addNeighbor( ni + nj*_dimension[0] );

	      // right bin
	      ni = i+1;
	      if ( ni < _dimension[0] )
	        _bins[index].addNeighbor( ni + nj*_dimension[0] );
        
	      //--------------- top row ----------------//
	      nj = j+1;
	      if ( nj < _dimension[1] )
	      {
	        // left bin
	        ni = i-1;
	        if ( ni >= 0 )
		        _bins[index].addNeighbor( ni + nj*_dimension[0] );

	        // center bin
	        ni = i;
	        _bins[index].addNeighbor( ni + nj*_dimension[0] );

	        // right bin
	        ni = i+1;
	        if ( ni < _dimension[0] )
		        _bins[index].addNeighbor( ni + nj*_dimension[0] );
	      }	      
	    }
  }
  
  /**********************************************************************/
  //                              3D                                    //
  /**********************************************************************/
  else if ( dim == 3 )
  {
    // set the neighboring bins
    int ni, nj, nk;
    int index;
    for ( int i = 0; i < _dimension[0]; i++ )
      for ( int j = 0; j < _dimension[1]; j++ )
        for ( int k = 0; k < _dimension[2]; k++ )
	      {
	        index = i + j*_dimension[0] + k*_dimension[0]*_dimension[1];
                     
	        //------------------------------------------//
	        //               front slab                 //
	        //------------------------------------------//
	        nk = k-1;
	        if ( nk >= 0 )
	        {
	          //------------ bottom row ---------------//
	          nj = j-1;
	          if ( nj >= 0 )
		        {
		          // left bin
		          ni = i-1;
		          if ( ni >= 0 )
		            _bins[index].addNeighbor( ni +
					              nj*_dimension[0] +
					              nk*_dimension[0]*_dimension[1] );

		          // center bin
		          ni = i;
		          _bins[index].addNeighbor( ni +
					            nj*_dimension[0] +
					            nk*_dimension[0]*_dimension[1] );


		          // right bin
		          ni = i+1;
		          if ( ni < _dimension[0] )
		            _bins[index].addNeighbor( ni +
					              nj*_dimension[0] +
					              nk*_dimension[0]*_dimension[1] );

		        }
            
	          //------------ middle row --------------//
	          nj = j;

	          // left bin
	          ni = i-1;
	          if ( ni >= 0 )
		          _bins[index].addNeighbor( ni +
					            nj*_dimension[0] +
					            nk*_dimension[0]*_dimension[1] );

	          // center bin
	          ni = i;
	          _bins[index].addNeighbor( ni +
					    nj*_dimension[0] +
					    nk*_dimension[0]*_dimension[1] );

	          // right bin
	          ni = i+1;
	          if ( ni < _dimension[0] )
		          _bins[index].addNeighbor( ni +
					            nj*_dimension[0] +
					            nk*_dimension[0]*_dimension[1] );
            
	          //-------------- top row ---------------//
	          nj = j+1;
	          if ( nj < _dimension[1] )
		        {
		          // left bin
		          ni = i-1;
		          if ( ni >= 0 )
		            _bins[index].addNeighbor( ni +
					              nj*_dimension[0] +
					              nk*_dimension[0]*_dimension[1] );

		          // center bin
		          ni = i;
		          _bins[index].addNeighbor( ni +
					            nj*_dimension[0] +
					            nk*_dimension[0]*_dimension[1] );

		          // right bin
		          ni = i+1;
		          if ( ni < _dimension[0] )
		            _bins[index].addNeighbor( ni +
					              nj*_dimension[0] +
					              nk*_dimension[0]*_dimension[1] );
		        }
	        }

	        //------------------------------------------//
	        //               center slab                //
	        //------------------------------------------//
	        nk = k;

	        //-------------- bottom row ----------------//
	        nj = j-1;
	        if ( nj >= 0 )
	        {
	          // left bin
	          ni = i-1;
	          if ( ni >= 0 )
		          _bins[index].addNeighbor( ni +
					            nj*_dimension[0] +
					            nk*_dimension[0]*_dimension[1] );

	          // center bin
	          ni = i;
	          _bins[index].addNeighbor( ni +
					    nj*_dimension[0] +
					    nk*_dimension[0]*_dimension[1] );

	          // right bin
	          ni = i+1;
	          if ( ni < _dimension[0] )
		          _bins[index].addNeighbor( ni +
					            nj*_dimension[0] +
					            nk*_dimension[0]*_dimension[1] );
			    }
  	      

	        //-------------- middle row ---------------//
	        nj = j;

	        // left bin
	        ni = i-1;
	        if ( ni >= 0 )
	          _bins[index].addNeighbor( ni +
				            nj*_dimension[0] +
				            nk*_dimension[0]*_dimension[1] );

	        // right bin
	        ni = i+1;
	        if ( ni < _dimension[0] )
	          _bins[index].addNeighbor( ni +
				            nj*_dimension[0] +
				            nk*_dimension[0]*_dimension[1] );
          
	        //--------------- top row ----------------//
	        nj = j+1;
	        if ( nj < _dimension[1] )
	        {
	          // left bin
	          ni = i-1;
	          if ( ni >= 0 )
		          _bins[index].addNeighbor( ni +
					            nj*_dimension[0] +
					            nk*_dimension[0]*_dimension[1] );

	          // center bin
	          ni = i;
	          _bins[index].addNeighbor( ni +
					      nj*_dimension[0] +
					      nk*_dimension[0]*_dimension[1] );

	          // right bin
	          ni = i+1;
	          if ( ni < _dimension[0] )
		          _bins[index].addNeighbor( ni +
					            nj*_dimension[0] +
					            nk*_dimension[0]*_dimension[1] );
	        }

	        //------------------------------------------//
	        //               back slab                  //
	        //------------------------------------------//
	        nk = k+1;
	        if ( nk < _dimension[2] )
	        {
	          //------------- bottom row ---------------//
	          nj = j-1;
	          if ( nj >= 0 )
		        {
		          // left bin
		          ni = i-1;
		          if ( ni >= 0 )
		            _bins[index].addNeighbor( ni +
					              nj*_dimension[0] +
					              nk*_dimension[0]*_dimension[1] );

		          // center bin
		          ni = i;
		          _bins[index].addNeighbor( ni +
					            nj*_dimension[0] +
					            nk*_dimension[0]*_dimension[1] );

		          // right bin
		          ni = i+1;
		          if ( ni < _dimension[0] )
		            _bins[index].addNeighbor( ni +
					              nj*_dimension[0] +
					              nk*_dimension[0]*_dimension[1] );
		        }

	          //------------- middle row --------------//
	          nj = j;

	          // left bin
	          ni = i-1;
	          if ( ni >= 0 )
		          _bins[index].addNeighbor( ni +
					            nj*_dimension[0] +
					            nk*_dimension[0]*_dimension[1] );

	            // center bin
	            ni = i;
	            _bins[index].addNeighbor( ni +
					      nj*_dimension[0] +
					      nk*_dimension[0]*_dimension[1] );

	            // right bin
	            ni = i+1;
	            if ( ni < _dimension[0] )
		            _bins[index].addNeighbor( ni +
					              nj*_dimension[0] +
					              nk*_dimension[0]*_dimension[1] );
           
	          //-------------- top row ---------------//
	          nj = j+1;
	          if ( nj < _dimension[1] )
		        {
		          // left bin
		          ni = i-1;
		          if ( ni >= 0 )
		            _bins[index].addNeighbor( ni +
					              nj*_dimension[0] +
					              nk*_dimension[0]*_dimension[1] );

		          // center bin
		          ni = i;
		          _bins[index].addNeighbor( ni +
					            nj*_dimension[0] +
					            nk*_dimension[0]*_dimension[1] );

		          // right bin
		          ni = i+1;
		          if ( ni < _dimension[0] )
		            _bins[index].addNeighbor( ni +
					              nj*_dimension[0] +
					              nk*_dimension[0]*_dimension[1] );
		        }
	        }
	      }
  }
  
  /**********************************************************************/
  //                              nD                                    //
  /**********************************************************************/
  else {}
  
}

















