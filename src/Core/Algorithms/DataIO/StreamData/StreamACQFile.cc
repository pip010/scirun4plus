/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2009 Scientific Computing and Imaging Institute,
   University of Utah.

   
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/

#include <Core/Algorithms/DataIO/StreamData/StreamACQFile.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Datatypes/ColumnMatrix.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/SparseRowMatrix.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

#include <string>
#include <vector>

#ifdef _WIN32
  #include <fcntl.h>
  #include <stdio.h>
  #include <ctype.h>
  #include <io.h>



#else
  #ifndef __USE_LARGEFILE64
    #define __USE_LARGEFILE64
  #endif
  #include <unistd.h>
  #include <fcntl.h>
  #include <stdio.h>
  #include <ctype.h>
#endif
  
#ifndef O_LARGEFILE
  #define O_LARGEFILE 0
#endif

#include <Core/Algorithms/DataIO/share.h>


namespace SCIRunAlgo {

using namespace SCIRun;

class SCISHARE StreamACQFileReaderAlgo {
  public:
    StreamACQFileReaderAlgo(AlgoBase* algo) : algo_(algo) {}
    
    bool scan_header();
    
    bool get_matrix_from_frame_indices(SCIRun::MatrixHandle indices,
                                       SCIRun::MatrixHandle& matrix);

    bool get_matrix_from_frame_indices(SCIRun::index_type start,
                                       SCIRun::index_type end,
                                       SCIRun::MatrixHandle& matrix);

    bool get_matrix_from_lead_indices(SCIRun::MatrixHandle indices,
                                       SCIRun::MatrixHandle& matrix);

    bool get_matrix_from_lead_indices(SCIRun::index_type start,
                                       SCIRun::index_type end,
                                       SCIRun::MatrixHandle& matrix);

    void swap_bytes(void *vbuffer,size_t elsize,size_t size);
    
  private:
    AlgoBase* algo_;
    
    short       num_leads_;
    int         num_frames_;
    std::string time_stamp_;
    std::string label_;
    
    
    std::string filename_;
    bool        swap_bytes_;
};

void 
StreamACQFileReaderAlgo::
swap_bytes(void *vbuffer,size_t elsize,size_t size)
{
   char temp;
   char *buffer = static_cast<char *>(vbuffer);

   size *= elsize;

   switch(elsize)
   {
      case 0:
      case 1:
         // Do nothing. Element size is 1 byte, so there is nothing to swap
         break;
      case 2:  
		// Do a 2 bytes element byte swap. 
		for(size_t p=0;p<size;p+=2)
		  { temp = buffer[p]; buffer[p] = buffer[p+1]; buffer[p+1] = temp; }
		break;
      case 4:
		// Do a 4 bytes element byte swap.
		for(size_t p=0;p<size;p+=4)
		  { temp = buffer[p]; buffer[p] = buffer[p+3]; buffer[p+3] = temp; 
			temp = buffer[p+1]; buffer[p+1] = buffer[p+2]; buffer[p+2] = temp; }
		break;
      case 8:
		// Do a 8 bytes element byte swap.
		for(size_t p=0;p<size;p+=8)
		  { temp = buffer[p]; buffer[p] = buffer[p+7]; buffer[p+7] = temp; 
			temp = buffer[p+1]; buffer[p+1] = buffer[p+6]; buffer[p+6] = temp; 
			temp = buffer[p+2]; buffer[p+2] = buffer[p+5]; buffer[p+5] = temp; 
			temp = buffer[p+3]; buffer[p+3] = buffer[p+4]; buffer[p+4] = temp; }
   	    break;
      default:
       return;    
   }  
}

bool
StreamACQFileReaderAlgo::
scan_header()
{
  // Set the filename  
  filename_ = algo_->get_filename("filename");
  
  // Determine whether we need to swap bytes
  // ACQ files are always big endian
  swap_bytes_ = false;
  short test = 0x00FF;
  char *testptr = reinterpret_cast<char *>(&test);
  if (testptr[0]) swap_bytes_ = true; 
  
  // Use Unix open/read/write so we can use files larger than 2Gb
  int acqfile = 0;
  acqfile = ::open(filename_.c_str(),O_RDONLY|O_LARGEFILE,0);
  if (acqfile < 0)
  {
    algo_->error("Could not find/open "+filename_);
    return (false);
  }

  // Read number of leads
  if (::lseek(acqfile,606,SEEK_SET)<0)
  {
    ::close(acqfile);
    algo_->error("Could not read file "+filename_);
    return (false);
  }  
  if (::read(acqfile,&(num_leads_),sizeof(short))!= sizeof(short))
  {
    ::close(acqfile);
    algo_->error("Could not read file "+filename_);
    return (false);
  }  
  if (swap_bytes_) swap_bytes(&(num_leads_),sizeof(short),1);
  
  // Read number of frames
  if (::lseek(acqfile,608,SEEK_SET)<0)
  {
    ::close(acqfile);
    algo_->error("Could not read file "+filename_);
    return (false);
  }  
  if (::read(acqfile,&(num_frames_),sizeof(int)) != sizeof(int))
  {
    ::close(acqfile);
    algo_->error("Could not read file "+filename_);
    return (false);
  }  
  if (swap_bytes_) swap_bytes(&(num_frames_),sizeof(int),1);

  // Read Time stamp
  char buffer[82];
  
  if (::lseek(acqfile,580,SEEK_SET)<0)
  {
    ::close(acqfile);
    algo_->error("Could not read file "+filename_);
    return (false);
  }

  if (::read(acqfile,buffer,sizeof(char)*12) != 12*sizeof(char))
  {
    ::close(acqfile);
    algo_->error("Could not read file "+filename_);
    return (false);
  } 
  for (size_t j=0; j<12;j++)
    if (buffer[j] == '\n' || buffer[j] == '\t' || buffer[j] == '\r') buffer[j] = ' ';
  
  buffer[12] = 0;
  time_stamp_ = std::string(buffer);
  
  // Read label
  if (::lseek(acqfile,122,SEEK_SET)<0)
  {
    ::close(acqfile);
    algo_->error("Could not read file "+filename_);
    return (false);
  }

  if (::read(acqfile,buffer,sizeof(char)*80) != 80*sizeof(char))
  {
    ::close(acqfile);
    algo_->error("Could not read file "+filename_);
    return (false);
  }
  
  for (size_t j=0; j<80;j++)
    if (buffer[j] == '\n' || buffer[j] == '\t' || buffer[j] == '\r') buffer[j] = ' ';

    
  buffer[80] = 0;
  label_ = std::string(buffer);  

  ::close(acqfile);
  
  algo_->set_int("num_leads",num_leads_);
  algo_->set_int("num_frames",num_frames_);
  algo_->set_string("time_stamp",time_stamp_);
  algo_->set_string("label",label_);
  
  return(true);
}

bool
StreamACQFileReaderAlgo::
get_matrix_from_frame_indices(SCIRun::MatrixHandle indices,
                              SCIRun::MatrixHandle& matrix)
{
  if(!(scan_header())) return (false);
  
  // Check whether we have an index or multiple indices
  if (indices.get_rep() == 0)
  {
    algo_->error("No indices given for extraction of colums");
    return (false);
  }
  
  // Copy the indices into a more accessible vector
  // Get the number of elements (sparse or dense it does not matter)
  std::vector<index_type> idx(indices->get_data_size());
  double *idataptr =  indices->get_data_pointer();

  if (idataptr == 0)
  {
    algo_->error("Could not read index matrix");
    return (false);
  }
    
  // Copy and cast the indices to integers
  for (size_t p=0; p<idx.size(); p++) 
  {
    idx[p] = static_cast<int>(idataptr[p]);
    if (idx[p] < 0) idx[p] = 0;
    if (idx[p] >= num_frames_) idx[p] = num_frames_-1; 
  }
  
  // Create the output matrix
  MatrixHandle temp_matrix = new SCIRun::DenseMatrix(idx.size(),num_leads_);
  if (temp_matrix.get_rep() == 0)
  {
    algo_->error("Could not allocate output matrix");
    return (false);
  }   
  
  // Get the data pointer where we can store the data
  char* buffer = reinterpret_cast<char *>(temp_matrix->get_data_pointer());

  if (buffer == 0)
  {
    algo_->error("Internal error: encountered an matrix with no data");
    return (false);
  }

  int acqfile = 0;
  acqfile = ::open(filename_.c_str(),O_RDONLY|O_LARGEFILE,0);
  if (acqfile < 0)
  {
    algo_->error("Could not find/open "+filename_);
    return (false);
  }


  size_type num_indices = static_cast<size_type>(idx.size());
  // Loop over all the indices
  for (index_type k=0; k< num_indices ; k++)
  {
    index_type coffset = idx[k];
    size_type cnt = 1;
    while ((k+1 < num_indices) && (idx[k+1] == idx[k]+1)) { cnt++; k++; }
    
    if (::lseek(acqfile,static_cast<off_t>(2)*static_cast<off_t>(num_leads_)*
          static_cast<off_t>(coffset)+static_cast<off_t>(1024),SEEK_SET)<0)
    {
      ::close(acqfile);
      algo_->error("Improper ACQ file, check number of columns and rows in header file");
      return (false);
    }    

    size_t read_size = static_cast<size_t>(num_leads_)*sizeof(short)*cnt;
    size_t ret = ::read(acqfile,reinterpret_cast<void*>(buffer),read_size);
    if (read_size != ret)
    {
      ::close(acqfile);
      algo_->error("Improper ACQ file, check number of columns and rows in header file");
      return (false);
    }
    
    // Move pointer to next column to be read
    buffer += (read_size);
  }
  
  ::close(acqfile);
  
  // Get the starting pointer again.
  buffer = reinterpret_cast<char *>(temp_matrix->get_data_pointer());
  
  // Do the byte swapping and then  unpack the data into doubles
  if (swap_bytes_) swap_bytes(reinterpret_cast<void*>(buffer),sizeof(short),idx.size()*num_leads_);
  short *fbuffer = reinterpret_cast<short *>(buffer); 
  double *dbuffer = reinterpret_cast<double *>(buffer);  
  MatrixHandle GainTable;
  algo_->get_handle("gain_table",GainTable);
  if (GainTable.get_rep())
  {
    size_type m = GainTable->nrows();
    size_type n = GainTable->ncols();
    double* gains = GainTable->get_data_pointer();
    
    if ((m != 1)&&(m != num_leads_))
    {
      algo_->error("GainTable needs to have the number of rows equal to the number of leads");
      return (false);
    }
    
    if (m == 1 && n == 1)
    {
      for (index_type j=(idx.size()*num_leads_-1);j>=0;j--) 
        { dbuffer[j] = gains[0]*static_cast<double>((0x0FFF&fbuffer[j])-2048); }
    }
    else if (m == num_leads_ && n == 1)
    {
      index_type i = num_leads_;
      for (index_type j=(idx.size()*num_leads_-1);j>=0;j--) 
      {
        i--;
        dbuffer[j] = gains[i]*static_cast<double>((0x0FFF&fbuffer[j])-2048); 
        if (i == 0) i = num_leads_;
      }    
    }
    else if (m == 1 && n == 8)
    {
      for (index_type j=(idx.size()*num_leads_-1);j>=0;j--) 
      {
        index_type gain  = (0x7000&fbuffer[j])>>12;
        dbuffer[j] = gains[gain]*static_cast<double>((0x0FFF&fbuffer[j])-2048); 
      }
    }
    else if (m == num_leads_ && n == 8)
    {
      index_type i = num_leads_;
      for (index_type j=(idx.size()*num_leads_-1);j>=0;j--) 
      { 
        i--;
        index_type gain  = (0x7000&fbuffer[j])>>12;
        dbuffer[j] = gains[gain+(i<<3)]*static_cast<double>((0x0FFF&fbuffer[j])-2048); 
        if (i == 0) i = num_leads_;        
      }    
    }
    else
    {
      algo_->error("Internal error: algorithm got to a point it should not get");
      return (false);
    }
  }
  else
  {    
    for (index_type j=(idx.size()*num_leads_-1);j>=0;j--) 
      { dbuffer[j] = static_cast<double>((0x0FFF&fbuffer[j])-2048); }
  }
  
  // From Fortran to C indexing
  matrix = temp_matrix->make_transpose();
  
  if (matrix.get_rep()  == 0)
  {
    algo_->error("Could not reorder data");
    return (false);  
  }

  return (true);
}


bool
StreamACQFileReaderAlgo::
get_matrix_from_lead_indices(SCIRun::MatrixHandle indices,
                             SCIRun::MatrixHandle& matrix)
{
  if(!(scan_header())) return (false);
  
  // Check whether we have an index or multiple indices
  if (indices.get_rep() == 0)
  {
    algo_->error("No indices given for extraction of colums");
    return (false);
  }
  
  // Copy the indices into a more accessible vector
  // Get the number of elements (sparse or dense it does not matter)
  std::vector<index_type> idx(indices->get_data_size());
  double *idataptr =  indices->get_data_pointer();

  if (idataptr == 0)
  {
    algo_->error("Could not read index matrix");
    return (false);
  }
    
  // Copy and cast the indices to integers
  for (size_t p=0; p<idx.size(); p++) 
  {
    idx[p] = static_cast<int>(idataptr[p]);
    if (idx[p] < 0) idx[p] = 0;
    if (idx[p] >= num_leads_) idx[p] = num_leads_-1; 
  }
  
  // Create the output matrix
  MatrixHandle temp_matrix = new SCIRun::DenseMatrix(num_frames_,idx.size());
  if (temp_matrix.get_rep() == 0)
  {
    algo_->error("Could not allocate output matrix");
    return (false);
  }   
  
  // Get the data pointer where we can store the data
  char* buffer = reinterpret_cast<char *>(temp_matrix->get_data_pointer());

  if (buffer == 0)
  {
    algo_->error("Internal error: encountered an matrix with no data");
    return (false);
  }

  int acqfile = 0;
  acqfile = ::open(filename_.c_str(),O_RDONLY|O_LARGEFILE,0);
  if (acqfile < 0)
  {
    algo_->error("Could not find/open "+filename_);
    return (false);
  }

  size_type num_indices = static_cast<size_type>(idx.size());
  
  // Loop over all the indices
  for (index_type j=0; j< num_frames_; j++)
  {
    for (index_type k=0; k< num_indices ; k++)
    {
      size_type cnt = 1;
      while ((k+1 < num_indices) && (idx[k+1] == idx[k]+1)) { cnt++; k++; }
      
      if (::lseek(acqfile,static_cast<off_t>(2)*static_cast<off_t>(num_leads_)*
            static_cast<off_t>(j)+static_cast<off_t>(1024)+
            static_cast<off_t>(2)*static_cast<off_t>(idx[k]),SEEK_SET)<0)
      {
        ::close(acqfile);
        algo_->error("Improper ACQ file, check number of columns and rows in header file");
        return (false);
      }    

      size_t read_size = sizeof(short)*cnt;
      size_t ret = ::read(acqfile,reinterpret_cast<void*>(&buffer),read_size);
      if (read_size != ret)
      {
        ::close(acqfile);
        algo_->error("Improper ACQ file, check number of columns and rows in header file");
        return (false);
      }
      
      buffer += read_size;
      // Move pointer to next column to be read
    }
  }
  
  ::close(acqfile);
  
  // Get the starting pointer again.
  buffer = reinterpret_cast<char *>(temp_matrix->get_data_pointer());
  
  // Do the byte swapping and then  unpack the data into doubles
  if (swap_bytes_) swap_bytes(reinterpret_cast<void*>(buffer),sizeof(short),idx.size()*num_frames_);
  short *fbuffer = reinterpret_cast<short *>(buffer); 
  double *dbuffer = reinterpret_cast<double *>(buffer);  
  
    MatrixHandle GainTable;
  algo_->get_handle("gain_table",GainTable);
  if (GainTable.get_rep())
  {
    size_type m = GainTable->nrows();
    size_type n = GainTable->ncols();
    double* gains = GainTable->get_data_pointer();

    if ((m != 1)&&(m != num_leads_))
    {
      algo_->error("GainTable needs to have the number of rows equal to the number of leads");
      return (false);
    }
        
    if (m == 1 && n == 1)
    {
      for (index_type j=(idx.size()*num_frames_-1);j>=0;j--) 
        { dbuffer[j] = gains[0]*static_cast<double>((0x0FFF&fbuffer[j])-2048); }
    }
    else if (m == num_leads_ && n == 1)
    {
      index_type i=idx.size();
      for (index_type j=(idx.size()*num_frames_-1);j>=0;j--) 
      { 
        i--;
        dbuffer[j] = gains[idx[i]]*static_cast<double>((0x0FFF&fbuffer[j])-2048); 
        if (i==0) i = idx.size();
      }    
    }
    else if (m == 1 && n == 8)
    {
      for (index_type j=(idx.size()*num_frames_-1);j>=0;j--) 
      {
        index_type gain  = (0x7000&fbuffer[j])>>12;
        dbuffer[j] = gains[gain]*static_cast<double>((0x0FFF&fbuffer[j])-2048); 
      }
    }
    else if (m == num_leads_ && n == 8)
    {
      index_type i = idx.size();
      for (index_type j=(idx.size()*num_frames_-1);j>=0;j--) 
      { 
        i--;
        index_type gain  = (0x7000&fbuffer[j])>>12;
        dbuffer[j] = gains[gain+((idx[i])<<3)]*static_cast<double>((0x0FFF&fbuffer[j])-2048); 
        if (i == 0) i = idx.size();
      }    
    }
    else
    {
      algo_->error("Internal error: algorithm got to a point it should not get");
      return (false);
    }
  }
  else
  {      
    for (index_type j=(idx.size()*num_frames_-1);j>=0;j--) 
      { dbuffer[j] = static_cast<double>((0x0FFF&fbuffer[j])-2048); }
  }
  // From Fortran to C indexing
  matrix = temp_matrix->make_transpose();
  
  if (matrix.get_rep()  == 0)
  {
    algo_->error("Could not reorder data");
    return (false);  
  }
  
  return (true);
}


bool
StreamACQFileReaderAlgo::
get_matrix_from_frame_indices(SCIRun::index_type start,
                              SCIRun::index_type end,
                              SCIRun::MatrixHandle& matrix)
{
  if(!(scan_header())) return (false);
  
  if (start<0) start = 0;
  if (end<0) end = 0;
  if (start >= num_frames_) start = num_frames_-1;
  if (end >= num_frames_) end = num_frames_-1;

  if (start >= end)
  {
    algo_->error("Start index is equal or larger than end index");
    return (false);
  }
  
  end++;
  SCIRun::size_type num_indices = (end-start);
  
  // Create the output matrix
  MatrixHandle temp_matrix = new SCIRun::DenseMatrix(num_indices,num_leads_);
  if (temp_matrix.get_rep() == 0)
  {
    algo_->error("Could not allocate output matrix");
    return (false);
  }   
  
  // Get the data pointer where we can store the data
  char* buffer = reinterpret_cast<char *>(temp_matrix->get_data_pointer());

  if (buffer == 0)
  {
    algo_->error("Internal error: encountered an matrix with no data");
    return (false);
  }

  int acqfile = 0;
  acqfile = ::open(filename_.c_str(),O_RDONLY|O_LARGEFILE,0);
  if (acqfile < 0)
  {
    algo_->error("Could not find/open "+filename_);
    return (false);
  }


  if (::lseek(acqfile,static_cast<off_t>(2)*static_cast<off_t>(num_leads_)*
        static_cast<off_t>(start)+static_cast<off_t>(1024),SEEK_SET)<0)
  {
    ::close(acqfile);
    algo_->error("Improper ACQ file, check number of columns and rows in header file");
    return (false);
  }    

  size_t read_size = static_cast<size_t>(num_leads_)*sizeof(short)*num_indices;
  size_t ret = ::read(acqfile,reinterpret_cast<void*>(buffer),read_size);
  if (read_size != ret)
  {
    ::close(acqfile);
    algo_->error("Improper ACQ file, check number of columns and rows in header file");
    return (false);
  }
  
  ::close(acqfile);
  
  // Get the starting pointer again.
  buffer = reinterpret_cast<char *>(temp_matrix->get_data_pointer());
  
  // Do the byte swapping and then  unpack the data into doubles
  if (swap_bytes_) swap_bytes(reinterpret_cast<void*>(buffer),sizeof(short),num_indices*num_leads_);
  short *fbuffer = reinterpret_cast<short *>(buffer); 
  double *dbuffer = reinterpret_cast<double *>(buffer); 
  
  MatrixHandle GainTable;
  algo_->get_handle("gain_table",GainTable);
  if (GainTable.get_rep())
  {
    size_type m = GainTable->nrows();
    size_type n = GainTable->ncols();
    double* gains = GainTable->get_data_pointer();
    
    if ((m != 1)&&(m != num_leads_))
    {
      algo_->error("GainTable needs to have the number of rows equal to the number of leads");
      return (false);
    }
    
    if (m == 1 && n == 1)
    {
      for (index_type j=(num_indices*num_leads_-1);j>=0;j--) 
        { dbuffer[j] = gains[0]*static_cast<double>((0x0FFF&fbuffer[j])-2048); }
    }
    else if (m == num_leads_ && n == 1)
    {
      index_type i = num_leads_;
      for (index_type j=(num_indices*num_leads_-1);j>=0;j--) 
      {
        i--;
        dbuffer[j] = gains[i]*static_cast<double>((0x0FFF&fbuffer[j])-2048); 
        if (i == 0) i = num_leads_;
      }    
    }
    else if (m == 1 && n == 8)
    {
      for (index_type j=(num_indices*num_leads_-1);j>=0;j--) 
      {
        index_type gain  = (0x7000&fbuffer[j])>>12;
        dbuffer[j] = gains[gain]*static_cast<double>((0x0FFF&fbuffer[j])-2048); 
      }
    }
    else if (m == num_leads_ && n == 8)
    {
      index_type i = num_leads_;
      for (index_type j=(num_indices*num_leads_-1);j>=0;j--) 
      { 
        i--;
        index_type gain  = (0x7000&fbuffer[j])>>12;
        dbuffer[j] = gains[gain+(i<<3)]*static_cast<double>((0x0FFF&fbuffer[j])-2048); 
        if (i == 0) i = num_leads_;        
      }    
    }
    else
    {
      algo_->error("Internal error: algorithm got to a point it should not get");
      return (false);
    }
  }
  else
  {  
    for (index_type j=(num_indices*num_leads_-1);j>=0;j--) 
      { dbuffer[j] = static_cast<double>((0x0FFF&fbuffer[j])-2048); }
  }

  // From Fortran to C indexing
  matrix = temp_matrix->make_transpose();
  
  if (matrix.get_rep()  == 0)
  {
    algo_->error("Could not reorder data");
    return (false);  
  }
  
  return (true);
}



bool
StreamACQFileReaderAlgo::
get_matrix_from_lead_indices(SCIRun::index_type start,
                             SCIRun::index_type end,
                             SCIRun::MatrixHandle& matrix)
{
  if(!(scan_header())) return (false);

  if (start < 0) start = 0;
  if (end < 0) end = 0;
  if (start >= num_leads_) start = num_leads_-1;
  if (end >= num_leads_) end = num_leads_-1;
  
  end++;
  
  size_type num_indices = end-start;
  // Create the output matrix
  MatrixHandle temp_matrix = new SCIRun::DenseMatrix(num_frames_,num_indices);
  if (temp_matrix.get_rep() == 0)
  {
    algo_->error("Could not allocate output matrix");
    return (false);
  }   
  
  // Get the data pointer where we can store the data
  char* buffer = reinterpret_cast<char *>(temp_matrix->get_data_pointer());

  if (buffer == 0)
  {
    algo_->error("Internal error: encountered an matrix with no data");
    return (false);
  }

  int acqfile = 0;
  acqfile = ::open(filename_.c_str(),O_RDONLY|O_LARGEFILE,0);
  if (acqfile < 0)
  {
    algo_->error("Could not find/open "+filename_);
    return (false);
  }
  
  // Loop over all the indices
  for (index_type j=0; j< num_frames_; j++)
  {
    for (index_type k=0; k< num_indices ; k++)
    {      
      if (::lseek(acqfile,static_cast<off_t>(2)*static_cast<off_t>(num_leads_)*
            static_cast<off_t>(j)+static_cast<off_t>(1024)+
            static_cast<off_t>(start)*static_cast<off_t>(2),SEEK_SET)<0)
      {
        ::close(acqfile);
        algo_->error("Improper ACQ file, check number of columns and rows in header file");
        return (false);
      }    

      size_t read_size = sizeof(short)*num_indices;
      size_t ret = ::read(acqfile,reinterpret_cast<void*>(&buffer),read_size);
      if (read_size != ret)
      {
        ::close(acqfile);
        algo_->error("Improper ACQ file, check number of columns and rows in header file");
        return (false);
      }
      
      buffer += read_size;
      // Move pointer to next column to be read
    }
  }
  
  ::close(acqfile);
  
  // Get the starting pointer again.
  buffer = reinterpret_cast<char *>(temp_matrix->get_data_pointer());
  
  // Do the byte swapping and then  unpack the data into doubles
  if (swap_bytes_) swap_bytes(reinterpret_cast<void*>(buffer),sizeof(short),num_indices*num_frames_);
  short *fbuffer = reinterpret_cast<short *>(buffer); 
  double *dbuffer = reinterpret_cast<double *>(buffer);
  
  MatrixHandle GainTable;
  algo_->get_handle("gain_table",GainTable);
  if (GainTable.get_rep())
  {
    size_type m = GainTable->nrows();
    size_type n = GainTable->ncols();
    double* gains = GainTable->get_data_pointer();

    if ((m != 1)&&(m != num_leads_))
    {
      algo_->error("GainTable needs to have the number of rows equal to the number of leads");
      return (false);
    }
        
    if (m == 1 && n == 1)
    {
      for (index_type j=(num_indices*num_frames_-1);j>=0;j--) 
        { dbuffer[j] = gains[0]*static_cast<double>((0x0FFF&fbuffer[j])-2048); }
    }
    else if (m == num_leads_ && n == 1)
    {
      index_type i=num_indices;
      for (index_type j=(num_indices*num_frames_-1);j>=0;j--) 
      { 
        i--;
        dbuffer[j] = gains[(start+i)]*static_cast<double>((0x0FFF&fbuffer[j])-2048); 
        if (i==0) i = num_indices;
      }    
    }
    else if (m == 1 && n == 8)
    {
      for (index_type j=(num_indices*num_frames_-1);j>=0;j--) 
      {
        index_type gain  = (0x7000&fbuffer[j])>>12;
        dbuffer[j] = gains[gain]*static_cast<double>((0x0FFF&fbuffer[j])-2048); 
      }
    }
    else if (m == num_leads_ && n == 8)
    {
      index_type i = num_indices;
      for (index_type j=(num_indices*num_frames_-1);j>=0;j--) 
      { 
        i--;
        index_type gain  = (0x7000&fbuffer[j])>>12;
        dbuffer[j] = gains[gain+((start+i)<<3)]*static_cast<double>((0x0FFF&fbuffer[j])-2048); 
        if (i == 0) i = num_indices;
      }    
    }
    else
    {
      algo_->error("Internal error: algorithm got to a point it should not get");
      return (false);
    }
  }
  else
  {      
    for (index_type j=(num_indices*num_frames_-1);j>=0;j--) 
      { dbuffer[j] = static_cast<double>((0x0FFF&fbuffer[j])-2048); }
  }
  
  // From Fortran to C indexing
  matrix = temp_matrix->make_transpose();
  
  if (matrix.get_rep()  == 0)
  {
    algo_->error("Could not reorder data");
    return (false);  
  }
  
  return (true);
}

bool
StreamACQFileAlgo::read_header()
{
  StreamACQFileReaderAlgo algo_(this);
  return(algo_.scan_header());
}


bool
StreamACQFileAlgo::run(MatrixHandle indices,
             MatrixHandle& output)
{
  algo_start("StreamACQFile");
  
  std::string method = get_option("method");

  // CHeck format of the GainTable if one is given
  MatrixHandle GainTable;
  get_handle("gain_table",GainTable);
  if (GainTable.get_rep())
  {
    if ((!(matrix_is::dense(GainTable)))&&(!(matrix_is::column(GainTable))))
    {
      error("GainTable needs to be a dense matrix or a column matrix");
      algo_end(); return (false);
    }
    
    if ((GainTable->ncols() != 8)&&(GainTable->ncols() != 1))
    {
      error("GainTable needs to have 1 or 8 columns");
      algo_end(); return (false);
    }
  }


  bool ret;
  if (method == "frame_indices")
  {
    StreamACQFileReaderAlgo algo_(this);
    ret = algo_.get_matrix_from_frame_indices(indices,output);
  }
  else if (method == "lead_indices")
  {
    StreamACQFileReaderAlgo algo_(this);
    ret = algo_.get_matrix_from_lead_indices(indices,output);  
  }
  else
  {
    error("Method has not yet been implemented");
    algo_end();
  }
    
  algo_end();
  return (true);
}


bool
StreamACQFileAlgo::run(index_type start,
                   index_type end,
                   MatrixHandle& output)
{
  algo_start("StreamACQFile");
  
  std::string method = get_option("method");
  bool ret;
  if (method == "frame_indices")
  {
    StreamACQFileReaderAlgo algo_(this);
    ret = algo_.get_matrix_from_frame_indices(start,end,output);
  }
  else if (method == "lead_indices")
  {
    StreamACQFileReaderAlgo algo_(this);
    ret = algo_.get_matrix_from_lead_indices(start,end,output);  
  }
  else
  {
    error("Method has not yet been implemented");
    algo_end();
  }
  
  algo_end();
  return (true);
}

} // End namespace SCIRunAlgo
