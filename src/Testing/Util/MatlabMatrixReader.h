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

#ifndef MATLAB_MATRIX_READER_H__
#define MATLAB_MATRIX_READER_H__ 1

#include <stdexcept>
#include <boost/filesystem.hpp>
#include <Core/Datatypes/Matrix.h>
#include <Core/Matlab/matlabarray.h>
#include <Core/Matlab/matlabconverter.h>

namespace SCIRun
{
  namespace TestUtils
  {
    class EasyMatlabMatrixReader
    {
    public:
      explicit EasyMatlabMatrixReader(const boost::filesystem::path& filename)
      {
        init(filename.string());
      }
      explicit EasyMatlabMatrixReader(const std::string& filename)
      {
        init(filename);
      }

      ~EasyMatlabMatrixReader()
      {
        file_.close();
      }

      MatrixHandle read_matrix(const std::string& name)
      {
        try
        {
          MatlabIO::matlabarray marray = file_.getmatlabarray(name);
          MatrixHandle matrix;
          converter_.mlArrayTOsciMatrix(marray, matrix);
          return matrix;
        }
        catch (...)
        {
          return MatrixHandle();
        }
      }

    private:
      MatlabIO::matlabfile file_;
      MatlabIO::matlabconverter converter_;

      void init(const std::string& filename)
      {
        using namespace boost::filesystem;
        if (!exists(filename))
        {
          std::string fullFileName = system_complete(path(filename)).string();
          throw std::invalid_argument(fullFileName + " does not exist.");
        }
        file_.open(filename, "r");
      }
    };

    typedef boost::shared_ptr<EasyMatlabMatrixReader> MatlabMatrixReaderHandle;

    inline MatlabMatrixReaderHandle open_matlab_file(const std::string& filename)
    {
      return MatlabMatrixReaderHandle(new EasyMatlabMatrixReader(filename));
    }
  }
}



#endif