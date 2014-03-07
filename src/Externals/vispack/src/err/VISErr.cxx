// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: VISErr.cxx,v 1.1.1.1 2003/02/12 16:51:52 whitaker Exp $
// Orig author: Josh Cates (jecates@utk.edu), Feb 2000
// Implements heirarchy of exception classes for Vispack

#include <ioiostream>
#include <string>
#include "VISErr.h"

namespace VISErr
{
  const std::string Err::package_name = "VISPack";
  const std::string VISWarningStack::DELIM = "\n";
  VISWarningStack Warnings;

  void *New(const Err_t &t)
	{
	  switch (t)
		{
		case (out_of_range_t):
		 return new OutOfRange();
		 break;
		case (undef_t):
		  return new Err();
		  break;
		case (io_t):
		  return new IO();
		  break;
		case (file_not_found_t):
		  return new FileNotFound();
		  break;
		case (read_out_of_data_t):
		  return new ReadOutOfData();
		  break;
		case (matrix_t):
		  return new Matrix();
		  break;
		case (non_square_matrix_t):
		  return new NonSquareMatrix();
		  break;
		case (zero_size_matrix_t):
		  return new ZeroSizeMatrix();
		  break;
		case (viserr_debug_t):
		  return new VISErrDebug();
		  break;
		case (incompatible_operands_t):
		  return new IncompatibleOperands();
		  break;
		case (operation_failed_t):
		  return new OperationFailed();
		  break;
		case (bad_operand_t):
		  return new BadOperand();
		  break;
		case (missing_package_t):
		  return new MissingPackage();
		  break;
		default:
		  return new Err();
		  break;
		}  
	  return 0;
	}

  Err *Err::makeCopy()
	{
	  Err *temp = (Err *)New(this->getType());
	  temp->setLocation(this->getLocation());
	  temp->setMemo(this->getMemo());
	  return temp;
	}

  Err *VISErrDebug::makeCopy()
	{
	  VISErrDebug *temp = (VISErrDebug *)New(this->getType());
	  temp->setLocation(this->getLocation());
	  temp->setMemo(this->getMemo());
	  temp->setRef_cnt(this->getRef_cnt());
	  this->setRef_cnt((this->getRef_cnt())+1);
	  return (Err *)temp;
	}

  //
  // TEST CODE
  //
  void VISWarningStack::copyStack(const std::stack<Err *>& orig)
	{
	  int i=0;
	  Err **errArr = new Err *[orig.size()];
	  std::stack<Err *> temp_stack = orig;
	  Err *temp_err;
	  
	  // Copy stack in reverse order
	  while ( !temp_stack.empty() )
		{
		  temp_err = (temp_stack.top())->makeCopy();
		  errArr[i++] = temp_err;
		  temp_stack.pop();
		}
	  
	  // Place in proper order in the new stack
	  for (i = orig.size()-1; i>=0; i--) { w_stack.push(errArr[i]); }
	  
	  delete [] errArr;
	  
	}
  
  void VISWarningStack::write(ostream &stream)
	{
	  // Copy the stack so we can pop elements off without disturbing the
	  // original, then write, pop, each element of the copied stack.
	  VISWarningStack temp=*this;

	  while (! temp.empty() )
		{
		  (temp.top())->write(stream);
		  stream << this->DELIM;
		  temp.pop();
		}

	} 

    void VISWarningStack::handle()
	{
	  // Copy the stack so we can pop elements off without disturbing the
	  // original, then handle, pop, each element of the copied stack.
	  VISWarningStack temp=*this;

	  while (! temp.empty() )
		{
		  (temp.top())->handle();
		  temp.pop();
		}

	} 

  
}
