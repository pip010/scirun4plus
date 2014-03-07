// $Id: VISErr.h,v 1.1.1.1 2003/02/12 16:51:52 whitaker Exp $
// Orig author: Josh Cates,  Feb 2000
// This file defines the VISErr namespace used for exception handling in
// Vispack. 
//
// WarningStack

// Err
//   |--OutOfRange
//   |
//   |--IncompatibleOperands
//   |
//   |--OperationFailed
//   |
//   |--Matrix
//   |     |--NonSquareMatrix
//   |     |--ZeroSize
//   |
//   |--FileIO
//   |     |--FileNotFound
//   |     |--ReadOutOfData

#ifndef __VISErr
#define __VISErr

#include <string>
#include <stack>
#include <ioiostream>

namespace VISErr
{
  const short unsigned int ErrType_sz = 50;
  enum Err_t {undef_t=0,
			  out_of_range_t=1,
			  io_t=2,
			  file_not_found_t=3,
			  read_out_of_data_t=4,
			  matrix_t=5,
			  non_square_matrix_t=6,
			  zero_size_matrix_t=7,
			  viserr_debug_t=8,
			  incompatible_operands_t=9,
			  operation_failed_t=10,
			  bad_operand_t,
			  missing_package_t
			  
  };

  const char ErrType[][ErrType_sz] = { "Undefined",
									   "Out of range",
									   "Input/Output",
									   "File not found",
									   "Read out of data",
									   "Matrix",
									   "Non-square matrix",
									   "Matrix of zero size",
									   "VISErr debug type",
									   "Incompatible operands",
									   "Operation failed",
									   "Bad operand",
									   "Missing package"
  };

  // Some useful error strings
  const char UNDEFINED_RESULT[] = "Result is undefined on supplied operands";
  const char UNDEFINED_RESULT_DIMENSIONS_UNEQUAL[] = "Result is undefined on operands of unequal dimensionality";
  const char OUT_OF_RANGE[] = "Index is out of range";
  const char NON_SQUARE_MATRIX[] = "Non-square matrix";
  const char MISSING_LAPACK[] = "Lapack package required for this operation";
  const char CONCAT_ZERO_SIZED_OPERANDS[] = "Concatenated zero sized operands";
  
  // Base class for all Vispack error classes
  class Err
  {
	static const std::string package_name;
  
	std::string location;
	std::string memo;
	Err_t type;
	
  protected:
	void setType(const Err_t &t) { type = t; }
	void assign(const Err& orig)
	{
	  this->location = orig.location;
	  this->memo     = orig.memo;
	  this->type     = orig.type;
	}
	
  public:
	Err() { setType(undef_t); }
	Err(const char *l) { setType(undef_t); setLocation(l); }
	Err(const char *l, const char *m)
	  {setType(undef_t); setMemo(m); setLocation(l); }
	~Err() { }
	Err(const Err& orig) { assign(orig); }

	void setLocation(const char *l) { location = l; }
	void setMemo(const char *m)     { memo = m;     }
	void setLocation(const std::string l) { location = l; }
	void setMemo(const std::string m)     { memo = m;     }

	Err_t getType() const { return type; }
	std::string getLocation() const { return location; }
	std::string getMemo() const { return memo; }

	Err& operator=(const Err& orig) { assign(orig); return *this; }
	
	virtual void write(ostream &stream)
	  {
		stream << package_name + " detected exception in \""  + location +
		  "\" of type \"" + ErrType[type] + "\" --> " + memo;
	  }
	virtual void handle() { this->write(cerr); exit (1); }

	virtual Err *makeCopy();
  };

  // Generic inserter for VISErr errors
  inline ostream& operator<<(ostream& stream, Err &e)
	{
	  (&e)->write(stream);
	  return stream;
	}

  // VISErr debug type for testing
  class VISErrDebug : public Err
  {
	int ref_cnt;
  private:
	void assign(const VISErrDebug& orig)
	  {
		Err::assign(orig);
		ref_cnt =  orig.ref_cnt;
		ref_cnt++;
	  }	
  public:
	VISErrDebug() : Err() { setType(viserr_debug_t); ref_cnt=1;}
	~VISErrDebug() {};
	VISErrDebug(const char *l) : Err(l)
	  {setType(viserr_debug_t); ref_cnt=1; }
	VISErrDebug(const char *l, const char *m) : Err(l, m)
	  { setType(viserr_debug_t); ref_cnt=1; }
	VISErrDebug(const VISErrDebug& orig) { assign(orig); }
	VISErrDebug& operator=(const VISErrDebug& orig) {
	  assign(orig); return *this;}
	
	void write(ostream &stream) {
	  Err::write(stream);
	  stream << " debug ref_cnt=" << ref_cnt; }
	Err *makeCopy();
	int getRef_cnt() { return ref_cnt; };
	void setRef_cnt(const int &i) { ref_cnt=i; }
	
  };	  

  class MissingPackage : public Err
  {
  public:
	MissingPackage() : Err() { setType(missing_package_t); }
	~MissingPackage() {}
	MissingPackage(const char *m) : Err(m)
	  { setType(missing_package_t); }
	MissingPackage(const char *m, const char *l) : Err(m,l)
	  { setType(missing_package_t); }
	MissingPackage(const MissingPackage& orig) { assign(orig); }
	MissingPackage& operator=(const MissingPackage& orig)
	  { assign(orig); return *this; }
  };
  
  // Matrix error
  class Matrix : public Err
  {
  public:
	Matrix() : Err() { setType(matrix_t); }
	~Matrix() {}
	Matrix(const char *m) : Err(m) {setType(matrix_t); }
	Matrix(const char *m, const char *l) : Err(m,l) {setType(matrix_t);}
	Matrix(const Matrix& orig) { assign(orig); }
	Matrix& operator=(const Matrix& orig) { assign(orig); return *this; }
  };

  ///////// PROBABLY DON"T NEED THIS-- ISN"T IT JUST A BAD OPERANDS?  
  // Non square matrix
  class NonSquareMatrix : public Matrix
  {
  public:
	NonSquareMatrix() : Matrix() { setType(non_square_matrix_t); }
	~NonSquareMatrix() {}
	NonSquareMatrix(const char *m) : Matrix(m)
	  {setType(non_square_matrix_t);}
	NonSquareMatrix(const char *m, const char *l) : Matrix(m,l)
	  { setType(non_square_matrix_t); }
	NonSquareMatrix(const NonSquareMatrix& orig) { assign(orig); }
	NonSquareMatrix& operator=(const NonSquareMatrix& orig)
	  { assign(orig); return *this; }
  };
  ////////////////////////////////////////////

  // Zero size matrix
  class ZeroSizeMatrix : public Matrix
  {
  public:
	ZeroSizeMatrix() : Matrix() { setType(zero_size_matrix_t); }
	~ZeroSizeMatrix() {}
	ZeroSizeMatrix(const char *m) : Matrix(m)
	  {setType(zero_size_matrix_t); }
	ZeroSizeMatrix(const char *m, const char *l) : Matrix(m,l)
	  {setType(zero_size_matrix_t); }
	ZeroSizeMatrix(const ZeroSizeMatrix& orig) { assign(orig); }
	ZeroSizeMatrix& operator=(const ZeroSizeMatrix& orig) { assign(orig); return *this; }
  };

  // BadOperand
  class BadOperand : public Err
  {
  public:
	BadOperand() : Err() { setType(bad_operand_t); }
	~BadOperand() {};
    BadOperand(const char *l) : Err(l) { setType(bad_operand_t);}
    BadOperand(const char *l, const char *m) : Err(l,m)
	  { setType(bad_operand_t);}
    BadOperand(const BadOperand &orig) { assign(orig); }
	BadOperand& operator=( const BadOperand &orig)
	  { assign(orig); return *this; }
  };
  
  // Incompatible operands
  class IncompatibleOperands : public Err
  {
  public:
	IncompatibleOperands() : Err() {setType(incompatible_operands_t); };
	~IncompatibleOperands() {};
    IncompatibleOperands(const char *l) : Err(l)
	  { setType(incompatible_operands_t);}
    IncompatibleOperands(const char *l, const char *m) : Err(l,m)
	  { setType(incompatible_operands_t);}
	IncompatibleOperands(const IncompatibleOperands &orig) { assign(orig); }
	IncompatibleOperands& operator=( const IncompatibleOperands &orig)
	  {assign(orig); return *this; }
  };
	
  // Input/Output (IO) error
  class IO : public Err
  {
  public:
	IO() : Err() { setType(io_t); }
	~IO() {};
	IO(const char *m) : Err(m) { setType(io_t); }
	IO(const char *m, const char *l) : Err(m,l) { setType(io_t); }
	IO(const IO& orig) {assign(orig);}
	IO& operator=(const IO& orig) { assign(orig); return *this; }
  };	

    // File not found  (IO) error
  class FileNotFound : public IO
  {
  public:
	FileNotFound() : IO() { setType(file_not_found_t); }
	~FileNotFound() {};
	FileNotFound(const char *m) : IO(m)
	  { setType(file_not_found_t);}
	FileNotFound(const char *m, const char *l) : IO(m,l)
	  { setType(file_not_found_t); }
	FileNotFound(const FileNotFound& orig) {assign(orig);}
	FileNotFound& operator=(const FileNotFound& orig) { assign(orig); return *this; }
  };	

  // Read out of data error
  class ReadOutOfData : public IO
  {
  public:
	ReadOutOfData() : IO() { setType(read_out_of_data_t); }
	~ReadOutOfData() {};
	ReadOutOfData(const char *m) : IO(m)
	  { setType(read_out_of_data_t); }
	ReadOutOfData(const char *m, const char *l) : IO(m,l)
	  { setType(read_out_of_data_t); }
	ReadOutOfData(const ReadOutOfData& orig) { assign(orig); }
	ReadOutOfData& operator=(const ReadOutOfData& orig)
	  { assign(orig); return *this; }
  };

  // Operation failed
  class OperationFailed : public Err
  {
  public:
	OperationFailed() : Err() { setType(operation_failed_t); }
	~OperationFailed() {}
	OperationFailed(const char *m) : Err(m)
	  {setType(operation_failed_t); }
	OperationFailed(const char *m, const char *l) : Err(m,l)
	  {setType(operation_failed_t); }
	OperationFailed(const OperationFailed& orig) { assign(orig); }
	OperationFailed& operator=(const OperationFailed& orig) { assign(orig); return *this; }
  };

  // Out of range error
  class OutOfRange : public Err
  {
  private:
	void assign(const OutOfRange& orig) { Err::assign(orig);}
  public:
	OutOfRange() : Err() { setType(out_of_range_t);}
	~OutOfRange() {};
	OutOfRange(const char *m) : Err(m) {setType(out_of_range_t); }
	OutOfRange(const char *m, const char *l) : Err(m,l)
	  { setType(out_of_range_t); }
	OutOfRange(const OutOfRange& orig) { assign(orig); }
	OutOfRange& operator=(const OutOfRange& orig)
	  {assign(orig); return *this;	}
  };
  

  
  // Object factory method
  extern void *New(const Err_t &);

  // Warning stack for holding Err's
  class VISWarningStack
  {
    std::stack<Err *> w_stack;
	static const std::string DELIM;
	void copyStack(const std::stack<Err *> &) ;
	void assign(const VISWarningStack &orig) { copyStack(orig.w_stack); }
 	
  public:
	VISWarningStack() {};
	~VISWarningStack() { this->purge(); }

	VISWarningStack(const VISWarningStack& orig) { assign(orig);}

	VISWarningStack &operator=(const VISWarningStack& orig)
  	  { assign(orig); return *this;}
	
	void push(Err *e) { w_stack.push(e); }
    void pop() { w_stack.pop(); }
	Err *top() { return w_stack.top(); }  // Precondition: empty() is false
    bool empty()  { return w_stack.empty(); }
    unsigned int size() { return w_stack.size(); }
	
	void write(ostream &);
	void handle();
	void purge() { while (!w_stack.empty()) { delete w_stack.top();  w_stack.pop();} };

  };
  
  // Inserter for VISWarningStack
  inline ostream& operator<<(ostream& stream, VISWarningStack &e)
	{e.write(stream); return stream;}
  
  extern VISWarningStack Warnings;  // VISErr::WarningStack
  


  // Some predefined inline functions for easy error handling

  inline void WARN(const Err_t type) { Warnings.push((Err *)New(type)); }
  inline void WARN(const Err_t type, const char *l)
	{
	  Err *e = (Err *)New(type);
	  e->setLocation(l);
	  Warnings.push(e);
	}
  inline void WARN(const Err_t type, const char *l, const char *m)
	{
	  Err *e = (Err *)New(type);
	  e->setLocation(l);
	  e->setMemo(m);
	  Warnings.push(e);
	}

}
#endif
