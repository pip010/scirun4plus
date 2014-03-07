// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: param.cxx,v 1.1.1.1 2003/02/12 16:51:55 whitaker Exp $
// Josh Cates 03/07/00 jecates@utk.edu
// University of Tennessee, Knoxville
//
// This file defines the parameter file syntax tree.

#include "param.h"

namespace VISParameterFile
{
  Parameter ParameterList::NullParameter;
  Parameter ParameterFile::NullParameter;

  Value *Value::Copy(const Value &o)
	{
	  switch ( (&o)->getType() )
		{
		case integer_t:
		  return new IntegerVal((IntegerVal *)(&o));
		  break;
		case decimal_t:
		  return new DecimalVal((DecimalVal *)(&o));
		  break;
		case string_t:
		  return new StringVal((StringVal *)(&o));
		  break;
		default:
		  return 0;
		  break;
		}
	}

  void ValueList::print()
	{
	  cout << "      mem = " << this << " list size="
		   << this->size() << endl;
	  cout << "      contents: " << endl;
	  for (int i = 0; i< this->size() ; i++)
		{
		  cout << "      [" << i << "] ";
		  cout << " mem = " << List[i];
		  (*List[i]).print();
		}
	}

  void ValueList::clear()
    {
  	  // This function deletes all Value structures on the pointer list (cleans
	  // up all memory) and then clears the pointer list
	  //cout << "About to clear Valuelist of size " << this->size() << endl;
	  for (int i=0; i<List.size(); i++) delete List[i];
	  List.clear();
	}
  Parameter::Parameter(const Parameter &o)
	{
	  this->List = o.List;
	  this->Name = o.Name;
	  this->validity = o.validity;
	}
  Parameter &Parameter::operator=(const Parameter &o)
	{

	  this->clear();
	  this->List = o.List;
	  this->Name = o.Name;
	  this->validity = o.validity;

	  return *this;
	}
  
  void ParameterList::clear()
    {
  	  // This function deletes all Value structures on the pointer list (cleans
	  // up all memory) and then clears the pointer list
	  for (int i=0; i<List.size(); i++) delete List[i];
	  List.clear();
	}
  std::vector<Value *> ValueList::deepCopyList(const std::vector<Value *>& origList)
	{
	  // This function returns a a deep copy of a vector(Value *> list
	  std::vector<Value *> newList;
	  for (int i=0; i<origList.size(); i++)
		newList.push_back(Value::Copy(*origList[i]));
	  return newList;
	}

  std::vector<Parameter *>
	ParameterList::deepCopyList (const std::vector<Parameter *>& origList)
    {
	  //cout << "in ParameterList::deepCopyList" << endl;
	  //cout << "copying mem= " << &origList << " to mem= " << this << endl;
	  // This function returns a a deep copy of a vector(Parameter *> list
	  std::vector<Parameter *> newList;
	  for (int i=0; i<origList.size(); i++)
		newList.push_back(new Parameter(*origList[i]));
	  return newList;
	}	
  
  void Parameter::print()
	{
	  cout << "  Parameter: mem = " << this << endl;
	  cout << "    Name=";  Name.print();
	  cout << "    ValueList: ";  List.print();
	}
  void ParameterList::print()
	{
	  cout << "ParameterList:  mem = " << this << " size="
		   << this->size() << endl;
	  cout << "List contains:"<<endl;
	  for (int i = 0; i< this->size() ; i++) {
		cout << "[" << i << "] ";
		List[i]->print();
	  }
	}

  Parameter &ParameterList::getParameter(const char *str)
    {
	  // Returns the first occurrence of a parameter in the List whose Name
	  // value equals "str".  If such a parameter is not found, return the
	  // ParameterList's static NullParameter object.  Just in case something
	  // fishy has gone on, we reset the NullParameter object's validity to
	  // false before beginning.
	  NullParameter.setValid(false);
	  //std::string s(str);
	  Jstring s(str);
	  for (int i = 0; i< List.size(); i++)
		{
		  if (s == List[i]->getName()) return *List[i];
		}
	  //	  cout << "about to return null param from getParameter\n";
	  return NullParameter;
	}
  
  ParameterFile::ParameterFile(const ParameterFile &orig)
  {
	//cout << "IN PARAMETERFILE COPY CONSTRUCTOR\n";
	//cout << "  copying mem= " << &orig << " to mem =" << this <<endl;
	  // Point List pointer to a new ParameterList object that is a deep
	  // copy of the original ParameterList.
	  this->parameterlist_ptr = new ParameterList(*(orig.parameterlist_ptr));
	  if (this->parameterlist_ptr == 0)
		{
		  cerr << "ParameterFile: failed to allocate memory" << endl;
		  exit(3);
		}
	  this->validity = orig.validity;
	}
  
  ParameterFile &ParameterFile::operator=(const ParameterFile &orig)
	{
	  if (parameterlist_ptr)
		{
		  this->parameterlist_ptr->clear(); // Clean up current memory, clear list
	      *(this->parameterlist_ptr) = *(orig.parameterlist_ptr);
		}
	  else  parameterlist_ptr = new ParameterList(*(orig.parameterlist_ptr));
	  this->validity = orig.validity;
	  return *this;
	}
  Parameter &ParameterFile::operator[](ptrdiff_t i)
	{
	  NullParameter.setValid(false);
	  if (parameterlist_ptr) return (*parameterlist_ptr)[i];
	  else return NullParameter;
	}
  Parameter &ParameterFile::operator[](const char *c)
	{
	  NullParameter.setValid(false);
	  if (parameterlist_ptr) return (*parameterlist_ptr)[c];
	  else
		{
		  return NullParameter;
		}
	}
  
  
  bool set(int &i, const Value &v)
	{
	  if ((&v)->getType() != integer_t) return false;
	  i = (int) ((IntegerVal*)(&v))->getValue();
	  return true;
	}
  bool set(long int &i, const Value &v)
	{
	  if ((&v)->getType() != integer_t) return false;
	  i = (long int) ((IntegerVal*)(&v))->getValue();
	  return true;
	}
  bool set(float &i, const Value &v)
	{
	  if ((&v)->getType() != decimal_t) return false;
	  i = (float) ((DecimalVal*)(&v))->getValue();
	  return true;
	}
  //Copies 1st chars into the memory pointed to by str
  bool set(char *str, const Value &v)
	{
	  if ((&v)->getType() != string_t) return false;
	  strcpy(str, ((StringVal *)(&v))->getValue());
	  return true;
	}
  //  bool set(char *str, const Value &v, const int len)
  //	{
	  //std::string temp;
  //	  Jstring temp;
  //	  if ((&v)->getType() != string_t) return false;
  //	  temp = ((StringVal *)(&v))->getString();

  //	}
	  
  //  bool set(std::string &s, const Value&v)
  //	{
  //	  if ((&v)->getType() != string_t) return false;
  //	  s = ((StringVal *)(&v))->getString();
  //	  return true;
  //	}
  
}
