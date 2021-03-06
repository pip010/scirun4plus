/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMetaBlobConverter.txx,v $
  Language:  C++
  Date:      $Date: 2008-01-07 21:48:41 $
  Version:   $Revision: 1.8 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMetaBlobConverter_txx
#define __itkMetaBlobConverter_txx

#include "itkMetaBlobConverter.h"

namespace itk  
{

/** Constructor */ 
template <unsigned int NDimensions>
MetaBlobConverter<NDimensions>
::MetaBlobConverter()
{
  
}

/** Convert a metaBlob into an Blob SpatialObject  */
template <unsigned int NDimensions>
typename MetaBlobConverter<NDimensions>::SpatialObjectPointer
MetaBlobConverter<NDimensions>
::MetaBlobToBlobSpatialObject(MetaBlob * Blob)
{ 

  typedef itk::BlobSpatialObject<NDimensions> BlobSpatialObjectType;
  typename BlobSpatialObjectType::Pointer blob = BlobSpatialObjectType::New();
  
  //typedef BlobSpatialObjectType::VectorType VectorType;
  typedef vnl_vector<double> VectorType;

  unsigned int ndims = Blob->NDims();
  double spacing[NDimensions];
  for(unsigned int ii=0;ii<ndims;ii++)
    {
    spacing[ii]=Blob->ElementSpacing()[ii];
    }

  blob->GetIndexToObjectTransform()->SetScaleComponent(spacing);
  blob->GetProperty()->SetName(Blob->Name());
  blob->SetId(Blob->ID());
  blob->SetParentId(Blob->ParentID());
  blob->GetProperty()->SetRed(Blob->Color()[0]);
  blob->GetProperty()->SetGreen(Blob->Color()[1]);
  blob->GetProperty()->SetBlue(Blob->Color()[2]);
  blob->GetProperty()->SetAlpha(Blob->Color()[3]);

  typedef itk::SpatialObjectPoint<NDimensions> BlobPointType;
  typedef BlobPointType*                       BlobPointPointer;

  
  typedef MetaBlob::PointListType ListType;
  ListType::iterator it2 = Blob->GetPoints().begin();
    
  vnl_vector<double> v(ndims);
  
  for(unsigned int identifier=0;identifier< Blob->GetPoints().size();identifier++)
    {
    BlobPointType pnt;
    
    typedef typename BlobSpatialObjectType::PointType PointType;
    PointType point;

    for(unsigned int ii=0;ii<ndims;ii++)
      {
      point[ii]=(*it2)->m_X[ii];
      }

    pnt.SetPosition(point);

    pnt.SetRed((*it2)->m_Color[0]);
    pnt.SetGreen((*it2)->m_Color[1]);
    pnt.SetBlue((*it2)->m_Color[2]);
    pnt.SetAlpha((*it2)->m_Color[3]);

    blob->GetPoints().push_back(pnt);
    it2++;
    }
 
  return blob;
}

/** Convert an Blob SpatialObject into a metaBlob */
template <unsigned int NDimensions>
MetaBlob*
MetaBlobConverter<NDimensions>
::BlobSpatialObjectToMetaBlob(SpatialObjectType * spatialObject)
{ 
  MetaBlob* Blob = new MetaBlob(NDimensions);

  // fill in the Blob information
   
  typename SpatialObjectType::PointListType::const_iterator i;
  for(i = dynamic_cast<SpatialObjectType*>(spatialObject)->GetPoints().begin();
      i != dynamic_cast<SpatialObjectType*>(spatialObject)->GetPoints().end(); 
      i++)
    {
    BlobPnt* pnt = new BlobPnt(NDimensions);

    for(unsigned int d=0;d<NDimensions;d++)
      {
      pnt->m_X[d]=(*i).GetPosition()[d];
      }
     
    pnt->m_Color[0] = (*i).GetRed();
    pnt->m_Color[1] = (*i).GetGreen();
    pnt->m_Color[2] = (*i).GetBlue();
    pnt->m_Color[3] = (*i).GetAlpha();

    Blob->GetPoints().push_back(pnt); 
    }
    
  if(NDimensions == 2)
    {
    Blob->PointDim("x y red green blue alpha");
    }
  else
    {
    Blob->PointDim("x y z red green blue alpha");
    }

  float color[4];
  for(unsigned int ii=0;ii<4;ii++)
    {
    color[ii]=spatialObject->GetProperty()->GetColor()[ii];
    }

  Blob->Color(color);
  Blob->ID( spatialObject->GetId());
  if(spatialObject->GetParent())
    {
    Blob->ParentID(spatialObject->GetParent()->GetId());
    }
  Blob->NPoints(Blob->GetPoints().size());

  
  for(unsigned int ii=0;ii<NDimensions;ii++)
    {
    Blob->ElementSpacing(ii, spatialObject->GetIndexToObjectTransform()
                                         ->GetScaleComponent()[ii]);
    }

  return Blob;
}


/** Read a meta file give the type */
template <unsigned int NDimensions>
typename MetaBlobConverter<NDimensions>::SpatialObjectPointer
MetaBlobConverter<NDimensions>
::ReadMeta(const char* name)
{
  SpatialObjectPointer spatialObject;
  MetaBlob* Blob = new MetaBlob();
  Blob->Read(name);
  spatialObject = MetaBlobToBlobSpatialObject(Blob);

  return spatialObject;
}


/** Write a meta Blob file */
template <unsigned int NDimensions>
bool
MetaBlobConverter<NDimensions>
::WriteMeta(SpatialObjectType* spatialObject,const char* name)
{
  MetaBlob* Blob = BlobSpatialObjectToMetaBlob(spatialObject);
  Blob->BinaryData(true);
  Blob->Write(name);
  return true;
}

} // end namespace itk 

#endif
