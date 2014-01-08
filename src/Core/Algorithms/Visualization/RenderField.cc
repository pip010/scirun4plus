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

//    File   : RenderField.cc
//    Author : Martin Cole
//    Date   : Tue May 22 10:57:12 2001

#include <Core/Algorithms/Visualization/RenderField.h>
#include <Core/Thread/Time.h>

namespace SCIRun {



RenderFieldBase::RenderFieldBase() :
  UsedWithLockingHandle<Mutex>("RenderFieldBase"),
  node_switch_(0),
  edge_switch_(0),
  face_switch_(0)
{}

RenderFieldBase::~RenderFieldBase()
{}

GeomHandle
RenderFieldV::render_mesh_nodes(FieldHandle field_handle, 
				const std::string &display_mode,
				double scale,
				int resolution,
				unsigned int render_state)
{
  return render_nodes(field_handle, display_mode,
		      scale, resolution, render_state );
}


GeomHandle
RenderFieldV::render_mesh_edges(FieldHandle field_handle,
				const std::string &display_mode,
				double scale,
				int resolution,
				unsigned int render_state,
				unsigned int approx_div)
{
  VField *fld = field_handle->vfield();
  VMesh  *mesh = field_handle->vmesh();
  bool do_linear = ( fld->basis_order() < 2 &&
		     mesh->basis_order() < 2 &&
		     approx_div == 1);

  if (do_linear) 
  {
    return render_edges_linear(field_handle, display_mode,
			       scale, resolution, render_state);
  } 
  else 
  {
    return render_edges(field_handle, display_mode,
			scale, resolution, render_state, approx_div);
  }
}


GeomHandle
RenderFieldV::render_mesh_faces(FieldHandle field_handle, 
				ColorMapHandle colormap_handle,
				unsigned int render_state,
				unsigned int approx_div)
{
  VField *fld = field_handle->vfield();
  VMesh  *mesh = field_handle->vmesh();
  
  if(mesh->is_regularmesh() && mesh->is_surface() &&
     get_flag(render_state, USE_TEXTURE))
  {
    return render_faces_texture(field_handle, colormap_handle, render_state);
  }

  else
  {
    bool do_linear = (fld->basis_order() < 2 && mesh->basis_order() < 2 &&
		      approx_div == 1);
    
    if (do_linear)
      return render_faces_linear(field_handle, render_state );
    else
      return render_faces(field_handle, render_state, approx_div);
  }
}


GeomHandle
RenderFieldV::render_nodes(FieldHandle field_handle, 
			   const std::string &display_mode,
			   double scale,
			   int resolution,
			   unsigned int render_state )
{
  VField *fld = field_handle->vfield();
  VMesh  *mesh = field_handle->vmesh();

  GeomPoints    *points;
  GeomGlyphBase *glyphs;
  GeomHandle data_switch;

  // 0 Points 1 Spheres
  bool points_p  = (display_mode != "Spheres");
  bool spheres_p = (display_mode == "Spheres");

  if (points_p) // Points
  {
    if (get_flag(render_state, USE_TRANSPARENCY))
      points = new GeomTranspPoints();
    else
      points = new GeomPoints();

    data_switch = new GeomDL(points);

    points->setPointSize(resolution/5.0);
  }
  else if (spheres_p) // Spheres
  {
    if (get_flag(render_state, USE_TRANSPARENCY))
      glyphs = new GeomTranspGlyph();
    else
      glyphs = new GeomGlyph();

    data_switch = new GeomDL(glyphs->getObj());
  }

  double sval;
  Vector vval;
  Tensor tval;

  unsigned int color_scheme = 0;
  double scol;
  MaterialHandle vcol(0);

  if (fld->basis_order() < 0 ||
      (fld->basis_order() == 0 && mesh->dimensionality() != 0) ||
      get_flag(render_state, USE_DEFAULT_COLOR))
  {
    color_scheme = 0; // Default color
  }
  else if (get_flag(render_state, USE_COLORMAP))
  {
    // Color std::map lookup using either a scalar value
    // the vector magnitude, the tensor.
    color_scheme = 1;
  }
  else // if (fld->basis_order() >= 1)
  {
    color_scheme = 2; // Value become RGB
    vcol = new Material(Color(1.0, 1.0, 1.0));
    if (get_flag(render_state, USE_TRANSPARENCY))
      vcol->transparency = 0.75;
    else
      vcol->transparency = 1.0;
  }

  mesh->synchronize(Mesh::NODES_E);
  
  VMesh::Node::iterator iter, iter_end;
  mesh->begin(iter);
  mesh->end(iter_end);

  while (iter != iter_end) 
  {
    Point p;
    mesh->get_point(p, *iter);

    if( color_scheme ) 
    {
      if (fld->is_scalar())
      { 
        fld->get_value(sval,*iter); 
        value_to_color( color_scheme, sval, scol, vcol );
      }
      if (fld->is_vector())
      { 
        fld->get_value(vval,*iter); 
        value_to_color( color_scheme, vval, scol, vcol );
      }
      if (fld->is_tensor())
      { 
        fld->get_value(tval,*iter); 
        value_to_color( color_scheme, tval, scol, vcol );
      }
    }

    if (points_p)
    {  
      if (color_scheme == 0)
        points->add(p);
      else if (color_scheme == 1)
        points->add(p, scol);
      else //if (color_scheme == 2)
        points->add(p, vcol);
    }
    else if (spheres_p)
    {
      if (color_scheme == 0)
        glyphs->add_sphere(p, scale, resolution, resolution);
      else if (color_scheme == 1)
        glyphs->add_sphere(p, scale, scol, resolution, resolution);
      else //if (color_scheme == 2)
        glyphs->add_sphere(p, scale, vcol, resolution, resolution);
    }

    ++iter;
  }

  return data_switch;
}

void
RenderFieldV::add_edge_geom(GeomLines *lines,
			    GeomGlyphBase *glyphs,
			    const Point &p0, const Point &p1,
			    const double radius,
			    const int resolution,
			    bool cylinders_p,
			    unsigned int color_scheme,
			    double scol0,
			    double scol1,
			    MaterialHandle& vcol0,
			    MaterialHandle& vcol1)
{
  Vector vec = (Vector) (p1 - p0);

  if (color_scheme == 0)
  {
    if (cylinders_p)
      glyphs->add_cylinder(p0, vec, radius, radius, vec.length(), resolution );
    else
      lines->add(p0, p1);
  }
  else if (color_scheme == 1)
  {
    if (cylinders_p)
      glyphs->add_cylinder(p0, vec, radius, radius, vec.length(),
			   scol0, scol1, resolution);
    else
      lines->add(p0, scol0, p1, scol1);
  }
  else if (color_scheme == 2)
  {
    if ( cylinders_p )
      glyphs->add_cylinder(p0, vec, radius, radius, vec.length(),
			   vcol0, vcol1, resolution);
    else
      lines->add(p0, vcol0, p1, vcol1);
  }
}


GeomHandle
RenderFieldV::render_edges(FieldHandle field_handle,
			   const std::string &edge_display_mode,
			   double scale,
			   int resolution,
			   unsigned int render_state,
			   unsigned div)
{
  VField *fld = field_handle->vfield();
  VMesh  *mesh = field_handle->vmesh();

  GeomLines*     lines;
  GeomGlyphBase* glyphs;
  GeomHandle data_switch;

  // 0 Lines 1 Cylinders
  bool lines_p = (edge_display_mode == "Lines");
  bool cylinders_p = (edge_display_mode == "Cylinders");
  
  if (lines_p) // Lines
  {
    if (get_flag(render_state, USE_TRANSPARENCY))
      lines = new GeomTranspLines;
    else
      lines = new GeomLines;

    data_switch = new GeomDL(lines);

    lines->setLineWidth(resolution/5.0);
  }
  else if (cylinders_p) // Cylinders
  {
    if (get_flag(render_state, USE_TRANSPARENCY))
      glyphs = new GeomTranspGlyph();
    else
      glyphs = new GeomGlyph();

    data_switch = new GeomDL(glyphs->getObj());
  }

  double sval0, sval1;
  Vector vval0, vval1;
  Tensor tval0, tval1;

  unsigned int color_scheme = 0;
  double scol0 = 0.0, scol1 = 0.0;
  MaterialHandle vcol0(0);
  MaterialHandle vcol1(0);

  if (fld->basis_order() < 0 ||
      (fld->basis_order() == 0 && mesh->dimensionality() != 1) ||
      get_flag(render_state, USE_DEFAULT_COLOR))
  {
    color_scheme = 0; // Default color
  }
  else if (get_flag(render_state, USE_COLORMAP))
  {
    // Color map lookup using either a scalar value
    // the vector magnitude, the tensor.
    color_scheme = 1;
  }
  else // if (fld->basis_order() >= 1)
  {
    color_scheme = 2; // Values become RGB
    vcol0 = new Material(Color(1.0, 1.0, 1.0));
    vcol1 = new Material(Color(1.0, 1.0, 1.0));
    if (get_flag(render_state, USE_TRANSPARENCY))
      vcol0->transparency = vcol1->transparency = 0.75;
    else
      vcol0->transparency = vcol1->transparency = 1.0;
  }

#if defined(_MSC_VER) || defined(__ECC)
  typedef hash_set<std::string> edge_ht_t;
#else
  typedef hash_set<std::string, str_hasher, equal_str> edge_ht_t;
#endif
  edge_ht_t rendered_edges; 

  mesh->synchronize(Mesh::EDGES_E | Mesh::FACES_E | Mesh::CELLS_E);
  VMesh::Elem::iterator eiter, eiter_end;
  VMesh::Edge::array_type edges;
  
  mesh->begin(eiter);
  mesh->end(eiter_end);

  while (eiter != eiter_end) 
  {  
    mesh->get_edges(edges, *eiter);

    VMesh::Edge::array_type::iterator edge_iter;
    edge_iter = edges.begin();
    int ecount = 0;
    while (edge_iter != edges.end()) 
    {

      VMesh::Node::array_type nodes;
      VMesh::Edge::index_type eidx = *edge_iter++;

      Point cntr;
      mesh->get_center(cntr, eidx);
      std::ostringstream pstr;  
      pstr << setiosflags(std::ios::scientific);
      pstr << std::setprecision(7); 
      pstr << cntr.x() << cntr.y() << cntr.z();
      
      edge_ht_t::const_iterator it = rendered_edges.find(pstr.str());

      if (it != rendered_edges.end()) 
      {
        ++ecount;
        continue;
      } 
      else 
      {
        rendered_edges.insert(pstr.str());
      }
      // following print is useful for debugging edge ordering
      //      cout << "elem: " << *eiter << " count " << ecount 
      //	   << " edge" << eidx << std::endl;
      VMesh::coords_array_type coords;
      mesh->pwl_approx_edge(coords, *eiter, ecount++, div);
      VMesh::coords_array_type::iterator coord_iter = coords.begin();
      do 
      {
        VMesh::coords_type &c0 = *coord_iter++;
        if (coord_iter == coords.end()) break;
        VMesh::coords_type &c1 = *coord_iter;
        Point p0, p1;      

        // get the geometry at the approx.
        mesh->interpolate(p0, c0, *eiter);
        mesh->interpolate(p1, c1, *eiter);

        // get the field variables values at the approx (if they exist)
        if (color_scheme) 
        {
          if (fld->is_scalar())
          {
            fld->interpolate(sval0, c0, *eiter);
            fld->interpolate(sval1, c1, *eiter);

            value_to_color( color_scheme, sval0, scol0, vcol0 );
            value_to_color( color_scheme, sval1, scol1, vcol1 );
          }
          if (fld->is_vector())
          {
            fld->interpolate(vval0, c0, *eiter);
            fld->interpolate(vval1, c1, *eiter);

            value_to_color( color_scheme, vval0, scol0, vcol0 );
            value_to_color( color_scheme, vval1, scol1, vcol1 );
          }
          if (fld->is_tensor())
          {
            fld->interpolate(tval0, c0, *eiter);
            fld->interpolate(tval1, c1, *eiter);

            value_to_color( color_scheme, tval0, scol0, vcol0 );
            value_to_color( color_scheme, tval1, scol1, vcol1 );
          }
        }

        add_edge_geom(lines, glyphs, p0, p1,
		      scale, resolution, cylinders_p,
		      color_scheme, scol0, scol1, vcol0, vcol1);
	
      } 
      while (coords.size() > 1 && coord_iter != coords.end()); 
    }
    
    ++eiter;
  }
  
  return data_switch;
}


GeomHandle
RenderFieldV::render_edges_linear(FieldHandle field_handle,
				  const std::string &display_mode,
				  double scale,
				  int resolution,
				  unsigned int render_state) 
{
  VField *fld = field_handle->vfield();
  VMesh  *mesh = field_handle->vmesh();

  GeomLines*     lines = 0;
  GeomGlyphBase* glyphs = 0;
  GeomHandle data_switch;

  // 0 Lines 1 Cylinders
  bool lines_p     = (display_mode == "Lines");
  bool cylinders_p = (display_mode == "Cylinders");

  if (lines_p) // Lines
  {
    if (get_flag(render_state, USE_TRANSPARENCY))
    {
      lines = new GeomTranspLines;
      data_switch = lines;
    }
    else
    {
      lines = new GeomLines;
      data_switch = new GeomDL(lines);
    }

    lines->setLineWidth(resolution/5.0);
  }
  else if (cylinders_p) // Cylinders
  {
    if (get_flag(render_state, USE_TRANSPARENCY))
    {
      glyphs = new GeomTranspGlyph();
      data_switch = glyphs->getObj();
    }
    else
    {
      glyphs = new GeomGlyph();
      data_switch = new GeomDL(glyphs->getObj());
    }
  }

  double sval0, sval1;
  Vector vval0, vval1;
  Tensor tval0, tval1;

  unsigned int color_scheme = 0;
  double scol0 = 0.0, scol1 = 0.0;
  MaterialHandle vcol0(0);
  MaterialHandle vcol1(0);

  if (fld->basis_order() < 0 ||
      (fld->basis_order() == 0 && mesh->dimensionality() != 1) ||
      get_flag(render_state, USE_DEFAULT_COLOR))
  {
    color_scheme = 0; // Default color
  }
  else if (get_flag(render_state, USE_COLORMAP))
  {
    // Color map lookup using either a scalar value
    // the vector magnitude, the tensor.
    color_scheme = 1;
  }
  else // if (fld->basis_order() >= 1)
  {
    color_scheme = 2; // Vector values become RGB

    vcol0 = new Material(Color(1.0, 1.0, 1.0));
    vcol1 = new Material(Color(1.0, 1.0, 1.0));

    if (get_flag(render_state, USE_TRANSPARENCY))
      vcol0->transparency = vcol1->transparency = 0.75;
    else
      vcol0->transparency = vcol1->transparency = 1.0;
  }

  mesh->synchronize(Mesh::EDGES_E);
  VMesh::Edge::iterator eiter, eiter_end;

  mesh->begin(eiter);
  mesh->end(eiter_end);
 
  while (eiter != eiter_end) 
  {  
    VMesh::Node::array_type nodes;
    mesh->get_nodes(nodes, *eiter);
      
    Point p0, p1;
    mesh->get_point(p0, nodes[0]);
    mesh->get_point(p1, nodes[1]);

    if( color_scheme )
    {
      if (fld->is_scalar())
      {
        if (fld->basis_order() == 1)
        {
          fld->get_value(sval0, nodes[0]);
          fld->get_value(sval1, nodes[1]);
        }
        else //if (mesh->dimensionality() == 1)
        {
          fld->get_value(sval0, *eiter);
          
          sval1 = sval0;
        }
        
        value_to_color( color_scheme, sval0, scol0, vcol0 );
        value_to_color( color_scheme, sval1, scol1, vcol1 );
      }
      else if (fld->is_vector())
      {
        if (fld->basis_order() == 1)
        {
          fld->get_value(vval0, nodes[0]);
          fld->get_value(vval1, nodes[1]);
        }
        else //if (mesh->dimensionality() == 1)
        {
          fld->get_value(vval0, *eiter);
          
          vval1 = vval0;
        }

        value_to_color( color_scheme, vval0, scol0, vcol0 );
        value_to_color( color_scheme, vval1, scol1, vcol1 );
      }
      else if (fld->is_tensor())
      {
        if (fld->basis_order() == 1)
        {
          fld->get_value(tval0, nodes[0]);
          fld->get_value(tval1, nodes[1]);
        }
        else //if (mesh->dimensionality() == 1)
        {
          fld->get_value(tval0, *eiter);
          
          tval1 = tval0;
        }
        
        value_to_color( color_scheme, tval0, scol0, vcol0 );
        value_to_color( color_scheme, tval1, scol1, vcol1 );
      }
    }
    
    add_edge_geom(lines, glyphs, p0, p1,
		  scale, resolution, cylinders_p,
		  color_scheme, scol0, scol1, vcol0, vcol1);

    ++eiter;
  }

  return data_switch;
}


void
RenderFieldV::add_face_geom(GeomFastTriangles *faces,
			    GeomFastQuads *qfaces,
			    const std::vector<Point>  &points,
			    const std::vector<Vector> &normals,
			    bool with_normals,
			    unsigned int color_scheme,
			    std::vector<double> &scols,
			    std::vector<MaterialHandle> &vcols )
{
  if (color_scheme == 0)
  {
    if (points.size() == 4)
    {
      if (with_normals)
      {
        qfaces->add(points[0], normals[0],
		    points[1], normals[1],
		    points[2], normals[2],
		    points[3], normals[3]);
      }
      else
      {
        qfaces->add(points[0], points[1], points[2], points[3]);
      }
    }
    else
    {
      for (size_t i = 2; i < points.size(); i++)
      {
        if (with_normals)
        {
          faces->add(points[0],   normals[0],
		     points[i-1], normals[i-1],
		     points[i],   normals[i]);
        }
        else
        {
          faces->add(points[0], points[i-1], points[i]);
        }
      }
    }    
  }
  else if (color_scheme == 1)
  {
    if (points.size() == 4)
    {
      if (with_normals)
      {
        qfaces->add(points[0], normals[0], scols[0],
		    points[1], normals[1], scols[1],
		    points[2], normals[2], scols[2],
		    points[3], normals[3], scols[3]);
      }
      else
      {
        qfaces->add(points[0], scols[0],
		    points[1], scols[1],
		    points[2], scols[2],
		    points[3], scols[3]);
      }
    }
    else
    {
      for (size_t i = 2; i < points.size(); i++)
      {
        if (with_normals)
        {
          faces->add(points[0],   normals[0],   scols[0],
		     points[i-1], normals[i-1], scols[i-1],
		     points[i],   normals[i],   scols[i]);
        }
        else
        {
          faces->add(points[0],   scols[0],
		     points[i-1], scols[i-1],
		     points[i],   scols[i]);
        }
      }
    }
  }
  else if (color_scheme == 2)
  {
    if (points.size() == 4)
    {
      if (with_normals)
      {
        qfaces->add(points[0], normals[0], vcols[0],
		    points[1], normals[1], vcols[1],
		    points[2], normals[2], vcols[2],
		    points[3], normals[3], vcols[3]);
      }
      else
      {
        qfaces->add(points[0], vcols[0],
		    points[1], vcols[1],
		    points[2], vcols[2],
		    points[3], vcols[3]);
      }
    }
    else
    {
      for (size_t i = 2; i < points.size(); i++)
      {
        if (with_normals)
        {
          faces->add(points[0],   normals[0],   vcols[0],
		     points[i-1], normals[i-1], vcols[i-1],
		     points[i],   normals[i],   vcols[i]);
        }
        else
        {
          faces->add(points[0],   vcols[0],
		     points[i-1], vcols[i-1],
		     points[i],   vcols[i]);
        }
      }
    }
  }
}


GeomHandle 
RenderFieldV::render_faces(FieldHandle field_handle,
			   unsigned int render_state,
			   unsigned int div)
{
  VField *fld = field_handle->vfield();
  VMesh  *mesh = field_handle->vmesh();

  const bool with_normals = (get_flag(render_state, USE_NORMALS)
			     && mesh->has_normals());

  GeomFastTriangles* tfaces;
  GeomFastQuads* qfaces;
  GeomGroup *grp = new GeomGroup();
  GeomHandle face_switch = new GeomDL(grp);

  if (get_flag(render_state, USE_TRANSPARENCY))
  {
    tfaces = new GeomTranspTriangles;
    qfaces = new GeomTranspQuads;
  }
  else
  {
    tfaces = new GeomFastTriangles;
    qfaces = new GeomFastQuads;
  }

  grp->add(tfaces);
  grp->add(qfaces);

  unsigned int color_scheme = 0;

  std::vector<double> scols(10);
  std::vector<MaterialHandle> vcols(10, (Material *) 0);

  if (fld->basis_order() < 0 ||
      get_flag(render_state, USE_DEFAULT_COLOR))
  {
    color_scheme = 0; // Default color
  }
  else if (get_flag(render_state, USE_COLORMAP))
  {
    // Color map lookup using either a scalar value
    // the vector magnitude, the tensor.
    color_scheme = 1;
  }
  else // if (fld->basis_order() >= 0)
  {
    color_scheme = 2; // Values become RGB

    for (unsigned int i=0; i<10; i++)
    {
      vcols[i] = new Material(Color(1.0, 1.0, 1.0));

      if (get_flag(render_state, USE_TRANSPARENCY))
        vcols[i]->transparency = 0.75;
      else
        vcols[i]->transparency = 1.0;
    }
  }

#if defined(_MSC_VER) || defined(__ECC)
  typedef hash_set<std::string> face_ht_t;
#else
  typedef hash_set<std::string, str_hasher, equal_str> face_ht_t;
#endif
  face_ht_t rendered_faces; 
  
  mesh->synchronize(Mesh::FACES_E | Mesh::EDGES_E | Mesh::CELLS_E);
  VMesh::Elem::iterator eiter, eiter_end;
  mesh->begin(eiter);
  mesh->end(eiter_end);

  while (eiter != eiter_end) 
  {  
    VMesh::Face::array_type face_indecies;
    mesh->get_faces(face_indecies, *eiter);

    VMesh::Face::array_type::iterator face_iter = face_indecies.begin();
    int fcount = 0;
    while (face_iter != face_indecies.end()) 
    {
      VMesh::Node::array_type nodes;
      VMesh::Face::index_type fidx = *face_iter++;

      Point cntr;
      mesh->get_center(cntr, fidx);
      std::ostringstream pstr;
      pstr << setiosflags(std::ios::scientific);
      pstr << std::setprecision(7); 
      pstr << cntr.x() << cntr.y() << cntr.z();
      
      face_ht_t::const_iterator it = rendered_faces.find(pstr.str());

      if (it != rendered_faces.end()) 
      {
        ++fcount;
        continue;
      } 
      else 
      {
        rendered_faces.insert(pstr.str());
      }

      //coords organized as scanlines of quad/tri strips.
      VMesh::coords_array2_type coords;
      mesh->pwl_approx_face(coords, *eiter, fcount, div);

      const int face_sz = mesh->num_nodes_per_face();

      std::vector<Point> points(face_sz);
      std::vector<Vector> normals(face_sz);
      std::vector<double> svals(face_sz);
      std::vector<Vector> vvals(face_sz);
      std::vector<Tensor> tvals(face_sz);

      VMesh::coords_array2_type::iterator coord_iter = coords.begin();

      // TRI STRIPS
      if (face_sz == 3)
      {
        while (coord_iter != coords.end()) 
        {
          VMesh::coords_array_type &sl = *coord_iter++;
          VMesh::coords_array_type::iterator sliter = sl.begin();

          for (size_t i=0; i<sl.size()-2; i++) 
          {
            VMesh::coords_array_type::iterator it0, it1;
            
            VMesh::coords_type &c0 = !(i%2) ? sl[i] : sl[i+1];
            VMesh::coords_type &c1 = !(i%2) ? sl[i+1] : sl[i];
            VMesh::coords_type &c2 = sl[i+2];
                  
            // get the geometry at the approx.
            mesh->interpolate(points[0], c0, *eiter);
            mesh->interpolate(points[1], c1, *eiter);
            mesh->interpolate(points[2], c2, *eiter);

            if (color_scheme) 
            {
              if (fld->is_scalar())
              {
                // get the field variables values at the approx (if they exist)
                fld->interpolate(svals[0], c0, *eiter);
                fld->interpolate(svals[1], c1, *eiter);
                fld->interpolate(svals[2], c2, *eiter);

                value_to_color( color_scheme, svals[0], scols[0], vcols[0] );
                value_to_color( color_scheme, svals[1], scols[1], vcols[1] );
                value_to_color( color_scheme, svals[2], scols[2], vcols[2] );
              }
              else if (fld->is_vector())
              {
                // get the field variables values at the approx (if they exist)
                fld->interpolate(vvals[0], c0, *eiter);
                fld->interpolate(vvals[1], c1, *eiter);
                fld->interpolate(vvals[2], c2, *eiter);

                value_to_color( color_scheme, vvals[0], scols[0], vcols[0] );
                value_to_color( color_scheme, vvals[1], scols[1], vcols[1] );
                value_to_color( color_scheme, vvals[2], scols[2], vcols[2] );
              }
              else if (fld->is_tensor())
              {      
                // get the field variables values at the approx (if they exist)
                fld->interpolate(tvals[0], c0, *eiter);
                fld->interpolate(tvals[1], c1, *eiter);
                fld->interpolate(tvals[2], c2, *eiter);

                value_to_color( color_scheme, tvals[0], scols[0], vcols[0] );
                value_to_color( color_scheme, tvals[1], scols[1], vcols[1] );
                value_to_color( color_scheme, tvals[2], scols[2], vcols[2] );
              }
            }
  
            if (with_normals) 
            {	      
              mesh->get_normal(normals[0], c0, *eiter, fcount);
              mesh->get_normal(normals[1], c1, *eiter, fcount);
              mesh->get_normal(normals[2], c2, *eiter, fcount);
            }

            add_face_geom(tfaces, qfaces, points, normals, with_normals,
              color_scheme, scols, vcols);                        
          }
        }
      }
	
      // QUADS
      else 
      {
        while (coord_iter != coords.end()) 
        {
          VMesh::coords_array_type &sl = *coord_iter++;
          VMesh::coords_array_type::iterator sliter = sl.begin();

          for (size_t i=0; i<sl.size()-2; i++) 
          {
            VMesh::coords_type &c0 = *sliter++;
            if (sliter == sl.end()) break;
            VMesh::coords_type &c1 = *sliter++;
            if (sliter == sl.end()) break;
            VMesh::coords_type &c2 = *sliter;
            VMesh::coords_type &c3 = *(sliter + 1);

            // get the geometry at the approx.
            mesh->interpolate(points[0], c2, *eiter);
            mesh->interpolate(points[1], c3, *eiter);
            mesh->interpolate(points[2], c1, *eiter);
            mesh->interpolate(points[3], c0, *eiter);

            if (color_scheme) 
            {
            if (fld->is_scalar())
              {
                // get the field variables values at the approx (if they exist)
                fld->interpolate(svals[0], c2, *eiter);
                fld->interpolate(svals[1], c3, *eiter);
                fld->interpolate(svals[2], c1, *eiter);
                fld->interpolate(svals[3], c0, *eiter);

                value_to_color( color_scheme, svals[0], scols[0], vcols[0] );
                value_to_color( color_scheme, svals[1], scols[1], vcols[1] );
                value_to_color( color_scheme, svals[2], scols[2], vcols[2] );
                value_to_color( color_scheme, svals[3], scols[3], vcols[3] );
              }
              else if (fld->is_vector())
              {
                // get the field variables values at the approx (if they exist)
                fld->interpolate(vvals[0], c2, *eiter);
                fld->interpolate(vvals[1], c3, *eiter);
                fld->interpolate(vvals[2], c1, *eiter);
                fld->interpolate(vvals[3], c0, *eiter);

                value_to_color( color_scheme, vvals[0], scols[0], vcols[0] );
                value_to_color( color_scheme, vvals[1], scols[1], vcols[1] );
                value_to_color( color_scheme, vvals[2], scols[2], vcols[2] );
                value_to_color( color_scheme, vvals[3], scols[3], vcols[3] );
              }
              else if (fld->is_tensor())
              {
                // get the field variables values at the approx (if they exist)
                fld->interpolate(tvals[0], c2, *eiter);
                fld->interpolate(tvals[1], c3, *eiter);
                fld->interpolate(tvals[2], c1, *eiter);
                fld->interpolate(tvals[3], c0, *eiter);

                value_to_color( color_scheme, tvals[0], scols[0], vcols[0] );
                value_to_color( color_scheme, tvals[1], scols[1], vcols[1] );
                value_to_color( color_scheme, tvals[2], scols[2], vcols[2] );
                value_to_color( color_scheme, tvals[3], scols[3], vcols[3] );
              }
            }
                   
            if (with_normals) 
            {	      
              mesh->get_normal(normals[0], c2, *eiter, fcount);
              mesh->get_normal(normals[1], c3, *eiter, fcount);
              mesh->get_normal(normals[2], c1, *eiter, fcount);
              mesh->get_normal(normals[3], c0, *eiter, fcount);
            }
            
            add_face_geom(tfaces, qfaces, points, normals, with_normals,
              color_scheme, scols, vcols);
          }
        }
      }
      ++fcount;
    }
    ++eiter;
  }
  return face_switch;
}

GeomHandle 
RenderFieldV::render_faces_linear(FieldHandle field_handle,
				  unsigned int render_state)
{
  VField *fld = field_handle->vfield();
  VMesh  *mesh = field_handle->vmesh();

  const bool with_normals = (get_flag(render_state, USE_NORMALS)
			     && mesh->has_normals());

  GeomFastTriangles* tfaces = 0;
  GeomFastQuads* qfaces = 0;
  GeomFastTrianglesTwoSided* ttfaces = 0;
  GeomFastQuadsTwoSided* tqfaces = 0;

  GeomGroup *grp = new GeomGroup();
  GeomHandle face_switch = new GeomDL(grp);

  unsigned int color_scheme = 0;
  std::vector<double> svals(10);
  std::vector<Vector> vvals(10);
  std::vector<Tensor> tvals(10);

  std::vector<MaterialHandle> vcols(10, (Material *)NULL);
  std::vector<double> scols(10);

  if (fld->basis_order() < 0 ||
      get_flag(render_state, USE_DEFAULT_COLOR))
  {
    color_scheme = 0; // Default color
  }
  else if (get_flag(render_state, USE_COLORMAP))
  {
    // Color map lookup using either a scalar value,
    // the vector magnitude, the tensor.
    color_scheme = 1;
  }
  else // if (fld->basis_order() >= 0)
  {
    color_scheme = 2; // Values become RGB

    for (unsigned int i=0; i<10; i++)
    {
      vcols[i] = new Material(Color(1.0, 1.0, 1.0));

      if (get_flag(render_state, USE_TRANSPARENCY))
        vcols[i]->transparency = 0.75;
      else
        vcols[i]->transparency = 1.0;
    }
  }

  // Special case for cell centered data
  if ((fld->basis_order() == 0) && (mesh->dimensionality() == 3) &&
      (color_scheme > 0))
  {
    if (get_flag(render_state, USE_TRANSPARENCY))
    {
      ttfaces = new GeomTranspTrianglesTwoSided;
      tqfaces = new GeomTranspQuadsTwoSided;
      grp->add(ttfaces);
      grp->add(tqfaces);
    }
    else
    {
      ttfaces = new GeomFastTrianglesTwoSided;
      tqfaces = new GeomFastQuadsTwoSided;
      grp->add(ttfaces);
      grp->add(tqfaces);
    }
  }
  else
  {
    if (get_flag(render_state, USE_TRANSPARENCY))
    {
      tfaces = new GeomTranspTriangles;
      qfaces = new GeomTranspQuads;
      
      grp->add(tfaces);
      grp->add(qfaces);
    }
    else
    {
      tfaces = new GeomFastTriangles;
      qfaces = new GeomFastQuads;
    
      grp->add(tfaces);
      grp->add(qfaces);
    }
  }

  if (with_normals) mesh->synchronize(Mesh::NORMALS_E);

  mesh->synchronize(Mesh::FACES_E);
  VMesh::Face::iterator fiter, fiter_end;
  VMesh::Node::array_type nodes;

  mesh->begin(fiter);
  mesh->end(fiter_end);

  VMesh::Face::size_type f;
  VMesh::Cell::size_type c;

  mesh->size(f);
  mesh->size(c);

  while (fiter != fiter_end) 
  {
    mesh->get_nodes(nodes, *fiter);
 
    std::vector<Point> points(nodes.size());
    std::vector<Vector> normals(nodes.size());

    for (size_t i=0; i<nodes.size(); i++)
      mesh->get_point(points[i], nodes[i]);

    if (with_normals) 
    {
      for (size_t i=0; i<nodes.size(); i++)
        mesh->get_normal(normals[i], nodes[i]);
    }

    // Default color single face no matter the element data.
    if (color_scheme == 0)
    {
      add_face_geom(tfaces, qfaces, points, normals, with_normals,
		    color_scheme, scols, vcols);                        
    }
    // Element data (Cells) so two sided faces.
    else if (fld->basis_order() == 0 && mesh->dimensionality() == 3)
    {
      VMesh::Elem::array_type cells;
      mesh->get_elems(cells, *fiter);
      
      if (fld->is_scalar())
      {
        fld->get_value(svals[0], cells[0]);
        
        if (cells.size() > 1)
          fld->get_value(svals[1], cells[1]);
        else
          svals[1] = svals[0];
        
        value_to_color( color_scheme, svals[0], scols[0], vcols[0] );
        value_to_color( color_scheme, svals[1], scols[1], vcols[1] );
      }
      else if (fld->is_vector())
      {
        fld->get_value(vvals[0], cells[0]);
        
        if (cells.size() > 1)
          fld->get_value(vvals[1], cells[1]);
        else
          svals[1] = svals[0];
        
        value_to_color( color_scheme, vvals[0], scols[0], vcols[0] );
        value_to_color( color_scheme, vvals[1], scols[1], vcols[1] );
      }
      else if (fld->is_tensor())
      {
        fld->get_value(tvals[0], cells[0]);
        
        if (cells.size() > 1)
          fld->get_value(tvals[1], cells[1]);
        else
          svals[1] = svals[0];
        
        value_to_color( color_scheme, tvals[0], scols[0], vcols[0] );
        value_to_color( color_scheme, tvals[1], scols[1], vcols[1] );
      }

      if (color_scheme == 1)
      {
        if (nodes.size() == 4)
        {
          if (with_normals)
            tqfaces->add(points[0], normals[0], scols[0], scols[1],
                         points[1], normals[1], scols[0], scols[1],
                         points[2], normals[2], scols[0], scols[1],
                         points[3], normals[3], scols[0], scols[1]);
          else
            tqfaces->add(points[0], scols[0], scols[1],
                         points[1], scols[0], scols[1],
                         points[2], scols[0], scols[1],
                         points[3], scols[0], scols[1]);
        }
        else
        {
          for (size_t i=2; i<nodes.size(); i++)
          {
            if (with_normals)
              ttfaces->add(points[0],   normals[0],   scols[0], scols[1],
                           points[i-1], normals[i-1], scols[0], scols[1],
                           points[i],   normals[i],   scols[0], scols[1]);
            else
              ttfaces->add(points[0],   scols[0], scols[1],
                           points[i-1], scols[0], scols[1],
                           points[i],   scols[0], scols[1]);
          }
        }
      }
      else //if (color_scheme == 2)
      {
        if (nodes.size() == 4)
        {
          if (with_normals)
            tqfaces->add(points[0], normals[0], vcols[0], vcols[1],
                         points[1], normals[1], vcols[0], vcols[1],
                         points[2], normals[2], vcols[0], vcols[1],
                         points[3], normals[3], vcols[0], vcols[1]);
          else
            tqfaces->add(points[0], vcols[0], vcols[1],
                         points[1], vcols[0], vcols[1],
                         points[2], vcols[0], vcols[1],
                         points[3], vcols[0], vcols[1]);
        }
        else
        {
          for (size_t i=2; i<nodes.size(); i++)
          {
            if (with_normals)
              ttfaces->add(points[0],   normals[0],   vcols[0], vcols[1],
                           points[i-1], normals[i-1], vcols[0], vcols[1],
                           points[i],   normals[i],   vcols[0], vcols[1]);
            else
              ttfaces->add(points[0],   vcols[0], vcols[1],
                           points[i-1], vcols[0], vcols[1],
                           points[i],   vcols[0], vcols[1]);
          }
        }
      }
    }
    
    // Element data (faces)
    else if (fld->basis_order() == 0 && mesh->dimensionality() == 2)
    {
      if (fld->is_scalar())
      {
        fld->get_value(svals[0], *fiter);
        value_to_color( color_scheme, svals[0], scols[0], vcols[0] );
      }
      else if (fld->is_vector())
      {
        fld->get_value(vvals[0], *fiter);
        value_to_color( color_scheme, vvals[0], scols[0], vcols[0] );
      }
      else if (fld->is_tensor())
      {
        fld->get_value(tvals[0], *fiter);
        value_to_color( color_scheme, tvals[0], scols[0], vcols[0] );
      }

      // Same color at all corners.
      for(size_t i=0; i<nodes.size(); ++i)
      {
        scols[i] = scols[0];
        vcols[i] = vcols[0];
      }
      
      add_face_geom(tfaces, qfaces, points, normals, with_normals,
		    color_scheme, scols, vcols);                        
    }

    // Data at nodes
    else if (fld->basis_order() == 1)
    {
      if (fld->is_scalar())
      {
        for (size_t i=0; i<nodes.size(); i++)
        {
          fld->get_value(svals[i], nodes[i]);
          value_to_color( color_scheme, svals[i], scols[i], vcols[i] );
        }
      }
      else if (fld->is_vector())
      {
        for (size_t i=0; i<nodes.size(); i++)
        {
          fld->get_value(vvals[i], nodes[i]);
          value_to_color( color_scheme, vvals[i], scols[i], vcols[i] );
        }
      }      
      else if (fld->is_tensor())
      {
        for (size_t i=0; i<nodes.size(); i++)
        {
          fld->get_value(tvals[i], nodes[i]);
          value_to_color( color_scheme, tvals[i], scols[i], vcols[i] );
        }
      }
      
      add_face_geom(tfaces, qfaces, points, normals, with_normals,
		    color_scheme, scols, vcols);
    }

    ++fiter;     
  }

  return face_switch;
}


GeomHandle
RenderFieldV::render_faces_texture(FieldHandle field_handle,
					ColorMapHandle colormap_handle,
					unsigned int render_state)
{
  VField *fld = field_handle->vfield();
  VMesh  *mesh = field_handle->vmesh();

  GeomHandle texture_face;
  float tex_coords[8];
  float pos_coords[12];
  const int colorbytes = 4;

  GeomTexRectangle *tr = new GeomTexRectangle();
  texture_face = tr;

  VMesh::size_type ni = mesh->get_ni();
  VMesh::size_type nj = mesh->get_nj();

  // Set up the texture parameters, power of 2 dimensions.
  int width = Pow2(ni);
  int height = Pow2(nj);

  // Use for the texture coordinates 
  double tmin_x, tmax_x, tmin_y, tmax_y;

  // Create texture array 
  unsigned char * texture = new unsigned char[colorbytes*width*height];

  //***************************************************
  // we need to find the corners of the square in space
  // use the node indices to grab the corner points
  
  Point p1, p2, p3, p4;
  mesh->get_center(p1, VMesh::Node::index_type(0));
  pos_coords[0] = p1.x();
  pos_coords[1] = p1.y();
  pos_coords[2] = p1.z();

  mesh->get_center(p2, VMesh::Node::index_type(ni-1));
  pos_coords[3] = p2.x();
  pos_coords[4] = p2.y();
  pos_coords[5] = p2.z();

  mesh->get_center(p3, VMesh::Node::index_type(ni*(nj-1)));
  pos_coords[6] = p3.x();
  pos_coords[7] = p3.y();
  pos_coords[8] = p3.z();

  mesh->get_center(p4, VMesh::Node::index_type((ni*nj) -1));
  pos_coords[9] = p4.x();
  pos_coords[10] = p4.y();
  pos_coords[11] = p4.z();

  double dval;

  double sval;
  Vector vval;
  Tensor tval;

  if (fld->basis_order() == 1)
  {
    tr->interpolate(true);

    tmin_x = 0.5/static_cast<double>(width);
    tmax_x = (static_cast<double>(ni)- 0.5)/static_cast<double>(width);
    tmin_y = 0.5/static_cast<double>(height);
    tmax_y = (static_cast<double>(nj)-0.5)/static_cast<double>(height);

    tex_coords[0] = tmin_x; tex_coords[1] = tmin_y;
    tex_coords[2] = tmax_x; tex_coords[3] = tmin_y;
    tex_coords[4] = tmax_x; tex_coords[5] = tmax_y;
    tex_coords[6] = tmin_x; tex_coords[7] = tmax_y;

    VMesh::Node::array_type nodes;

    VMesh::Node::index_type idx = 0;
    for (index_type i=0;i<ni;i++)
    {
      for (index_type j=0;j<nj;j++)
      {
        mesh->to_index(idx,i,j);
        
        // Convert data values to double.
        if (fld->is_scalar())
        {
          fld->get_value(sval, idx);
          to_double(sval, dval);
        }
        if (fld->is_vector())
        {
          fld->get_value(vval, idx);
          to_double(vval, dval);
        }
        if (fld->is_tensor())
        {
          fld->get_value(tval, idx);
          to_double(tval, dval);
        }
                
        // Compute the ColorMap index and retreive the color.
        const double cmin = colormap_handle->getMin();
        const double cmax = colormap_handle->getMax();
        const double index = Clamp((dval - cmin)/(cmax - cmin), 0.0, 1.0);
        float r,g,b,a;
        colormap_handle->get_color(index, r, g, b, a);
        const float zro = 0.0;
        const float one = 1.0;

        // Fill the texture.
        index_type k = (i+width*j)*colorbytes;
        texture[k]   = (unsigned char)(Clamp(r, zro, one)*255);
        texture[k+1] = (unsigned char)(Clamp(g, zro, one)*255);
        texture[k+2] = (unsigned char)(Clamp(b, zro, one)*255);
        texture[k+3] = (unsigned char)(Clamp(a, zro, one)*255);
      }
    }
  }
  else if (fld->basis_order() == 0)
  {
    tr->interpolate( false );
    tmin_x = 0.0;
    tmax_x = (ni-1)/static_cast<double>(width);
    tmin_y = 0.0;
    tmax_y = (nj-1)/static_cast<double>(height);
    tex_coords[0] = tmin_x; tex_coords[1] = tmin_y;
    tex_coords[2] = tmax_x; tex_coords[3] = tmin_y;
    tex_coords[4] = tmax_x; tex_coords[5] = tmax_y;
    tex_coords[6] = tmin_x; tex_coords[7] = tmax_y;

    VMesh::Elem::index_type idx = 0;
    for (index_type i=0;i<ni;i++)
    {
      for (index_type j=0;j<nj;j++)
      {
        mesh->to_index(idx,i,j);
        // Convert data values to double.
        if (fld->is_scalar())
        {
          fld->value(sval, idx);
          to_double(sval, dval);
        }
        if (fld->is_vector())
        {
          fld->value(vval, idx);
          to_double(vval, dval);
        }
        if (fld->is_tensor())
        {
          fld->value(tval, idx);
          to_double(tval, dval);
        }
     
        // Compute the ColorMap index and retreive the color.
        const double cmin = colormap_handle->getMin();
        const double cmax = colormap_handle->getMax();
        const double index = Clamp((dval - cmin)/(cmax - cmin), 0.0, 1.0);
        float r,g,b,a;
        colormap_handle->get_color(index, r, g, b, a);
        const float zro = 0.0;
        const float one = 1.0;

        // Fill the texture.
        index_type k = (i+width*j)*colorbytes;
        texture[k]   = (unsigned char)(Clamp(r, zro, one)*255);
        texture[k+1] = (unsigned char)(Clamp(g, zro, one)*255);
        texture[k+2] = (unsigned char)(Clamp(b, zro, one)*255);
        texture[k+3] = (unsigned char)(Clamp(a, zro, one)*255);
      }
    }
  }

  // Set normal for lighting.
  Vector normal = Cross( p2 - p1, p4 - p1 );
  normal.normalize();
  float n[3];
  n[0] = normal.x(); n[1] = normal.y(); n[2] = normal.z();
  tr->set_normal( n );

  tr->set_transparency( get_flag(render_state,
				 RenderStateBase::USE_TRANSPARENCY) );
  tr->set_coords(tex_coords, pos_coords);
  tr->set_texture( texture, colorbytes, width, height );

  delete [] texture;
  
std::cerr << "texture="<<  texture_face.get_rep() << "\n";
  return texture_face;
}


GeomHandle 
RenderFieldV::render_text(FieldHandle field_handle,
			  ColorMapHandle colormap_handle,
			  bool use_colormap,
			  bool use_default_color,
			  bool backface_cull_p,
			  int fontsize,
			  int precision,
			  bool render_locations,
			  bool render_data,
			  bool render_nodes,
			  bool render_edges,
			  bool render_faces,
			  bool render_cells,
			  bool always_visible)
{
  GeomGroup *texts = new GeomGroup;
  GeomHandle text_switch = new GeomSwitch(texts);

  if (render_data)
  {
    texts->add(render_text_data(field_handle, colormap_handle, 
				use_colormap,
				use_default_color,
				backface_cull_p,
				fontsize, precision, always_visible));
  }
  if (render_nodes)
  {
    texts->add(render_text_nodes(field_handle, colormap_handle, 
				 use_colormap,
				 use_default_color,
				 backface_cull_p,
				 fontsize, precision, 
				 render_locations, always_visible));
  }
  if (render_edges)
  {
    texts->add(render_text_edges(field_handle, colormap_handle, 
				 use_colormap,
				 use_default_color,
				 fontsize, precision, 
				 render_locations, always_visible));
  }
  if (render_faces)
  {
    texts->add(render_text_faces(field_handle, colormap_handle, 
				 use_colormap,
				 use_default_color,
				 fontsize, precision, 
				 render_locations, always_visible));
  }
  if (render_cells)
  {
    texts->add(render_text_cells(field_handle, colormap_handle, 
				 use_colormap,
				 use_default_color,
				 fontsize, precision, 
				 render_locations, always_visible));
  }
  return text_switch;
}


GeomHandle 
RenderFieldV::render_text_data(FieldHandle field_handle,
			       ColorMapHandle colormap_handle, 
			       bool use_color_map,
			       bool use_default_color,
			       bool backface_cull_p,
			       int fontsize,
			       int precision,
			       bool always_visible)
{
  if (backface_cull_p && field_handle->basis_order() == 1)
  {
    return render_text_data_nodes(field_handle, colormap_handle, use_color_map,
				  use_default_color,
				  backface_cull_p, fontsize,
				  precision, always_visible);
  }

  VField* fld = field_handle->vfield();
  VMesh* mesh = field_handle->vmesh();

  GeomTexts *texts = new GeomTexts();
  if (always_visible) texts->set_always_visible();
  GeomHandle text_switch = new GeomSwitch(new GeomDL(texts));
  texts->set_font_index(fontsize);

  std::string buffer;

  double sval;
  Vector vval;
  Tensor tval;

  unsigned int color_scheme = 0;
  double scol;
  MaterialHandle vcol(0);

  if( use_default_color )
    color_scheme = 0;
  else if( use_color_map )
    color_scheme = 1;
  else
  {
    color_scheme = 2;
    vcol = new Material(Color(1.0, 1.0, 1.0));
  }

  if (fld->basis_order() == 0)
  {
    VMesh::Elem::iterator iter, end;
    mesh->begin(iter);
    mesh->end(end);
    Point p;
    
    while (iter != end)
    {
      mesh->get_center(p, *iter);

      if (fld->is_scalar())
      {
        fld->get_value(sval, *iter);
        buffer = to_string(sval, precision);
        value_to_color( color_scheme, sval, scol, vcol );
      }
      else if (fld->is_vector())
      {
        fld->get_value(vval, *iter);
        buffer = to_string(vval, precision);
        value_to_color( color_scheme, vval, scol, vcol );
      }
      else if (fld->is_tensor())
      {
        fld->get_value(tval, *iter);
        buffer = to_string(tval, precision);
        value_to_color( color_scheme, tval, scol, vcol );
      }

      if (color_scheme == 0)
      {
        texts->add(buffer, p);
      }
      else if (color_scheme == 1)
      {
        // Compute the ColorMap index and retreive the color.
        const double cmin = colormap_handle->getMin();
        const double cmax = colormap_handle->getMax();
        const double index = Clamp((scol - cmin)/(cmax - cmin), 0.0, 1.0);
        float r,g,b,a;
        colormap_handle->get_color(index, r, g, b, a);
        const Color c(r, g, b);
              
        texts->add(buffer, p, c);
      }
      else // if (color_scheme == 2)
      {
        texts->add(buffer, p, vcol->diffuse);
      }

      ++iter;
    }
  }
  else if (fld->basis_order() == 1)
  {
    VMesh::Node::iterator iter, end;
    mesh->begin(iter);
    mesh->end(end);
    Point p;
    
    while (iter != end)
    {
      mesh->get_center(p, *iter);

      if (fld->is_scalar())
      {
        fld->get_value(sval, *iter);
        buffer = to_string(sval, precision);
        value_to_color( color_scheme, sval, scol, vcol );
      }
      else if (fld->is_vector())
      {
        fld->get_value(vval, *iter);
        buffer = to_string(vval, precision);
        value_to_color( color_scheme, vval, scol, vcol );
      }
      else if (fld->is_tensor())
      {
        fld->get_value(tval, *iter);
        buffer = to_string(tval, precision);
        value_to_color( color_scheme, tval, scol, vcol );
      }

      if (color_scheme == 0)
      {
        texts->add(buffer, p);
      }
      else if (color_scheme == 1)
      {
        // Compute the ColorMap index and retreive the color.
        const double cmin = colormap_handle->getMin();
        const double cmax = colormap_handle->getMax();
        const double index = Clamp((scol - cmin)/(cmax - cmin), 0.0, 1.0);
        float r,g,b,a;
        colormap_handle->get_color(index, r, g, b, a);
        const Color c(r, g, b);
              
        texts->add(buffer, p, c);
      }
      else // if (color_scheme == 2)
      {
        texts->add(buffer, p, vcol->diffuse);
      }

      ++iter;
    }
  }



  
  return text_switch;
}


GeomHandle 
RenderFieldV::render_text_data_nodes(FieldHandle field_handle,
				     ColorMapHandle colormap_handle,
				     bool use_color_map,
				     bool use_default_color,
				     bool backface_cull_p,
				     int fontsize,
				     int precision,
				     bool always_visible)
{
  VField *fld = field_handle->vfield();
  VMesh* mesh = field_handle->vmesh();

  GeomTexts *texts = 0;
  GeomTextsCulled *ctexts = 0;
  GeomHandle text_switch = 0;

  const bool culling_p = backface_cull_p && mesh->has_normals();

  if (culling_p)
  {
    mesh->synchronize(Mesh::NORMALS_E);
    ctexts = new GeomTextsCulled();
    text_switch = new GeomSwitch(ctexts);
    ctexts->set_font_index(fontsize);
  }
  else
  {
    texts = new GeomTexts();
    if (always_visible) texts->set_always_visible();
    text_switch = new GeomSwitch(new GeomDL(texts));
    texts->set_font_index(fontsize);
  }

  std::string buffer;

  double sval;
  Vector vval;
  Tensor tval;

  unsigned int color_scheme = 0;
  double scol;
  MaterialHandle vcol(0);

  if( fld->basis_order() < 0 ||
      (fld->basis_order() == 0 && mesh->dimensionality() != 0) ||
      use_default_color )
    color_scheme = 0;
  else if( use_color_map )
    color_scheme = 1;
  else
  {
    color_scheme = 2;
    vcol = new Material(Color(1.0, 1.0, 1.0));
  }

  VMesh::Node::iterator iter, end;
  mesh->begin(iter);
  mesh->end(end);
  Point p;
  Vector n;
  
  while (iter != end) 
  {
    mesh->get_center(p, *iter);

    if (fld->is_scalar())
    {
      fld->get_value(sval, *iter);
      buffer = to_string(sval, precision);
      value_to_color( color_scheme, sval, scol, vcol );
    }
    else if (fld->is_vector())
    {
      fld->get_value(vval, *iter);
      buffer = to_string(vval, precision);
      value_to_color( color_scheme, vval, scol, vcol );
    }
    else if (fld->is_tensor())
    {
      fld->get_value(tval, *iter);
      buffer = to_string(tval, precision);
      value_to_color( color_scheme, tval, scol, vcol );
    }
    
    if (color_scheme == 0)
    {
      if (culling_p)
      {
        mesh->get_normal(n, *iter);
        ctexts->add(buffer, p, n);
      }
      else
      {
        texts->add(buffer, p);
      }
    }
    else if (color_scheme == 1)
    {
      // Compute the ColorMap index and retreive the color.
      const double cmin = colormap_handle->getMin();
      const double cmax = colormap_handle->getMax();
      const double index = Clamp((scol - cmin)/(cmax - cmin), 0.0, 1.0);
      float r,g,b,a;
      colormap_handle->get_color(index, r, g, b, a);
      const Color c(r, g, b);
      
      if (culling_p)
      {
        mesh->get_normal(n, *iter);
        ctexts->add(buffer, p, n, c);
      }
      else
      {
        texts->add(buffer, p, c);
      }
    }
    else // if (color_scheme == 2)
    {
      if (culling_p)
      {
        mesh->get_normal(n, *iter);
        ctexts->add(buffer, p, n, vcol->diffuse);
      }
      else
      {
        texts->add(buffer, p, vcol->diffuse);
      }
    }
            
    ++iter;
  }

  return text_switch;
}


GeomHandle 
RenderFieldV::render_text_nodes(FieldHandle field_handle,
				ColorMapHandle colormap_handle,
				bool use_color_map,
				bool use_default_color,
				bool backface_cull_p,
				int fontsize,
				int precision,
				bool render_locations,
				bool always_visible)
{
  VField *fld = field_handle->vfield();
  VMesh* mesh = field_handle->vmesh();

  GeomTexts *texts = 0;
  GeomTextsCulled *ctexts = 0;
  GeomHandle text_switch = 0;

  const bool culling_p = backface_cull_p && mesh->has_normals();

  if (culling_p)
  {
    mesh->synchronize(Mesh::NORMALS_E);
    ctexts = new GeomTextsCulled();
    text_switch = new GeomSwitch(ctexts);
    ctexts->set_font_index(fontsize);
  }
  else
  {
    texts = new GeomTexts();
    if (always_visible) texts->set_always_visible();
    text_switch = new GeomSwitch(new GeomDL(texts));
    texts->set_font_index(fontsize);
  }

  std::ostringstream buffer;
  buffer.precision(precision);

  double sval;
  Vector vval;
  Tensor tval;

  unsigned int color_scheme = 0;
  double scol;
  MaterialHandle vcol(0);

  if( (fld->basis_order() <  0) ||
      (fld->basis_order() == 0 && mesh->dimensionality() != 0) ||
      use_default_color )
    color_scheme = 0;
  else if( use_color_map )
    color_scheme = 1;
  else
  {
    color_scheme = 2;
    vcol = new Material(Color(1.0, 1.0, 1.0));
  }

  VMesh::Node::iterator iter, end;
  mesh->begin(iter);
  mesh->end(end);
  Point p;
  Vector n;

  while (iter != end)
  {
    mesh->get_center(p, *iter);

    buffer.str("");
    if (render_locations)
    {
      buffer << p;
    }
    else
    {
      (*iter).str_render(buffer);
    }

    if (color_scheme) 
    {
      if (fld->is_scalar())
      {
        fld->get_value(sval, *iter);
        value_to_color( color_scheme, sval, scol, vcol );
      }
      else if (fld->is_vector())
      {
        fld->get_value(vval, *iter);
        value_to_color( color_scheme, vval, scol, vcol );
      }
      else if (fld->is_tensor())
      {
        fld->get_value(tval, *iter);
        value_to_color( color_scheme, tval, scol, vcol );
      }
    }
    
    if (color_scheme == 0)
    {
      if (culling_p)
      {
        mesh->get_normal(n, *iter);
        ctexts->add(buffer.str(), p, n);
      }
      else
      {
        texts->add(buffer.str(), p);
      }
    }
    else if (color_scheme == 1)
    {
      // Compute the ColorMap index and retreive the color.
      const double cmin = colormap_handle->getMin();
      const double cmax = colormap_handle->getMax();
      const double index = Clamp((scol - cmin)/(cmax - cmin), 0.0, 1.0);
      float r,g,b,a;
      colormap_handle->get_color(index, r, g, b, a);
      const Color c(r, g, b);
      
      if (culling_p)
      {
        mesh->get_normal(n, *iter);
        ctexts->add(buffer.str(), p, n, c);
      }
      else
      {
        texts->add(buffer.str(), p, c);
      }
    }
    else // if (color_scheme == 2)
    {
      if (culling_p)
      {
        mesh->get_normal(n, *iter);
        ctexts->add(buffer.str(), p, n, vcol->diffuse);
      }
      else
      {
        texts->add(buffer.str(), p, vcol->diffuse);
      }
    }

    ++iter;
  }
  return text_switch;
}


GeomHandle 
RenderFieldV::render_text_edges(FieldHandle field_handle,
				ColorMapHandle colormap_handle,
				bool use_color_map,
				bool use_default_color,
				int fontsize,
				int precision,
				bool render_locations,
				bool always_visible)
{
  VField *fld = field_handle->vfield();
  VMesh* mesh = field_handle->vmesh();

  mesh->synchronize(Mesh::EDGES_E);

  GeomTexts *texts = new GeomTexts();
  if (always_visible) texts->set_always_visible();
  GeomHandle text_switch = new GeomSwitch(new GeomDL(texts));
  texts->set_font_index(fontsize);

  std::ostringstream buffer;
  buffer.precision(precision);

  double sval;
  Vector vval;
  Tensor tval;

  unsigned int color_scheme = 0;
  double scol;
  MaterialHandle vcol(0);

  if( fld->basis_order() != 0 || mesh->dimensionality() != 1 ||
      use_default_color )
    color_scheme = 0;
  else if( use_color_map )
    color_scheme = 1;
  else
  {
    color_scheme = 2;
    vcol = new Material(Color(1.0, 1.0, 1.0));
  }
 
  VMesh::Edge::iterator iter, end;
  mesh->begin(iter);
  mesh->end(end);
  Point p;

  while (iter != end)
  {
    mesh->get_center(p, *iter);

    buffer.str("");
    if (render_locations)
    {
      buffer << p;
    }
    else
    {
      buffer << (int)(*iter);
    }

    if (color_scheme) 
    {
      if (fld->is_scalar())
      {
        fld->get_value(sval, *iter);
        value_to_color( color_scheme, sval, scol, vcol );
      }
      else if (fld->is_vector())
      {
        fld->get_value(vval, *iter);
        value_to_color( color_scheme, vval, scol, vcol );
      }
      else if (fld->is_tensor())
      {
        fld->get_value(tval, *iter);
        value_to_color( color_scheme, tval, scol, vcol );
      }
    }
    
    if (color_scheme == 0)
    {
      texts->add(buffer.str(), p);
    }
    else if (color_scheme == 1)
    {
      // Compute the ColorMap index and retreive the color.
      const double cmin = colormap_handle->getMin();
      const double cmax = colormap_handle->getMax();
      const double index = Clamp((scol - cmin)/(cmax - cmin), 0.0, 1.0);
      float r,g,b,a;
      colormap_handle->get_color(index, r, g, b, a);
      const Color c(r, g, b);
      
      texts->add(buffer.str(), p, c);
    }
    else // if (color_scheme == 2)
    {
      texts->add(buffer.str(), p, vcol->diffuse);
    }

    ++iter;
  }
  return text_switch;
}

GeomHandle 
RenderFieldV::render_text_faces(FieldHandle field_handle,
				ColorMapHandle colormap_handle,
				bool use_color_map,
				bool use_default_color,
				int fontsize,
				int precision,
				bool render_locations,
				bool always_visible)
{
  VField *fld = field_handle->vfield();
  VMesh* mesh = field_handle->vmesh();

  mesh->synchronize(Mesh::FACES_E);

  GeomTexts *texts = new GeomTexts;
  if (always_visible) texts->set_always_visible();
  GeomHandle text_switch = new GeomSwitch(new GeomDL(texts));
  texts->set_font_index(fontsize);

  std::ostringstream buffer;
  buffer.precision(precision);

  double sval;
  Vector vval;
  Tensor tval;

  unsigned int color_scheme = 0;
  double scol;
  MaterialHandle vcol(0);

  if( fld->basis_order() != 0 || mesh->dimensionality() != 2 ||
      use_default_color )
    color_scheme = 0;
  else if( use_color_map )
    color_scheme = 1;
  else
  {
    color_scheme = 2;
    vcol = new Material(Color(1.0, 1.0, 1.0));
  }

  VMesh::Face::iterator iter, end;
  mesh->begin(iter);
  mesh->end(end);
  Point p;

  while (iter != end)
  {
    mesh->get_center(p, *iter);

    buffer.str("");
    if (render_locations)
    {
      buffer << p;
    }
    else
    {
      (*iter).str_render(buffer);
    }

    if (color_scheme) 
    {
      if (fld->is_scalar())
      {
        fld->get_value(sval, *iter);
        value_to_color( color_scheme, sval, scol, vcol );
      }
      else if (fld->is_vector())
      {
        fld->get_value(vval, *iter);
        value_to_color( color_scheme, vval, scol, vcol );
      }
      else if (fld->is_tensor())
      {
        fld->get_value(tval, *iter);
        value_to_color( color_scheme, tval, scol, vcol );
      }
    }
    
    if (color_scheme == 0)
    {
      texts->add(buffer.str(), p);
    }
    else if (color_scheme == 1)
    {
      // Compute the ColorMap index and retreive the color.
      const double cmin = colormap_handle->getMin();
      const double cmax = colormap_handle->getMax();
      const double index = Clamp((scol - cmin)/(cmax - cmin), 0.0, 1.0);
      float r,g,b,a;
      colormap_handle->get_color(index, r, g, b, a);
      const Color c(r, g, b);
      
      texts->add(buffer.str(), p, c);
    }
    else // if (color_scheme == 2)
    {
      texts->add(buffer.str(), p, vcol->diffuse);
    }

    ++iter;
  }
  return text_switch;
}


GeomHandle 
RenderFieldV::render_text_cells(FieldHandle field_handle,
				ColorMapHandle colormap_handle,
				bool use_color_map,
				bool use_default_color,
				int fontsize,
				int precision,
				bool render_locations,
				bool always_visible)
{
  VField *fld = field_handle->vfield();
  VMesh* mesh = field_handle->vmesh();

  mesh->synchronize(Mesh::CELLS_E);

  GeomTexts *texts = new GeomTexts;
  if (always_visible) texts->set_always_visible();
  GeomHandle text_switch = new GeomSwitch(new GeomDL(texts));
  texts->set_font_index(fontsize);

  std::ostringstream buffer;
  buffer.precision(precision);

  double sval;
  Vector vval;
  Tensor tval;

  unsigned int color_scheme = 0;
  double scol;
  MaterialHandle vcol(0);

  if( fld->basis_order() != 0 || use_default_color )
    color_scheme = 0;
  else if( use_color_map )
    color_scheme = 1;
  else
  {
    color_scheme = 2;
    vcol = new Material(Color(1.0, 1.0, 1.0));
  }

  VMesh::Cell::iterator iter, end;
  mesh->begin(iter);
  mesh->end(end);
  Point p;

  while (iter != end)
  {
    mesh->get_center(p, *iter);

    buffer.str("");
    if (render_locations)
    {
      buffer << p;
    }
    else
    {
      (*iter).str_render(buffer);
    }

    if (color_scheme) 
    {
      if (fld->is_scalar())
      {
        fld->get_value(sval, *iter);
        value_to_color( color_scheme, sval, scol, vcol );
      }
      else if (fld->is_vector())
      {
        fld->get_value(vval, *iter);
        value_to_color( color_scheme, vval, scol, vcol );
      }
      else if (fld->is_tensor())
      {
        fld->get_value(tval, *iter);
        value_to_color( color_scheme, tval, scol, vcol );
      }
    }
    
    if (color_scheme == 0)
    {
      texts->add(buffer.str(), p);
    }
    else if (color_scheme == 1)
    {
      // Compute the ColorMap index and retreive the color.
      const double cmin = colormap_handle->getMin();
      const double cmax = colormap_handle->getMax();
      const double index = Clamp((scol - cmin)/(cmax - cmin), 0.0, 1.0);
      float r,g,b,a;
      colormap_handle->get_color(index, r, g, b, a);
      const Color c(r, g, b);
      
      texts->add(buffer.str(), p, c);
    }
    else // if (color_scheme == 2)
    {
      texts->add(buffer.str(), p, vcol->diffuse);
    }
    
    ++iter;
  }
  return text_switch;
}

} // end namespace SCIRun
