#include "modelnode.h"
#include "util/mathutil.h"
#include "util/array.h"
#include "vol/volume.h"
#include "image/imagefile.h"


boolean checkIndex(int i, int j, int k)
{
    return(((i == 14)&&(j == 11)&&(k == 6)));
}

OctModel::OctModel(const VISVolume<float>& vol)
{
    setup();
    _values->assignVol(vol, -1.0f);
}

OctModel::OctModel(unsigned levels)
{
    setup();
    ModelNode node;
    node.state(UNKNOWN_STATUS);
    node.value(-1.0f);
    values(new OctNode(node, levels));
}

void OctModel::construct_lists(float scale_factor)
{
/* construct a list of active pixels */
    if (_active_list != NULL)
	delete _active_list;
    
    _active_list = _values->makeActiveVISList();
    makeLayers(scale_factor);

    OctMeshNode<ModelNode>* leaf;
    int h, d, w;
    
    h = d = w = _values->size();

    for (int i = 0; i < w; i++)
	for (int j = 0; j < h; j++)
	    for (int k = 0; k < d; k++)
		{
		    leaf = _values->leafRef(i, j, k);
		    if ((leaf->value())->state() == UNKNOWN_STATUS)
			if ((leaf->value())->value() < 0.0f)
			    (leaf->valueRef())->value(-1.0f);
			else
			    (leaf->valueRef())->value(1.0f);
		}

}

void OctModel::makeLayers(float scale_factor)
{
    int N[6][3];
    N[0][0] = 1;
    N[0][1] = 0;
    N[0][2] = 0;
    N[1][0] = -1;
    N[1][1] = 0;
    N[1][2] = 0;
    N[2][0] = 0;
    N[2][1] = 1;
    N[2][2] = 0;
    N[3][0] = 0;
    N[3][1] = -1;
    N[3][2] = 0;
    N[4][0] = 0;
    N[4][1] = 0;
    N[4][2] = 1;
    N[5][0] = 0;
    N[5][1] = 0;
    N[5][2] = -1;

/* construct a list of active pixels */

    VISVolIndex index;    
    VolIndexValue v_index;    
    int x, y, z;
    int k;
    unsigned h = _values->height(), 
	w = _values->width(), d = _values->depth();
    
    _inside_list[0].clean();
    _outside_list[0].clean();
	    

// set a layer pixels on the inside and the outside to sandwitch
// the active pixels 6-15-95

    OctLeafNode *neighbor_tmp, *oct_leaf;
    ModelNode new_value;
    ModelNode *the_node;
    float value;

    _active_list->reset();
    while (_active_list->valid())
	{
	    v_index = _active_list->itemAtCurrent();
	    x = v_index.a();
	    y = v_index.b();
	    z = v_index.c();
	    if ((x < (w - 1))&&(x > 0)&&(y < (h - 1))&&(y > 0)
		&&(z < (d - 1))&&(z > 0))
		{

		    oct_leaf = (OctLeafNode*)
			(_values->leafRefCreate(x, y, z, 0));
		    the_node = oct_leaf->valueRef();
		    value = the_node->value();
		    the_node->value(value*scale_factor);
		    
		    for (k = 0; k < 6; k++)
			{
// if the mesh can get you there, then use that
			    if (!(neighbor_tmp = (OctLeafNode*)
				  oct_leaf->neighbor(k)))
				neighbor_tmp = (OctLeafNode*)
				    _values->leafRefCreate
				    (x + N[k][0],
				     y + N[k][1],
				     z + N[k][2], 0);

			    the_node = neighbor_tmp->valueRef();
			    
			    if (the_node->state() == UNKNOWN_STATUS)
				{
				    value = the_node->value();
				    if (value > 0.0f)
					{
					    _inside_list[0].appendItem
						(VISVolIndex
						 (x + N[k][0],
						  y + N[k][1],
						  z + N[k][2]));
					    
					    the_node->value
						(value*scale_factor);
					    
					    the_node->state
						(INSIDE_STATUS);
					}
				    else		 
					{
					    _outside_list[0]
						.appendItem
						(VISVolIndex
						 (x + N[k][0],
						  y + N[k][1],
						  z + N[k][2]));
					    the_node->value
						(value*scale_factor);
					    the_node->state
						(OUTSIDE_STATUS); 
					}
				} // if status not known
			} // for each neighbor direction
		} // if in bounds
	    _active_list->stepForward();
	} // while deal with active list

//    set another layer outside that 

    _inside_list[1].clean();
    _outside_list[1].clean();
	    
    _inside_list[0].reset();
    _outside_list[0].reset();
    
    while (_inside_list[0].valid())
	{
	    index = _inside_list[0].itemAtCurrent();
	    x = index.a();
	    y = index.b();
	    z = index.c();
	    if ((x < (w - 1))&&(x > 0)&&(y < (h - 1))&&(y > 0)
		&&(z < (d - 1))&&(z > 0))
		{
		    oct_leaf = (OctLeafNode*)
			(_values->leafRefCreate(x, y, z, 0));
		    
		    for (k = 0; k < 6; k++)
			{
// if the mesh can get you there, then use that
			    if (!(neighbor_tmp = (OctLeafNode*)
				  oct_leaf->neighbor(k)))
				neighbor_tmp = (OctLeafNode*)
				    _values->leafRefCreate
				    (x + N[k][0],
				     y + N[k][1],
				     z + N[k][2], 0);

			    the_node = neighbor_tmp->valueRef();
			    value = the_node->value();

			    if (the_node->state() == UNKNOWN_STATUS)
				{
				    _inside_list[1].appendItem
					(VISVolIndex(x + N[k][0],
						      y + N[k][1],
						      z + N[k][2]));
				    
				    the_node->value
					(value*scale_factor);
				    the_node->state
					(INSIDE_STATUS_2);
				}
			}

		}
	    _inside_list[0].stepForward();
	}
    
    while (_outside_list[0].valid())
	{
	    index = _outside_list[0].itemAtCurrent();
	    x = index.a();
	    y = index.b();
	    z = index.c();


	    if ((x < (w - 1))&&(x > 0)&&(y < (h - 1))&&(y > 0)
		&&(z < (d - 1))&&(z > 0))
		{
		    oct_leaf = (OctLeafNode*)
			(_values->leafRefCreate(x, y, z, 0));
			
		    
		    for (k = 0; k < 6; k++)
			{
// if the mesh can get you there, then use that

			    if (!(neighbor_tmp = (OctLeafNode*)
				  oct_leaf->neighbor(k)))
				neighbor_tmp = (OctLeafNode*)
				    _values->leafRefCreate
				    (x + N[k][0],
				     y + N[k][1],
				     z + N[k][2], 0);
		    
			    the_node = neighbor_tmp->valueRef();
			    value = the_node->value();
			    
			    if (the_node->state() == UNKNOWN_STATUS)
				{
				    _outside_list[1].appendItem
					(VISVolIndex(x + N[k][0],
						      y + N[k][1],
						      z + N[k][2]));
				
				    the_node->value
					(value*scale_factor);

				    the_node->state
					(OUTSIDE_STATUS_2);
				}
			}
		}
	    _outside_list[0].stepForward();
	}
}


void OctModel::setup()
{
    _values = new OctNode();
    _active_list = new VolIndexValueVISList();
    _status_up_list = new VISVolIndexVISList[2];
    _status_down_list = new VISVolIndexVISList[2];
    _inside_list = new VolIndexValueVISList[2];
    _outside_list = new VolIndexValueVISList[2];

    _wave_dt = WAVE_DT;
    _curve_dt = CURVE_DT;

}


void OctNode::setBorder(ModelNode value, unsigned w)
{

    unsigned width_tmp;
    int i, j, k;

    width_tmp = MIN(w, depth());
    for (i = 0; i < width_tmp; i++)
	for (j = 0; j < height(); j++)	
	    for (k = 0; k < width(); k++)	
		{
		    setValue(value, k, j, i);
		    setValue(value, k, j, (depth() - 1) - i);
		}

    width_tmp = MIN(w, height());
    for (i = 0; i < depth(); i++)
	for (j = 0; j < width_tmp; j++)	
	    for (k = 0; k < width(); k++)	
		{
		    setValue(value, k, j, i);
		    setValue(value, k, (height() - 1) - j, i);
		}

    width_tmp = MIN(w, width());
    for (i = 0; i < depth(); i++)
	for (j = 0; j < height(); j++)	
	    for (k = 0; k < width_tmp; k++)	
		{
		    setValue(value, k, j, i);
		    setValue(value, (width() - 1) - k, j, i);
		}
}


const OctNode& OctNode::assignVol(const VISVolume<float>& vol, float padding)
{

    unsigned max_size;
    max_size = VISmax(vol.depth(), VISmax(vol.height(), vol.width()));

    int levels = ceil(log((float)max_size)/log(2.0));

    deleteChildren();
    _value = new ModelNode(padding);
    level(levels);

    setValue(padding, 0, 0, 0, levels);
    
    int i, j, k;

    for (k = 0; k < vol.depth(); k++)
	for (j = 0; j < vol.depth(); j++)
	    for (i = 0; i < vol.depth(); i++)
		setValue(vol.itemAt(i, j, k), i, j, k);

    return(*this);

}


VolIndexValueVISList* OctNode::makeActiveVISList()
{
    VolIndexValueVISList* r = new VolIndexValueVISList();
    
    float this_value, that_value;
    ModelNode *this_node, *that_node;
    unsigned int w, h, d;

    h = d = w = size();
    
    int i, j, k;

    r->clear();
    r->reset();

    for (k = 0; k < d; k++)
	for (i = 0; i < h; i++)
	    {
		for (j = 0; j < w - 1; j++)
		    {
			this_node = (leafRef(j, i, k))->valueRef();
			that_node = (leafRef(j + 1, i, k))->valueRef();
			this_value = this_node->value();
			that_value = that_node->value();
			if (((this_value < 0.0f)&&(that_value > 0.0f))
			    ||((this_value > 0.0f)&&(that_value < 0.0f))
			    ||((this_value == 0.0f)&&(that_value != 0.0f))
			    ||((this_value != 0.0f)&&(that_value == 0.0f))
			    )
			    {
				if (fabs(this_value) > fabs(that_value))
				    {
					if (that_node->state() != ACTIVE_STATUS)
					    {
						that_node
						    ->state(ACTIVE_STATUS);
						r->appendItem
						    (VolIndexValue(j + 1, 
								   i, k));
					    }
				    }
				
				else
				    {
					if (this_node->state() != ACTIVE_STATUS)
					    {
						this_node
						    ->state(ACTIVE_STATUS);
						r->appendItem
						    (VolIndexValue
						     (j, i, k));
					    }
				    }
			    }
		    }
	    }

    for (k = 0; k < d; k++)
	for (j = 0; j < w; j++)
	    {
		for (i = 0; i < h - 1; i++)
		    {
			this_node = (leafRef(j, i, k))->valueRef();
			that_node = (leafRef(j, i + 1, k))->valueRef();
			this_value = this_node->value();
			that_value = that_node->value();
			if (((this_value < 0.0f)&&(that_value > 0.0f))
			    ||((this_value > 0.0f)&&(that_value < 0.0f))
			    ||((this_value == 0.0f)&&(that_value != 0.0f))
			    ||((this_value != 0.0f)&&(that_value == 0.0f))
			    )
			    {
				if (fabs(this_value) > fabs(that_value))
				    {
					if (that_node->state() 
					    != ACTIVE_STATUS)
					    {
						that_node
						    ->state(ACTIVE_STATUS);
						r->appendItem
						    (VolIndexValue
						     (j, i + 1, k));
					    }
				    }
						
				else
				    {
					if (this_node->state() 
					    != ACTIVE_STATUS)
					    {
						this_node
						    ->state(ACTIVE_STATUS);
						r->appendItem
						    (VolIndexValue
						     (j, i, k));
					    }
				    }
			    }
		    }
	    }

    for (j = 0; j < w; j++)
	for (i = 0; i < h; i++)
	    {
		for (k = 0; k < d - 1; k++)
		    {
			this_node = (leafRef(j, i, k))->valueRef();
			that_node = (leafRef(j, i, k + 1))->valueRef();
			this_value = this_node->value();
			that_value = that_node->value();
			if (((this_value < 0.0f)&&(that_value > 0.0f))
			    ||((this_value > 0.0f)&&(that_value < 0.0f))
			    ||((this_value == 0.0f)&&(that_value != 0.0f))
			    ||((this_value != 0.0f)&&(that_value == 0.0f))
			    )
			    {
				if (fabs(this_value) > fabs(that_value))
				    {
					if (that_node
					    ->state() != ACTIVE_STATUS)
					    {
						that_node
						    ->state(ACTIVE_STATUS);
						r->appendItem
						    (VolIndexValue(j, i, 
								   k + 1));
					    }
				    }
						
				else
				    {
					if (this_node
					    ->state() != ACTIVE_STATUS)
					    {
						this_node
						    ->state(ACTIVE_STATUS);
						r->appendItem
						    (VolIndexValue
						     (j, i, k));
					    }
				    }
			    }
		    }
	    }

    return(r);
}


VISVolume<float> OctNode::values() const
{
    int w, h, d;

    w = width();
    h = height();
    d = depth();

    VISVolume<float> r(w, h, d);

    int i, j, k;
    
    for (i = 0; i < w; i++)
	for (j = 0; j < h; j++)
	    for (k = 0; k < d; k++)
		{
		    r.at(i, j, k) = value(i, j, k).value();
		}
    return(r);
}



VISVolume<byte> OctNode::state() const
{
    int w, h, d;

    w = width();
    h = height();
    d = depth();

    VISVolume<float> r(w, h, d);

    int i, j, k;
    
    for (i = 0; i < w; i++)
	for (j = 0; j < h; j++)
	    for (k = 0; k < d; k++)
		{
		    r.at(i, j, k) = value(i, j, k).state();
		}
    return(r);
}


float OctModel::calculate_change()
{
    VISVolIndex index;
    _active_list->reset();
//    int h, i, j;
    while (_active_list->valid())
	{
	    index = _active_list->itemAtCurrent();
	    (_active_list->atCurrent()).value(1.0f);
	    _active_list->stepForward();
	}
    return(WAVE_DT);
}

float OctModel::calculate_change_curve()
{
    VISVolIndex index;
    _active_list->reset();
    VolHood hood;
    VISArray<float> derivatives;
    float n_x, n_y, n_z, d_x, d_y, d_z, d_xx, d_yy, d_zz, d_xy, d_xz, d_yz;
    float grad_mag, mean_curve, grad_mag_sq, L_ww;
    
    int h, i, j;
    
    while (_active_list->valid())
	{
	    index = _active_list->itemAtCurrent();
	    j = index.a();
	    i = index.b();
	    h = index.c();
// calculate the level-set curvature
//---------------------------------
	    hood = _values->Neighborhood18(j, i, h);
	    derivatives = hood.derivatives();

	    d_x = derivatives.itemAt(D_X);
	    d_y = derivatives.itemAt(D_Y);
	    d_z = derivatives.itemAt(D_Z);	    

	    d_xx = derivatives.itemAt(D_XX);
	    d_yy = derivatives.itemAt(D_YY);
	    d_zz = derivatives.itemAt(D_ZZ);	    

	    d_xy = derivatives.itemAt(D_XY);
	    d_xz = derivatives.itemAt(D_XZ);
	    d_yz = derivatives.itemAt(D_YZ);	    
	    
	    grad_mag_sq = d_x*d_x + d_y*d_y + d_z*d_z;
		    
	    grad_mag = (float)sqrt((double)grad_mag_sq);
	    
	    if (grad_mag > 0.0)
		{
		    n_x = d_x/grad_mag;
		    n_y = d_y/grad_mag;
		    n_z = d_z/grad_mag;
		}
	    else
		n_x = n_y = n_z = 0.0f;

	    if (grad_mag_sq > 0.0f)
		{
// would it be faster to do the f_x, f_y, f_z, f_xx, f_yy, f_zz, 
// f_xy, f_yz, f_xz
		    L_ww = 
			(n_x*n_x*d_xx + n_y*n_y*d_yy
			 + n_z*n_z*d_zz)
			+ 2.0f*(n_x*n_y*d_xy + 
				n_x*n_z*d_xz + 
				n_y*n_z*d_yz);
		    mean_curve = d_xx + d_yy + d_zz - L_ww;
		}
	    else
		{
		    mean_curve = d_xx + d_yy + d_zz;
		}

	    (_active_list->atCurrent()).value(mean_curve);

//---------------------------------

	    _active_list->stepForward();
	}
    return(CURVE_DT);
}



float OctModel::calculate_change_all()
{
    VISVolIndex index;
    _active_list->reset();
    VolHood hood;
    VISArray<float> derivatives;
    float n_x, n_y, n_z, d_x, d_y, d_z, d_xx, d_yy, d_zz, d_xy, d_xz, d_yz;
    float grad_mag, mean_curve, grad_mag_sq, L_ww;
    
    int i, j, k;
    
    boolean calc_curve, calc_force, calc_grow;
    float curve_change, force_change, grow_change;
    float grow_factor;
    float max_force_change = 0, force_change_abs, max_grow_change = 0;

    calc_curve = (_curve_weight != 0.0);
    calc_force = (_force_weight != 0.0f)||(_normal_force_weight != 0.0f);
    calc_grow = ((_grow_weight != 0.0f)||
		 (_const_grow_weight != 0.0f));

    Point3 the_force, p3, d_F;
    float this_value;
	
    
    while (_active_list->valid())
	{
	    index = _active_list->itemAtCurrent();
	    i = index.a();
	    j = index.b();
	    k = index.c();
// calculate the level-set curvature
//---------------------------------
	    hood = _values->Neighborhood18(i, j, k);
	    this_value = hood.itemAt(1, 1, 1);
	    derivatives = hood.derivatives();

	    d_x = derivatives.itemAt(D_X);
	    d_y = derivatives.itemAt(D_Y);
	    d_z = derivatives.itemAt(D_Z);	    

	    d_xx = derivatives.itemAt(D_XX);
	    d_yy = derivatives.itemAt(D_YY);
	    d_zz = derivatives.itemAt(D_ZZ);	    

	    d_xy = derivatives.itemAt(D_XY);
	    d_xz = derivatives.itemAt(D_XZ);
	    d_yz = derivatives.itemAt(D_YZ);	    
	    
	    grad_mag_sq = d_x*d_x + d_y*d_y + d_z*d_z;
		    
	    grad_mag = (float)sqrt((double)grad_mag_sq);
	    
	    if (grad_mag > 0.0)
		{
		    n_x = d_x/grad_mag;
		    n_y = d_y/grad_mag;
		    n_z = d_z/grad_mag;
		}
	    else
		n_x = n_y = n_z = 0.0f;

	    if (grad_mag_sq > 0.0f)
		{
// would it be faster to do the f_x, f_y, f_z, f_xx, f_yy, f_zz, 
// f_xy, f_yz, f_xz
		    L_ww = 
			(n_x*n_x*d_xx + n_y*n_y*d_yy
			 + n_z*n_z*d_zz)
			+ 2.0f*(n_x*n_y*d_xy + 
				n_x*n_z*d_xz + 
				n_y*n_z*d_yz);
		    mean_curve = d_xx + d_yy + d_zz - L_ww;
		}
	    else
		{
		    mean_curve = d_xx + d_yy + d_zz;
		}
	    

	    curve_change = _curve_weight*mean_curve;
//--- calculate the adjusted position of the voxel

	    float len;
	    float grad_x, grad_y, grad_z, grad;

	    float grad_before_x, grad_after_x,
		grad_before_y, grad_after_y,
		grad_before_z, grad_after_z;

	    
	    if ((calc_force)||(calc_grow))
		{

		    derivatives = hood.derivsHalf();
		    
		    grad_before_x = derivatives.itemAt(D_XF);
		    grad_after_x = derivatives.itemAt(D_XB);

		    grad_before_y = derivatives.itemAt(D_YF);
		    grad_after_y = derivatives.itemAt(D_YB);

		    grad_before_z = derivatives.itemAt(D_ZF);
		    grad_after_z = derivatives.itemAt(D_ZB);

		    if (this_value > 0.0f)
			{
			    grad_x = 
				VISmin(grad_before_x, 0.0f)
				+ VISmax(grad_after_x, 0.0f);
			    grad_y = 
				VISmin(grad_before_y, 0.0f)
				+ VISmax(grad_after_y, 0.0f);
			    grad_z = 
				VISmin(grad_before_z, 0.0f)
				+ VISmax(grad_after_z, 0.0f);
			}
		    else if (this_value < 0.0f)
			{
			    grad_x = 
				VISmax(grad_before_x, 0.0f)
				+ VISmin(grad_after_x, 0.0f);
			    grad_y = 
				VISmax(grad_before_y, 0.0f)
				+ VISmin(grad_after_y, 0.0f);
			    grad_z = 
				VISmax(grad_before_z, 0.0f)
				+ VISmin(grad_after_z, 0.0f);
			}
		    else 
			grad_x = grad_y = grad_z = 0;

		    

#ifdef not_for_now		    
		    if (fabs(grad_before_x) > fabs(grad_after_x))
			grad_x = grad_before_x;
		    else
			grad_x = grad_after_x;
	    
		    if (fabs(grad_before_y) > fabs(grad_after_y))
			grad_y = grad_before_y;
		    else
			grad_y = grad_after_y;

		    if (fabs(grad_before_z) > fabs(grad_after_z))
			grad_z = grad_before_z;
		    else
			grad_z = grad_after_z;

#endif

		    d_F = Point3(grad_x, grad_y, grad_z);
		    len = sqrt(d_F*d_F);

		    p3 = Point3((float)i, (float)j, (float)k);


		    if ((len) > 0.0)
			{
			    d_F = d_F/len;
			    p3 = p3 - (d_F)
				*(this_value/(len));
//			    normal_half = d_F/sqrt(len);
			}
//		    else
//			normal_half = Point3(0.0f, 0.0f, 0.0f);



		}

	    if (calc_force)
		{
		    if (_force_weight != 0.0f)
			the_force = _force_weight
			    *force(p3.x(), p3.y(), p3.z());
		    else 
			the_force = Point3(0.0f, 0.0f, 0.0f);

// this is for a force that depend on the normal direction
		    if (_normal_force_weight != 0.0f)
			the_force += _normal_force_weight
			    *normal_force(p3.x(), p3.y(), p3.z(), 
					  d_F.x(), d_F.y(), d_F.z());

		    if ((the_force.x()) > 0.0)
			grad_x = grad_before_x;
		    else 
			grad_x = grad_after_x;
		
		    if ((the_force.y()) > 0.0)
			grad_y = grad_before_y;
		    else 
			grad_y = grad_after_y;


		    if ((the_force.z()) > 0.0)
			grad_z = grad_before_z;
		    else 
			grad_z = grad_after_z;

		    force_change = (grad_x*the_force.x() + 
				     grad_y*the_force.y() + 
				     grad_z*the_force.z());
	    
		    force_change_abs 
			= fabs(the_force.x())
			+ fabs(the_force.y())
			+ fabs(the_force.z());

		    max_force_change = VISmax(max_force_change, 
					    force_change_abs);
		} //if calc_force
	    else 
		force_change = 0;
	    
	    if (calc_grow)
		{
		    grow_factor = _const_grow_weight + 
			_grow_weight*grow(p3.x(), p3.y(), p3.z(), 
					  d_F.x(), d_F.y(), d_F.z()); 


	    if ((checkIndex(i, j, k)||checkIndex(i - 1, j, k)))
		{
		    printf("change %d %d %d\n Neighborhood\n", i, j, k);
		    hood.printData();
		    printf("point p3 %f %f %f\n", p3.x(), p3.y(), 
			   p3.z());
		    printf("grow %f\n", grow(p3.x(), p3.y(), p3.z(), 
					     d_F.x(), d_F.y(), d_F.z()));
		}


		    if (grow_factor > 0.0f)
			{
// 
// Changed for debugging --- should change back Ross 7-10-97
//
// CHANGE THIS!!!
//
			    grad_before_x = VISmax(grad_before_x, 0.0f);
			    grad_after_x = VISmax(-grad_after_x, 0.0f);	      
			    grad_x = grad_after_x*grad_after_x 
				+ grad_before_x*grad_before_x;

			    grad_before_y = VISmax(grad_before_y, 0.0f);
			    grad_after_y = VISmax(-grad_after_y, 0.0f); 
			    grad_y = grad_after_y*grad_after_y 
				+ grad_before_y*grad_before_y;

			    grad_before_z = VISmax(grad_before_z, 0.0f);
			    grad_after_z = VISmax(-grad_after_z, 0.0f);		       
			    grad_z = grad_after_z*grad_after_z 
				+ grad_before_z*grad_before_z;
			}
		    else if (grow_factor < 0.0f)
			{

			    grad_before_x = VISmin(grad_before_x, 0.0f);
			    grad_after_x = VISmin(-grad_after_x, 0.0f);	      
			    grad_x = grad_after_x*grad_after_x 
				+ grad_before_x*grad_before_x;

			    grad_before_y = VISmin(grad_before_y, 0.0f);
			    grad_after_y = VISmin(-grad_after_y, 0.0f); 
			    grad_y = grad_after_y*grad_after_y 
				+ grad_before_y*grad_before_y;

			    grad_before_z = VISmin(grad_before_z, 0.0f);
			    grad_after_z = VISmin(-grad_after_z, 0.0f);		       
			    grad_z = grad_after_z*grad_after_z 
				+ grad_before_z*grad_before_z;
			}
		    else 
			grad_x = grad_y = grad_z = 0.0f;

		    
		    grad = sqrt(grad_x + grad_y + grad_z);

//		    if (fabs(grow_factor) > max_grow_change)
//			{
//			    max_grow_change = fabs(grow_factor);
//			    printf("max change %3.2f\n", max_grow_change);
//			    hood.printData();
//			}

		    max_grow_change = VISmax(max_grow_change,
					  (float)fabs(grow_factor));

		    grow_change = grad*grow_factor;
		}


//---------------------------------

	    (_active_list->atCurrent()).value(curve_change
					      + grow_change
					      + force_change);

	    _active_list->stepForward();
	}

    max_force_change += max_grow_change; 

    float dt;

    if (_curve_weight > 0.0)
	{
	    if (max_force_change > 0.0)
		dt = VISmin((_wave_dt/max_force_change), 
			 (_curve_dt/_curve_weight));
	    else
		dt = _curve_dt/_curve_weight;
	}
    else
	{
	    if (max_force_change > 0.0)
		dt = _wave_dt/max_force_change;
	    else 
		dt = 0.0;
	}
    
   return(dt);

}


float OctModel::update(float dt)
{
//    int w, d, he;

    int N[6][3];
    N[0][0] = 1;
    N[0][1] = 0;
    N[0][2] = 0;
    N[1][0] = -1;
    N[1][1] = 0;
    N[1][2] = 0;
    N[2][0] = 0;
    N[2][1] = 1;
    N[2][2] = 0;
    N[3][0] = 0;
    N[3][1] = -1;
    N[3][2] = 0;
    N[4][0] = 0;
    N[4][1] = 0;
    N[4][2] = 1;
    N[5][0] = 0;
    N[5][1] = 0;
    N[5][2] = -1;
    
    int h, i, j, k;

    int width, depth, height;
    width = _values->width();
    height = _values->height();    
    depth = _values->depth();    

    VISVolIndex index;
    VolIndexValue v_index;

    OctLeafNode *neighbor_tmp, *oct_leaf, 
	*neighbor_tmp_next;
    ModelNode *the_node;
    int hh, ii, jj, kk;
    float new_value;

    ii = 0;
    
// 2) Pass through the active pixels again, update their values and create
// a list of pixels that are moving to and from the various lists

    _active_list->reset();

    float total_change, max_change, change;

    total_change = 0.0; 
    max_change = 0.0;

    float value_tmp;

    while (_active_list->valid())
	{
	    ii++;
	    v_index = _active_list->itemAtCurrent();
	    j = v_index.a();
	    i = v_index.b();
	    h = v_index.c();

	    change = _active_list->itemAtCurrent().value();

	    total_change += change*change;
	    max_change = VISmax(max_change, (float)fabs(change));


	    oct_leaf = (OctLeafNode*)(_values->leafRefCreate(j, i, h, 0));
	    new_value = (the_node = oct_leaf->valueRef())->value()
		+ dt*change;


	    value_tmp = (oct_leaf->valueRef())->value();
// debugging
//	    if (fabs(value_tmp) > CHANGE_FACTOR)
//		printf("this is up %d %d %d\n", j, i, h);

	    //
	    // now you have to look at the inside and outside pixels as well 
	    // for the NEW sparse-field algorithm 6-15-95
	    //

// sign change does not appear to be necessary with the new algorithm.
	    the_node->value(new_value);

	    if (new_value > CHANGE_FACTOR*1.0)
		{
		    if ((j > 3)&&(j < (width - 4))&&
			(i > 3)&&(i < (height - 4))&&
			(h > 3)&&(h < (depth - 4)))
			{
			    _status_up_list[0]
				.appendItem(VISVolIndex(j, i, h));
			    _active_list->removeCurrent();
//			    if (checkIndex(j, i, h))
//				printf("removed %d %d %d from active list\n",
//				       j, i, h);
			}
		    else
			{
			    the_node->value(CHANGE_FACTOR*1.0);
			    _active_list->stepForward();

			    if ((j < 3)&&(j > (width - 3))&&
				(i < 3)&&(i > (height - 3))&&
				(h < 3)&&(h > (depth - 3)))
				printf("got one outside\n");

			}
		}
	    else if ((new_value < CHANGE_FACTOR*(-1.0)))
		{
		    _status_down_list[0].appendItem(VISVolIndex(j, i, h));
		    _active_list->removeCurrent();
//		    if (checkIndex(j, i, h))
//			printf("removed %d %d %d from active list\n",
//			       j, i, h);
		}
	    else
		{
		    _active_list->stepForward();
		}
	}
//	printf("there are %d active pixels\n", ii);

    total_change = sqrt(total_change/(float)ii);

	
	// at this point you have the new values for the active list,
	// you have removed those pixels from the active list that
	// are changing state, and you have added those pixels to either
	// the "status_up" or "status_down" lists.
	// 3) Pass through the inside and outside pixels and update their 
	// values
    updateLayers();

// you need this tmp variable ???
//    ModelNode *modelnode;

	
//    printf("there are %d inside 1 pixels\n", ii);

// 4) Pass through the all of the lists that are changing status, and 
// change their status.

    ii = 0;
    _status_up_list[0].reset();    
    while (_status_up_list[0].valid())
	{
	    ii++;
	    index = _status_up_list[0].itemAtCurrent();
	    j = index.a();
	    i = index.b();
	    h = index.c();

	    _inside_list[0].appendItem(VISVolIndex(j, i, h));

	    oct_leaf = (OctLeafNode*)(_values->leafRefCreate(j, i, h, 0));
	    oct_leaf->valueRef()->state(INSIDE_STATUS);

	    _status_up_list[0].stepForward();


	}

//    printf("there are %d status up 0 pixels\n", ii);

    ii = 0;
    _status_down_list[0].reset();
    while (_status_down_list[0].valid())
	{
	    ii++;
	    index = _status_down_list[0].itemAtCurrent();
	    j = index.a();
	    i = index.b();
	    h = index.c();

	    _outside_list[0].appendItem(VISVolIndex(j, i, h));

	    oct_leaf = (OctLeafNode*)(_values->leafRefCreate(j, i, h, 0));
	    oct_leaf->valueRef()->state(OUTSIDE_STATUS);

	    _status_down_list[0].stepForward();
	}

//    printf("there are %d status down 0 pixels\n", ii);

    ii = 0;
    _status_up_list[1].reset();    
    while (_status_up_list[1].valid())
	{
	    ii++;
	    index = _status_up_list[1].itemAtCurrent();
	    j = index.a();
	    i = index.b();
	    h = index.c();
	    _inside_list[1].appendItem(VISVolIndex(j, i, h));
	    _status_up_list[1].removeCurrent();

	    oct_leaf = (OctLeafNode*)(_values->leafRefCreate(j, i, h, 0));
	    oct_leaf->valueRef()->state(INSIDE_STATUS_2);
	}
    
//    printf("there are %d status up 1 pixels\n", ii);

    ii = 0;
    _status_down_list[1].reset();
    while (_status_down_list[1].valid())
	{
	    ii++;
	    index = _status_down_list[1].itemAtCurrent();
	    j = index.a();
	    i = index.b();
	    h = index.c();
	    _outside_list[1].appendItem(VISVolIndex(j, i, h));
	    _status_down_list[1].removeCurrent();

	    oct_leaf = (OctLeafNode*)(_values->leafRefCreate(j, i, h, 0));
	    oct_leaf->valueRef()->state(OUTSIDE_STATUS_2);
	}


//    printf("there are %d status down 1 pixels\n", ii);
	   
    ii = 0;
    _status_up_list[0].reset();
    while (_status_up_list[0].valid())
	{
	    ii++;
	    index = _status_up_list[0].itemAtCurrent();
	    j = index.a();
	    i = index.b();
	    h = index.c();

	    oct_leaf = (OctLeafNode*)(_values->leafRefCreate(j, i, h, 0));

	    for (k = 0; k < 6; k++)
		{
		    if (!(neighbor_tmp = (OctLeafNode*)oct_leaf->neighbor(k)))
			neighbor_tmp = (OctLeafNode*)_values->leafRefCreate
			    (j + N[k][0], 
			     i + N[k][1],  
			     h + N[k][2], 0);
		    
		    the_node = neighbor_tmp->valueRef();
		    
		    if (the_node->state() == OUTSIDE_STATUS)
			{
			    the_node->value(VISmax(-CHANGE_FACTOR, 
						the_node->value()));
			    the_node->state(ACTIVE_STATUS);
			    

			    _active_list->appendItem
				(VolIndexValue(j + N[k][0], 
					       i + N[k][1],  
					       h + N[k][2], 0.0));

			    if (checkIndex(j + N[k][0], 
					   i + N[k][1],  
					   h + N[k][2]))
				printf("added %d %d %d to active list\n",
				       j + N[k][0], 
				       i + N[k][1],  
				       h + N[k][2]);

				    _status_up_list[1].appendItem
				(VISVolIndex(j + N[k][0], 
					      i + N[k][1], 
					      h + N[k][2]));  
			}
		    
		}

	    _status_up_list[0].removeCurrent();
	}

//	printf("there are %d status up 0 pixels\n", ii);

    ii = 0;
    _status_down_list[0].reset();
    while (_status_down_list[0].valid())
	{
	    ii++;
	    index = _status_down_list[0].itemAtCurrent();
	    j = index.a();
	    i = index.b();
	    h = index.c();

	    oct_leaf = (OctLeafNode*)_values->leafRefCreate(j, i, h, 0);

	    for (k = 0; k < 6; k++)
		{

		    if (!(neighbor_tmp = (OctLeafNode*)oct_leaf->neighbor(k)))
			neighbor_tmp = (OctLeafNode*)_values->leafRefCreate
			    (j + N[k][0], 
			     i + N[k][1],  
			     h + N[k][2], 0);
		    
		    the_node = neighbor_tmp->valueRef();
		    
		    
		    if (the_node->state() == INSIDE_STATUS)
			{
			    the_node->value(VISmin(CHANGE_FACTOR, 
						the_node->value()));
			    the_node->state(ACTIVE_STATUS);
			    
			    _active_list->appendItem
				(VolIndexValue(j + N[k][0], 
					       i + N[k][1],  
					       h + N[k][2], 0.0));

			    if (checkIndex(j + N[k][0], 
					   i + N[k][1],  
					   h + N[k][2]))
				printf("added %d %d %d from active list\n",
				       j + N[k][0], 
				       i + N[k][1],  
				       h + N[k][2]);

			    _status_down_list[1].appendItem
				(VISVolIndex(j + N[k][0], 
					      i + N[k][1], 
					      h + N[k][2]));  
			}
		    
		}

	    _status_down_list[0].removeCurrent();
	}


//    printf("there are %d status down 0 pixels\n", ii);
	   
    ii = 0;
    _status_up_list[1].reset();
    while (_status_up_list[1].valid())
	{
	    ii++;
	    index = _status_up_list[1].itemAtCurrent();
	    j = index.a();
	    i = index.b();
	    h = index.c();

	    oct_leaf = (OctLeafNode*)_values->leafRefCreate(j, i, h, 0);

// what happens when things move out from under these pointers
	    for (k = 0; k < 6; k++)
		{
		    if (!(neighbor_tmp = (OctLeafNode*)oct_leaf->neighbor(k)))
			neighbor_tmp = (OctLeafNode*)_values->leafRefCreate
			    (j + N[k][0], 
			     i + N[k][1],  
			     h + N[k][2], 0);
		    
		    the_node = neighbor_tmp->valueRef();

		    if (the_node->state() == OUTSIDE_STATUS_2)
			{

			    _outside_list[0].appendItem
				(VolIndexValue(j + N[k][0], 
					       i + N[k][1],  
					       h + N[k][2], 0.0));
			    
			    the_node->state(OUTSIDE_STATUS);
			    
			    hh = h + N[k][2];
			    ii = i + N[k][1];
			    jj = j + N[k][0];
		    

			    for (kk = 0; kk < 6; kk++)
				{


				    if (!(neighbor_tmp_next 
					  = neighbor_tmp->neighbor(kk)))
					neighbor_tmp_next 
					    = (OctLeafNode*)
					    _values->leafRefCreate
					    (jj + N[kk][0], 
					     ii + N[kk][1],  
					     hh + N[kk][2], 0);
		    
				    the_node = neighbor_tmp_next->valueRef();
				    
				    if (the_node->state() == UNKNOWN_STATUS)
					{
					    the_node->state(OUTSIDE_STATUS_2);
					    
					    _outside_list[1].appendItem
						(VISVolIndex(jj + N[kk][0], 
							      ii + N[kk][1], 
							      hh + N[kk][2]));
					}
				}
			}
		}
	    _status_up_list[1].removeCurrent();
	}

//    printf("there are %d status up 1 pixels\n", ii);

	   
    ii = 0;
    _status_down_list[1].reset();
    while (_status_down_list[1].valid())
	{
	    ii++;
	    index = _status_down_list[1].itemAtCurrent();
	    j = index.a();
	    i = index.b();
	    h = index.c();

	    oct_leaf = (OctLeafNode*)_values->leafRefCreate(j, i, h, 0);

// what happens when things move out from under these pointers
	    for (k = 0; k < 6; k++)
		{
		    if (!(neighbor_tmp = (OctLeafNode*)oct_leaf->neighbor(k)))
			neighbor_tmp = (OctLeafNode*)_values->leafRefCreate
			    (j + N[k][0], 
			     i + N[k][1],  
			     h + N[k][2], 0);
		    
		    the_node = neighbor_tmp->valueRef();

		    if (the_node->state() == INSIDE_STATUS_2)
			{

			    _inside_list[0].appendItem
				(VolIndexValue(j + N[k][0], 
					       i + N[k][1],  
					       h + N[k][2], 0.0));
			    
			    the_node->state(INSIDE_STATUS);
			    
			    hh = h + N[k][2];
			    ii = i + N[k][1];
			    jj = j + N[k][0];
		    

			    for (kk = 0; kk < 6; kk++)
				{


				    if (!(neighbor_tmp_next 
					  = neighbor_tmp->neighbor(kk)))
					neighbor_tmp_next 
					    = (OctLeafNode*)
					    _values->leafRefCreate
					    (jj + N[kk][0], 
					     ii + N[kk][1],  
					     hh + N[kk][2], 0);
		    
				    the_node = neighbor_tmp_next->valueRef();
				    
				    if (the_node->state() == UNKNOWN_STATUS)
					{
					    the_node->state(INSIDE_STATUS_2);
					    _inside_list[1].appendItem
						(VISVolIndex(jj + N[kk][0], 
							      ii + N[kk][1], 
							      hh + N[kk][2]));
					}
				}
			}
		}
	    _status_down_list[1].removeCurrent();
	}


//    printf("there are %d status down 1 pixels\n", ii);
//
//    ii = 0;
//    _status_up_list[0].reset();
//    while (_status_up_list[0].valid())
//	{
//	    _status_up_list[0].stepForward();	
//	    ii++;
//	}
//    printf("there are %d status up 0 pixels\n", ii);
//
//    ii = 0;
//    _status_up_list[1].reset();
//    while (_status_up_list[1].valid())
//	{
//	    _status_up_list[1].stepForward();	
	    ii++;
//	}
//    printf("there are %d status up 1 pixels\n", ii);
//
//    ii = 0;
//    _status_down_list[0].reset();
//    while (_status_down_list[0].valid())
//	{
//	    _status_down_list[0].stepForward();	
//	    ii++;
//	}
//    printf("there are %d status down 0 pixels\n", ii);
//
//    ii = 0;
//    _status_down_list[1].reset();
//    while (_status_down_list[1].valid())
//	{
//	    _status_down_list[1].stepForward();	
//	    ii++;
//	}
//    printf("there are %d status down 1 pixels\n\n", ii);



//    for (i = 0; i < NUM_DISTANCE_UPDATES; i++)
//	updateDistance();

    return(total_change);

}


void OctModel::updateLayers()
{
    int N[6][3];
    N[0][0] = 1;
    N[0][1] = 0;
    N[0][2] = 0;
    N[1][0] = -1;
    N[1][1] = 0;
    N[1][2] = 0;
    N[2][0] = 0;
    N[2][1] = 1;
    N[2][2] = 0;
    N[3][0] = 0;
    N[3][1] = -1;
    N[3][2] = 0;
    N[4][0] = 0;
    N[4][1] = 0;
    N[4][2] = 1;
    N[5][0] = 0;
    N[5][1] = 0;
    N[5][2] = -1;
	
	// at this point you have the new values for the active list,
	// you have removed those pixels from the active list that
	// are changing state, and you have added those pixels to either
	// the "status_up" or "status_down" lists.
	// 3) Pass through the inside and outside pixels and update their 
	// values
	
    float this_value;
    int total_active_neighbors;
    float neighbor_value;
    VISVolIndex index;
    int i, j, k, h;

    OctLeafNode *oct_leaf, *neighbor_tmp;
    ModelNode *the_node, *the_node_next;
    
    
    
    int ii = 0;		
    _outside_list[0].reset();
    while (_outside_list[0].valid())
	{
	    ii++;
	    index = _outside_list[0].itemAtCurrent();
	    j = index.a();
	    i = index.b();
	    h = index.c();

// If the status has been changed by some other action, then 
// just remove it from the list (in the "else" clause)
//	    if ((_values.checkBounds(j, i, h))&&
//		(_value.value(j, i, h)) == OUTSIDE_STATUS))

	    oct_leaf = (OctLeafNode*)_values->leafRefCreate(j, i, h, 0);
	    the_node = oct_leaf->valueRef();
	    
	    if (the_node->state() == OUTSIDE_STATUS)
		{
		    total_active_neighbors = 0;
		    neighbor_value = -FLT_MAX;

		    for (k = 0; k < 6; k++)
			{
			    if (!(neighbor_tmp
				  = oct_leaf->neighbor(k)))
				neighbor_tmp
				    = (OctLeafNode*)_values->leafRefCreate
				    (j + N[k][0], 
				     i + N[k][1],  
				     h + N[k][2], 0);
		    
			    the_node_next = neighbor_tmp->valueRef();

			    if (the_node_next->state() == ACTIVE_STATUS)
				{
				    neighbor_value 
					= VISmax(neighbor_value, 
					      the_node_next->value());
				    total_active_neighbors++;
				}
			}

		    if (total_active_neighbors == 0)
			{
			    _outside_list[0].removeCurrent();
			    _status_down_list[1]
				.appendItem(VISVolIndex(j, i, h));
			}
		    else  
			{
			    _outside_list[0].stepForward();
			    the_node->value(neighbor_value 
					    - DIFFERENCE_FACTOR);
			}
		}
	    else
		_outside_list[0].removeCurrent();
	}
//    printf("there are %d outside pixels\n", ii);

    
    _inside_list[0].reset();
    while (_inside_list[0].valid())
	{
	    ii++;
	    index = _inside_list[0].itemAtCurrent();
	    j = index.a();
	    i = index.b();
	    h = index.c();

// If the status has been changed by some other action, then 
// just remove it from the list (in the "else" clause)
//	    if ((_values->checkBounds(j, i, h))&&
//		(_value.value(j, i, h)) == OUTSIDE_STATUS))

	    oct_leaf = (OctLeafNode*)_values->leafRefCreate(j, i, h, 0);
	    the_node = oct_leaf->valueRef();
	    
	    if (the_node->state() == INSIDE_STATUS)
		{
		    total_active_neighbors = 0;
		    neighbor_value = FLT_MAX;

		    for (k = 0; k < 6; k++)
			{
			    if (!(neighbor_tmp
				  = oct_leaf->neighbor(k)))
				neighbor_tmp
				    = (OctLeafNode*)_values->leafRefCreate
				    (j + N[k][0], 
				     i + N[k][1],  
				     h + N[k][2], 0);
		    
			    the_node_next = 
				neighbor_tmp->valueRef();

			    if (the_node_next->state() == ACTIVE_STATUS)
				{
				    neighbor_value 
					= VISmin(neighbor_value, 
					      the_node_next->value());
				    total_active_neighbors++;
				}
			}

		    if (total_active_neighbors == 0)
			{
			    _inside_list[0].removeCurrent();
			    _status_up_list[1]
				.appendItem(VISVolIndex(j, i, h));
			}
		    else  
			{
			    _inside_list[0].stepForward();
			    the_node->value(neighbor_value 
					    + DIFFERENCE_FACTOR);
			}
		}
	    else
		_inside_list[0].removeCurrent();
	}



    _outside_list[1].reset();
    while (_outside_list[1].valid())
	{
	    ii++;
	    index = _outside_list[1].itemAtCurrent();
	    j = index.a();
	    i = index.b();
	    h = index.c();

// If the status has been changed by some other action, then 
// just remove it from the list (in the "else" clause)
// if ((_values->checkBounds(j, i, h))&&
// (_value.value(j, i, h)) == OUTSIDE_STATUS))

	    oct_leaf = (OctLeafNode*)_values->leafRefCreate(j, i, h, 0);
	    the_node = oct_leaf->valueRef();
	    
	    if (the_node->state() == OUTSIDE_STATUS_2)
		{
		    total_active_neighbors = 0;
		    neighbor_value = -FLT_MAX;

		    for (k = 0; k < 6; k++)
			{
			    if (!(neighbor_tmp
				  = oct_leaf->neighbor(k)))
				neighbor_tmp
				    = (OctLeafNode*)_values->leafRefCreate
				    (j + N[k][0], 
				     i + N[k][1],  
				     h + N[k][2], 0);
		    
			    the_node_next = neighbor_tmp->valueRef();

			    if (the_node_next->state() == OUTSIDE_STATUS)
				{
				    neighbor_value 
					= VISmax(neighbor_value, 
					      the_node_next->value());
				    total_active_neighbors++;
				}
			}

		    if (total_active_neighbors == 0)
			{
			    _outside_list[1].removeCurrent();
			    the_node->value(-1.0f);
			    the_node->state(UNKNOWN_STATUS);
			    oct_leaf->resolveChange();
			}
		    else  
			{
			    _outside_list[1].stepForward();
			    the_node->value(VISmax(neighbor_value 
					    - DIFFERENCE_FACTOR, -1.0f));
			}
		}
	    else
		_outside_list[1].removeCurrent();
	}
//    printf("there are %d outside pixels\n", ii);


    _inside_list[1].reset();
    while (_inside_list[1].valid())
	{
	    ii++;
	    index = _inside_list[1].itemAtCurrent();
	    j = index.a();
	    i = index.b();
	    h = index.c();

// If the status has been changed by some other action, then 
// just remove it from the list (in the "else" clause)
//	    if ((_values->checkBounds(j, i, h))&&
//		(_value.value(j, i, h)) == OUTSIDE_STATUS))

	    oct_leaf = (OctLeafNode*)_values->leafRefCreate(j, i, h, 0);
	    the_node = oct_leaf->valueRef();
	    
	    if (the_node->state() == INSIDE_STATUS_2)
		{
		    total_active_neighbors = 0;
		    neighbor_value = FLT_MAX;

		    for (k = 0; k < 6; k++)
			{
			    if (!(neighbor_tmp
				  = oct_leaf->neighbor(k)))
				neighbor_tmp
				    = (OctLeafNode*)_values->leafRefCreate
				    (j + N[k][0], 
				     i + N[k][1],  
				     h + N[k][2], 0);
		    
			    the_node_next = neighbor_tmp->valueRef();

			    if (the_node_next->state() == INSIDE_STATUS)
				{
				    neighbor_value = VISmin(neighbor_value, 
							 the_node_next
							 ->value());
				    total_active_neighbors++;
				}
			}

		    if (total_active_neighbors == 0)
			{
			    _inside_list[1].removeCurrent();
			    the_node->value(1.0f);
			    the_node->state(UNKNOWN_STATUS);
			    oct_leaf->resolveChange();
			}
		    else  
			{
			    _inside_list[1].stepForward();
			    the_node->value(VISmin(neighbor_value 
						+ DIFFERENCE_FACTOR, 1.0f));
			}
		}
	    else
		_inside_list[1].removeCurrent();
	}
//    printf("there are %d outside pixels\n", ii);

}



VISArray<float> VolHood::derivatives()
{
    VISArray<float> r(9);
    r = 0.0f;
    const float* buf = rep()->buffer();
    
    r.at(D_X) = 0.5f*(buf[14] - buf[12]);
    r.at(D_Y) = 0.5f*(buf[16] - buf[10]);
    r.at(D_Z) = 0.5f*(buf[22] - buf[4]);

    r.at(D_XX) = buf[14] + buf[12] - 2.0f*buf[13];
    r.at(D_YY) = buf[16] + buf[10] - 2.0f*buf[13];
    r.at(D_ZZ) = buf[22] + buf[4] - 2.0f*buf[13];

    r.at(D_XY) = 0.25f*(buf[17] + buf[9] - buf[11] - buf[15]);
    r.at(D_XZ) = 0.25f*(buf[23] + buf[3] - buf[21] - buf[5]);
    r.at(D_YZ) = 0.25f*(buf[25] + buf[1] - buf[19] - buf[7]);

    return(r);
}


VISArray<float> VolHood::derivsHalf()
{

    VISArray<float> r(6);
    const float* buf = rep()->buffer();
    
    r.at(D_XF) = buf[14] - buf[13];
    r.at(D_XB) = buf[13] - buf[12];

    r.at(D_YF) = buf[16] - buf[13];
    r.at(D_YB) = buf[13] - buf[10];

    r.at(D_ZF) = buf[22] - buf[13];
    r.at(D_ZB) = buf[13] - buf[4];

    return(r);
}

VolHood OctNode::Neighborhood26(unsigned x, unsigned y, unsigned z)
{
    VolHood r;

    r = 0.0f;
    
    OctMeshNode<ModelNode> *the_node, *the_neighbor, *the_next_neighbor;

    the_node = leafRef(x, y, z);
    if (the_node->level() == 0)
	{
// middle
	    r.at(1,1,1) = (the_node->value())->value();

	    the_neighbor = the_node->getNeighborNode(X_NEG);
// faces
	    r.at(0,1,1) = (the_neighbor->value())->value();
	    the_next_neighbor 
		= the_neighbor->getNeighborNode(Y_NEG);
// edge
	    r.at(0,0,1) = (the_next_neighbor->value())->value();

// corners
	    r.at(0,0,0) = ((the_next_neighbor->getNeighborNode(Z_NEG))->value())
		->value();
	    r.at(0,0,2) = ((the_next_neighbor->getNeighborNode(Z_POS))->value())
		->value();

	    the_next_neighbor = the_neighbor->getNeighborNode(Y_POS);
// edge
	    r.at(0,2,1) = (the_next_neighbor->value())->value();

// corners
	    r.at(0,2,0) = ((the_next_neighbor->getNeighborNode(Z_NEG))->value())
		->value();
	    r.at(0,2,2) = ((the_next_neighbor->getNeighborNode(Z_POS))->value())
		->value();

// faces
	    r.at(0,1,0) = ((the_neighbor->getNeighborNode(Z_NEG))
			   ->value())->value();

	    r.at(0,1,2) = ((the_neighbor->getNeighborNode(Z_POS))
			   ->value())->value();

	    the_neighbor = the_node->getNeighborNode(X_POS);
	    r.at(2,1,1) = (the_neighbor->value())->value();

	    the_next_neighbor = the_neighbor->getNeighborNode(Y_NEG);
// edge
	    r.at(2,0,1) = (the_next_neighbor->value())->value();

// corners
	    r.at(2,0,0) = ((the_next_neighbor->getNeighborNode(Z_NEG))->value())
		->value();
	    r.at(2,0,2) = ((the_next_neighbor->getNeighborNode(Z_POS))->value())
		->value();

	    the_next_neighbor = the_neighbor->getNeighborNode(Y_POS);
// edge
	    r.at(2,2,1) = (the_next_neighbor->value())->value();

// corners
	    r.at(2,2,0) = ((the_next_neighbor->getNeighborNode(Z_NEG))->value())
		->value();
	    r.at(2,2,2) = ((the_next_neighbor->getNeighborNode(Z_POS))->value())
		->value();

// edges
	    r.at(2,1,0) = ((the_neighbor->getNeighborNode(Z_NEG))
			   ->value())->value();

// edges
	    r.at(2,1,2) = ((the_neighbor->getNeighborNode(Z_POS))
			   ->value())->value();


	    the_neighbor = the_node->getNeighborNode(Y_NEG);
// face
	    r.at(1,0,1) = (the_neighbor->value())->value();

	    r.at(1,0,0) = ((the_neighbor->getNeighborNode(Z_NEG))
			   ->value())->value();

	    r.at(1,0,2) = ((the_neighbor->getNeighborNode(Z_POS))
			   ->value())->value();


	    the_neighbor = the_node->getNeighborNode(Y_POS);
	    r.at(1,2,1) = (the_neighbor->value())->value();

	    r.at(1,2,0) = ((the_neighbor->getNeighborNode(Z_NEG))
			   ->value())->value();

	    r.at(1,2,2) = ((the_neighbor->getNeighborNode(Z_POS))
			   ->value())->value();


	    r.at(1,1,0) = ((the_node->getNeighborNode(Z_NEG))->value())->value();

	    r.at(1,1,2) = ((the_node->getNeighborNode(Z_POS))->value())->value();
	}
    else
	{
// middle
	    r.at(1,1,1) = value(x, y, z).value();
	    r.at(0,1,1) = value(x - 1, y, z).value();
	    r.at(0,0,1) = value(x - 1, y - 1, z).value();
	    r.at(0,0,0) = value(x - 1, y - 1, z - 1).value();
	    r.at(0,0,2) = value(x - 1, y - 1, z + 1).value();
	    r.at(0,2,1) = value(x - 1, y + 1, z).value();
	    r.at(0,2,0) = value(x - 1, y + 1, z - 1).value();
	    r.at(0,2,2) = value(x - 1, y + 1, z + 1).value();
	    r.at(0,1,0) = value(x - 1, y, z - 1).value();
	    r.at(0,1,2) = value(x - 1, y, z + 1).value();
	    r.at(2,1,1) = value(x + 1, y, z).value();
	    r.at(2,0,1) = value(x + 1, y - 1, z).value();
	    r.at(2,0,0) = value(x + 1, y - 1, z - 1).value();
	    r.at(2,0,2) = value(x + 1, y - 1, z + 1).value();
	    r.at(2,2,1) = value(x + 1, y + 1, z).value();
	    r.at(2,2,0) = value(x + 1, y + 1, z - 1).value();
	    r.at(2,2,2) = value(x + 1, y + 1, z + 1).value();
	    r.at(2,1,0) = value(x + 1, y, z - 1).value();
	    r.at(2,1,2) = value(x + 1, y, z + 1).value();
	    r.at(1,0,1) = value(x, y - 1, z).value();
	    r.at(1,0,0) = value(x, y - 1, z - 1).value();
	    r.at(1,0,2) = value(x, y - 1, z + 1).value();
	    r.at(1,2,1) = value(x, y + 1, z).value();
	    r.at(1,2,0) = value(x, y + 1, z - 1).value();
	    r.at(1,2,2) = value(x, y + 1, z + 1).value();
	    r.at(1,1,0) = value(x, y, z - 1).value();
	    r.at(1,1,2) = value(x, y, z + 1).value();
	}
	
    
    return(r);
}

VolHood OctNode::Neighborhood18(unsigned x, unsigned y, unsigned z)
{
    VolHood r;

    r = 0.0f;
    
    OctMeshNode<ModelNode> *the_node, *the_neighbor;

    the_node = leafRef(x, y, z);

    if (the_node->level() == 0)
	{
	    r.at(1,1,1) = (the_node->value())->value();

	    the_neighbor = the_node->getNeighborNode(X_NEG);
	    r.at(0,1,1) = (the_neighbor->value())->value();

	    r.at(0,0,1) = ((the_neighbor->getNeighborNode(Y_NEG))
			   ->value())->value();

	    r.at(0,2,1) = ((the_neighbor->getNeighborNode(Y_POS))
			   ->value())->value();

	    r.at(0,1,0) = ((the_neighbor->getNeighborNode(Z_NEG))
			   ->value())->value();

	    r.at(0,1,2) = ((the_neighbor->getNeighborNode(Z_POS))
			   ->value())->value();


	    the_neighbor = the_node->getNeighborNode(X_POS);
	    r.at(2,1,1) = (the_neighbor->value())->value();

	    r.at(2,0,1) = ((the_neighbor->getNeighborNode(Y_NEG))
			   ->value())->value();

	    r.at(2,2,1) = ((the_neighbor->getNeighborNode(Y_POS))
			   ->value())->value();

	    r.at(2,1,0) = ((the_neighbor->getNeighborNode(Z_NEG))
			   ->value())->value();

	    r.at(2,1,2) = ((the_neighbor->getNeighborNode(Z_POS))
			   ->value())->value();


	    the_neighbor = the_node->getNeighborNode(Y_NEG);
	    r.at(1,0,1) = (the_neighbor->value())->value();

	    r.at(1,0,0) = ((the_neighbor->getNeighborNode(Z_NEG))
			   ->value())->value();

	    r.at(1,0,2) = ((the_neighbor->getNeighborNode(Z_POS))
			   ->value())->value();


	    the_neighbor = the_node->getNeighborNode(Y_POS);
	    r.at(1,2,1) = (the_neighbor->value())->value();

	    r.at(1,2,0) = ((the_neighbor->getNeighborNode(Z_NEG))
			   ->value())->value();

	    r.at(1,2,2) = ((the_neighbor->getNeighborNode(Z_POS))
			   ->value())->value();


	    r.at(1,1,0) = ((the_node->getNeighborNode(Z_NEG))->value())->value();

	    r.at(1,1,2) = ((the_node->getNeighborNode(Z_POS))->value())->value();
	}
    else
	{
// middle
	    r.at(1,1,1) = value(x, y, z).value();
	    r.at(0,1,1) = value(x - 1, y, z).value();
	    r.at(0,0,1) = value(x - 1, y - 1, z).value();
	    r.at(0,2,1) = value(x - 1, y + 1, z).value();
	    r.at(0,1,0) = value(x - 1, y, z - 1).value();
	    r.at(0,1,2) = value(x - 1, y, z + 1).value();
	    r.at(2,1,1) = value(x + 1, y, z).value();
	    r.at(2,0,1) = value(x + 1, y - 1, z).value();
	    r.at(2,2,1) = value(x + 1, y + 1, z).value();
	    r.at(2,1,0) = value(x + 1, y, z - 1).value();
	    r.at(2,1,2) = value(x + 1, y, z + 1).value();
	    r.at(1,0,1) = value(x, y - 1, z).value();
	    r.at(1,0,0) = value(x, y - 1, z - 1).value();
	    r.at(1,0,2) = value(x, y - 1, z + 1).value();
	    r.at(1,2,1) = value(x, y + 1, z).value();
	    r.at(1,2,0) = value(x, y + 1, z - 1).value();
	    r.at(1,2,2) = value(x, y + 1, z + 1).value();
	    r.at(1,1,0) = value(x, y, z - 1).value();
	    r.at(1,1,2) = value(x, y, z + 1).value();
	}

    
    return(r);
}


VolHood OctNode::Neighborhood6(unsigned x, unsigned y, unsigned z)
{
    VolHood r;

    r = 0.0f;
    
    OctMeshNode<ModelNode> *the_node, *the_neighbor;

    the_node = leafRef(x, y, z);

    if (the_node->level() == 0)
	{
//	    printf("got leaf neigh\n");
	    r.at(1,1,1) = (the_node->value())->value();

	    the_neighbor = the_node->getNeighborNode(X_NEG);
	    r.at(0,1,1) = (the_neighbor->value())->value();

	    the_neighbor = the_node->getNeighborNode(X_POS);
	    r.at(2,1,1) = (the_neighbor->value())->value();

	    the_neighbor = the_node->getNeighborNode(Y_NEG);
	    r.at(1,0,1) = (the_neighbor->value())->value();

	    the_neighbor = the_node->getNeighborNode(Y_POS);
	    r.at(1,2,1) = (the_neighbor->value())->value();

	    the_neighbor = the_node->getNeighborNode(Z_NEG);
	    r.at(1,1,0) = (the_neighbor->value())->value();

	    the_neighbor = the_node->getNeighborNode(Z_POS);
	    r.at(1,1,2) = (the_neighbor->value())->value();
	}
    else
	{
// middle
//	    printf("got non-leaf neigh\n");
	    r.at(1,1,1) = value(x, y, z).value();
	    r.at(0,1,1) = value(x - 1, y, z).value();
	    r.at(2,1,1) = value(x + 1, y, z).value();
	    r.at(1,0,1) = value(x, y - 1, z).value();
	    r.at(1,2,1) = value(x, y + 1, z).value();
	    r.at(1,1,0) = value(x, y, z - 1).value();
	    r.at(1,1,2) = value(x, y, z + 1).value();
	}

//    r.printData();

    return(r);
}



void OctMultiModel::higherResFunc()
{
    // you must
    //
    // (1) promote the tree
    //
    // (2) put all of the children of the active and 
    // first-layer nodes onto a temporary list
    //
    // (3) set the new values of these cells based on the field function
    //
    // (4) construct a set of zero crossing (active cells)
    // from this list
    // 
    // (5) construct the new neighborhoods
    // 

    VISVolIndexVISList *new_list = new VISVolIndexVISList();
    new_list->clear();

    // SOMEWHERE HERE !!!!
    // you must go through the second layer lists and set 
    // the values and states
    //

    VISVolIndex index;
    int x, y, z;

    ModelNode new_value;
    new_value.state(UNKNOWN_STATUS);
    new_value.value(1.0f);

    _inside_list[1].reset();
    while (_inside_list[1].valid())
	{
	    index = _inside_list[1].itemAtCurrent();
	    x = index.a();
	    y = index.b();
	    z = index.c();

	    _values->setValue(new_value, x, y, z);
	    _inside_list[1].stepForward();
	}

    new_value.value(-1.0f);
    _outside_list[1].reset();
    while (_outside_list[1].valid())
	{
	    index = _outside_list[1].itemAtCurrent();
	    x = index.a();
	    y = index.b();
	    z = index.c();

	    _values->setValue(new_value, x, y, z);
	    _outside_list[1].stepForward();
	}

    // (1) promote the tree

    int old_size = _values->size();

    _values->promote(1);

    int size = _values->size();
    
    _x_offset = 2.0f*_x_offset + 0.5;
    _y_offset = 2.0f*_y_offset + 0.5;
    _z_offset = 2.0f*_z_offset + 0.5;

    _scale = _scale/2.0f;    
//    _scale = 1.0f/(2.0f/_scale + 1.0f);

    // (2) put all of the children of the active and 
    // first-layer nodes onto a temporary list
    
    int i, j, k;
    int low_x, high_x, low_y, high_y, low_z, high_z; 

    _active_list->reset();
    while (_active_list->valid())
	{
	    index = _active_list->itemAtCurrent();
	    x = index.a();
	    y = index.b();
	    z = index.c();

	    if (x == BORDER_WIDTH)
		low_x = BORDER_WIDTH;
	    else
		low_x = 0;

	    if (x == ((old_size - BORDER_WIDTH) - 1))
		high_x = BORDER_WIDTH;
	    else
		high_x = 0;

	    if (y == BORDER_WIDTH)
		low_y = BORDER_WIDTH;
	    else
		low_y = 0;

	    if (y == ((old_size - BORDER_WIDTH) - 1))
		high_y = BORDER_WIDTH;
	    else
		high_y = 0;

	    if (z == BORDER_WIDTH)
		low_z = BORDER_WIDTH;
	    else
		low_z = 0;

	    if (z == ((old_size - BORDER_WIDTH) - 1))
		high_z = BORDER_WIDTH;
	    else
		high_z = 0;

	    for (i = -low_x; i < (2 + high_x); i++)
		for (j = -low_y; j < (2 + high_y); j++)
		    for (k = -low_z; k < (2 + high_z); k++)
			new_list->appendItem(VISVolIndex(2*x + i, 
							  2*y + j, 
							  2*z + k));
	    _active_list->stepForward();
	}

    _inside_list[0].reset();
    while (_inside_list[0].valid())
	{
	    index = _inside_list[0].itemAtCurrent();
	    x = index.a();
	    y = index.b();
	    z = index.c();

	    if (x == BORDER_WIDTH)
		low_x = BORDER_WIDTH;
	    else
		low_x = 0;

	    if (x == ((old_size - BORDER_WIDTH) - 1))
		high_x = BORDER_WIDTH;
	    else
		high_x = 0;

	    if (y == BORDER_WIDTH)
		low_y = BORDER_WIDTH;
	    else
		low_y = 0;

	    if (y == ((old_size - BORDER_WIDTH) - 1))
		high_y = BORDER_WIDTH;
	    else
		high_y = 0;

	    if (z == BORDER_WIDTH)
		low_z = BORDER_WIDTH;
	    else
		low_z = 0;

	    if (z == ((old_size - BORDER_WIDTH) - 1))
		high_z = BORDER_WIDTH;
	    else
		high_z = 0;

	    for (i = -low_x; i < (2 + high_x); i++)
		for (j = -low_y; j < (2 + high_y); j++)
		    for (k = -low_z; k < (2 + high_z); k++)
			new_list->appendItem(VISVolIndex(2*x + i, 
							  2*y + j, 
							  2*z + k));
	    _inside_list[0].stepForward();
	}

    _outside_list[0].reset();
    while (_outside_list[0].valid())
	{
	    index = _outside_list[0].itemAtCurrent();
	    x = index.a();
	    y = index.b();
	    z = index.c();

	    if (x == BORDER_WIDTH)
		low_x = BORDER_WIDTH;
	    else
		low_x = 0;

	    if (x == ((old_size - BORDER_WIDTH) - 1))
		high_x = BORDER_WIDTH;
	    else
		high_x = 0;

	    if (y == BORDER_WIDTH)
		low_y = BORDER_WIDTH;
	    else
		low_y = 0;

	    if (y == ((old_size - BORDER_WIDTH) - 1))
		high_y = BORDER_WIDTH;
	    else
		high_y = 0;

	    if (z == BORDER_WIDTH)
		low_z = BORDER_WIDTH;
	    else
		low_z = 0;

	    if (z == ((old_size - BORDER_WIDTH) - 1))
		high_z = BORDER_WIDTH;
	    else
		high_z = 0;

	    for (i = -low_x; i < (2 + high_x); i++)
		for (j = -low_y; j < (2 + high_y); j++)
		    for (k = -low_z; k < (2 + high_z); k++)
			new_list->appendItem(VISVolIndex(2*x + i, 
							  2*y + j, 
							  2*z + k));
	    _outside_list[0].stepForward();
	}


    //
    // (3) set the new values of these cells based on the field function
    //

    OctLeafNode *new_node;
    int vox_count = 0;

    new_list->reset();
    while (new_list->valid())
	{
	    index = new_list[0].itemAtCurrent();
	    x = index.a();
	    y = index.b();
	    z = index.c();

	    new_node = (OctLeafNode*)
		_values->leafRefCreate(x, y, z, 0);

	    new_value.value(fieldFunction(position((unsigned)x, 
						   (unsigned)y, 
						   (unsigned)z))/_scale);
	    
	    new_node->value(new_value);
	    
	    new_list->stepForward();
	    vox_count++;
	}

    printf("total # and percent of voxels visited is %d %f\n",
	   vox_count, (float)vox_count/power((float)_values->size(), 3));

//    VISVolume<byte> state_vol(_values->size(), 
//			       _values->size(), 
//			       _values->size());
//    state_vol = 0;

//    new_list->reset();
//    while (new_list->valid())
//	{
//	    index = new_list->itemAtCurrent();
//	    x = index.a();
//	    y = index.b();
//	    z = index.c();
//	    state_vol.at(x, y, z) = 100;
//	    new_list->stepForward();
//	}

//    VISImageFile im_file;
//    im_file.write(((state_vol.image()).becomeFlat()), "state_list2.tif");


    //
    // (4) construct a set of zero crossing (active cells)
    // from this list
    // 

    if (_active_list)
	delete _active_list;
    _active_list = makeActiveVISList(new_list);
    makeLayers(DIFFERENCE_FACTOR);


    new_list->reset();
    new_list->reset();
    while (new_list->valid())
	{
	    index = new_list->itemAtCurrent();
	    x = index.a();
	    y = index.b();
	    z = index.c();

	    new_value = _values->value(x, y, z);

	    if (new_value.state() == UNKNOWN_STATUS)
		{
		    if (new_value.value() > 0.0)
			new_value.value(1.0f);
		    else
			new_value.value(-1.0f);
		    
		    _values->setValue(new_value, x, y, z);
		}
	    
	    new_list->stepForward();
	}


    delete new_list;

    //
    // (5) construct the new neighborhoods
    // 
    
}

void OctMultiModel::higherResInterp()
{
    // you must
    //
    // (1) promote the tree
    //
    // (2) put all of the children of the active and 
    // first-layer nodes onto a temporary list
    //
    // (3) set the new values of these cells based on the field function
    //
    // (4) construct a set of zero crossing (active cells)
    // from this list
    // 
    // (5) construct the new neighborhoods
    // 

    VolIndexValueVISList *new_list = new VolIndexValueVISList();

    new_list->clear();

    // these are the weights for the averaging (with only the face neighbors)
    float w_1 = 1.0f/3.0f;
    float w_2 = 2.0f/9.0f;
    
    float n_x, n_y, n_z;

    // SOMEWHERE HERE !!!!
    // you must go through the second layer lists and set 
    // the values and states
    //

    VISVolIndex index;
    VolIndexValue value_index;
    int x, y, z;
    int i, j, k;

    VolHood hood;

    OctNode *this_leaf;

    _active_list->reset();
    while (_active_list->valid())
	{
	    index = _active_list->itemAtCurrent();
	    x = index.a();
	    y = index.b();
	    z = index.c();
	    
	    hood = _values->Neighborhood26(x, y, z);
	    this_leaf = (OctNode*)_values->leafRef(x, y, z);

	    for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
		    for (k = 0; k < 2; k++)
			new_list->appendItem
			    (VolIndexValue
			     (2*x + i, 
			      2*y + j, 
			      2*z + k, 
			      hood.interp(0.75 + 0.5*i, 
					  0.75 + 0.5*j, 
					  0.75 + 0.5*k)
			      )
			     );
	    _active_list->stepForward();
	}


    _inside_list[0].reset();
    while (_inside_list[0].valid())
	{
	    index = _inside_list[0].itemAtCurrent();
	    x = index.a();
	    y = index.b();
	    z = index.c();
	    
	    hood = _values->Neighborhood26(x, y, z);
	    this_leaf = (OctNode*)_values->leafRef(x, y, z);

	    for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
		    for (k = 0; k < 2; k++)
			new_list->appendItem
			    (VolIndexValue
			     (2*x + i, 
			      2*y + j, 
			      2*z + k, 
			      hood.interp(0.75 + 0.5*i, 
					  0.75 + 0.5*j, 
					  0.75 + 0.5*k)
			      )
			     );
	    _inside_list[0].stepForward();
	}

    _outside_list[0].reset();
    while (_outside_list[0].valid())
	{
	    index = _outside_list[0].itemAtCurrent();
	    x = index.a();
	    y = index.b();
	    z = index.c();
	    
	    hood = _values->Neighborhood26(x, y, z);
	    this_leaf = (OctNode*)_values->leafRef(x, y, z);

	    for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
		    for (k = 0; k < 2; k++)
			new_list->appendItem
			    (VolIndexValue
			     (2*x + i, 
			      2*y + j, 
			      2*z + k, 
			      hood.interp(0.75 + 0.5*i, 
					  0.75 + 0.5*j, 
					  0.75 + 0.5*k)
			      )
			     );
	    _outside_list[0].stepForward();
	}

// get rid of the outermost and innermost layers
    ModelNode new_value;
    new_value.state(UNKNOWN_STATUS);
    new_value.value(1.0f);

    _inside_list[1].reset();
    while (_inside_list[1].valid())
	{
	    index = _inside_list[1].itemAtCurrent();
	    x = index.a();
	    y = index.b();
	    z = index.c();

	    _values->setValue(new_value, x, y, z);
	    _inside_list[1].stepForward();
	}

    new_value.value(-1.0f);
    _outside_list[1].reset();
    while (_outside_list[1].valid())
	{
	    index = _outside_list[1].itemAtCurrent();
	    x = index.a();
	    y = index.b();
	    z = index.c();

	    _values->setValue(new_value, x, y, z);
	    _outside_list[1].stepForward();
	}

    // (1) promote the tree

    _values->promote(1);

    int size = _values->size();

    _x_offset = 2.0f*_x_offset + 0.5;
    _y_offset = 2.0f*_y_offset + 0.5;
    _z_offset = 2.0f*_z_offset + 0.5;

    
//    _scale = 1.0f/(2.0f/_scale + 1.0f);
    _scale = _scale/2.0f;
    
    // (2) put all of the children of the active and 
    // first-layer nodes onto a temporary list
    

    //
    // (3) set the new values of these cells based on the interpolation
    //

    VISVolIndexVISList *the_list = new VISVolIndexVISList();
    OctLeafNode *new_node;

    new_list->reset();
    the_list->reset();
    while (new_list->valid())
	{
	    value_index = new_list->itemAtCurrent();
	    x = value_index.a();
	    y = value_index.b();
	    z = value_index.c();

	    new_node = (OctLeafNode*)
		_values->leafRefCreate(x, y, z, 0);

	    new_value.value(value_index.value());
	    new_node->value(new_value);

	    the_list->appendItem(VISVolIndex(x, y, z));
	    new_list->removeCurrent();
	}


    // (4) construct a set of zero crossing (active cells)
    // from this list
    // 

    if (_active_list)
	delete _active_list;

    _active_list = makeActiveVISList(the_list);

    //
    // (5) construct the new neighborhoods
    // 

    makeLayers(1.0f);

    //
    // (6) get rid of stuff that is no longer needed
    //
 
    the_list->reset();
    while (the_list->valid())
	{
	    index = the_list->itemAtCurrent();
	    x = index.a();
	    y = index.b();
	    z = index.c();

	    new_value = _values->value(x, y, z);

	    if (new_value.state() == UNKNOWN_STATUS)
		{
		    if (new_value.value() > 0.0)
			new_value.value(1.0f);
		    else
			new_value.value(-1.0f);
		    
		    _values->setValue(new_value, x, y, z);
		}
	    
	    the_list->stepForward();
	}


    delete new_list;
    delete the_list;

    updateLayers();
    
}


VolIndexValueVISList* OctMultiModel::makeActiveVISList(VISVolIndexVISList* list)
{
    VolIndexValueVISList* r = new VolIndexValueVISList();
    
    float this_value, that_value;
    ModelNode *this_node, *that_node;
    unsigned int w, h, d;
    OctLeafNode *this_leaf, *that_leaf;
    VISVolIndex index;

    h = d = w = _values->size();
    
    int i, j, k;
    int x, y, z;

    r->clear();
    r->reset();

    list->reset();
    while (list->valid())
	{
	    index = list->itemAtCurrent();
	    x = index.a();
	    y = index.b();
	    z = index.c();

	    this_leaf = (OctLeafNode*)_values->leafRefCreate(x, y, z, 0);
	    this_node = (_values->leafRef(x, y, z))->valueRef();
	    this_value = this_node->value();
	    // look at the neighbors
	    for (k = 0; k < 6; k++)
		{
		    that_leaf = this_leaf->getNeighborCreate((Direct)k);
		    that_node = that_leaf->valueRef();
		    that_value = that_node->value();

		    if (((this_value < 0.0f)&&(that_value > 0.0f))
			||((this_value > 0.0f)&&(that_value < 0.0f))
			||((this_value == 0.0f)&&(that_value != 0.0f))
			||((this_value != 0.0f)&&(that_value == 0.0f))
			)
			{
			    if (fabs(this_value) > fabs(that_value))
				{
				    if (that_node->state() != ACTIVE_STATUS)
					{
					    that_node
						->state(ACTIVE_STATUS);
					    r->appendItem
						(VolIndexValue
						 (x + DirectAdjust[0][k], 
						  y + DirectAdjust[1][k], 
						  z + DirectAdjust[2][k]));
					}
				}
			    else
				{
				    if (this_node->state() != ACTIVE_STATUS)
					{
					    this_node
						->state(ACTIVE_STATUS);
					    r->appendItem
						(VolIndexValue
						 (x, y, z));
					}
				}
			}
		}
	    list->stepForward();
	}
	    
    return(r);
}


void OctMultiModel::setLevels(unsigned levels)
{
    if (_values) 
	delete _values;

    _values = new OctNode(ModelNode(-1.0f), levels);
}


void OctMultiModel::initialize(unsigned levels)
{
    int x, y, z;
    ModelNode node_value;

    if (_values) 
	delete _values;

    _values = new OctNode(ModelNode(-1.0f), levels);
    int size = _values->size();

//    _z_offset = _y_offset = _x_offset = size/2.0f - 0.5;
//    _scale = 1.0f/((float)size - 1.0f);

    node_value.state(UNKNOWN_STATUS);

    for (x = 3; x < (size - 3); x++)
	for (y = 3; y < (size - 3); y++)
	    for (z = 3; z < (size - 3); z++)
		{
		    node_value.value(
			fieldFunction(
			    position((unsigned)x, 
				     (unsigned)y, 
				     (unsigned)z)
			    )
			/_scale);
		    _values->setValue(node_value, x, y, z);
		}
    construct_lists();
}


void OctMultiModel::initialize()
{
    int x, y, z;
    ModelNode node_value;

    int size = _values->size();
//
//    _z_offset = _y_offset = _x_offset = size/2.0f - 0.5;
//    _scale = 1.0f/((float)size - 1.0f);

    node_value.state(UNKNOWN_STATUS);

    Point3 p;
    for (x = BORDER_WIDTH; x < (size - BORDER_WIDTH); x++)
	for (y = BORDER_WIDTH; y < (size - BORDER_WIDTH); y++)
	    for (z = BORDER_WIDTH; z < (size - BORDER_WIDTH); z++)
		{
//		    printf("test 1");
		    p = position((unsigned)x, 
				 (unsigned)y, 
				 (unsigned)z);
//		    printf("test 2");
//		    p.printData();

		    node_value.value(fieldFunction(p)/_scale);
		    _values->setValue(node_value, x, y, z);
		}

    node_value.value(-1.0f);
    _values->setBorder(node_value, BORDER_WIDTH);
    construct_lists();
}


Point3 OctMultiModel::position(unsigned i, unsigned j, unsigned k)
{
    Point3 r;
    
    r.x(((float)i - _x_offset)*_scale);
    r.y(((float)j - _y_offset)*_scale);
    r.z(((float)k - _z_offset)*_scale);

    return(r);
}

Point3 OctMultiModel::position(float x, float y, float z)
{
    Point3 r;

    r.x((x - _x_offset)*_scale);
    r.y((y - _y_offset)*_scale);
    r.z((z - _z_offset)*_scale);

    return(r);
}


