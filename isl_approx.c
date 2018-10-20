//#include <isl/ctx.h>
//#include <isl/map.h>
#include <isl/set.h>
#include <isl/map.h>
#include<isl/isl_approx.h>
#include <isl/constraint.h>
#include <math.h>
#include <string.h>

__isl_give isl_basic_map * generate_rotation_constraint(__isl_take isl_basic_map * rotation_map, int inp_coeff_1, int inp_coeff_2, int out_coeff_1, int out_coeff_2)
{
	isl_constraint *c = isl_constraint_alloc_equality(isl_local_space_from_space(isl_basic_map_get_space(rotation_map)));
	c = isl_constraint_set_coefficient_si(c,isl_dim_in, 0, inp_coeff_1);
	c = isl_constraint_set_coefficient_si(c,isl_dim_in, 1, inp_coeff_2);
	c = isl_constraint_set_coefficient_si(c,isl_dim_out, 0, out_coeff_1);
	c = isl_constraint_set_coefficient_si(c,isl_dim_out, 1, out_coeff_2);	
	c = isl_constraint_set_constant_si(c, 0);
	rotation_map = isl_basic_map_add_constraint(rotation_map, c);
	return rotation_map;
}

__isl_give isl_basic_map * get_rotation_map(__isl_take isl_basic_set* shadow)
{
	isl_basic_map * rotation_map = isl_basic_map_from_domain(isl_basic_set_copy(shadow));
	rotation_map = isl_basic_map_drop_constraints_involving_dims(rotation_map,isl_dim_in,0,2);
	rotation_map = isl_basic_map_add_dims(rotation_map, isl_dim_out,2);
	rotation_map = isl_basic_map_set_dim_name(rotation_map, isl_dim_out,0,"o0");
	rotation_map = isl_basic_map_set_dim_name(rotation_map, isl_dim_out,1,"o1");
	rotation_map = generate_rotation_constraint(rotation_map, 1, 1, -1, 0);
	rotation_map = generate_rotation_constraint(rotation_map, 1, -1, 0, -1);
	return rotation_map;
}
__isl_give isl_constraint * generate_utvpi_constraint(isl_local_space *ls, int nDim, int i, int j, isl_val* coeff_i, isl_val* coeff_j, isl_val* constant)
{
	isl_constraint * halfspace = isl_constraint_alloc_inequality(ls);	
	int k;	
	for(k=0;k<nDim;k++)
	{
		if(k==i)
		{
			halfspace = isl_constraint_set_coefficient_val(halfspace, isl_dim_set, k, coeff_i);
		}
		else if(k==j)
		{
			halfspace = isl_constraint_set_coefficient_val(halfspace, isl_dim_set, k, coeff_j);
		}
		else
		{
			halfspace = isl_constraint_set_coefficient_si(halfspace, isl_dim_set, k,0);
		}
	}
	halfspace = isl_constraint_set_constant_val(halfspace, constant);
	return halfspace;
}
__isl_give isl_basic_set * foreach_constraint(__isl_take isl_basic_set * approx_set, __isl_take isl_constraint *c, int dim1, int dim2)
{
	isl_local_space * ls = isl_local_space_from_space(isl_basic_set_get_space(approx_set));
	int nDim = isl_basic_set_n_dim(approx_set);
	isl_val * coeff_i = isl_constraint_get_coefficient_val(c,3,0);
	isl_val * constant = isl_constraint_get_constant_val(c);
	isl_constraint * new_constraint=NULL;
	if(dim2==-1)
	{
		new_constraint = generate_utvpi_constraint(ls, nDim, dim1, -1, coeff_i, NULL, constant);
	}
	else
	{
		isl_val * coeff_j = isl_constraint_get_coefficient_val(c,3,1);
		new_constraint = generate_utvpi_constraint(ls, nDim, dim1, dim2, coeff_i, coeff_j, constant);
	}
	isl_basic_set * new_approx_set = isl_basic_set_add_constraint(approx_set, new_constraint);			
	return new_approx_set;
	
	
	//printf("\n\t New set %s",isl_basic_set_to_str(new_approx_set));
	//fflush(stdout);
	//return new_approx_set;
}

__isl_give isl_basic_set * generate_octagonal_constraints(isl_basic_set* shadow, int dim1, int dim2,__isl_keep isl_basic_set* approx_set, isl_basic_map* rotated_map_inv)
{
	isl_basic_set * rotated_interval_1 = isl_basic_set_project_out(isl_basic_set_copy(shadow),3,1,1);
	isl_basic_set * rotated_interval_2 = isl_basic_set_project_out(shadow,3,0,1);
	rotated_interval_1 = isl_basic_set_apply(rotated_interval_1, isl_basic_map_remove_dims(isl_basic_map_copy(rotated_map_inv),2,1,1));
	rotated_interval_2 = isl_basic_set_apply(rotated_interval_2, isl_basic_map_remove_dims(rotated_map_inv,2,0,1));
	printf("\n\t rotated_Interval1_+ %s %d",isl_basic_set_to_str(rotated_interval_1), dim1);
	printf("\n\t rotated_Interval2_- %s %d",isl_basic_set_to_str(rotated_interval_2), dim2);
	isl_constraint_list * constraints_1 = isl_basic_set_get_constraint_list(rotated_interval_1);
	isl_constraint_list * constraints_2 = isl_basic_set_get_constraint_list(rotated_interval_2);
	int i;
	for(i=0;i<isl_constraint_list_n_constraint(constraints_1);i++)
	{
		isl_constraint * c = isl_constraint_list_get_constraint(constraints_1,i);
		approx_set =  foreach_constraint(approx_set,c,dim1, dim2);
		isl_constraint_free(c);
	}
	for(i=0;i<isl_constraint_list_n_constraint(constraints_2);i++)
	{
		isl_constraint * c = isl_constraint_list_get_constraint(constraints_2,i);
		approx_set =  foreach_constraint(approx_set,c,dim1, dim2);
		isl_constraint_free(c);
	}
	isl_basic_set_free(rotated_interval_1);
	isl_basic_set_free(rotated_interval_2);
	isl_constraint_list_free(constraints_1);
	isl_constraint_list_free(constraints_2);
	return approx_set;
}

__isl_give isl_basic_set * generate_interval_constraints(isl_basic_set* shadow, int dim1, int dim2,__isl_keep isl_basic_set* approx_set)
{
	isl_basic_set * interval_1 = isl_basic_set_project_out(isl_basic_set_copy(shadow),3,1,1);
	isl_basic_set * interval_2 = isl_basic_set_project_out(shadow,3,0,1);
	printf("\n\t Interval1 %s %d",isl_basic_set_to_str(interval_1), dim1);
	printf("\n\t Interval2 %s %d",isl_basic_set_to_str(interval_2), dim2);
	isl_constraint_list * constraints_1 = isl_basic_set_get_constraint_list(interval_1);
	isl_constraint_list * constraints_2 = isl_basic_set_get_constraint_list(interval_2);
	int i;
	for(i=0;i<isl_constraint_list_n_constraint(constraints_1);i++)
	{
		isl_constraint * c = isl_constraint_list_get_constraint(constraints_1,i);
		approx_set =  foreach_constraint(approx_set,c,dim1, -1);
		isl_constraint_free(c);
	}
	for(i=0;i<isl_constraint_list_n_constraint(constraints_2);i++)
	{
		isl_constraint * c = isl_constraint_list_get_constraint(constraints_2,i);
		approx_set =  foreach_constraint(approx_set,c,dim2, -1);
		isl_constraint_free(c);
	}
	isl_basic_set_free(interval_1);
	isl_basic_set_free(interval_2);
	isl_constraint_list_free(constraints_1);
	isl_constraint_list_free(constraints_2);
	return approx_set;
}

__isl_give isl_basic_set * get_shadow(__isl_take isl_basic_set * org_set, int dim1, int dim2)
{
	int i=0;
	int nDim = isl_basic_set_n_dim(org_set);	
	isl_basic_set * shadow =  org_set;
	for(i=0;i<nDim;i++)
	{
		if(i!= dim1 && i!=dim2)
		{
			shadow = isl_basic_set_project_out(shadow,3,i,1);
		}
	}
	return shadow;
}

__isl_give isl_basic_set * utvpi_over_approximation_FM(__isl_keep isl_basic_set * org_set)
{
	printf("\n\t Original polyhedra: %s",isl_basic_set_to_str(org_set));
	isl_space * space = isl_basic_set_get_space(org_set);
	isl_ctx * ctx = isl_basic_set_get_ctx(org_set);
	int i=0,j;
	isl_basic_set * approx_set =  isl_basic_set_universe(isl_space_copy(space));
	int nDim = isl_basic_set_n_dim(org_set);
	isl_basic_set * shadow, * rotated_shadow = NULL;
	for(i=0;i<nDim;i++)
	{
		for(j=i+1;j<nDim;j++)
		{
			printf("\n\t\t Dimensions considering: %s with %s ",isl_basic_set_get_dim_name(org_set,3,i), isl_basic_set_get_dim_name(org_set,3,j));
			shadow = get_shadow(isl_basic_set_copy(org_set),i,j);			
			printf("\n\t\t shadow: %s",isl_basic_set_to_str(shadow));
			approx_set = generate_interval_constraints(isl_basic_set_copy(shadow), i, j, approx_set);
			isl_basic_map * rotation_map = get_rotation_map(isl_basic_set_copy(shadow));
			printf("\n\t rotation_map %s",isl_basic_map_to_str(rotation_map));
			printf("\n\t Original shadow %s",isl_basic_set_to_str(shadow));
			rotated_shadow = isl_basic_set_apply(shadow, isl_basic_map_copy(rotation_map));
			printf("\n\t rotated shadow %s",isl_basic_set_to_str(rotated_shadow));
			fflush(stdout);
			approx_set = generate_octagonal_constraints(rotated_shadow, i, j, approx_set, isl_basic_map_reverse(rotation_map));
			isl_basic_set_free(shadow);
		}
	}
	isl_space_free(space);
	return approx_set;

}

/*
Note that Apron interface is not yet fixed.
#include</home/abhishek/Approximate_computing/apron/octagons/oct.h>
-lapron -loctMPQ -lgmp -lmpfr -lm
__isl_give isl_basic_set * interface_apron(__isl_take isl_basic_set * bset)
{
	print_poly("p",p); 
	print_oct("o2",o2); 
	p = poly_of_oct(o);
	o = oct_of_poly(p);

//isl_vertices *isl_basic_set_compute_vertices

	//ap_manager_t* man = oct_manager_alloc();
	//ap_abstract0_t* poly;
	//ap_generator0_array_t genrep =  ap_abstract0_to_generator_array(man, poly);
// ap_abstract0_t* ap_abstract0_oct_of_generator_array(ap_manager_t* man, size_t intdim, size_t realdim,ap_generator0_array_t* array)	


	//ap_abstract0_oct_of_generator_array(man,NULL ,NULL , &genrep);
	//ap_abstract0_to_generator_array()
}
*/

isl_stat build_const_search_space(__isl_take isl_constraint *c, __isl_keep void *user)
{
	//isl_basic_set * const_search_space = (isl_basic_set *) user;
	//if(isl_constraint_is_equality(c))
	//	return isl_stat_ok;
	//printf("\n\n OK here");
	return isl_stat_ok;
}

__isl_give isl_point * find_hyperplane(__isl_keep isl_basic_set * search_space_1, __isl_keep isl_basic_set * search_space_2)
{
	isl_space * ss_space_1 = isl_basic_set_get_space(search_space_1);
	isl_space * ss_space_2 = isl_basic_set_get_space(search_space_2);

	int c_cst_1_index = isl_space_find_dim_by_name(ss_space_1, 3, "c_cst");
	int c_cst_2_index = isl_space_find_dim_by_name(ss_space_2, 3, "c_cst");

	isl_space_free(ss_space_1);
	isl_space_free(ss_space_2);

	int ndim_1 = isl_basic_set_n_dim(search_space_1);
	int ndim_2 = isl_basic_set_n_dim(search_space_2);
	if(c_cst_1_index==0)
	{
		search_space_1 = isl_basic_set_project_out(search_space_1,3,1,ndim_1-1);
		isl_basic_set_set_dim_name(search_space_1, 3, c_cst_1_index, "c1" );
	printf("\n\t polyhedra: %s",isl_basic_set_to_str(search_space_1));
	}
	else
	{
		printf("\n\t unhandled case");
		//TODO: Call ISL_DIE
	}
	if(c_cst_2_index==0)
	{
		search_space_2 = isl_basic_set_project_out(search_space_2,3,1,ndim_2-1);
		isl_basic_set_set_dim_name(search_space_2, 3, c_cst_2_index, "c2" );
	printf("\n\t polyhedra: %s",isl_basic_set_to_str(search_space_2));
	}
	else
	{
		printf("\n\t unhandled case");
		//TODO: Call ISL_DIE
	}
	isl_set * ss1 = isl_set_from_basic_set(isl_set_polyhedral_hull(isl_set_from_basic_set(search_space_1)));
	isl_set * ss2 = isl_set_from_basic_set(isl_set_polyhedral_hull(isl_set_from_basic_set(search_space_2)));
	isl_set* const_search_space;
	isl_point * sol_pt;
	if(!isl_set_is_empty(ss1) && !isl_set_is_empty(ss2))
	{
		const_search_space = isl_set_product(ss1, ss2);
		const_search_space = isl_set_set_tuple_name(const_search_space,"s12");	
	}
	else if(!isl_set_is_empty(ss2))
	{
		const_search_space = ss2;
		isl_set_free(ss1);
		const_search_space = isl_set_set_tuple_name(const_search_space,"s2");
	}
	else if(!isl_set_is_empty(ss1))
	{
		const_search_space = ss1;
		isl_set_free(ss2);
		const_search_space = isl_set_set_tuple_name(const_search_space,"s1");
	}
	else
	{
		printf("\n\t No search space. Returning null point as lexmin.");
		isl_set_free(const_search_space);
		
		isl_point * sol_void = NULL;
		isl_set_free(ss2);
		isl_set_free(ss1);
		//printf("\n Here1");
		//printf(" %d %d", isl_set_is_empty(ss1), isl_set_is_empty(ss2) );
		return sol_void;
	}
	printf("\n\t Search space %s",isl_set_to_str(const_search_space));
	isl_set* solution = isl_set_lexmin(const_search_space);
	printf("\n\t Solution %s : %d",isl_set_to_str(solution), isl_set_is_singleton(solution));
	fflush(stdout);
	if(isl_set_is_singleton(solution))
	{
		sol_pt = isl_set_sample_point(solution);
		printf("\n Returning %s",isl_point_to_str(sol_pt));	
	}
	else
	{
		// TODO Unhandled case
		// Call ISL_DIE
		return NULL;
	}	
	return sol_pt;


	//isl_constraint_list_free(ss_1);
	//isl_constraint_list_free(ss_2);
	
	//if(isl_set_is_singleton(solution))
	//{
		//sol_pt = isl_set_sample_point(solution);
	//}
	//else
	//{
	//	sol_pt = NULL;
	//	printf("\n ERROR======== Solution point is NULL");
		// TODO call ISL_DIE
	//}
}



__isl_give isl_basic_set * update_OA(isl_basic_set * approx_set, isl_point * hplanes_coeff, isl_val * fix_coeff_i, isl_val * fix_coeff_j, int i, int j)
{
	// Remember fix_coeff_i are for space 2. We negate it to obtain values for space 1
	isl_space * space = isl_point_get_space(hplanes_coeff);
	const char * flag = isl_space_get_tuple_name(space, 3);
	int ptDim = isl_space_dim(space, 3);
	isl_space_free(space);
	printf("\n Check point: %s \t %s \t %d",isl_point_to_str(hplanes_coeff), flag, ptDim);
	isl_local_space * ls = isl_local_space_from_space(isl_basic_set_get_space(approx_set));
	int nDim = isl_local_space_dim(ls,3);
	isl_val * c_2 = NULL, * c_1 = NULL;
	i=i-1;
	j=j-1;
	
	if(hplanes_coeff == NULL)
	{
		//c_2 = isl_val_infty(isl_basic_set_get_ctx(approx_set));
		//c_1 = isl_val_infty(isl_basic_set_get_ctx(approx_set));	
		//printf("\n UTVPI Constraint: %s",isl_constraint_to_str(hspace_1));
		//printf("\n UTVPI Constraint: %s",isl_constraint_to_str(hspace_2));
	}
	else
	{
		isl_val * temp = isl_point_get_coordinate_val(hplanes_coeff, 3, 0);
		if(ptDim==1)
		{
			// only one UTVPI constraint possible
			if(strstr(flag, "2") != NULL)
			{
				c_2 = temp;
				c_1 = NULL;	
			}
			else if(strstr(flag, "1") != NULL)
			{
				c_1 = temp;
				c_2 = NULL;
			}
		}
		else if(ptDim==2)
		{
			c_1 = temp;
			c_2 = isl_point_get_coordinate_val(hplanes_coeff, 3, 1);
	
		}

		if(c_2)
		{
			isl_constraint * hspace_2 = generate_utvpi_constraint(isl_local_space_copy(ls), nDim, i, j, isl_val_copy(fix_coeff_i), isl_val_copy(fix_coeff_j), c_2);


			approx_set = isl_basic_set_add_constraint(approx_set,hspace_2);	
		}
		if(c_1)
		{
			fix_coeff_i = isl_val_neg(fix_coeff_i);
			fix_coeff_j = isl_val_neg(fix_coeff_j);
	
			isl_constraint * hspace_1 = generate_utvpi_constraint(isl_local_space_copy(ls), nDim, i, j, isl_val_copy(fix_coeff_i), isl_val_copy(fix_coeff_j), c_1);
	
			approx_set = isl_basic_set_add_constraint(approx_set,hspace_1);	
		}	
 	}
	isl_local_space_free(ls);
	return approx_set;
//	printf("\n coeffs: %ld %ld nDim %d %d %d",c_1, c_2, nDim, i, j);
	
}

__isl_give isl_basic_set * specialize_dual(isl_basic_set * dual_set, int i, int j)
{
	isl_basic_set * dual_specialized = isl_basic_set_copy(dual_set);
	isl_ctx * ctx = isl_basic_set_get_ctx(dual_set);
	int nDim = isl_basic_set_n_dim(dual_set);
	isl_val * fix_coeff = isl_val_int_from_si(ctx, 0);
	int k;
	for(k=1;k<nDim; k++)
	{
		if(k!=i && k!=j)
		{
			dual_specialized = isl_basic_set_fix_val(dual_specialized,3,k,isl_val_copy(fix_coeff));
		}
	}
	isl_val_free(fix_coeff);
	return dual_specialized;
}

__isl_give isl_basic_set * utvpi_over_approximation_farkas(__isl_keep isl_basic_set * org_set)
{
	printf("\n\t Original polyhedra: %s",isl_basic_set_to_str(org_set));
	isl_space * space = isl_basic_set_get_space(org_set);
	isl_ctx * ctx = isl_basic_set_get_ctx(org_set);
	int i=0,j;
	long vi, vj;
	isl_basic_set * approx_set =  isl_basic_set_universe(isl_space_copy(space));
	isl_basic_set * search_space_1, * search_space_2;
	isl_basic_set * dual_set = isl_basic_set_coefficients(isl_basic_set_copy(org_set)); 
	int nDim = isl_basic_set_n_dim(dual_set);
	// Loop will start from 1 as first dimension is constant dimension.	
	for(i=1;i<nDim;i++)
	{
		for(j=i+1;j<nDim;j++)
		{
			printf("\n\t\t Dimensions considering: %s with %s ",isl_basic_set_get_dim_name(dual_set,3,i), isl_basic_set_get_dim_name(dual_set,3,j));
			for(vi=-1;vi<1;vi++)
			{
				for(vj=-1;vj<2;vj++)
				{
					if(!((vi==0 && vj==0) || (vi==0 && vj==1)))
					{
						isl_val * fix_coeff_i = isl_val_int_from_si(ctx, vi);
						isl_val * fix_coeff_j = isl_val_int_from_si(ctx, vj);
						isl_basic_set * dual_specialized = specialize_dual(dual_set, i, j);
						search_space_1 = isl_basic_set_fix_val(isl_basic_set_copy(dual_specialized),3,i,isl_val_copy(fix_coeff_i));
						search_space_1 = isl_basic_set_fix_val(search_space_1,3,j,isl_val_copy(fix_coeff_j));	

			printf("\n\t\t Space1: %s \n\n",isl_basic_set_to_str(search_space_1));

						fix_coeff_i = isl_val_neg(fix_coeff_i);
						fix_coeff_j = isl_val_neg(fix_coeff_j);

						search_space_2 = isl_basic_set_fix_val(isl_basic_set_copy(dual_specialized),3,i,isl_val_copy(fix_coeff_i));
						search_space_2 = isl_basic_set_fix_val(search_space_2,3,j,isl_val_copy(fix_coeff_j));	
			printf("\n\t\t Space2: %s \n\n",isl_basic_set_to_str(search_space_2));

						isl_point* hplane_coeff = find_hyperplane(search_space_1, search_space_2);		
						approx_set = update_OA(approx_set,hplane_coeff, fix_coeff_i, fix_coeff_j, i,j);
						isl_point_free(hplane_coeff);
						isl_val_free(fix_coeff_i);
						isl_val_free(fix_coeff_j);
						isl_basic_set_free(dual_specialized);
					}
				}
			}			
	


			//printf("\n\t\t approx1:  %s",isl_basic_set_to_str(approx_set));
			

		}		
	}
	printf("\n\t\t Dual polyhedra: %s",isl_basic_set_to_str(dual_set));
	isl_basic_set_free(dual_set);
	isl_space_free(space);
	return approx_set;

	//isl_pw_aff_to_str("\n\t DIM:Max %s",isl_set_dim_max(isl_set_from_basic_set(isl_basic_set_copy(bset)), i));
}

__isl_give isl_basic_set * isl_basic_set_overapproximate(__isl_take isl_basic_set * bset)
{
	isl_basic_set * approxset ; //= isl_basic_set_copy(bset);
	//approxset = utvpi_over_approximation_farkas(bset);
	approxset = utvpi_over_approximation_FM(bset);	
	//approxset = interface_apron(approxset);	
	//approxset = bset;
	printf("\n\t\t Overapproximated farkas polyhedra: %s, %d",isl_basic_set_to_str(approxset), isl_basic_set_is_empty(approxset));
	return approxset;
}

