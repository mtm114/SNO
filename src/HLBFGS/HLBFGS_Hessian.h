#ifndef HLBFGS_HESSIAN_H
#define HLBFGS_HESSIAN_H

#include <vector>
#include "Lite_Sparse_Matrix.h"

//! ICFS_INFO stores ICFS's working arrays
class ICFS_INFO
{
public:
	ICFS_INFO()
	{
		p = 15;
	}
	~ICFS_INFO()
	{
	}
	void allocate_mem(const size_t N)
	{
		if (N > 0)
		{
			lcol_ptr.resize(N + 1);
			ldiag.resize(N);
			iwa.resize(3 * N);
			wa1.resize(N);
			wa2.resize(N);
			r.resize(N);
			p = 15;
			CGP.resize(N);
			CGR.resize(N);
			CGQ.resize(N);
			CGZ.resize(N);
		}
	}
	size_t * get_lcol_ptr()
	{
		return &lcol_ptr[0];
	}
	size_t * get_lrow_ind()
	{
		return &lrow_ind[0];
	}
	double * get_ldiag()
	{
		return &ldiag[0];
	}
	double * get_l()
	{
		return &l[0];
	}
	size_t * get_iwa()
	{
		return &iwa[0];
	}
	double * get_wa1()
	{
		return &wa1[0];
	}
	double * get_wa2()
	{
		return &wa2[0];
	}
	size_t & get_p()
	{
		return p;
	}
	double * get_r()
	{
		return &r[0];
	}
	double * get_CGP()
	{
		return &CGP[0];
	}
	double * get_CGQ()
	{
		return &CGQ[0];
	}
	double * get_CGR()
	{
		return &CGR[0];
	}
	double * get_CGZ()
	{
		return &CGZ[0];
	}
	double & get_icfs_alpha()
	{
		return icfs_alpha;
	}
	void set_lrow_ind_size(const size_t size)
	{
		lrow_ind.resize(size);
	}
	void set_l_size(const size_t size)
	{
		l.resize(size);
	}
private:
	std::vector<size_t> lcol_ptr;
	std::vector<size_t> lrow_ind;
	std::vector<double> ldiag;
	std::vector<double> l;
	std::vector<size_t> iwa;
	std::vector<double> wa1;
	std::vector<double> wa2;
	size_t p;
	std::vector<double> r;
	std::vector<double> CGP;
	std::vector<double> CGQ;
	std::vector<double> CGR;
	std::vector<double> CGZ;
	double icfs_alpha;
};
//////////////////////////////////////////////////////////////////////////
//! Stores the pointers of hessian matrix
class HESSIAN_MATRIX
{
public:
	HESSIAN_MATRIX()
	{
		mat = 0;
	}
	~HESSIAN_MATRIX()
	{
		if (mat)
		{
			delete mat; mat = 0;
		}
	}
	size_t get_dimension()
	{
		return mat->rows();
	}
	size_t get_nonzeros()
	{
		return mat->get_nonzero();
	}
	double * get_values()
	{
		return mat->get_values();
	}
	size_t * get_rowind()
	{
		return mat->get_rowind();
	}
	size_t * get_colptr()
	{
		return mat->get_colptr();
	}
	double * get_diag()
	{
		return mat->get_diag();
	}
	ICFS_INFO& get_icfs_info()
	{
		return l_info;
	}
	Lite_Sparse_Matrix * create_mat(const size_t dim)
	{
		if (mat)
		{
			delete mat; mat = 0;
		}
		mat = new Lite_Sparse_Matrix(dim, dim, SPARSE_SYM_LOWER, SPARSE_CCS, SPARSE_FORTRAN_TYPE, true);
		return mat;
	}
private:
	ICFS_INFO l_info;
	Lite_Sparse_Matrix *mat;
};

#endif
