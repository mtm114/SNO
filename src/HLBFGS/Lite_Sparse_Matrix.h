#ifndef LITE_SPARSE_MATRIX_H
#define LITE_SPARSE_MATRIX_H

#include <vector>
#include <algorithm>
#include <cassert>
#include <ostream>
#include <cstring>
#include "Lite_Sparse_Entry.h"

// \addtogroup MathSuite
//@{
//! Symmetric status
enum Lite_SPARSE_SYMMETRIC_STATE
{
	SPARSE_NOSYM, /*!< general case */
	SPARSE_SYM_UPPER, /*!< symmetric (store upper triangular part) */
	SPARSE_SYM_LOWER, /*!< symmetric (store lower triangular part) */
	SPARSE_SYM_BOTH   /*!< symmetric (store both upper and lower triangular part) */
};

//!     Storage
enum Lite_SPARSE_STORAGE_TYPE
{
	SPARSE_CCS, /*!< compress column format */
	SPARSE_CRS, /*!< compress row format */
	SPARSE_TRIPLE /*!< row-wise coordinate format */
};

//! Array type
enum Lite_SPARSE_ARRAY_TYPE
{
	SPARSE_FORTRAN_TYPE, /*!< the index starts from 1 */
	SPARSE_C_TYPE /*!< the index starts from 0 */
};

//////////////////////////////////////////////////////////////////////////
//! Lite Sparse Matrix Class
class Lite_Sparse_Matrix
{
private:

	//!     Status for creating sparse solver
	enum STORE_STATE
	{
		ENABLE, DISABLE, LOCK
	};

	STORE_STATE state_fill_entry;
	Lite_SPARSE_SYMMETRIC_STATE sym_state;
	Lite_SPARSE_STORAGE_TYPE s_store;
	Lite_SPARSE_ARRAY_TYPE arraytype;

	size_t nrows; //!< number of rows
	size_t ncols; //!< number of columns
	size_t nonzero; //!< number of nonzeros
	//! pointers to where columns begin in rowind and values 0-based, length is (col+1)
	/*!
	 * When s_store is CRS, colptr stores column indices;
	 */
	std::vector<size_t> colptr;
	//! row indices, 0-based
	/*!
	 * When s_store is CRS, rowind stores row-pointers
	 */
	std::vector<size_t> rowind;
	std::vector<double> values; //!< nonzero values of the sparse matrix

	std::vector<std::vector<Lite_Sparse_Entry> > entryset; //!< store temporary sparse entries

	std::vector<double> diag; //! special usage for some libraries

	bool save_diag_separetely;

public:

	//! Sparse matrix constructor
	/*!
	 * \param m row dimension
	 * \param n column dimension
	 * \param symmetric_state
	 * \param m_store the storage format
	 * \param atype Fortran or C type of array
	 */
	Lite_Sparse_Matrix(size_t m, size_t n, Lite_SPARSE_SYMMETRIC_STATE symmetric_state =
		SPARSE_NOSYM, Lite_SPARSE_STORAGE_TYPE m_store = SPARSE_CCS, Lite_SPARSE_ARRAY_TYPE atype = SPARSE_C_TYPE,
		bool save_diag = false) :
		state_fill_entry(DISABLE), sym_state(symmetric_state), s_store(m_store), arraytype(atype),
		nrows(m), ncols(n), nonzero(0), save_diag_separetely(save_diag)
	{
		if (m != n)
		{
			symmetric_state = SPARSE_NOSYM;
		}

		size_t nn = (m_store == SPARSE_CCS ? ncols : nrows);
		entryset.resize(nn);
		if (save_diag_separetely)
		{
			diag.resize(nrows < ncols ? nrows : ncols);
			//std::fill(diag.begin(), diag.end(), 0.0);
			memset(&diag[0], 0, sizeof(double) * diag.size());
		}
	}

	//! Sparse matrix destructor
	~Lite_Sparse_Matrix()
	{
		clear_mem();
	}
	//! Start to build sparse matrix pattern
	inline void begin_fill_entry()
	{
		state_fill_entry = ENABLE;
	}
	//! Construct sparse pattern
	void end_fill_entry()
	{
		assert(state_fill_entry == ENABLE);

		clear_mem();

		state_fill_entry = LOCK;

		int inc = (arraytype == SPARSE_FORTRAN_TYPE ? 1 : 0);

		if (s_store == SPARSE_CCS)
		{
			//construct map and ccs matrix
			size_t i, j, k = 0;
			colptr.resize(ncols + 1);
			colptr[0] = inc;
			for (j = 1; j < ncols + 1; j++)
			{
				colptr[j] = (int)entryset[j - 1].size() + colptr[j - 1];
			}

			nonzero = colptr[ncols];

			if (nonzero > 0)
			{
				rowind.resize(nonzero);
				values.resize(nonzero);

				for (j = 0; j < ncols; j++)
				{
					for (i = 0; i < colptr[j + 1] - colptr[j]; i++)
					{
						rowind[k] = entryset[j][i].index + inc;
						values[k] = entryset[j][i].value;
						k++;
					}
				}
			}
		}
		else if (s_store == SPARSE_CRS)
		{
			//construct map and crs matrix
			size_t i, j, k = 0;
			rowind.resize(nrows + 1);
			rowind[0] = inc;
			for (j = 1; j < nrows + 1; j++)
			{
				rowind[j] = (int)entryset[j - 1].size() + rowind[j - 1];
			}
			nonzero = rowind[nrows];
			if (nonzero > 0)
			{
				colptr.resize(nonzero);
				values.resize(nonzero);

				for (j = 0; j < nrows; j++)
				{
					for (i = 0; i < rowind[j + 1] - rowind[j]; i++)
					{
						colptr[k] = entryset[j][i].index + inc;
						values[k] = entryset[j][i].value;
						k++;
					}
				}
			}
		}
		else if (s_store == SPARSE_TRIPLE)
		{
			size_t i, j, k = 0;
			nonzero = 0;
			for (i = 0; i < nrows; i++)
			{
				nonzero += (int)entryset[i].size();
			}

			if (nonzero > 0)
			{
				rowind.resize(nonzero);
				colptr.resize(nonzero);
				values.resize(nonzero);

				for (i = 0; i < nrows; i++)
				{
					size_t jsize = entryset[i].size();
					for (j = 0; j < jsize; j++)
					{
						rowind[k] = i + inc;
						colptr[k] = entryset[i][j].index + inc;
						values[k] = entryset[i][j].value;
						k++;
					}
				}
			}
		}
		entryset.clear();
	}

	//! Fill matrix entry \f$  Mat_{row_index, col_index} += val \f$
	void fill_entry(size_t row_index, size_t col_index, double val = 0)
	{
		if (row_index >= nrows || col_index >= ncols)
			return;

		if (save_diag_separetely && row_index == col_index)
		{
			diag[row_index] += val;
			return;
		}

		if (sym_state == SPARSE_NOSYM)
		{
			fill_entry_internal(row_index, col_index, val);
		}
		else if (sym_state == SPARSE_SYM_UPPER)
		{
			if (row_index <= col_index)
			{
				fill_entry_internal(row_index, col_index, val);
			}
			else
			{
				fill_entry_internal(col_index, row_index, val);
			}
		}
		else if (sym_state == SPARSE_SYM_LOWER)
		{
			if (row_index <= col_index)
			{
				fill_entry_internal(col_index, row_index, val);
			}
			else
			{
				fill_entry_internal(row_index, col_index, val);
			}
		}
		else if (sym_state == SPARSE_SYM_BOTH)
		{
			fill_entry_internal(row_index, col_index, val);

			if (row_index != col_index)
			{
				fill_entry_internal(col_index, row_index, val);
			}
		}
	}

	//fill the diagonal entry
	inline void fill_diag(size_t diagid, double v = 0)
	{
		if (save_diag_separetely)
		{
			if (diag.size() == 0)
			{
				diag.resize(nrows < ncols ? nrows : ncols);
				//std::fill(diag.begin(), diag.end(), 0.0);
				memset(&diag[0], 0, sizeof(double) * diag.size());
			}
			diag[diagid] += v;
		}
	}

	//! get the number of nonzeros
	inline size_t get_nonzero()
	{
		return nonzero;
	}
	//! get the row dimension
	inline size_t rows()
	{
		return nrows;
	}
	//! get the column dimension
	inline size_t cols()
	{
		return ncols;
	}
	//! return the symmetric state
	inline bool issymmetric()
	{
		return sym_state != SPARSE_NOSYM;
	}

	//! tell whether the matrix is upper or lower symmetric
	inline bool issym_store_upper_or_lower()
	{
		return (sym_state == SPARSE_SYM_LOWER) || (sym_state == SPARSE_SYM_UPPER);
	}

	//! return symmetric state
	inline Lite_SPARSE_SYMMETRIC_STATE symmetric_state()
	{
		return sym_state;
	}

	//! tell whether the matrix is square
	inline bool issquare()
	{
		return nrows == ncols;
	}

	//! return the storage format
	inline Lite_SPARSE_STORAGE_TYPE storage()
	{
		return s_store;
	}

	//! return array type
	inline Lite_SPARSE_ARRAY_TYPE get_arraytype()
	{
		return arraytype;
	}

	//! get rowind
	inline size_t *get_rowind()
	{
		return rowind.empty() ? nullptr : &rowind[0];
	}

	inline const size_t *get_rowind() const
	{
		return rowind.empty() ? nullptr : &rowind[0];
	}

	//! get colptr
	inline size_t *get_colptr()
	{
		return colptr.empty() ? nullptr : &colptr[0];
	}

	inline const size_t *get_colptr() const
	{
		return colptr.empty() ? nullptr : &colptr[0];
	}

	//! get the values array
	inline double *get_values()
	{
		return values.empty() ? nullptr : &values[0];
	}

	inline const double *get_values() const
	{
		return values.empty() ? nullptr : &values[0];
	}

	//! get the diagonal array
	inline double *get_diag()
	{
		return diag.empty() ? nullptr : &diag[0];
	}

	inline const double *get_diag() const
	{
		return diag.empty() ? nullptr : &diag[0];
	}

	inline const bool is_diag_saved() const
	{
		return save_diag_separetely;
	}

	//////////////////////////////////////////////////////////////////////////
private:
	//! Clear memory
	inline void clear_mem()
	{
		colptr.clear();
		rowind.clear();
		values.clear();
	}

	//! fill matrix entry (internal) \f$ Mat[rowid][colid] += val \f$
	bool fill_entry_internal(size_t row_index, size_t col_index, double val = 0)
	{
		assert(state_fill_entry == ENABLE);

		size_t search_index = (s_store == SPARSE_CCS ? row_index : col_index);
		size_t pos_index = (s_store == SPARSE_CCS ? col_index : row_index);

		Lite_Sparse_Entry forcompare(search_index);

		std::vector<Lite_Sparse_Entry>::iterator iter = std::lower_bound(
			entryset[pos_index].begin(), entryset[pos_index].end(),
			forcompare);
		if (iter != entryset[pos_index].end())
		{
			if (iter->index == search_index)
			{
				iter->value += val;
			}
			else
				entryset[pos_index].insert(iter,
				Lite_Sparse_Entry(search_index, val));
		}
		else
		{
			entryset[pos_index].push_back(Lite_Sparse_Entry(search_index, val));
		}
		return true;
	}
	//////////////////////////////////////////////////////////////////////////
};

//! print sparse matrix
std::ostream & operator<<(std::ostream & s, const Lite_Sparse_Matrix * A);

//@}

#endif //Lite_Sparse_Matrix_H
