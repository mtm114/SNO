#include "Lite_Sparse_Matrix.h"

std::ostream & operator<<(std::ostream & s, Lite_Sparse_Matrix * A)
{
	s.precision(16);
	if (A == nullptr)
	{
		s << "the matrix does not exist !\n ";
	}

	const size_t row = A->rows();
	const size_t col = A->cols();
	const size_t nonzero = A->get_nonzero();
	const size_t *rowind = A->get_rowind();
	const size_t *colptr = A->get_colptr();
	const double *values = A->get_values();

	s << "row :" << row << " col :" << col << " Nonzero: " << nonzero << "\n\n";

	s << "matrix --- (i, j, value)\n\n";

	Lite_SPARSE_STORAGE_TYPE s_store = A->storage();
	int inc = (A->get_arraytype() == SPARSE_FORTRAN_TYPE ? -1 : 0);
	if (s_store == SPARSE_CCS)
	{
		size_t k = 0;
		for (size_t i = 1; i < col + 1; i++)
		{
			for (size_t j = 0; j < colptr[i] - colptr[i - 1]; j++)
			{
				s << rowind[k] + inc << " " << i - 1 << " " << std::scientific
						<< values[k] << "\n";
				k++;
			}
		}
	}
	else if (s_store == SPARSE_CRS)
	{
		size_t k = 0;
		for (size_t i = 1; i < row + 1; i++)
		{
			for (size_t j = 0; j < rowind[i] - rowind[i - 1]; j++)
			{
				s << i - 1 << " " << colptr[k] + inc << " " << std::scientific
						<< values[k] << "\n";
				k++;
			}
		}
	}
	else if (s_store == SPARSE_TRIPLE)
	{
		for (size_t k = 0; k < nonzero; k++)
		{
			s << rowind[k] + inc << " " << colptr[k] + inc << " "
					<< std::scientific << values[k] << "\n";
		}
	}

	if (A->is_diag_saved())
	{
		double *diag = A->get_diag();
		const size_t diag_size = std::min(row, col);
		for (size_t k = 0; k < diag_size; k++)
		{
			s << k << " " << k <<" " << diag[k] << "\n";
		}
	}

	return s;
}
