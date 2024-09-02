#ifndef Lite_Sparse_Entry_H
#define Lite_Sparse_Entry_H

//! Lite Sparse Entry class \ingroup MathSuite
class Lite_Sparse_Entry
{
public:

	//! Index ID
	size_t index;

	//! Real value
	double value;

public:
	//! constructor
	Lite_Sparse_Entry(size_t ind, double v = 0) :
		index(ind), value(v)
	{
	}

	//! destructor
	~Lite_Sparse_Entry()
	{
	}

	//! The compare function for sorting
	inline bool operator<(const Lite_Sparse_Entry & m_r) const
	{
		return index < m_r.index;
	}
};

#endif //Lite_Sparse_Entry_H
