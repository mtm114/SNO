#ifndef HLBFGS_H
#define HLBFGS_H

//#define USE_OPENMP

#include <vector>
#include <map>
#include <cstring>

#include "HLBFGS_Hessian.h"

typedef void(*eval_funcgrad_fp)(const size_t n_vars, const std::vector<double>& variables,
	double &func_value, std::vector<double>& gradient, void *user_pointer);
typedef void(*eval_hessian_fp)(const size_t n_vars, const std::vector<double>& variables,
	double &func_value, std::vector<double>& gradient, HESSIAN_MATRIX& hessian,
	void *user_pointer);
typedef void(*eval_constraints_fp)(const size_t n_eqns, const size_t n_ieqns,
	const std::vector<double>& variables, std::vector<double>& func_values, std::vector< std::vector<std::pair<size_t, double> > >& constraint_jacobian,
	void *user_pointer);
//typedef void(*newiteration_fp)(const size_t niter, const size_t call_iter,
//	const size_t n_vars, const std::vector<double>& variables, const double &func_value,
//	const std::vector<double>& gradient, const double &gnorm, void *user_pointer);
typedef void(*newiteration_constraints_fp)(const size_t niter, const size_t n_vars, const std::vector<double>& variables,
	const double f, const std::vector<double> &constraint_values, void *user_pointer);
typedef void(*newiteration_fp)(const size_t n_vars, const std::vector<double>& variables,
	double &func_value, std::vector<double>& gradient, void *user_pointer);

class HLBFGS
{
public:
	//////////////////////////////////////////////////////////////////////////
	HLBFGS();
	~HLBFGS();
	//////////////////////////////////////////////////////////////////////////
	inline void set_number_of_variables(size_t N)
	{
		n_N_ = N;
	}
	inline size_t get_num_of_variables()
	{
		return n_N_;
	}
	inline void set_M(size_t M)
	{
		n_M_ = M;
	}
	inline size_t get_M()
	{
		return n_M_;
	}
	inline void set_T(size_t T)
	{
		info[6] = T;
	}
	inline size_t get_T()
	{
		return info[6];
	}
	inline void set_max_constraint(double val)
	{
		parameters[9] = val;
	}
	inline void set_init_sigma(double val)
	{
		parameters[10] = val;
	}
	inline void set_number_of_equalities(size_t Ne)
	{
		n_E_ = Ne;
	}
	inline size_t get_number_of_equalities()
	{
		return n_E_;
	}
	inline void set_number_of_inequalities(size_t Ne)
	{
		n_IE_ = Ne;
	}
	inline size_t get_number_of_inequalities()
	{
		return n_IE_;
	}
	//////////////////////////////////////////////////////////////////////////
	inline void set_func_callback(eval_funcgrad_fp fp, eval_hessian_fp hfp = 0,
		eval_constraints_fp cfp = 0, newiteration_fp nfp = 0, newiteration_constraints_fp cnfp = 0)
	{
		evalfuncgrad_callback = fp;
		evalhessian_callback = hfp;
		eval_constraints_callback = cfp;
		newiteration_callback = nfp;
		newiteration_constraints_callback = cnfp;
	}
	//////////////////////////////////////////////////////////////////////////
	inline void get_advanced_setting(double hlbfgs_parameters[], size_t hlbfgs_info[])
	{
		memcpy(&hlbfgs_parameters[0], &parameters[0], sizeof(double) * 20);
		memcpy(&hlbfgs_info[0], &info[0], sizeof(size_t) * 20);
	}
	//////////////////////////////////////////////////////////////////////////
	inline void set_advanced_setting(double hlbfgs_parameters[], size_t hlbfgs_info[])
	{
		memcpy(&parameters[0], &hlbfgs_parameters[0], sizeof(double) * 20);
		memcpy(&info[0], &hlbfgs_info[0], sizeof(size_t) * 20);
	}
	//////////////////////////////////////////////////////////////////////////
	inline void set_verbose(bool verbose)
	{
		verbose_ = verbose;
	}
	//////////////////////////////////////////////////////////////////////////
	void optimize_without_constraints(double *init_sol, size_t max_iter, void *user_pointer = 0, int *it = 0);
	void optimize_with_constraints(double *init_sol, size_t max_iter, size_t local_iter = 30, void *user_pointer = 0, bool warm_start = false, double *init_lambda = 0, double *init_sigma = 0);
	inline void set_weak_wolf(bool weak_wolf) { weak_wolf_ = weak_wolf; }
	inline void set_qp_solver(bool use_qp) { qp_solver_ = use_qp; }
	void copy_lambda(double *lambda_copy) { if (lambda_copy) memcpy(lambda_copy, &lambda[0], sizeof(double) * lambda.size()); }
	void copy_sigma(double *sigma_copy) { if (sigma_copy) *sigma_copy = sigma; }
	void set_Info(int it, int ini) { info[5] = ini; info[16] = it; }
protected:
	void mem_allocation();
	bool self_check();
	void initialize_setting();
	void print_message(bool print, size_t id);

	void update_first_step(double *q, double *s, double *y, double *rho,
		double *alpha, size_t bound, size_t curpos, size_t iter);
	void update_hessian(double *q, double *s, double *y, size_t curpos,
		double *diag);
	void update_second_step(double *q, double *s, double *y, double *rho,
		double *alpha, size_t bound, size_t curpos, size_t iter);
	void conjugate_gradient_update(double *q, double *prev_q_update,
		double *prev_q_first_stage);

	void update_q_by_inverse_hessian();

	void build_hessian_info();

	void my_eval_funcgrad(double &f, double *guess_sigma, void *user_pointer);

	double max_absolute_value();

protected:
	size_t n_N_, n_M_, n_E_, n_IE_;

	eval_funcgrad_fp evalfuncgrad_callback;
	eval_hessian_fp evalhessian_callback;
	eval_constraints_fp eval_constraints_callback;
	newiteration_fp newiteration_callback;
	newiteration_constraints_fp newiteration_constraints_callback;
	double parameters[20];
	size_t info[20];
	bool verbose_;
	bool weak_wolf_;
	bool qp_solver_;

	std::vector<double> variables_;
	std::vector<double> gradient_;

	size_t cur_pos;
	std::vector<double> q_vec, alpha_vec, rho_vec, s_vec, y_vec,
		prev_x_vec, prev_g_vec, diag_vec, wa_vec;

	HESSIAN_MATRIX m_hessian;

	std::vector<double> prev_q_first_stage_vec, prev_q_update_vec;

	std::vector<double> lambda, d;
	double sigma;
	std::vector<double> dynamic_sigma, prev_K;
	std::vector<double> constraint_values_;
	std::vector< std::vector<std::pair<size_t, double> > > constraint_jacobian_;
	double object_func_;
};

#endif
