//
//  dimacs.c
//

#include <assert.h>
#include <stdio.h>
#include <stdbool.h>

#include "cnf.h"
#include "dimacs.h"
#include "wdsat_utils.h"
#include "xorset.h" //needed for counting sat

static bool _dimacs_header_read = false;
static bool _dimacs_read = false;

/// Number of atoms
static int_t dimacs_nb_of_vars = 0LL;
/// Number of rounds
static int_t dimacs_nb_of_rounds = 0LL;
/// Number of models
static int_t dimacs_nb_of_models = 0LL;
/// Number of atoms
static int_t dimacs_nb_of_atoms_be = 0LL;
static int_t dimacs_nb_of_atoms_xe = 0LL;

/// from meaning
static char * dimacs_meaning[__ID_SIZE__];
static byte_t dimacs_truth_value[__ID_SIZE__];

/// number of boolean equations
static int_t dimacs_nb_of_boolean_equations;
static int_t dimacs_boolean_equation[__MAX_EQ__][__MAX_EQ_SIZE__];
static int_t dimacs_size_of_boolean_equation[__MAX_EQ__];
static int_t dimacs_size_max_among_boolean_equations;
static int_t dimacs_size_min_among_boolean_equations;

/// number of xor equations
static int_t dimacs_nb_of_xor_equations;
static int_t dimacs_xor_equation[__MAX_XEQ__][__MAX_XEQ_SIZE__];
static int_t dimacs_size_of_xor_equation[__MAX_XEQ__];
static int_t dimacs_size_max_among_xor_equations;
static int_t dimacs_size_min_among_xor_equations;
static bool dimacs_xor_equation_constant[__MAX_XEQ__];

extern byte_t thread_digits;

static int_t dimacs_nb_of_unary_vars = 0LL;
inline const int_t dimacs_nb_unary_vars() { return dimacs_nb_of_unary_vars; }

#ifdef __XG_ENHANCED__
static int_t dimacs_nb_of_eq = 0LL;
inline const int_t _dimacs_nb_of_eq() { return dimacs_nb_of_eq; } //real count of ANF clauses
//monomials degree > 1
uint_t (*monomials_to_column)[__ID_SIZE__][__MAX_DEGREE__ - 1]; //uint_t monomials_to_column[__MAX_ANF_ID__][__ID_SIZE__][__MAX_DEGREE__ - 1];

boolean_t dimacs_current_degree[__ID_SIZE__];
#endif

/* ----------------- Getters, Setters */
inline const int_t dimacs_nb_vars() { return dimacs_nb_of_vars; }
inline const int_t dimacs_nb_equations() { return dimacs_nb_of_boolean_equations; }
inline const int_t dimacs_nb_atoms() { return dimacs_nb_of_atoms_be; }
inline const int_t dimacs_nb_xor_equations() { return dimacs_nb_of_xor_equations; }
inline const int_t dimacs_nb_xor_atoms() { return dimacs_nb_of_atoms_xe; }
inline int_t * dimacs_equation(const int_t i) { return dimacs_boolean_equation[i]; }
inline int_t * dimacs_size_of_equations() { return dimacs_size_of_boolean_equation; }
inline const int_t dimacs_size_of_equation(const int_t i) { return dimacs_size_of_boolean_equation[i]; }
inline int_t * dimacs_xor(const int_t i) { return dimacs_xor_equation[i]; }
inline int_t * dimacs_size_of_xor_equations() { return dimacs_size_of_xor_equation; }
inline const int_t dimacs_size_of_xor(const int_t i) { return dimacs_size_of_xor_equation[i]; }
inline const bool dimacs_is_header_read() { return(_dimacs_header_read); }
inline const bool dimacs_is_read() { return(_dimacs_read); }
inline int_t ** get_dimacs_xor_equation() { return dimacs_xor_equation; }
inline bool dimacs_xor_constant(const int_t i) { return dimacs_xor_equation_constant[i]; }
#ifdef __XG_ENHANCED__
inline boolean_t dimacs_get_current_degree(const int_t i) { return dimacs_current_degree[i]; }
#endif
/* ----------------- */

void dimacs_generate_meaning() {
	assert(dimacs_nb_of_vars > 0);
	int_t i, n_vsz = fast_int_log10(dimacs_nb_of_vars);
	char * tmp_str;
	dimacs_meaning[0] = NULL;
	for(i = 1; i <= dimacs_nb_of_vars; ++i) {
		tmp_str = NULL;
		_get_mem(char, tmp_str, (n_vsz + 2));
		sprintf(tmp_str, "x%0*lld", (int) n_vsz, i);
		dimacs_meaning[i] = tmp_str;
		dimacs_truth_value[i] = __UNDEF__;
	}
}

void dimacs_read_formula(FILE *f) {
	char str_clause[__STATIC_CLAUSE_STRING_SIZE__] = {0};
	int_t l, i, j, deg, a_temp, d, k;
	int_t a[__MAX_DEGREE__];
	char *str_l;
	
	if(f == NULL) return;
	
	//set all values to 0. 0 will be used to delimit a list
	//this can cause a problem if __MAX_DEGREE__ is not set properly
	memset(monomials_to_column, 0, sizeof(uint_t)*__MAX_ANF_ID__*__ID_SIZE__*(__MAX_DEGREE__ - 1));
	
	fgets(str_clause, sizeof(str_clause), f);
	str_l = strtok (str_clause, " ");
	str_l = strtok (NULL, " ");
	str_l = strtok (NULL, " ");
	dimacs_nb_of_unary_vars = atoi(str_l);
	str_l = strtok (NULL, " ");
	dimacs_nb_of_eq = atoi(str_l);
	
	assert(dimacs_nb_of_unary_vars <= __MAX_ANF_ID__);
	
	/// start reading the file
	dimacs_nb_of_xor_equations = 0ULL;
	memset((void *) dimacs_size_of_xor_equation, 0, __MAX_XEQ__ * sizeof(int_t));
	
	dimacs_nb_of_vars = dimacs_nb_of_unary_vars;
	while (fgets(str_clause, sizeof(str_clause), f) != NULL )
	{
		assert(dimacs_nb_of_xor_equations < dimacs_nb_of_eq);
		dimacs_xor_equation_constant[dimacs_nb_of_xor_equations] = false;
		str_l = strtok (str_clause, " ");
		if(str_l[0] == 'x')
		{
			str_l = strtok (NULL, " ");
			while(str_l != NULL)
			{
				if(str_l[0] == '0')
					break;
				if(strcmp(str_l, "T") == 0)
				{
					dimacs_xor_equation_constant[dimacs_nb_of_xor_equations] = (dimacs_xor_equation_constant[dimacs_nb_of_xor_equations] + 1) % 2;
					str_l = strtok (NULL, " ");
					continue;
				}
				if(str_l[0] == '.')
				{
					deg = str_l[1] - '0'; //att if __MAX_DEGREE__ > 9
					assert(deg < __MAX_DEGREE__);
					for(i = 0; i < deg; i++)
					{
						str_l = strtok (NULL, " ");
						assert(str_l != NULL);
						a_temp = atoi(str_l);
						assert(a_temp <= dimacs_nb_of_unary_vars);
						for(j = i; j > 0; j--)
						{
							if(a[j - 1] > a_temp)
							{
								a[j] = a[j - 1];
							}
							else
							{
								break;
							}
						}
						a[j] = a_temp;
					}
					
					//check if monomial already occured
					l = 0;
					i = 0;
					while(monomials_to_column[a[0]][i][0])
					{
						for(j = 1; j < deg; j++)
						{
							if(monomials_to_column[a[0]][i][j] != a[j])
							{
								break;
							}
						}
						if(j == deg)
						{
							if(j == __MAX_DEGREE__ - 1 || monomials_to_column[a[0]][i][j] == 0) //check if right degree monomial
							{
								l = monomials_to_column[a[0]][i][0];
								break;
							}
						}
						i++;
					}
					
					if(!l)
					{
						l = ++dimacs_nb_of_vars;
						assert(dimacs_nb_of_vars <= __ID_SIZE__);
						for(j = 0; j < deg; j++)
						{
							i = 0;
							while(monomials_to_column[a[j]][i][0] > 0) i++;
							monomials_to_column[a[j]][i][0] = l;
							k = 1;
							for(d = 0; d < deg; d++)
							{
								if(d != j)
								{
									monomials_to_column[a[j]][i][k] = a[d];
									k++;
								}
							}
							//cnf part
							dimacs_boolean_equation[dimacs_nb_of_boolean_equations + j][0] = -l;
							++dimacs_nb_of_atoms_be;
							dimacs_boolean_equation[dimacs_nb_of_boolean_equations + j][1] = a[j];
							++dimacs_nb_of_atoms_be;
							dimacs_size_of_boolean_equation[dimacs_nb_of_boolean_equations + j] = 2;
							
							dimacs_boolean_equation[dimacs_nb_of_boolean_equations + deg][j] = -a[j];
							++dimacs_nb_of_atoms_be;
						}
						dimacs_size_of_boolean_equation[dimacs_nb_of_boolean_equations + deg] = deg + 1;
						dimacs_boolean_equation[dimacs_nb_of_boolean_equations + deg][deg] = l;
						++dimacs_nb_of_atoms_be;
						assert(dimacs_nb_of_atoms_be < __MAX_BUFFER_SIZE__);
						dimacs_nb_of_boolean_equations += (deg + 1);
						assert(dimacs_nb_of_boolean_equations < __MAX_EQ__);
						assert(deg + 1 < __MAX_EQ_SIZE__);
						if(deg + 1 > dimacs_size_max_among_boolean_equations) dimacs_size_max_among_boolean_equations = deg + 1;
						if((!dimacs_size_min_among_boolean_equations) || (2 < dimacs_size_min_among_boolean_equations))
							dimacs_size_min_among_boolean_equations = 2;
						dimacs_current_degree[l] = deg;
						
						//transform to dimacs for other solvers
						/*
						for(j = 0; j < deg; j++)
						{
							fprintf(stderr,"%lld %lld 0\n", -l, a[j]);
						}
						fprintf(stderr,"%lld ", l);
						for(j = 0; j < deg; j++)
						{
							fprintf(stderr,"%lld ", -a[j]);
						}
						fprintf(stderr,"0\n");
						*/
					}
				}
				else
				{
					l = atoi(str_l);
					assert(l <= dimacs_nb_of_unary_vars);
				}
				dimacs_xor_equation[dimacs_nb_of_xor_equations][dimacs_size_of_xor_equation[dimacs_nb_of_xor_equations]++] = l;
				++dimacs_nb_of_atoms_xe;
				assert(dimacs_nb_of_atoms_xe < __MAX_BUFFER_SIZE__);
				assert(dimacs_size_of_xor_equation[dimacs_nb_of_xor_equations] < __MAX_XEQ_SIZE__);
				
				str_l = strtok (NULL, " ");
			}
			if(dimacs_size_of_xor_equation[dimacs_nb_of_xor_equations] > dimacs_size_max_among_xor_equations) dimacs_size_max_among_xor_equations = dimacs_size_of_xor_equation[dimacs_nb_of_xor_equations];
			if((!dimacs_size_min_among_xor_equations) || (dimacs_size_of_xor_equation[dimacs_nb_of_xor_equations] < dimacs_size_min_among_xor_equations)) dimacs_size_min_among_xor_equations = dimacs_size_of_xor_equation[dimacs_nb_of_xor_equations];
			if(dimacs_xor_constant(dimacs_nb_of_xor_equations))
			{
				dimacs_xor_equation[dimacs_nb_of_xor_equations][0] = -dimacs_xor_equation[dimacs_nb_of_xor_equations][0];
			}
			assert(dimacs_nb_of_xor_equations < __MAX_XEQ__);
			dimacs_nb_of_xor_equations++;
		}
		else
		{
			dimacs_size_of_boolean_equation[dimacs_nb_of_boolean_equations] = 1;
			dimacs_boolean_equation[dimacs_nb_of_boolean_equations++][0] = atoi(str_l);
			++dimacs_nb_of_atoms_be;
			dimacs_size_min_among_boolean_equations = 1;
			assert(dimacs_nb_of_atoms_be < __MAX_BUFFER_SIZE__);
			assert(dimacs_nb_of_boolean_equations < __MAX_EQ__);
			str_l = strtok (NULL, " ");
			assert(str_l[0] == '0');
		}
	}
	assert(dimacs_nb_of_vars <= __ID_SIZE__);
}

void dimacs_read_header(FILE *f) {
	//dm: add read header here:
	char str_clause[__STATIC_CLAUSE_STRING_SIZE__] = {0};
	char *str_l;
	
	if(f == NULL) return;
	
	fgets(str_clause, sizeof(str_clause), f);
	str_l = strtok (str_clause, " ");
	str_l = strtok (NULL, " ");
	str_l = strtok (NULL, " ");
	dimacs_nb_of_unary_vars = atoi(str_l);
	str_l = strtok (NULL, " ");
	dimacs_nb_of_eq = atoi(str_l);
	printf("c dimacs_nb_of_unary_vars: %lld\n", dimacs_nb_of_unary_vars);
	printf("c dimacs_nb_of_eq: %lld\n", dimacs_nb_of_eq);
	// now we can set:  __MAX_XEQ__            = 648  c __MAX_ANF_ID__         = 329
	//D_MAX_ANF_ID = dimacs_nb_of_unary_vars + 1;
	//D_MAX_EQ     = dimacs_nb_of_eq;
	return;

}


void dimacs_print_header() {
	int_t i;
	int_t n_vsz = fast_int_log10(dimacs_nb_of_vars);
	_cid_cout("#atoms : %d\n", dimacs_nb_of_vars);
	_cid_cout("#rounds: %d\n", dimacs_nb_of_rounds);
	_cid_cout("#models: %d\n", dimacs_nb_of_models);
	_cid_cout("#IDs:%s\n", "");
	for(i = 1LL; i <= dimacs_nb_of_vars; ++i) {
		if(dimacs_meaning[i] == NULL) continue;
		_cid_cout("   %0*d [%s]: %s\n", n_vsz, i, dimacs_meaning[i], (dimacs_truth_value[i] == __UNDEF__) ? "undef" :  (dimacs_truth_value[i] == __ON__) ? "true" : (dimacs_truth_value[i] == __OFF__) ? "false" : "n/a");
	}
}

void dimacs_print_formula() {
	int_t i, j;
	int_t n_besz = fast_int_log10(dimacs_nb_of_boolean_equations);
	int_t n_xesz = fast_int_log10(dimacs_nb_of_xor_equations);
	int_t mx_xesz = fast_int_log10(dimacs_size_max_among_xor_equations);
	int_t mx_besz = fast_int_log10(dimacs_size_max_among_boolean_equations);
	
	_cid_cout("#atoms         : %d\n", dimacs_nb_of_vars);
	_cid_cout("#OR equations  : %d\n", dimacs_nb_of_boolean_equations);
	_cid_cout("#XOR equations : %d\n", dimacs_nb_of_xor_equations);
	_cid_cout("#EQ: [%llu]\n", dimacs_nb_of_boolean_equations);
	
	printf("p cnf %lld %lld\n",dimacs_nb_of_vars,dimacs_nb_of_boolean_equations+dimacs_nb_of_xor_equations);
	for(i = 0LL; i < dimacs_nb_of_boolean_equations; ++i) {
		_cid_cout("   [%0*lld][%0*lld]:", n_besz, i, mx_besz, dimacs_size_of_boolean_equation[i]);
		for(j = 0LL; j < dimacs_size_of_boolean_equation[i]; ++j) printf(" %lld", dimacs_boolean_equation[i][j]);
		printf("\n");
	}
	_cid_cout("#XEQ: [%llu]\n", dimacs_nb_of_xor_equations);
	for(i = 0LL; i < dimacs_nb_of_xor_equations; ++i) {
		_cid_cout("   [%0*lld][%0*lld]:", n_xesz, i, mx_xesz, dimacs_size_of_xor_equation[i]);
		for(j = 0LL; j < dimacs_size_of_xor_equation[i]; ++j) printf(" %lld", dimacs_xor_equation[i][j]);
		printf("\n");
	}
}

// counts number of satisfied clauses: (dm)
int_t dimacs_count_sat(void) {
    
	//return 0; //speed

    int_t count_sat = 0;
    int_t i,j;
    int_t n_besz = fast_int_log10(dimacs_nb_of_boolean_equations);
	int_t n_xesz = fast_int_log10(dimacs_nb_of_xor_equations);
	int_t mx_xesz = fast_int_log10(dimacs_size_max_among_xor_equations);
	int_t mx_besz = fast_int_log10(dimacs_size_max_among_boolean_equations);

    for(i = 0LL; i < dimacs_nb_of_boolean_equations; ++i) {
    	for(j = 0LL; j < dimacs_size_of_boolean_equation[i]; ++j) {
    		int_t lit = dimacs_boolean_equation[i][j];
    		int lit_sign = lit>0? 1:0;
    		if (cnf_assignment[llabs(lit)]==lit_sign) {
    			count_sat++;
    			break;
    		}
    	}
    }
    
    return count_sat;
}
