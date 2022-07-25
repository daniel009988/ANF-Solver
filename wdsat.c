//
//  wdsat.c
//

#include <string.h>
#include<stdio.h>
#include<time.h>
#include <stdlib.h>

#include "wdsat.h"
#include "cnf.h"
#include "xorset.h"
#include "xorgauss.h"
#include "dimacs.h"

#define _LA_LEVEL_ 17

int output_cnt = 0;

int_t l;
int_t m;
int test[100] = {0};

/// @var int_t wdsat_cnf_up_stack[__ID_SIZE__];
/// @brief unit propagation cnf stack
static int_t wdsat_cnf_up_stack[__ID_SIZE__];

/// @var int_t wdsat_cnf_up_top_stack;
/// @brief unit propagation cnf stack top
static int_t wdsat_cnf_up_top_stack;

/// @var int_t wdsat_xorset_up_stack[__ID_SIZE__];
/// @brief unit propagation xorset stack
static int_t wdsat_xorset_up_stack[__ID_SIZE__];

/// @var int_t wdsat_cnf_up_top_stack;
/// @brief unit propagation cnf stack top
static int_t wdsat_xorset_up_top_stack;

static int_t set[__ID_SIZE__];

bool wdsat_set_true(const int_t l) {
    bool _next_loop;
    int_t _l;
    wdsat_cnf_up_top_stack = 0LL;
    wdsat_cnf_up_stack[wdsat_cnf_up_top_stack++] = l;
    wdsat_xorset_up_top_stack = 0LL;
    wdsat_xorset_up_stack[wdsat_xorset_up_top_stack++] = l;
    _next_loop = true;
    while(_next_loop) {
		_next_loop = false;
		while(wdsat_cnf_up_top_stack) {
			_l = wdsat_cnf_up_stack[--wdsat_cnf_up_top_stack];
            if(_cnf_is_undef(_l)) _next_loop = true;
			if(!cnf_set_true(_l)) {/*printf("ter contr %lld\n",_l);*/return false;}
		}
		while(wdsat_xorset_up_top_stack) {
			_l = wdsat_xorset_up_stack[--wdsat_xorset_up_top_stack];
			if(_xorset_is_undef(_l)) _next_loop = true;
			if(!xorset_set_true(_l)) {/*printf("xor contr %lld\n",_l);*/return false;}
		}
		wdsat_cnf_up_top_stack = xorset_last_assigned(wdsat_cnf_up_stack);
		wdsat_xorset_up_top_stack = cnf_last_assigned(wdsat_xorset_up_stack);
	}
    return true;
}

bool wdsat_set_unitary(void) {
	bool _next_loop;
	int_t _l;
	wdsat_cnf_up_top_stack = 0LL;
	wdsat_xorset_up_top_stack = 0LL;
	
	if(!cnf_set_unitary()) return false;
	if(!xorset_set_unitary()) return false;
	wdsat_cnf_up_top_stack = xorset_last_assigned(wdsat_cnf_up_stack);
	wdsat_xorset_up_top_stack = cnf_last_assigned(wdsat_xorset_up_stack);
	_next_loop = true;
	while(_next_loop) {
		_next_loop = false;
		while(wdsat_cnf_up_top_stack) {
			_l = wdsat_cnf_up_stack[--wdsat_cnf_up_top_stack];
			if(_cnf_is_undef(_l)) _next_loop = true;
			if(!cnf_set_true(_l)) return false;
		}
		while(wdsat_xorset_up_top_stack) {
			_l = wdsat_xorset_up_stack[--wdsat_xorset_up_top_stack];
			if(_xorset_is_undef(_l)) _next_loop = true;
			if(!xorset_set_true(_l)) return false;
		}
		wdsat_cnf_up_top_stack = xorset_last_assigned(wdsat_cnf_up_stack);
		wdsat_xorset_up_top_stack = cnf_last_assigned(wdsat_xorset_up_stack);
	}
	return true;
}

int global = 0;
int global_assigned = 0;
int assigned_cnt = 0;
int global_sat = 0;
int sat;

// recursive procedure, main procedure. 
// l = variable position/level (initially=0) 
// nb_min_vars = maximum variable position
// conf = conflict counter
bool wdsat_solve_rest_XG(int_t l, int_t nb_min_vars, int_t conf[]) {

	//timeout?
	wct_end = clock();
	double wct_current = (double)(wct_end - wct_start) / CLOCKS_PER_SEC;
	
	/*if (wct_current>60) {
		printf("c TIMEOUT.\n");
		return false;
	}*/
	

	//we are finished when variable level equals max:
	if(l > nb_min_vars) {
		return true;	
	} 

	//is the variable set[level] already set currently? then we go one level up recursively (aka not undef)
	if(!_cnf_is_undef(set[l])) return wdsat_solve_rest_XG(l + 1, nb_min_vars, conf);

	//dm output assigned vars and count them:
	assigned_cnt = 0;
	if (usevarmap==0) {
		printf(TEXT_YELLOW);
		printf("b ");
		for(int j = 1; j <= nb_min_vars; j++)
		{
			if (j<32) {
				//keyvar?
				bool keyvar = false;
				#pragma unroll
				for (int xx=0; xx<512; xx++) if (j==keyvars[xx]) keyvar = true;
				if (cnf_assignment[j]==0 || cnf_assignment[j]==1) {
					if (keyvar) {printf(TEXT_YELLOW);} else {printf(TEXT_BLUE);}
					printf("%d",cnf_assignment[j]);
					printf(TEXT_DEFAULT);
				} else {
					if (keyvar) {printf(TEXT_YELLOW);} else {printf(TEXT_BLUE);}
					printf("_");
				}
			}
			if (cnf_assignment[j]<2) ++assigned_cnt;
		}
		printf("_\n");
		printf(TEXT_DEFAULT);
		//printf(" Propagated = %d\n",assigned_cnt);
		// update highest number of assigned vars:
		if (assigned_cnt>global_assigned) global_assigned = assigned_cnt;
	}
	
	//update global (highest level):
	if (l>global) {
		global = l;
		//printf("c NEW HIGHEST ASSIGNMENT LEVEL %d\n",global);
	}

	//count best #sat clauses:
	//disabled: speed! 
	sat = dimacs_count_sat(); if (sat>global_sat) global_sat = sat;

	// output header every x conflicts
	if (output_cnt % 30 ==0) {
		//if (conf[0]>1) printf("\n");
		printf(TEXT_BLUE);
		printf("c           | CURRENT   | LEVEL               | PROPAGATED                     | #SAT     | #SAT     |                |\n");
		printf("c TIME (s)  | VAR       | CURRENT  | MAX      | CURRENT  | MAX      | TOTAL    | CLAUSES  | BEST     | CONFLICTS      |\n");
		printf("c ----------|-----------|----------|----------|----------|----------|----------|----------|----------|----------------|\n");
		printf(TEXT_DEFAULT);
		output_cnt++;
	}
	if (usevarmap==0)
		printf("\rc %9.2f | %9lld | %8lld | %8d | %8d | %8d | %8lld | %8d | %8d | %14lld |",wct_current, set[l], l, global, assigned_cnt, global_assigned,nb_min_vars,sat,global_sat,conf[0]);
			if (l<=2) printf("\n----------------------------------------------------------------------------------------------------------------------|\n");
	if (usevarmap==1)
		printf("\rc %9.2f | %s     | %8lld | %8d | %8d | %8d | %8lld | %8d | %8d | %14lld |",wct_current, varmap[set[l]], l, global, assigned_cnt, global_assigned,nb_min_vars,sat,global_sat,conf[0]);

	//show variable names of assigned vars:
	if (usevarmap==1) {
		//printf(TEXT_SILVER);
		printf("\nc ");
		for (int j=1; j<=nb_min_vars; j++) {
			if (cnf_assignment[j]==0 || cnf_assignment[j]==1 || set[l]==j) {
				int setvar=false; for (int x=0; x<=l; x++) if (j==set[x]) setvar=true;
				if (setvar) printf(TEXT_GREEN);
				if (set[l]==j) printf(TEXT_YELLOW);
				printf("%s=%d ",varmap[j],cnf_assignment[j]);
				if (setvar) printf(TEXT_DEFAULT);
			}
		}
		printf("\n");
		//printf(TEXT_DEFAULT);	
	}

	fflush(stdout);
	
	_cnf_breakpoint;
	_xorset_breakpoint;
	_xorgauss_breakpoint;
	conf[0]++;

	// try to set var to 0: ------------------------------------------------------------------------------
	int try_var_first = -set[l]; //standard: first assignment is negative
	if (usehess==1) try_var_first = set[l]*hess[set[l]]; //using hess dump - first assignment is from hess
	if (firstone==1) try_var_first = -try_var_first; //first assignment is positive
	if(!wdsat_infer(try_var_first)) { 
		/// propagation created conflict:
		cnf_undo();
		xorset_undo();
		xorgauss_undo();
		#ifdef __DEBUG__
				printf("lev:%d--undo on 0\n",set[l]);
				for(int i = 1; i <= dimacs_nb_unary_vars(); i++)
					printf("%d", xorgauss_assignment[i]);
				printf("\n");
		#endif
	}
	else { 
		/// propagation did not create conflict:
		/// move to next l:
		if(!wdsat_solve_rest_XG(l + 1, nb_min_vars, conf))
		{
			cnf_undo();
			xorset_undo();
			xorgauss_undo();
			#ifdef __DEBUG__
						printf("lev:%d--undo on 0 profond\n",set[l]);
						for(int i = 1; i <= dimacs_nb_unary_vars(); i++)
							printf("%d", xorgauss_assignment[i]);
						printf("\n");
			#endif
		}
		// also next l didn't create conflict:
		else
		{
			#ifdef __DEBUG__
						printf("lev:%d--ok on 0\n",set[l]);
			#endif
			_cnf_mergepoint;
			_xorset_mergepoint;
			_xorgauss_mergepoint;
			return true;
		}
	}
	// try to set var to 1: -------------------------------------------------------------------------------
	if(!wdsat_infer(-try_var_first)) { 
		/// propagation created conflict:
		#ifdef __DEBUG__
				printf("lev:%d--undo on 1\n",set[l]);
		#endif
		return false;
	}
	#ifdef __DEBUG__
		for(int i = 1; i <= dimacs_nb_unary_vars(); i++)
			printf("%d", xorgauss_assignment[i]);
		printf("\n");
	#endif
	/// propagation did not create conflict:
	/// move to next l:
	return wdsat_solve_rest_XG(l + 1, nb_min_vars, conf);
}

// the infer algorithm - propagates CNF, then XOR with Gauss Jordan: 
bool wdsat_infer(const int_t l) {
	bool _loop_pass = true;
	bool _continue;
	int_t cnf_history_it;
	int_t cnf_history_last = cnf_history_top;
	int_t xorgauss_history_it;
	int_t xorgauss_history_last = xorgauss_history_top;
	int_t _l;
	
	if(!wdsat_set_true(l)) return false;
	while(_loop_pass) {
		// finalyse with XORGAUSS
		_continue = false;
		cnf_history_it = cnf_history_top;
		while(cnf_history_it > cnf_history_last) {
			_l = cnf_history[--cnf_history_it];
			if(_xorgauss_is_undef(_l)) {
				if(!xorgauss_set_true(_l)) return false; // <-- slow
				_continue = true;
			}
		}
		cnf_history_last = cnf_history_top;
		_loop_pass = false;
		if(_continue) {
			// get list of literal set thanks to XORGAUSS
			xorgauss_history_it = xorgauss_history_top;
			while(xorgauss_history_it > xorgauss_history_last) {
				_l = xorgauss_history[--xorgauss_history_it];
				if(_cnf_is_false(_l)) return false;
				if(_cnf_is_undef(_l)) {
					_loop_pass = true;
					if(!wdsat_set_true(_l)) return false;
				}
			}
			xorgauss_history_last = xorgauss_history_top;
		}
	}
	return true;
}

// this function is executed once before the real solving process:
bool wdsat_infer_unitary() {
	bool _loop_pass = true;
	bool _continue;
	int_t cnf_history_it;
	int_t cnf_history_last = cnf_history_top;
	int_t xorgauss_history_it;
	int_t xorgauss_history_last = xorgauss_history_top;
	int_t _l;
	
	if(!wdsat_set_unitary()) return false;
	while(_loop_pass) {
		// finalyse with XORGAUSS
		_continue = false;
		cnf_history_it = cnf_history_top;
		while(cnf_history_it > cnf_history_last) {
			_l = cnf_history[--cnf_history_it];
			if(_xorgauss_is_undef(_l)) {
				if(!xorgauss_set_true(_l)) return false;
				_continue = true;
			}
		}
		cnf_history_last = cnf_history_top;
		_loop_pass = false;
		if(_continue) {
			// get list of literal set thanks to XORGAUSS
			xorgauss_history_it = xorgauss_history_top;
			while(xorgauss_history_it > xorgauss_history_last) {
				_l = xorgauss_history[--xorgauss_history_it];
				if(_cnf_is_false(_l)) return false;
				if(_cnf_is_undef(_l)) {
					_loop_pass = true;
					if(!wdsat_set_true(_l)) return false;
				}
			}
			xorgauss_history_last = xorgauss_history_top;
		}
	}
	return true;
}

/// @fn solve();
/// @return false if formula is unsatisfiable and true otherwise
/// mode = 0: use cnf_upstack[] for cubes (fast, false positives)
/// mode = 1: use wdsat_set_true() for cubes (slow, correct)
bool wdsat_solve(int_t n, char mvc_graph[100000], char thread[1000], char cubes[100000], int mode) {
	char cubes_rem[100000] = ""; strcat(cubes_rem, cubes);
	if (mode==0) printf("c SOLVING MODE=0 (FAST, NEEDS VERIFICATION OF RESULT)\n");
	if (mode==1) printf("c SOLVING MODE=1 (SLOW, ACCURATE)\n");
	int_t j;
	int_t nb_min_vars;
	int_t conf[1]={0};
	cnf_initiate_from_dimacs();
	xorset_initiate_from_dimacs();
	if(!xorgauss_initiate_from_dimacs())
	{
		printf("c UNSAT on XORGAUSS init\n");
		return false;
	}
	cpy_from_dimacs();
	
	if(strlen(thread) > 0)
	{
		char *str_l;
		int_t l;
		str_l = strtok (thread, ",");
		while(str_l != NULL)
		{
			l = atoi(str_l);
			cnf_up_stack[cnf_up_top_stack++] = l;
			assert(cnf_up_top_stack < __ID_SIZE__);
			str_l = strtok (NULL, ",");
		}
		printf("c [%s] PREFIX generated\n",thread);
	}
	// end code for multithread (this has to be done before wdsat_infer_unitary();

	

	//cubes provided? --------------------------------------------------------------------------------
	if(strlen(cubes) > 0){
		char *str_l;
		int set_cnt = 0;
		str_l = strtok(cubes,",");
		printf("c CUBES: ");
		while(str_l != NULL) {
			printf("%d ",atoi(str_l));
			if (usevarmap==1) printf("(%s) ",varmap[abs(atoi(str_l))]);
			
			if (mode==1) wdsat_set_true(atoi(str_l)); 

			if (mode==0) cnf_up_stack[cnf_up_top_stack++] = atoi(str_l); assert(cnf_up_top_stack < __ID_SIZE__);
			
			str_l = strtok (NULL, ",");
			set_cnt++;
		}
		printf("\n");
		printf("c %d CUBES set.\n",set_cnt);
	}
	
	wdsat_infer_unitary();
	
	//MVC graph provided? -----------------------------------------------------------------------------
	if(strlen(mvc_graph) > 0) {
		nb_min_vars = 0;
		char *str_l;
		j = 0;
		int set_cnt = 0;
		str_l = strtok (mvc_graph, ",");
		while(str_l != NULL) {
			set[j] = atoi(str_l);
			str_l = strtok (NULL, ",");
			nb_min_vars++;
			j++;
			set_cnt++;
		}
		//dm: now add the remaining vars:
		for (int i=1;i<=dimacs_nb_unary_vars();i++) {
			bool found = false;
			for (int q=0; q<set_cnt; q++) if (set[q]==i) found = true;
			if (!found) {
				set[j] = i;
				nb_min_vars++;
				j++;
			}
		}
		//output set:
		printf("c USING MVC-SET: ");
		for (int i=0; i<nb_min_vars; i++) printf("%lld ",set[i]);
		printf("\n");
	}
	//no MVC graph, set is just the sequence 1,2,3,...
	else {
		nb_min_vars = dimacs_nb_unary_vars(); //number of vars in ANF
		for(j = 1; j <= nb_min_vars; j++) {
			set[j - 1] = j;
		}
	}

	//reverse set? ----------------------------------------------------------------------------------
	if (reverse==1) {
		for (int i=0; i<nb_min_vars/2; i++) {
			int tmp = set[i];
			set[i] = set[nb_min_vars-i-1];
			set[nb_min_vars-i-1] = tmp;
		}
		printf("c REVERSE-SET: ");
		for (int i=0; i<nb_min_vars; i++) printf("%lld ",set[i]);
		printf("\n");	
	}
	
	/// main solve: -----------------------------------------------------------------------------------------
	conflict_block = 0;
	wct_start = clock();
	/// set[level] -> is the array of the variable sequence span the search tree
	/// if no MVC graph, set[] = 1,2,3,4,5,6...dimacs_nb_unary_vars
	if(!wdsat_solve_rest_XG(0, nb_min_vars - 1, conf)) {
		printf("\nc UNSAT\n");
		printf("c conflicts: %lld\n",conf[0]);
		total_conflicts = total_conflicts + conf[0];
		printf("c total conflicts: %lld\n",total_conflicts);
		return false;
	}
	/// return=TRUE 
	wct_end = clock();
	double time_spent = (double)(wct_end - wct_start_total) / CLOCKS_PER_SEC;
	//check if the result is SAT (if mode=0):
	if (mode==0) {
		printf(TEXT_YELLOW);
		printf("\nc SEARCH SPACE COMPLETED - VERIFICATION...\n");
		printf(TEXT_DEFAULT);
		// running in mode=1 to verify the result:
		if (!wdsat_solve(n , mvc_graph, thread, cubes_rem, 1)) {
			printf("\nc UNSAT\n");
			printf("c conflicts: %lld\n",conf[0]);
			printf("c total conflicts: %lld\n",total_conflicts);
			return false;
		}	
	}
	/// SAT:
	if (mode==1) {
		printf("\nc SAT\n");
		printf("v ");
		for(j = 1; j <= dimacs_nb_unary_vars(); j++)
		{
			if (usevarmap==1) {
				if (cnf_assignment[j]==0) printf("%lld (%s) ",j*-1,varmap[j]);
				if (cnf_assignment[j]==1) printf("%lld (%s) ",j,varmap[j]);
			} else {
				if (cnf_assignment[j]==0) printf("%lld ",j*-1);
				if (cnf_assignment[j]==1) printf("%lld ",j);
			}
			
		}
		printf("\n");
		
		printf("c conflicts: %lld\n",conf[0]);
		printf("c SAT clauses: %lld\n",dimacs_count_sat());
		total_conflicts = total_conflicts + conf[0];
		printf("c TOTAL CONFLICTS: %lld\n",total_conflicts);
		printf("c TOTAL TIME SPENT: %fs\n",time_spent);	
	}
	
	return (true);

}
