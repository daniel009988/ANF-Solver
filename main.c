// compile: make
// run: ./q_solver -i <instance.anf>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>

#include "wdsat_utils.h"
#include "dimacs.h"
#include "cnf.h"
#include "xorset.h"
#include "xorgauss.h"
#include "wdsat.h"

extern byte_t thread_digits;

char input_filename[__STATIC_STRING_SIZE__];
char input_varmap[__STATIC_STRING_SIZE__];
char cnf_input_varmap[__STATIC_STRING_SIZE__];
char cubes_filename[__STATIC_STRING_SIZE__];
char output_filename[__STATIC_STRING_SIZE__];
char hess_filename[__STATIC_STRING_SIZE__];
char irr[160] = "11001";
char X3[160] = "0001";
int action = 2;
int n = 4;
int _l = 0;
int _m = 3;
int occ = 0;
int s = 0;
int xg = 1;
int br_sym = 0;
char mvc_graph[100000] = "";
char thread[1000] = "";
char cubes[100000] = "";
int reverse = 0;
int firstone = 0;
int usevarmap = 0;
int usecubesfile = 0;
usehess = 0;
/**
 * @fn int scan_opt(int argc, char **argv, const char *opt)
 * @brief It scans and checks command line options
 */
int scan_opt(int argc, char **argv, const char *opt) {
	char c;
	while ((c = getopt (argc, argv, opt)) != -1)
		switch (c) {
			case 'a': if(!strncmp(optarg, "compile", 7)) action = 0; if(!strncmp(optarg, "time", 4)) action = 2; if(!strncmp(optarg, "delete", 6)) action = 3; if(!strncmp(optarg, "color", 5)) action = 4; break;
			case 'i': strcpy(input_filename, optarg); break;
			case 'v': strcpy(input_varmap, optarg); break;
			case 'o': strcpy(output_filename, optarg); break;
			case 'n': n = atoi(optarg); break;
			case 'l': _l = atoi(optarg); break;
			case 'm': _m = atoi(optarg); break;
			case 'x': xg = 1; break;
			case 'b': br_sym = 1; break;
			//case 'c': occ = 1; break;
			case 's': s = 1; break;
			case 'r': reverse = 1; break;
			case 'f': firstone = 1; break;
			case 'g': strcpy(mvc_graph, optarg); break;
			case 't': strcpy(thread, optarg); break;
			case 'c': strcpy(cubes, optarg); break;
			case 'w': strcpy(cnf_input_varmap, optarg); break;
			case 'q': strcpy(cubes_filename, optarg); break;
			case 'h': strcpy(hess_filename, optarg); break;
			default: return(-1);
		}
	return(0);
}

boolean_t generate(char *ifn, char *ofn) {
	FILE *d;
	
	if((d = fopen(ifn, "r")) == NULL) {
		_cid_cout("%s :: unkown filename", ifn);
		return(__OFF__);
	}
	/// read the header if exists
	dimacs_read_header(d); //get # vars, etc
	fclose(d);
	
	/// read the formula
	d = fopen(ifn, "r");
	dimacs_read_formula(d);
	fclose(d);
	
	// output settings: ---------------------------------------------------------------------------------------------------
	printf("c\n");
	printf("c CONFIG.h COMPILATION ADVICE:\n");
	printf("c __MAX_BUFFER_SIZE__    = %lld\n",(dimacs_nb_atoms()>dimacs_nb_xor_atoms())?dimacs_nb_atoms()+1:dimacs_nb_xor_atoms()+1);
	printf("c __MAX_EQ__             = %lld\n",dimacs_nb_equations()+1);
	printf("c __MAX_XEQ__            = %lld\n",dimacs_nb_xor_equations());
	printf("c __MAX_ANF_ID__         = %lld\n",dimacs_nb_unary_vars()+1);
	printf("c __MAX_ID__             = %lld\n",dimacs_nb_vars());
	printf("c\n");
	printf("c using file: %s\n",ifn);
	printf("c ANF vars (unary): %lld\n",dimacs_nb_unary_vars());
	printf("c ANF clauses (unary): %lld\n",_dimacs_nb_of_eq());
	printf("c CNF vars: %lld\n",dimacs_nb_vars());
	printf("c CNF clauses: %lld\n",dimacs_nb_equations());
	printf("c atoms [boolean equations]: %lld\n", dimacs_nb_atoms());
	printf("c atoms [xor equations]: %lld\n", dimacs_nb_xor_atoms());
	printf("c size of xor gauss: %lld\n", __SZ_GAUSS__);
	printf("c thread set to: %s\n",thread);
	printf("c reverse: %d\n",reverse);
	printf("c first bit always: %d\n",firstone);

	//loading hess dump file: -----------------------------------------------------------------------------------------------
	printf("c Loading hess-dump file: %s\n",hess_filename);
	FILE *fp3;
    fp3 = fopen(hess_filename, "r");
    if (fp3 == NULL) {
        fprintf(stdout, "c Error: Cannot open hess-dump file.\n");
    } else {
        // load assignment file (comes from cadical):
        int assignment = 0;
        int ivar;
        int cnt = 0;
        while (fscanf(fp3, "%i,", &ivar) > 0) {
            if (ivar<0) {assignment=-1;} else {assignment=1;} //adjust to 1/0
            hess[abs(ivar)] = assignment;
            cnt++;
        }
        fclose(fp3);
        fprintf(stdout, "c %i hess assignments loaded.\n",cnt);
        usehess = 1;
        //for (int j=1; j<=cnt; j++) printf("%d->%d ",j,hess[j]);
    }

	//loading varmap: -------------------------------------------------------------------------------------------------------
	//varmap[anf-varnum] -> returns anf-varname
	printf("c Loading varmap: %s\n",input_varmap);
	FILE *fp;
	int varmap_cnt = 0;
	fp = fopen(input_varmap,"r");
	if (fp == NULL) {
        fprintf(stdout, "c Error: Cannot open varmap.\n");
    } else {
    	int varnum = 0;
		while (fscanf(fp,"%i,%s",&varnum,varmap[varnum+1]) > 0) varmap_cnt++;
    	fclose(fp);
    	printf("c %d varmap vars loaded.\n",varmap_cnt);
    	usevarmap=1;
    }
    //output: 
    //for (int i=1;i<varmap_cnt+1;i++) printf("%s = %d ",varmap[i],i);

    //loading cnf varmap: ----------------------------------------------------------------------------------------------------
    //cnfvarmap[cnf-varnum] -> returns anf-varnum
    printf("c Loading cnf-varmap: %s\n",cnf_input_varmap);
    int cnf_varmap_cnt = 0;
    FILE *fp2;
	fp2 = fopen(cnf_input_varmap,"r");
	if (fp2 == NULL) {
        fprintf(stdout, "c Error: Cannot open cnf-varmap.\n");
    } else {
    	int cnf_varnum = 0;
    	int anf_varnum = 0;
    	while (fscanf(fp2,"%i,%i",&cnf_varnum,&anf_varnum) > 0) {
			cnfvarmap[anf_varnum] = cnf_varnum;
			cnf_varmap_cnt++;
		}
    	fclose(fp2);
    	printf("c %d cnf varmap vars loaded.\n",cnf_varmap_cnt);
    }
    //output:
    //for (int i=1;i<cnf_varmap_cnt+1;i++) printf("%d = %d ",cnfvarmap[i],i);

    //using cubes file (gerated with march_cu): ----------------------------------------------------------------------------
    printf("c Using cubes file: %s\n",cubes_filename);
    FILE *cubesfile;
    char buffer[32];
    int cubevar;
    int cubes_cnt = 0;
    total_conflicts = 0;
    cubesfile = fopen(cubes_filename,"r");
    if (cubesfile == NULL) {
    	printf("c Error: Cannot open cubes file.\n");
    } else {
    	bool solvecube = false;
    	wct_start_total = clock();
    	while (fscanf(cubesfile,"%s ",buffer) > 0) {
    		//new cube starting:
    		if (strcmp(buffer,"a")==0) {
    			//printf("a "); 
    			strcpy(cubes, ""); //empty cubes string
    			solvecube = true;
    		}
    		//end of cube:
    		if (strcmp(buffer,"0")==0) {
    			//printf("\n"); 
    			cubes[strlen(cubes)-1] = '\0'; //remove last "," from cubes-string
    			///START SOLVER +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    			//IF CUBES ARE TOO SHORT, STILL DOING? -> NO. WE ARE LOOSING POTENTIAL CORRECT PATHS
    			//MAYBE ADD TIMEOUT?
    			cubes_cnt++;
    			if (solvecube) {
    				printf(TEXT_CYAN);
    				printf("c SOLVING CUBE %d: %s\n",cubes_cnt, cubes);
	    			printf(TEXT_DEFAULT);
					if (wdsat_solve(n , mvc_graph, thread, cubes, 0)) return(__ON__);
	    			//wdsat_solve(n , mvc_graph, thread, cubes, 0); //test: run all cubes to see if we find the correct solution
					wct_end = clock();
					double time_spent = (double)(wct_end - wct_start) / CLOCKS_PER_SEC;
					printf("c TIME SPENT: %fs\n",time_spent);	
    			}
    			///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    			
    		}
    		//cube vars:
    		if (strcmp(buffer,"a")!=0 && strcmp(buffer,"0")!=0) {
    			//add to cube-string the ANF-varnum (converted from CNF-varnum):
    			int _cube_var = atoi(buffer);
    			int cube_var = 0;
    			for (int i=1; i<cnf_varmap_cnt;i++ ) if (cnfvarmap[i]==abs(_cube_var)) cube_var=i;
    			if (_cube_var<0) cube_var = cube_var*-1;
    			//printf("%d(CNF)>%d(ANF) ",_cube_var, cube_var);
    			if (cube_var!=0) { //SHOULD WE SKIP THIS CUBE THEN? OTHERWISE WE RUN CUBES TWICE
    				char cube_var_str[12];
	    			snprintf(cube_var_str,12,"%d",cube_var);
					strcat(cubes, cube_var_str); 
	    			strcat(cubes,","); //add "," to cube-string	
    			} else {
    				printf(TEXT_CYAN);
    				printf("c CUBE %d WARNING - VARIABLE %d OUTSIDE SCOPE\n",cubes_cnt+1,_cube_var);
    				printf(TEXT_DEFAULT);
    				//solvecube = false; //up to d=10 cubes we have FULL real cubes; wont find sol for higher
    			}
    		}
    	}
    	fclose(cubesfile);
    	wct_end = clock();
		double time_spent = (double)(wct_end - wct_start_total) / CLOCKS_PER_SEC;
		printf("c TOTAL TIME SPENT: %fs\n",time_spent);
    	return(__ON__);
    }
	
	if(_l == 0)_l = n;

	/// start the solver (without cube file): --------------------------------------------------------------------
	printf("c STARTING SOLVER...\n");
	wdsat_solve(n , mvc_graph, thread, cubes, 1);
	wct_end = clock();
	double time_spent = (double)(wct_end - wct_start) / CLOCKS_PER_SEC;
	printf("c TIME SPENT: %fs\n",time_spent);
	
	return(__ON__);
}

int main(int argc, char * argv[]) {
	byte_t exit_value = (byte_t) EXIT_SUCCESS;

	char *syntax =
	"c          -x : to enable Gaussian Elimination (default: enabled)\n"
	"c          -i file    : input file\n"
	"c          -v file    : input varmap (ANF_VARNAME->ANF_VARNUM mapping)\n"
	"c          -w file    : input cnfvarmap (CNF_VARNAME->ANF_VARNUM mapping)\n"
	"c          -q file    : input cube file (march_cu compatible)\n"
	"c          -h file    : input hess dump file (hessanf compatible)\n"
	"c          -r : reverse variable order to by tried\n"
	"c          -g mvc    : where mvc is a string of comma-separated variables that defines statically the branching order\n"
	"c          -c cubes  : where cubes is a string of comma-sperated variables that set the cube\n"
	"c          -h : help (shows the argument list)\n"
	;
	
	goto on_continue;
on_break:
	printf("c Syntax: %s <... Args ...>\n", argv[0]);
	printf("c Args:\n");
	printf("%s", syntax);
	printf("\n");
	exit_value = 1;
	goto end;
	
on_continue:
	if(scan_opt(argc, argv, "i:o:w:q:f:r:c:h:a:v:n:l:xbshm:g:t:")) goto on_break;
	
	printf("c ----------------------------------------------------------------\n");
	printf("c ANF SOLVER\n");
	printf("c ----------------------------------------------------------------\n");

	/// allocate memory: -todo: we should allocate all memory here -----------------------------------------------
	//xorgauss.h:
	xorgauss_equivalency_history = malloc(sizeof(uint_t[__MAX_ANF_ID__][__ID_SIZE__][__SZ_GAUSS__]));
	xorgauss_equivalent_history = malloc(sizeof(bool[__MAX_ANF_ID__][__ID_SIZE__]));
	xorgauss_assignment_buffer_history = malloc(sizeof(boolean_t[__MAX_ANF_ID__][__SIGNED_ID_SIZE__]));
	xorgauss_current_degree_history = malloc(sizeof(boolean_t[__MAX_ANF_ID__][__ID_SIZE__]));
	//dimacs.h:
	monomials_to_column = malloc(sizeof(uint_t[__MAX_ANF_ID__][__ID_SIZE__][__MAX_DEGREE__ - 1]));

	/// -----------------------------------------------------------------------------------------------------------
	
    /// start:
	if(!generate(input_filename, output_filename)) {
		exit_value = (byte_t) EXIT_FAILURE;
	}
	
end:
	//printf("\n");
	return((int) exit_value);
}
