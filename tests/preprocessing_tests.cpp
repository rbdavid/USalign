/******************************************************************************* 
 * code to implement the preprocessing functions as written natively in USalign
 ******************************************************************************/

#include <stdio.h>
#include <iostream>
#include "SOIalign.h"

int main(int argc, char* argv[])
{
    // only need to read in a single structure file and a closeK_opt value
    if (argc > 3) 
    {
	std::cerr << "too many args! should be:\n ./preprocessing_tests pdb_file closeK_opt\n where `pdb_file` is a path to a PDB structure and `closeK_opt` is a positive integer value.";
    }

    std::string structure_file = argv[1];
    int closeK_opt = atoi(argv[2]);

    // prep the variables for reading the structure file; most of these do not
    // matter for the purposes of the test, but needed to use the 
    // USalign-native pdb parser
    int chainnum = 0;
    std::vector<std::vector<std::string> >PDB_lines;
    std::vector<std::string> chainID_list;
    std::vector<std::string> chain2parse1;
    std::vector<std::string> model2parse1;
    std::vector<std::string> resi_vec;
    std::vector<int> mol_vec;
    int read_resi =  0;
    int ter_opt   =  2; 
    int infmt_opt = -1; 
    int split_opt =  2; 
    int het_opt   =  0; 
    std::string atom_opt = "auto"; 
    bool autojustify =  true; 
    int len;
    double **xa;
    double **xk;
    int **sec_bond;
    char *seq;
    char *sec;
    int i;

    // read the structure file
    chainnum = get_PDB_lines(structure_file, 
		    	     PDB_lines, 
			     chainID_list, 
			     mol_vec,
			     ter_opt, 
			     infmt_opt, 
			     atom_opt, 
			     autojustify, 
			     split_opt, 
			     het_opt,
			     chain2parse1, 
			     model2parse1);
    if (!chainnum)
    {
        std::cerr<<"Warning! Cannot parse file: "<<structure_file
            <<". Chain number 0."<<std::endl;
    }
    len = PDB_lines[0].size();
    mol_vec[0] = -1;
    if (!len)
    {
        std::cerr<<"Warning! Cannot parse file: "<<structure_file
            <<". Chain length 0."<<std::endl;
    }
    else if (len<3)
    {
        std::cerr<<"Sequence is too short <3!: "<<structure_file<<std::endl;
    }
    NewArray(&xa, len, 3);
    seq = new char[len + 1];
    sec = new char[len + 1];
    len = read_PDB(PDB_lines[0], xa, seq, resi_vec, read_resi);
    make_sec(xa, len, sec);
    if (closeK_opt >= 3) 
    {
	NewArray(&xk, len*closeK_opt, 3);
	getCloseK(xa, len, closeK_opt, xk); 
    }

    NewArray(&sec_bond, len, 2);
    assign_sec_bond(sec_bond, sec, len);

    // output coords, 
    printf("First atom's coordinates: %6.3f, %6.3f, %6.3f\n",xa[0][0],xa[0][1],xa[0][2]);
    printf("Last atom's coordinates:  %6.3f, %6.3f, %6.3f\n",xa[len-1][0],xa[len-1][1],xa[len-1][2]);
    
    // seq, 
    printf("Sequence: %s\n",seq);
    
    // sec, 
    printf("2ndary structure string: %s\n",sec);
    
    // sec_bond, 
    printf("secondary struct boundaries array:\n[");
    for (i = 0; i < len; i++)
    {
	printf("[%d,%d]",sec_bond[i][0], sec_bond[i][1]);
	if (i != len-1) printf(",");
    }
    printf("]\n");

    // and xk array (k nearest neighbors' coords)
    if (closeK_opt >= 3) 
    {
        printf("closeK_array (nearest neighbors' coordinates):\n[");
        for (i = 0; i < len*closeK_opt; i++)
	{
	    printf("[%6.3f,%6.3f,%6.3f]",xk[i][0], xk[i][1], xk[i][2]);
	    if (i != len*closeK_opt-1) printf(",");
	}
    }

    // clean up
    vector<string>().swap(chain2parse1);
    vector<string>().swap(model2parse1);
    PDB_lines.clear();
    chainID_list.clear();
    mol_vec.clear();
    resi_vec.clear();
    
    delete [] seq;
    delete [] sec;

    DeleteArray(&xa, len);
    DeleteArray(&sec_bond, len);
    if (closeK_opt >= 3) DeleteArray(&xk, len*closeK_opt);
    return 0;
}
