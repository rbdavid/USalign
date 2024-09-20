/******************************************************************************* 
 * code to implement the preprocessing functions as written natively in USalign
 ******************************************************************************/

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include "SOIalign.h"

int main(int argc, char* argv[])
{
    // only need to read in a single structure file and a closeK_opt value
    if (argc != 4) 
    {
	std::cerr << "\nIncorrect number of input args!\n\n Usage: ./preprocessing_tests pdb_file closeK_opt out_stem\n where `pdb_file` is a path to a PDB structure, `closeK_opt` is a positive integer value, and `out_stem` is a string to be used as the stem to output file names.\n";
    }

    std::string structure_file = argv[1];
    int closeK_opt = atoi(argv[2]);
    std::string out_stem = argv[3];

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

    /*******************************
      CALCULATE THE INPUT FEATURES
    *******************************/

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
    
    // do the source's check of the pdb file...
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
    // call to fill xa array and seq char
    NewArray(&xa, len, 3);
    seq = new char[len + 1];
    len = read_PDB(PDB_lines[0], xa, seq, resi_vec, read_resi);
    // call to fill sec char
    sec = new char[len + 1];
    make_sec(xa, len, sec);
    // check if the closeK_opt parameter is an expected value for the 
    // neighborlist to be filled
    if (closeK_opt >= 3) 
    {
	NewArray(&xk, len*closeK_opt, 3);
	getCloseK(xa, len, closeK_opt, xk); 
    }
    // otherwise, just fill an empty array of zeros
    else 
    {
	NewArray(&xk, 1, 3);
	xk[0][0] = xk[0][1] = xk[0][2] = 0.0;
    }

    // call to fill the secondary struct boundaries list of lists
    NewArray(&sec_bond, len, 2);
    assign_sec_bond(sec_bond, sec, len);

    /***********************
      SAVE RESULTS TO FILES
    ***********************/
    
    // output coords, 
    // printf("First atom's coordinates: %6.3f, %6.3f, %6.3f\n",xa[0][0],xa[0][1],xa[0][2]);
    // printf("Last atom's coordinates:  %6.3f, %6.3f, %6.3f\n",xa[len-1][0],xa[len-1][1],xa[len-1][2]);
    ofstream coordFile;
    coordFile.open (out_stem + "_coords.dat");
    coordFile << "[";
    for (i = 0; i < len; i++)
    {
	coordFile << std::fixed << std::setprecision(4) <<"[" << xa[i][0] << "," << xa[i][1] << "," << xa[i][2] << "]";
	if (i != len-1) coordFile << ",";
    }
    coordFile << "]";
    coordFile.close();
    
    // output seq and sec strings
    ofstream stringsFile;
    stringsFile.open (out_stem + "_strings.dat");
    stringsFile << seq << "\n";
    stringsFile << sec << "\n";
    stringsFile.close();

    //printf("Sequence: %s\n",seq);
    //
    //printf("2ndary structure string: %s\n",sec);
    
    // output sec_bond list of lists in a single line, 
    ofstream bondFile;
    bondFile.open (out_stem + "_2ndaryboundaries.dat");
    bondFile << "[";
    for (i = 0; i < len; i++)
    {
	bondFile << "[" << sec_bond[i][0] << "," << sec_bond[i][1] << "]";
	if (i != len-1) bondFile << ",";
    }
    bondFile << "]";
    bondFile.close();

    //printf("secondary struct boundaries array:\n[");
    //for (i = 0; i < len; i++)
    //{
    //    printf("[%d,%d]",sec_bond[i][0], sec_bond[i][1]);
    //    if (i != len-1) printf(",");
    //}
    //printf("]\n");

    // and xk array (k nearest neighbors' coords)
    ofstream neighborFile;
    neighborFile.open (out_stem + "_neighborlist.dat");
    neighborFile << "[";
    if (closeK_opt >= 3) 
    {
        for (i = 0; i < len*closeK_opt; i++)
        {
            neighborFile << std::fixed << std::setprecision(4) << "[" << xk[i][0] << "," << xk[i][1] << "," << xk[i][2] << "]";
            if (i != len*closeK_opt-1) neighborFile << ",";
        }
    }
    else
    {
        neighborFile << std::fixed << std::setprecision(4) << "[" << xk[0][0] << "," << xk[0][1] << "," << xk[0][2] << "]";
    }
    neighborFile << "]";
    neighborFile.close();

    //if (closeK_opt >= 3) 
    //{
    //    printf("closeK_array (nearest neighbors' coordinates):\n[");
    //    for (i = 0; i < len*closeK_opt; i++)
    //    {
    //        printf("[%6.3f,%6.3f,%6.3f]",xk[i][0], xk[i][1], xk[i][2]);
    //        if (i != len*closeK_opt-1) printf(",");
    //    }
    //    printf("]\n");
    //}

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
    else DeleteArray(&xk, 1);
    return 0;
}
