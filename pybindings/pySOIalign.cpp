//essential preamble for pybind11 casting
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

#include "SOIalign.h"	// 

namespace py = pybind11;


/****************************************************************************
 * FUNCTIONS NOT TO PORT TO PYTHON
 * these are still important to the the ports
 * they are called within the ported code
 ****************************************************************************/

// input arrays from the python-side of things are always flattened. to convert
// to USalign-like array object (array of arrays), need to map the c-type flat 
// array to the 2d.
template <class A> void fill_input_array(py::array_t<double> flat_array,
		      const int rows,
		      const int columns,
		      A *** array)
{
    // create the to-be-filled array object
    NewArray(array, rows, columns);
    
    // get the buffer regions for the input array object
    py::buffer_info flat_array_info = flat_array.request();
    // create array filled with the pointers for the array elements
    auto inp_ptr  = static_cast <double *>(flat_array_info.ptr);
   
    // loop over input array dimensions and assign values to the USalign-like
    // array object
    int ndims = rows * columns; // number of elements in the flat_array
    int i, j, k; // index from ndims, row index, column index
    for (i = 0; i < ndims; i++) 
    {
        j = i / columns;
	k = i % columns;
	// fill USalign-like array
	(*array)[j][k] = inp_ptr[i];
    }
}


// input arrays from the python-side of things are always flattened. to convert
// to USalign-like array object (array of arrays), need to map the c-type flat 
// array to the 2d.
template <class A> py::array_t<double> fill_output_array(A *** array,
		    		      			 const int rows,
		    		      			 const int columns)
{
    // define the result array object to be filled
    auto result = py::array_t<double>(rows*columns);
    // get the buffer regions for the array object
    py::buffer_info res_info = result.request();
    // create array filled with the pointers for the array elements
    auto out_ptr = static_cast <double *>(res_info.ptr);
    
    // loop over numpy array dimensions and assign values from the USalign-like
    // array object to the numpy array elements
    int ndims = rows * columns; // number of elements in the flat_array
    int i, j, k; // index from ndims, row index, column index

    for (i = 0; i < ndims; i++) 
    {
        j = i / columns;
	k = i % columns;
	// fill USalign-like array
	out_ptr[i] = (*array)[j][k];
    }
    return result;
}


/****************************************************************************
 * FUNCTIONS TO PORT TO PYTHON
 * these are still important to the the ports
 ****************************************************************************/


/*
 * a wrapper function around the original USalign make_sec() function
 * specifically, TMalign.h lines 765 to 792
 */
std::string wrap_make_sec(py::array_t<double> coords,
		   	  int len)
{
    // fill a USalign-like array, xa, from the input array, coords
    double **xa;
    fill_input_array(coords, len, 3, &xa);

    // prep the secondary structure string variable
    char *sec;
    sec = new char[len + 1];
    // run the TMalign.h make_sec() (lines 765-792)
    make_sec(xa, len, sec);
    // convert the char variable to a string, ready for the python env
    std::string str(sec);

    DeleteArray(&xa, len);
    return sec;
}


/*
 * a wrapper function around the original SOIalign assign_sec_bond() function
 * specifically, SOIalign.h lines 17 to 57
 */
py::array_t<int> wrap_assign_sec_bond(const std::string sec, const int len)
{
    // declare the USalign-like array of shape (len x 2)
    int **sec_bond;
    NewArray(&sec_bond, len, 2);
    
    // run SOIalign.h assign_sec_bond() (lines 17-57)
    assign_sec_bond(sec_bond, sec.c_str(), len);

    // sec_bond is now filled with (len x 2) ints associated with start,stop 
    // residues; 
    // convert it to a py::array_t<int>
    py::array_t<double> results = fill_output_array(&sec_bond, len, 2);

    // clean up
    DeleteArray(&sec_bond, len);

    return results;
}


/*
 * a wrapper function around the original SOIalign getCloseK() function
 * specifically, SOIalign.h lines 58 to 90
 */
py::array_t<double> wrap_getCloseK(py::array_t<double> coords,
				   const int len,
				   const int closeK_opt)
{
    // fill a USalign-like array, xa, from the input array, coords
    double **xa;
    fill_input_array(coords, len, 3, &xa);
   
    // create the nearest-neighbor array to be filled in getCloseK() 
    double **xk;
    NewArray(&xk, len*closeK_opt, 3);

    // fill the xk array with coordinates of all nearest-neighbors
    getCloseK(xa, len, closeK_opt, xk);

    // fill the results array with the nearest-neighbors data
    py::array_t<double> results = fill_output_array(&xk, len*closeK_opt, 3);
    
    DeleteArray(&xa, len);
    DeleteArray(&xk, len*closeK_opt);
    return results;
}


/*******************************************************************************
 * define the classes for user-input and parameter inputs for SOIalign_main()
 ******************************************************************************/

/*
 * define the input class to hold data associated with each structure, will be
 * passed as arguments into the wrapper function for the SOIalign_main function
 * important variables for structures: 
 *     coordinates, 
 *     sequence,
 *     2ndary structures string,
 *     nResidues length,
 *     k-nearest neighbors' coordinates array,
 *     2ndary structure features boundary array
*/

//class to handle the input parameters that are structure-specific
class alnStruct {
    public:
	// flattened 2d array nAtoms x 3
	py::array_t<double> coords; 
	
	// flattened 2d k-nearest-neighbors array; expected shape is 
	// (len * closeK_opt * 3); only created if closeK_opt >= 3
	py::array_t<double> k_nearest;
	
	// flattened 2d 2ndary struct boundaries array; expected shape is
	// (len x 2); only needed if mm_opt == 6 (sNS method)
	py::array_t<int> sec_bond; 
	std::string seq; // string, length len; holds the sequence
	std::string sec; // string, length len; holds the 2ndary struct labels
	int len;  // number of atoms to be aligned from structure
};


// class to handle the input parameters that are not structure-specific;
// sets default values for _all_ parameters that USalign passes :(
class alnParameters {
    public:
        // default set to -1 and then to 0 or 5 depending on mm_opt
	int closeK_opt;
        int molec_types; // -1*(xlen+ylen); a protein-protein alignment
	int mm_opt; // SOI alignment switch; 5 for fNS, 6 for sNS
};


/*******************************************************************************
 * define the class for SOIalign output data
 ******************************************************************************/

class outputResults {
    public:
        py::array_t<double> translation_vector; // 1d array vector for translation
        py::array_t<double> rotation_matrix;    // flattened 3x3 matrix for rotation 
        std::string seqM;	// mapping string between aligned sequences
        std::string seqA_mobile;// aligned seq; ordered by alnment to target sequence
        std::string seqA_target;// aligned seq
        double TM1, TM2; // normed by mobile and target lens, respectively
	double TM3, d0a; // only used if a_opt
       	double TM4, d0u; // only used if u_opt
	double TM5, d0_scale; // only used if d_opt
	double rmsd0; // final rmsd from Kabsch algo
        double d0_out = 5.0; // norm d0 value; also signify mapping in sequence alignment
	double Liden = 0;    // number of mapped residues w/ identical res type
        int n_ali8; // number of residues w/in d0_out
	int L_ali;   // n residues aligned, can include d > d0 pairs...
	double d0_mobile, d0_target; // d0 values calc'd when norming by 1 struct
	//double TM_ali, rmsd_ali; // looks like i_opt uses these as output vars
	//double TM_0, d0_0; // not used in the original output_results() func
};


/*******************************************************************************
 * TO REMIND MYSELF:
 * ! -> denotes an important func argument to be INPUT
 * | -> denotes an important func argument to be OUTPUT
 * all other arguments will be ignored for these python bindings
 *
 *USalign.cpp main() calls: 
 *    SOIalign(xname, yname, fname_super, fname_lign, fname_matrix, sequence, 
 *             Lnorm_ass, d0_scale, m_opt, i_opt, o_opt, a_opt, u_opt, d_opt, 
 *             TMcut, infmt1_opt, infmt2_opt, ter_opt, split_opt, outfmt_opt, 
 *             fast_opt, cp_opt, mirror_opt, het_opt, atom_opt, autojustify, 
 *             mol_opt, dir_opt, dirpair_opt, dir1_opt, dir2_opt, chain2parse1, 
 *             chain2parse2, model2parse1, model2parse2, chain1_list, 
 *             chain2_list, se_opt, closeK_opt ! , mm_opt ! );
 *
 *
 *
 *Then, in USalign.cpp SOIalign(), SOIalign_main() is called:
 *    SOIalign_main(xa ! , ya ! , xk ! , yk ! , closeK_opt ! , seqx ! , seqy ! ,
 *    		    secx ! , secy ! , t0 | , u0 | , TM1 | , TM2 | , TM3, TM4, 
 *    		    TM5, d0_0 | , TM_0 | , d0A | , d0B | , d0u, d0a, d0_out, 
 *    		    seqM | , seqxA | , seqyA | , invmap | , rmsd0 | , L_ali | , 
 *    		    Liden | , 
 *    		    TM_ali, rmsd_ali, n_ali, n_ali8, xlen, ylen, sequence, 
 *                  Lnorm_ass, d0_scale, i_opt, a_opt, u_opt, d_opt, 
 *                  force_fast_opt, mol_vec1[chain_i]+mol_vec2[chain_j], 
 *                  dist_list, secx_bond, secy_bond, mm_opt);
 *
 *
 *
 *Results are output in the USalign.cpp SOIalign() call via a 
 *output_results() call:
 *    output_results(xname.substr(dir1_opt.size()+dir_opt.size()+dirpair_opt.size()),
 *                   yname.substr(dir2_opt.size()+dir_opt.size()+dirpair_opt.size()),
 *                   chainID_list1[chain_i], chainID_list2[chain_j], xlen, ylen,
 *                   t0, u0, TM1, TM2, TM3, TM4, TM5, rmsd0, d0_out, 
 *                   seqM.c_str(), seqxA.c_str(), seqyA.c_str(), Liden, n_ali8,
 *                   L_ali, TM_ali, rmsd_ali, TM_0, d0_0, d0A, d0B, Lnorm_ass, 
 *                   d0_scale, d0a, d0u, (m_opt?fname_matrix:"").c_str(),
 *                   outfmt_opt, ter_opt, false, split_opt, o_opt, fname_super,
 *                   i_opt, a_opt, u_opt, d_opt, mirror_opt, resi_vec1, resi_vec2); 
 *
 * only need t0, u0, TM1, TM2, TM3, TM4, TM5, rmsd0, d0_out, 
 *           seqM, seqxA, seqyA, 
 *           Liden, n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0, d0A, 
 *           d0B, d0a, d0u 
 *
 *
 ******************************************************************************/

/*
 *
 * *** =  important to return from python bindings; check on these in the 
 *        SOIalign() function.
 *
 * Other parameters used in SOIalign_main() function:
 * t0, u0, ***
 * TM1, TM2, *** 
 * TM3, TM4, TM5,
 * d0_0, ***
 * TM_0,
 * d0A, d0B, d0u, d0a, d0_out,
 * seqM, seqxA, seqyA, invmap, ***
 * rmsd0, L_ali, Liden, ***
 * TM_ali, rmsd_ali, n_ali, n_ali8,
 * sequence, Lnorm_ass, 
 * d0_scale,
 * i_opt, a_opt, u_opt, d_opt, fast_opt, 
 * mol_type, 
 * dist_list ***
 *
 * THESE PARAMETERS ARE PASSED TO SOIalign_main FOR SCOPE PURPOSES AND/OR FOR 
 * MEMORY CONTROL. 
 *
 * some of these parameters are hard-coded in the parent SOIalign() function in 
 * USalign.cpp. So I think these should be hardcoded in an input parameter 
 * class or in the wrapper function.
 */

/*******************************************************************************
 * define the wrapper function that reads input, runs SOIalign_main, and 
 * returns the output data.
 ******************************************************************************/

outputResults runSOIalign( alnStruct& mobile_data, 
			   alnStruct& target_data, 
			   alnParameters& parameters)
{
    // create variables to be used/filled w/in the SOIalign_main() call
    double t0[3], u0[3][3];
    double TM1, TM2;
    double TM3, TM4, TM5; 
    double d0_0, TM_0;
    double d0A, d0B, d0u, d0a;
    double d0_out=5.0;
    std::string seqM, seqxA, seqyA;	// for output alignment
    double rmsd0;
    int L_ali;                		// Aligned length in standard_TMscore
    double Liden;
    double TM_ali, rmsd_ali;  		// TMscore and rmsd in standard_TMscore
    int n_ali;
    int n_ali8;
    

    // hardcode some default parameters to avoid their implementation...
    // !!! design choice to still input these variables into the 
    //     SOIalign_main() call.
    int i_opt = 0;
    std::vector<std::string> sequence;
    int a_opt = 0; 		// > 0 if wish normalize TMscore by avg len of structs
    bool u_opt = false; 	// true if wish normalize TMscore by Lnorm_ass value
    double Lnorm_ass;
    bool d_opt = false; 	// true if wish to use d0_scale for calc of pair weights
    double d0_scale;
    bool fast_opt = false; 	// true if wish for fast but inaccurate alnment

    // convert input arrays to their USalign-like versions...
    int i,j,k;

    // mobile coords
    double **mobile_coords;
    fill_input_array(mobile_data.coords, mobile_data.len, 3, &mobile_coords);
    
    // mobile k_nearest array
    int **mobile_sec_bond;
    fill_input_array(mobile_data.sec_bond, mobile_data.len, 2, &mobile_sec_bond);
    
    // mobile k_nearest array
    double **mobile_k_nearest;
    fill_input_array(mobile_data.k_nearest, mobile_data.len * parameters.closeK_opt, 3, &mobile_k_nearest);
  
    // convert std::string to c-like char pointers
    char *mobile_seq = mobile_data.seq.data();
    char *mobile_sec = mobile_data.sec.data();

    // target coords
    double **target_coords;
    fill_input_array(target_data.coords, target_data.len, 3, &target_coords);
    
    // target k_nearest array
    int **target_sec_bond;
    fill_input_array(target_data.sec_bond, target_data.len, 2, &target_sec_bond);
    
    // target k_nearest array
    double **target_k_nearest;
    fill_input_array(target_data.k_nearest, target_data.len * parameters.closeK_opt, 3, &target_k_nearest);
    
    // convert std::string to c-like char pointers
    char *target_seq = target_data.seq.data();
    char *target_sec = target_data.sec.data();


    // prep other parameters

    // if min of len values are too large (>1500), then do fast alignment 
    // method...
    bool force_fast_opt=(std::min(mobile_data.len,target_data.len)>1500)?true:fast_opt;
    
    // these lines aren't gonna work...
    int *invmap = new int[target_data.len+1];
    double *dist_list = new double[target_data.len+1];
   
    // run the SOIalign_main function as defined in the SOIalign.cpp file,
    // for all its faults
    SOIalign_main(mobile_coords, target_coords, 
		  mobile_k_nearest, target_k_nearest,
		  parameters.closeK_opt,
		  mobile_seq, target_seq, mobile_sec, target_sec,
		  t0, u0,
		  TM1, TM2, TM3, TM4, TM5,
		  d0_0, TM_0,
		  d0A, d0B,
		  d0u, d0a, d0_out,
		  seqM, seqxA, seqyA,
		  invmap,
		  rmsd0,
		  L_ali, Liden,
		  TM_ali, rmsd_ali, n_ali,
		  n_ali8,
		  mobile_data.len, target_data.len, 
		  sequence,
		  Lnorm_ass, d0_scale,
		  i_opt, 
		  a_opt, u_opt, d_opt, 
		  force_fast_opt,
		  parameters.molec_types, dist_list,
		  mobile_sec_bond, target_sec_bond, parameters.mm_opt);

    // gather only the necessary output and return them in the output class as
    // python accessible data types

    outputResults out; // instantiate the output object

    // define the translation array object to be filled
    auto trans_array = py::array_t<double>(3);
    // get the buffer regions for the array object
    py::buffer_info trans_info = trans_array.request();
    // create array filled with the pointers for the array elements
    auto trans_ptr = static_cast <double *>(trans_info.ptr);
    // fill those elements 
    for (i = 0; i<3; i++)
    {
	trans_ptr[i] = t0[i];
    }
    out.translation_vector = trans_array;

    // gather the rotation array into a flat numpy array
    auto rot_array = py::array_t<double>(9);
    // get the buffer regions for the array object
    py::buffer_info rot_info = rot_array.request();
    // create array filled with the pointers for the array elements
    auto rot_ptr = static_cast <double *>(rot_info.ptr);
    // fill those elements 
    for (i = 0; i<3; i++)
    {
	for (int j = 0; j<3; j++)
	{
	    k = i*3 + j; 
            rot_ptr[k] = u0[i][j];
	}
    }
    out.rotation_matrix = rot_array;

    // seq results are char* so need to convert it to a std::string
    out.seqM = pybind11::str(seqM);
    out.seqA_mobile = pybind11::str(seqxA);
    out.seqA_target = pybind11::str(seqyA);

    // handling tm score values
    out.TM1 = TM1;
    out.TM2 = TM2;
    if (a_opt > 0) 
    {
	out.TM3 = TM3;
	out.d0a = d0a;
    }
    else 
    {
        out.TM3 = -1;
        out.d0a = -1;
    }

    if (u_opt > 0) 
    {
	out.TM4 = TM4;
        out.d0u = d0u;
    }
    else 
    {
	out.TM4 = -1;
        out.d0u = -1;
    }

    if (d_opt > 0) 
    {
	out.TM5 = TM5;
        out.d0_scale = d0_scale;
    }
    else 
    {
	out.TM5 = -1;
        out.d0_scale = -1;
    }

    // handling other  value
    out.rmsd0 = rmsd0;
    out.d0_out = d0_out;
    out.Liden = Liden;
    out.n_ali8 = n_ali8;
    out.L_ali = L_ali;
    out.d0_mobile = d0A;
    out.d0_target = d0B;

    // cleaning up
    DeleteArray(&mobile_coords, mobile_data.len);
    DeleteArray(&target_coords, target_data.len);
    DeleteArray(&mobile_sec_bond, mobile_data.len);
    DeleteArray(&target_sec_bond, target_data.len);
    DeleteArray(&mobile_k_nearest, mobile_data.len * parameters.closeK_opt);
    DeleteArray(&target_k_nearest, target_data.len * parameters.closeK_opt);

    return out;
}

/*******************************************************************************
 * defining the pybind11 wrapper for SOIalign_main
 ******************************************************************************/

PYBIND11_MODULE(pySOIalign, m) {
    m.doc() = "pybind11 port of SOIalign_main and related functions from USalign codebase"; 
    m.def("wrap_make_sec",
	  &wrap_make_sec,
	  "function to assign 2ndary structure character to each residue in the structure",
	  py::arg("coords"), 
	  py::arg("len"));

    m.def("wrap_getCloseK",
	  &wrap_getCloseK,
	  "function to determine nearest neighbors for each residue",
	  py::arg("coords"), 
	  py::arg("len"),
	  py::arg("closeK_opt"));
    
    m.def("wrap_assign_sec_bond",
	  &wrap_assign_sec_bond,
	  "function to identify the boundaries of large helix/sheet 2ndary structure elements",
	  py::arg("sec"), 
	  py::arg("len"));
    
    py::class_<alnStruct>(m, "alnStruct")
	.def(py::init<py::array_t<double>, 
		      py::array_t<double>, 
		      py::array_t<int>, 
		      std::string, 
		      std::string, 
		      int>())
	.def_readwrite("coords", &alnStruct::coords)
	.def_readwrite("k_nearest", &alnStruct::k_nearest)
	.def_readwrite("sec_bond", &alnStruct::sec_bond)
	.def_readwrite("seq", &alnStruct::seq)
	.def_readwrite("sec", &alnStruct::sec)
	.def_readwrite("len", &alnStruct::len);

    py::class_<alnParameters>(m, "alnParameters")
	.def(py::init<int, 
		      int, 
		      int>())
	.def_readwrite("closeK_opt", &alnParameters::closeK_opt)
	.def_readwrite("molec_types", &alnParameters::molec_types)
	.def_readwrite("mm_opt", &alnParameters::mm_opt);

    py::class_<outputResults>(m, "outputResults")
	.def(py::init<py::array_t<double>,
		      py::array_t<double>,
		      std::string,
		      std::string,
		      std::string,
		      double,
		      double,
		      double,
		      double,
		      double,
		      double,
		      double,
		      double,
		      double,
		      double,
		      double,
		      int,
		      int,
		      double,
		      double>())
	.def_readwrite("translation_vector", &outputResults::translation_vector)
	.def_readwrite("rotation_matrix", &outputResults::rotation_matrix)
	.def_readwrite("seqM", &outputResults::seqM)
	.def_readwrite("seqA_mobile", &outputResults::seqA_mobile)
	.def_readwrite("seqA_target", &outputResults::seqA_target)
	.def_readwrite("TM1", &outputResults::TM1)
	.def_readwrite("TM2", &outputResults::TM2)
	.def_readwrite("TM3", &outputResults::TM3)
	.def_readwrite("d0a", &outputResults::d0a)
	.def_readwrite("TM4", &outputResults::TM4)
	.def_readwrite("d0u", &outputResults::d0u)
	.def_readwrite("TM5", &outputResults::TM5)
	.def_readwrite("d0_scale", &outputResults::d0_scale)
	.def_readwrite("rmsd0", &outputResults::rmsd0)
	.def_readwrite("d0_out", &outputResults::d0_out)
	.def_readwrite("Liden", &outputResults::Liden)
	.def_readwrite("n_ali8", &outputResults::n_ali8)
	.def_readwrite("L_ali", &outputResults::L_ali)
	.def_readwrite("d0_mobile", &outputResults::d0_mobile)
	.def_readwrite("d0_target", &outputResults::d0_target);

    m.def("runSOIalign",
	  &runSOIalign,
	  "wrapper function for the SOIalign_main cpp function",
	  py::arg("mobile_data"),
	  py::arg("target_data"),
	  py::arg("parameters"),
	  py::return_value_policy::copy);
}
