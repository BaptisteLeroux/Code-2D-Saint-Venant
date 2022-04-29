#ifndef _FINITEVOLUME_CPP
#include "FiniteVolume.h"
#include <fstream>
#include <iostream>

using namespace std;
using namespace Eigen;

// Constructeur
FiniteVolume::FiniteVolume(Function* function, DataFile* data_file, Mesh2D* mesh) :
_fct(function), _df(data_file), _msh(mesh)
{
	std::cout << "Build finite volume class." << std::endl;
	std::cout << "-------------------------------------------------" << std::endl;
}

// Construit la matrice des flux
void FiniteVolume::Build_flux_mat_and_rhs(const double& t)
{
	// Matrix
	_mat_flux.resize(_msh->Get_triangles().size(),_msh->Get_triangles().size());
	// RHS
	_BC_RHS.resize(_msh->Get_triangles().size());
	_BC_RHS.setZero();
	vector<Triplet<double>> triplets;	triplets.clear();
	for (int i = 0; i < _msh->Get_edges().size(); i++)
	{
		int t1 = _msh -> Get_edges()[i].Get_T1();//ref du triangle 1 de l'arrête i
		int t2 = _msh -> Get_edges()[i].Get_T2();//Ref du triangle 2 de l'arrête i
		double pos_t1_abs = _msh -> Get_triangles_center()(t1,0);
		double pos_t1_ord = _msh -> Get_triangles_center()(t1,1);//vecteur pos du triangle T1	
		double vol_1 = _msh -> Get_triangles_area()(t1);//Volume/aire du triangle T1
		double lenght = _msh -> Get_edges_length()[i];//longueur de l'arrête i
		//centre de l'arrête i
		double pos_center_i_abs = _msh -> Get_edges_center()(i,0);
		double pos_center_i_ord = _msh -> Get_edges_center()(i,1);
		// vitesse au ventre de l'arête
		double v_x = _fct -> Velocity_x(pos_center_i_abs, pos_center_i_ord, t);
		double v_y = _fct -> Velocity_y(pos_center_i_abs, pos_center_i_ord, t);
		// normale de l'arête
		double n_e_x = _msh -> Get_edges_normal()(i,0);
		double n_e_y = _msh -> Get_edges_normal()(i,1);
		double vn = v_x*n_e_x + v_y*n_e_y;

		double alpha;
		double beta;
		double c;

		if (_df -> Get_numerical_flux_choice() == "rusanov") //Si schéma Rusanov
		{
			alpha=1/2*(v_x + v_y) - c;
			beta =1/2*(v_x + v_y) + c;
		}

		if (t2!=-1)
		{
			double pos_t2_abs = _msh -> Get_triangles_center()(t2,0);
			double pos_t2_ord = _msh -> Get_triangles_center()(t2,1);
			double delta = sqrt(pow(pos_t1_abs-pos_t2_abs,2)+pow(pos_t1_ord-pos_t2_ord,2));//distance entre le triangle T1 et le triangle T2
			double vol_2 = _msh -> Get_triangles_area()(t2);

			triplets.push_back({t1,t1,lenght*alpha/vol_1});//(i,i)
			triplets.push_back({t1,t2,lenght*beta/vol_1});//(i,k)
			triplets.push_back({t2,t1,-lenght*alpha/vol_2});//(k,i)
			triplets.push_back({t2,t2,-lenght*beta/vol_2});//(k,k)
		}
		else //Condition au bord
		{
			//distance triangle centre i et arrete i
			double delta = sqrt(pow(pos_t1_abs-pos_center_i_abs,2)+pow(pos_t1_ord-pos_center_i_ord,2));
			

			if (_msh -> Get_edges()[i].Get_BC() == "Neumann") // Neumann
			{

			}

			else if (_msh -> Get_edges()[i].Get_BC() == "Dirichlet") // condition de Dirichlet
			{

			}

	}
	_mat_flux.setFromTriplets(triplets.begin(), triplets.end());
}


// --- Déjà implémenté ---
// Construit la condition initiale au centre des triangles
VectorXd FiniteVolume::Initial_condition()
{
	VectorXd sol0(_msh->Get_triangles().size());

	for (int i = 0; i < _msh->Get_triangles().size(); i++)
	sol0(i) = _fct->Initial_condition(_msh->Get_triangles_center()(i,0),
	_msh->Get_triangles_center()(i,1));

	return sol0;
}

// Terme source au centre des triangles
VectorXd FiniteVolume::Source_term(double t)
{
	VectorXd sourceterm(_msh->Get_triangles().size());

	for (int i = 0; i < _msh->Get_triangles().size(); i++)
	{
		sourceterm(i) = _fct->Source_term(_msh->Get_triangles_center()(i,0),
		_msh->Get_triangles_center()(i,1), t);
	}
	return sourceterm;
}

// Solution exacte au centre des triangles
VectorXd FiniteVolume::Exact_solution(const double t)
{
	VectorXd exactsol(_msh->Get_triangles().size());

	for (int i = 0; i < _msh->Get_triangles().size(); i++)
	exactsol(i) = _fct->Exact_solution(_msh->Get_triangles_center()(i,0),
	_msh->Get_triangles_center()(i,1), t);

	return exactsol;
}

// Sauvegarde la solution
void FiniteVolume::Save_sol(const Eigen::VectorXd& sol, int n, std::string st)
{
	double norm = 0;
	for (int i = 0; i < sol.rows(); i++)
	norm += sol(i)*sol(i)*_msh->Get_triangles_area()[i];
	norm = sqrt(norm);

	if (st == "solution")
	{
		cout << "Norme de u = " << norm << endl;
	}

	string name_file = _df->Get_results() + "/" + st + "_" + std::to_string(n) + ".vtk";
	int nb_vert = _msh->Get_vertices().size();
	assert((sol.size() == _msh->Get_triangles().size())
	&& "The size of the solution vector is not the same than the number of _triangles !");

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

	solution << "# vtk DataFile Version 3.0 " << endl;
	solution << "2D Unstructured Grid" << endl;
	solution << "ASCII" << endl;
	solution << "DATASET UNSTRUCTURED_GRID" << endl;

	solution << "POINTS " << nb_vert << " float " << endl;
	for (int i = 0 ; i < nb_vert ; ++i)
	{
		solution << ((_msh->Get_vertices()[i]).Get_coor())[0] << " "
		<< ((_msh->Get_vertices()[i]).Get_coor())[1] << " 0." << endl;
	}
	solution << endl;

	solution << "CELLS " << _msh->Get_triangles().size() << " "
	<< _msh->Get_triangles().size()*4 << endl;
	for (int i = 0 ; i < _msh->Get_triangles().size() ; ++i)
	{
		solution << 3 << " " << ((_msh->Get_triangles()[i]).Get_vertices())[0]
		<< " " << ((_msh->Get_triangles()[i]).Get_vertices())[1]
		<< " " << ((_msh->Get_triangles()[i]).Get_vertices())[2] << endl;
	}
	solution << endl;

	solution << "CELL_TYPES " << _msh->Get_triangles().size() << endl;
	for (int i = 0 ; i < _msh->Get_triangles().size() ; ++i)
	{
		solution << 5 << endl;
	}
	solution << endl;

	solution << "CELL_DATA " << _msh->Get_triangles().size() << endl;
	solution << "SCALARS sol float 1" << endl;
	solution << "LOOKUP_TABLE default" << endl;
	// To avoid strange behaviour (which appear only with Apple)
	// with Paraview when we have very small data (e-35 for example)
	double eps = 1.0e-10;
	for (int i = 0 ; i < _msh->Get_triangles().size() ; ++i)
	{
		solution << max(eps,sol[i]) << endl;
	}
	solution << endl;

	//solution << "CELL_DATA " << _msh->Get_triangles().size() << endl;
	solution << "SCALARS CFL float 1" << endl;
	solution << "LOOKUP_TABLE default" << endl;
	// To avoid strange behaviour (which appear only with Apple)
	// with Paraview when we have very small data (e-35 for example)
	for (int i = 0 ; i < _msh->Get_triangles().size() ; ++i)
	{
		solution << max(eps,_df->Get_dt()*fabs(sol[i])/_msh->Get_triangles_length()(i)) << endl;
	}
	solution << endl;

	if (_df->Get_mu() > 1e-10)
	{
		solution << "SCALARS Pe float 1" << endl;
		solution << "LOOKUP_TABLE default" << endl;
		// To avoid strange behaviour (which appear only with Apple)
		// with Paraview when we have very small data (e-35 for example)
		for (int i = 0 ; i < _msh->Get_triangles().size() ; ++i)
		{
			solution << max(eps,_msh->Get_triangles_length()(i)*fabs(sol[i])/_df->Get_mu()) << endl;
		}
		solution << endl;
	}

	solution.close();
}

#define _FINITEVOLUME_CPP
#endif
