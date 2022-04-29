#ifndef _DATA_FILE_CPP

#include "DataFile.h"
#include "../../toml.hpp"
#include <fstream>
#include <iostream>

using namespace std;

DataFile::DataFile(std::string file_name)
: _file_name(file_name)
{
  // Lecture du fichier de données
  auto config = toml::parse(file_name);

  // Other
  const auto& other = toml::find(config, "other");
  _mesh_name = toml::find<std::string>(other, "mesh");
  _mu = toml::find<double>(other, "mu");
  _numerical_flux_choice = toml::find<std::string>(other, "numerical_flux");
  _results = toml::find<std::string>(other, "results");

  // Time
  const auto& time = toml::find(config, "time");
  _t0 = toml::find<double>(time, "t0");
  _tfinal = toml::find<double>(time, "tfinal");
  _dt = toml::find<double>(time, "dt");
  _scheme = toml::find<std::string>(time, "scheme");

  // Boundary conditions
  const auto& BC = toml::find(config, "BC");
  _BC_ref = toml::find<std::vector<int> >(BC, "ref");
  _BC_type = toml::find<std::vector<std::string> >(BC, "BC");

  // Scenarii
  const auto& scenarii = toml::find(config, "scenarii");
  _which_scenario = toml::find<std::string>(scenarii, "which_scenario");

  if ((_numerical_flux_choice != "centered") && (_numerical_flux_choice != "upwind"))
  {
    cout << "Only centered and upwind numerical flows are implemented." << endl;
    exit(0);
  }

  if ((_scheme != "ExplicitEuler") && (_scheme != "ImplicitEuler"))
  {
    cout << "Only Explicit Euler and Implicit Euler schemes are implemented." << endl;
    exit(0);
  }
  if (_scheme == "ExplicitEuler")
  {
    cout << "Beware to the CFL condition." << endl;
  }

  if ((_which_scenario == "diffusion_hom_neumann") || (_which_scenario == "diffusion_all_BC")
          || (_which_scenario == "advection_hom_neumann") ||  (_which_scenario == "advection_all_BC")
          || (_which_scenario == "diffusion_advection_all_BC") )
  {
    cout << "-------------------------------------------------" << endl;
    cout << "The test case: " << _which_scenario << " has been chosen." << endl;
    cout << "-------------------------------------------------" << endl;
    if (_mesh_name == "Meshes/square_mini.mesh")
    {
      _print_info = true;
    }
  }
  else if (_which_scenario == "none")
  {
    cout << "-------------------------------------------------" << endl;
    cout << "It is not a test case!" << endl;
    cout << "-------------------------------------------------" << endl;
  }
  else
  {
    cout << "-------------------------------------------------" << endl;
    cout << "A scenario has to be chosen. If it is not a test case, consider none." << endl;
    cout << "-------------------------------------------------" << endl;
    exit(0);
  }


  if ( (_which_scenario == "advection_hom_neumann") && (fabs(_mu) > 1e-6) )
  {
    cout << "Only advection: mu has been fixed at 0." << endl;
    _mu = 0;
  }

  if (_which_scenario == "advection_all_BC")
  {
    _which_scenario = "diffusion_advection_all_BC";
    if (fabs(_mu) > 1e-6)
    {
      cout << "Only advection: mu has been fixed at 0." << endl;
      _mu = 0;
    }
  }

  // Créer le dossier de résultats
  system(("mkdir -p ./" +_results).c_str());
  // Supprimer les anciens résultats
  system(("rm -f ./" +_results + "/*.vtk").c_str());
  // Copier le fichier de données dans le dossier résultats
  system(("cp -r ./" + _file_name + " ./"
          + _results + "/params.txt").c_str());
}

#define _DATA_FILE_CPP
#endif
