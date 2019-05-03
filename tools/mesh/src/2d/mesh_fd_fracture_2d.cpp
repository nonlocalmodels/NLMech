// Copyright (c)     2017
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

#include "mesh_twod.hpp"

static int counter=0;

namespace fd_fracture_2d {

struct FractureData{
   FractureData(){this->method = 0; this->fracture_enabled = false; this->along_x_axis = false; this->p1 = util::point(); this->p2 = util::point();};
   bool fracture_enabled;
   bool along_x_axis;
   int method;
   util::point p1;
   util::point p2;
};

//
// read data from mesh_fd_2d.yaml to generate initial configuration data
// 
void
readDataFile(std::string &geom_dir_file,
      std::string &bc_l_dir_file,
      std::string &bc_r_dir_file,
      std::string &bc_t_dir_file,
      std::string &bc_b_dir_file,
      std::string &bc_l_dir_wc_file,
      std::string &bc_r_dir_wc_file,
      std::string &u_ic_dir_file,
      std::string &v_ic_dir_file,
      std::string &fracture_flag_dir_file,
      std::string &geom_neu_file,
      std::string &bc_l_neu_file,
      std::string &bc_r_neu_file,
      std::string &bc_t_neu_file,
      std::string &bc_b_neu_file,
      std::string &bc_l_neu_wc_file,
      std::string &bc_r_neu_wc_file,
      std::string &u_ic_neu_file,
      std::string &v_ic_neu_file,
      std::string &fracture_flag_neu_file,
      bool & create_dirichlet_files,
      bool & create_neuman_files,      
      std::pair<std::vector<double>, std::vector<double>>   &Domain,
      double                            &Horizon,
      int                               &Ratio,
      double                            &meshSize,
      bool                              &vol_cor_flag,
      size_t &u_ic_flag,
      std::vector<double> &u_ic_params,
      size_t &v_ic_flag,
      std::vector<double> &v_ic_params,
      FractureData & fracture_data, 
      YAML::Node config) {

   // // YAML file
   // config = YAML::LoadFile(config_file);

   // read directory path where input files should be created
   std::string path; 

   std::string dummy;

   //
   // local variable
   //
   std::vector<double> d;

   // read
   if (config["OutputFile"]["Dirichlet_Filename"]) {
      create_dirichlet_files = true;
      
      util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename", "Path");
      path = config["OutputFile"]["Dirichlet_Filename"]["Path"].as<std::string>();

      util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename", "Geometry_File");
      dummy = config["OutputFile"]["Dirichlet_Filename"]["Geometry_File"].as<std::string>();
      geom_dir_file.append(path);
      geom_dir_file.append("/");
      geom_dir_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename", "BC_Displ_ID_Left_File");
      dummy = config["OutputFile"]["Dirichlet_Filename"]["BC_Displ_ID_Left_File"].as<std::string>();
      bc_l_dir_file.append(path);
      bc_l_dir_file.append("/");
      bc_l_dir_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename", "BC_Displ_ID_Right_File");
      dummy = config["OutputFile"]["Dirichlet_Filename"]["BC_Displ_ID_Right_File"].as<std::string>();
      bc_r_dir_file.append(path);
      bc_r_dir_file.append("/");
      bc_r_dir_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename", "BC_Displ_ID_Top_File");
      dummy = config["OutputFile"]["Dirichlet_Filename"]["BC_Displ_ID_Top_File"].as<std::string>();
      bc_t_dir_file.append(path);
      bc_t_dir_file.append("/");
      bc_t_dir_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename", "BC_Displ_ID_Bottom_File");
      dummy = config["OutputFile"]["Dirichlet_Filename"]["BC_Displ_ID_Bottom_File"].as<std::string>();
      bc_b_dir_file.append(path);
      bc_b_dir_file.append("/");
      bc_b_dir_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename", "IC_Displ_File");
      dummy = config["OutputFile"]["Dirichlet_Filename"]["IC_Displ_File"].as<std::string>();
      u_ic_dir_file.append(path);
      u_ic_dir_file.append("/");
      u_ic_dir_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename", "IC_Vel_File");
      dummy = config["OutputFile"]["Dirichlet_Filename"]["IC_Vel_File"].as<std::string>();
      v_ic_dir_file.append(path);
      v_ic_dir_file.append("/");
      v_ic_dir_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename", "IC_Vel_File");
      if (config["OutputFile"]["Dirichlet_Filename"]["Fracture_Flag_File"])
         dummy = config["OutputFile"]["Dirichlet_Filename"]["Fracture_Flag_File"].as<std::string>();
      else 
         dummy = "fracture_flag_dir.csv";

      fracture_flag_dir_file.append(path);
      fracture_flag_dir_file.append("/");
      fracture_flag_dir_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename", "BC_Displ_ID_Left_File");
      dummy = config["OutputFile"]["Dirichlet_Filename"]["BC_Displ_ID_Left_File"].as<std::string>();
      bc_l_dir_wc_file.append(path);
      bc_l_dir_wc_file.append("/");

      // remove csv from dummy
      std::string ext = ".csv";
      size_t found = dummy.find(ext);
      if (found != std::string::npos)
         dummy.erase(found, ext.length());

      // add wc to dummy
      dummy.append("_wc");

      // add csv back
      dummy.append(ext);

      bc_l_dir_wc_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename", "BC_Displ_ID_Right_File");
      dummy = config["OutputFile"]["Dirichlet_Filename"]["BC_Displ_ID_Right_File"].as<std::string>();
      bc_r_dir_wc_file.append(path);
      bc_r_dir_wc_file.append("/");

      // remove csv from dummy
      ext = ".csv";
      found = dummy.find(ext);
      if (found != std::string::npos)
         dummy.erase(found, ext.length());

      // add wc to dummy
      dummy.append("_wc");

      // add csv back
      dummy.append(ext);

      bc_r_dir_wc_file.append(dummy);
   }

   if (config["OutputFile"]["Neuman_Filename"]) {
      create_neuman_files = true;

      util::methods::checkForData(config, "OutputFile", "Neuman_Filename", "Path");
      path = config["OutputFile"]["Neuman_Filename"]["Path"].as<std::string>();

      util::methods::checkForData(config, "OutputFile", "Neuman_Filename", "Geometry_File");
      dummy = config["OutputFile"]["Neuman_Filename"]["Geometry_File"].as<std::string>();
      geom_neu_file.append(path);
      geom_neu_file.append("/");
      geom_neu_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Neuman_Filename", "BC_Displ_ID_Left_File");
      dummy = config["OutputFile"]["Neuman_Filename"]["BC_Displ_ID_Left_File"].as<std::string>();
      bc_l_neu_file.append(path);
      bc_l_neu_file.append("/");
      bc_l_neu_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Neuman_Filename", "BC_Displ_ID_Right_File");
      dummy = config["OutputFile"]["Neuman_Filename"]["BC_Displ_ID_Right_File"].as<std::string>();
      bc_r_neu_file.append(path);
      bc_r_neu_file.append("/");
      bc_r_neu_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Neuman_Filename", "BC_Displ_ID_Top_File");
      dummy = config["OutputFile"]["Neuman_Filename"]["BC_Displ_ID_Top_File"].as<std::string>();
      bc_t_neu_file.append(path);
      bc_t_neu_file.append("/");
      bc_t_neu_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Neuman_Filename", "BC_Displ_ID_Bottom_File");
      dummy = config["OutputFile"]["Neuman_Filename"]["BC_Displ_ID_Bottom_File"].as<std::string>();
      bc_b_neu_file.append(path);
      bc_b_neu_file.append("/");
      bc_b_neu_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Neuman_Filename", "IC_Displ_File");
      dummy = config["OutputFile"]["Neuman_Filename"]["IC_Displ_File"].as<std::string>();
      u_ic_neu_file.append(path);
      u_ic_neu_file.append("/");
      u_ic_neu_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Neuman_Filename", "IC_Vel_File");
      dummy = config["OutputFile"]["Neuman_Filename"]["IC_Vel_File"].as<std::string>();
      v_ic_neu_file.append(path);
      v_ic_neu_file.append("/");
      v_ic_neu_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Neuman_Filename", "IC_Vel_File");
      if (config["OutputFile"]["Neuman_Filename"]["Fracture_Flag_File"])
         dummy = config["OutputFile"]["Neuman_Filename"]["Fracture_Flag_File"].as<std::string>();
      else 
         dummy = "fracture_flag_neu.csv";
      
      fracture_flag_neu_file.append(path);
      fracture_flag_neu_file.append("/");
      fracture_flag_neu_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Neuman_Filename", "BC_Displ_ID_Left_File");
      dummy = config["OutputFile"]["Neuman_Filename"]["BC_Displ_ID_Left_File"].as<std::string>();
      bc_l_neu_wc_file.append(path);
      bc_l_neu_wc_file.append("/");

      // remove csv from dummy
      std::string ext = ".csv";
      size_t found = dummy.find(ext);
      if (found != std::string::npos)
         dummy.erase(found, ext.length());

      // add wc to dummy
      dummy.append("_wc");

      // add csv back
      dummy.append(ext);

      bc_l_neu_wc_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Neuman_Filename", "BC_Displ_ID_Right_File");
      dummy = config["OutputFile"]["Neuman_Filename"]["BC_Displ_ID_Right_File"].as<std::string>();
      bc_r_neu_wc_file.append(path);
      bc_r_neu_wc_file.append("/");

      // remove csv from dummy
      ext = ".csv";
      found = dummy.find(ext);
      if (found != std::string::npos)
         dummy.erase(found, ext.length());

      // add wc to dummy
      dummy.append("_wc");

      // add csv back
      dummy.append(ext);

      bc_r_neu_wc_file.append(dummy);
   }

   if (config["MeshData"]) {
      for (auto e : config["MeshData"]["Domain"]["Left_Bottom"]) 
         d.push_back(e.as<double>());

      Domain.first.resize(2);
      Domain.first[0] = d[0];
      Domain.first[1] = d[1];
      
      d.clear();
      for (auto e : config["MeshData"]["Domain"]["Right_Top"]) 
         d.push_back(e.as<double>());

      Domain.second.resize(2);
      Domain.second[0] = d[0];
      Domain.second[1] = d[1];

      Horizon = config["MeshData"]["Horizon"].as<double>();

      Ratio = config["MeshData"]["Horizon_Mesh_Ratio"].as<int>();

      // mesh size 
      meshSize = Horizon / (double(Ratio));

      // check if we need to output corrected volume
      vol_cor_flag = true;
      if (config["MeshData"]["Vol_Correct"])
         vol_cor_flag = config["MeshData"]["Vol_Correct"].as<bool>();
   }

   // no test for fracture mesh data
   bool testFlag = false;

   if (testFlag == false) {
      if (config["ICData"]["Displacement"]) {
    
         u_ic_flag = 
            config["ICData"]["Displacement"]["Flag"].as<size_t>();

         u_ic_params.clear();
         for (auto e : config["ICData"]["Displacement"]["Parameters"])
            u_ic_params.push_back(e.as<double>());
      }
      else
         u_ic_flag = 0;

      if (config["ICData"]["Velocity"]) {  

         v_ic_flag = 
            config["ICData"]["Velocity"]["Flag"].as<size_t>();

         v_ic_params.clear();
         for (auto e : config["ICData"]["Velocity"]["Parameters"])
            v_ic_params.push_back(e.as<double>());
      }
      else
         v_ic_flag = 0;
   }

   // read fracture data
   if (config["FractureData"]) {

      fracture_data.fracture_enabled = true;

   	// check if method is specified, else assign default value
   	if (config["FractureData"]["Method"])
   		fracture_data.method = config["FractureData"]["Method"].as<int>();
   	else
   		fracture_data.method = 0;

      // check if crack line is along x or y axis
      if (config["FractureData"]["Along_X_Axis"])
         fracture_data.along_x_axis = config["FractureData"]["Along_X_Axis"].as<bool>();
      else
         fracture_data.along_x_axis = false;

      if (config["FractureData"]["Loc_X"]) {

         double loc_x = config["FractureData"]["Loc_X"].as<double>();
         double loc_y = config["FractureData"]["Loc_Y"].as<double>();

         fracture_data.p1 = util::point(loc_x, 0.0, 0.0);
         fracture_data.p2 = util::point(loc_x, loc_y, 0.0);
      }

      if (config["FractureData"]["P1"]) {


         // read point 1
         auto e = config["FractureData"]["P1"];

         std::vector<double> locs;
         for (auto j : e)
            locs.push_back(j.as<double>());

         if (locs.size() < 2) {

            std::cerr<<"Error: Need atleast two data in P1.\n";
            exit(1);
         }

         fracture_data.p1 = util::point();
         for (size_t i = 0; i<locs.size(); i++)
            fracture_data.p1[i] = locs[i];

         // read point 2
         auto f = config["FractureData"]["P2"];

         locs.clear();
         for (auto j : f)
            locs.push_back(j.as<double>());

         if (locs.size() < 2) {

            std::cerr<<"Error: Need atleast two data in P2.\n";
            exit(1);
         }

         fracture_data.p2 = util::point();
         for (size_t i = 0; i<locs.size(); i++)
            fracture_data.p2[i] = locs[i];
      }      
   }
   else
      fracture_data.fracture_enabled = false;

   return;
} 

//
//	getGuassianField()
//
util::point getGuassianField(util::point x,
	util::point x_center,
	util::point Direction,
	double A,
	double Beta) {

	util::point v = util::point();
	
	util::point x_diff = x - x_center;
	double abs_x_diff = x_diff.length();

	v.x = A * std::exp(-abs_x_diff * abs_x_diff / Beta) * Direction.x;
	v.y = A * std::exp(-abs_x_diff * abs_x_diff / Beta) * Direction.y;

	return v;
}

//
// apply initial condition
//
void computeIC(util::point xi, 
   util::point &ui, 
   util::point &vi,
   std::pair<std::vector<double>, std::vector<double>> domain,
   bool test_flag,
   size_t u_ic_flag,
   std::vector<double> u_ic_params,
   size_t v_ic_flag,
   std::vector<double> v_ic_params) {

   // displacement
   if(u_ic_flag == 0) {
      ui.x = 0.0; 
      ui.y = 0.0;
   }
   else if(u_ic_flag == 1) {

   	util::point center = util::point();
   	util::point direction = util::point();
   	center.x = u_ic_params[0];
   	center.y = u_ic_params[1];

   	direction.x = u_ic_params[2];
   	direction.y = u_ic_params[3];

   	double beta = u_ic_params[4];
   	double a = u_ic_params[5];

   	ui = getGuassianField(xi, center, direction, a, beta);
   }
   else if(u_ic_flag == 2) {

   	util::point center = util::point();
   	util::point direction = util::point();

   	// first pulse
   	center.x = u_ic_params[0];
   	center.y = u_ic_params[1];

   	direction.x = u_ic_params[4];
   	direction.y = u_ic_params[5];

   	double beta = u_ic_params[8];
   	double a = u_ic_params[9];

   	ui = getGuassianField(xi, center, direction, a, beta);

   	// second pulse
   	center.x = u_ic_params[2];
   	center.y = u_ic_params[3];

   	direction.x = u_ic_params[6];
   	direction.y = u_ic_params[7];      	
   	
   	ui = ui + getGuassianField(xi, center, direction, a, beta);
   }            
   else {
      std::cerr <<"Mesh: Error in displacement ic flag = "<<u_ic_flag<<"\n";
      exit(EXIT_FAILURE);
   }

   // velocity
   if(v_ic_flag == 0) {
      vi.x = 0.0; 
      vi.y = 0.0; 
   }
   else if(v_ic_flag == 1) {

   	util::point center = util::point();
   	util::point direction = util::point();
   	center.x = v_ic_params[0];
   	center.y = v_ic_params[1];

   	direction.x = v_ic_params[2];
   	direction.y = v_ic_params[3];

   	double beta = v_ic_params[4];
   	double a = v_ic_params[5];

   	vi = getGuassianField(xi, center, direction, a, beta);
   }
   else if(v_ic_flag == 2) {

   	util::point center = util::point();
   	util::point direction = util::point();

   	// first pulse
   	center.x = v_ic_params[0];
   	center.y = v_ic_params[1];

   	direction.x = v_ic_params[4];
   	direction.y = v_ic_params[5];

   	double beta = v_ic_params[8];
   	double a = v_ic_params[9];

   	vi = getGuassianField(xi, center, direction, a, beta);

   	// second pulse
   	center.x = v_ic_params[2];
   	center.y = v_ic_params[3];

   	direction.x = v_ic_params[6];
   	direction.y = v_ic_params[7];      	
   	
   	vi = vi + getGuassianField(xi, center, direction, a, beta);
   }            
   else {
      std::cerr <<"Mesh: Error in velocity ic flag = "<<v_ic_flag<<"\n";
      exit(EXIT_FAILURE);
   }
   return;
}

} // namespace fd_fracture_2d

void mesh::fd_fracture_twod(YAML::Node config) {

   // //  local variables
   // YAML::Node config;
   // std::string config_filename = "mesh_fd_fracture_2d.yaml";

   std::pair<std::vector<double>, std::vector<double>> domain;
   double horizon;
   int ratio;
   double mesh_size;
   size_t displ_ic_flag;
   std::vector<double> displ_ic_params;
   size_t vel_ic_flag;
   std::vector<double> vel_ic_params;
   int center[3];

   bool test_flag = false; // no test for fracture problem

   // flag for volume correction
   bool vol_cor_flag = true;

   bool create_dirichlet_files = false;
   std::string geom_diri_filename;
   std::string displ_bc_left_diri_filename;
   std::string displ_bc_right_diri_filename;
   std::string displ_bc_top_diri_filename;
   std::string displ_bc_bottom_diri_filename;
   std::string displ_ic_diri_filename;
   std::string vel_ic_diri_filename;
   std::string fracture_flag_diri_filename;
   std::string displ_bc_left_diri_wc_filename;
   std::string displ_bc_right_diri_wc_filename;

   bool create_neuman_files = false;
   std::string geom_neum_filename;
   std::string displ_bc_left_neum_filename;
   std::string displ_bc_right_neum_filename;
   std::string displ_bc_top_neum_filename;
   std::string displ_bc_bottom_neum_filename;
   std::string displ_ic_neum_filename;
   std::string vel_ic_neum_filename;   
   std::string fracture_flag_neum_filename;
   std::string displ_bc_left_neum_wc_filename;
   std::string displ_bc_right_neum_wc_filename;

   // for fracture
   fd_fracture_2d::FractureData fracture_data;

   // read data from config file
   fd_fracture_2d::readDataFile(geom_diri_filename,
      displ_bc_left_diri_filename,
      displ_bc_right_diri_filename,
      displ_bc_top_diri_filename,
      displ_bc_bottom_diri_filename,
      displ_bc_left_diri_wc_filename,
      displ_bc_right_diri_wc_filename,
      displ_ic_diri_filename,
      vel_ic_diri_filename,
      fracture_flag_diri_filename,
      geom_neum_filename,
      displ_bc_left_neum_filename,
      displ_bc_right_neum_filename,
      displ_bc_top_neum_filename,
      displ_bc_bottom_neum_filename,
      displ_bc_left_neum_wc_filename,
      displ_bc_right_neum_wc_filename,
      displ_ic_neum_filename,
      vel_ic_neum_filename,
      fracture_flag_neum_filename,
      create_dirichlet_files,
      create_neuman_files,
      domain,
      horizon,
      ratio,
      mesh_size,
      vol_cor_flag,
      displ_ic_flag,
      displ_ic_params,
      vel_ic_flag,
      vel_ic_params,
      fracture_data,
      config);

   if (fracture_data.along_x_axis == true) {

      std::cout<<"Warning: Crack along x-axis has only been implemented for Neumann type with method = 0.\n";
   }

   //
   // Perform check locationComparedToCrackLine()
   //
   bool check_localcompare = false;
   if (check_localcompare) {

      // equation y = x + 1
      double x1 = -3.0, y1 = - 2.0; 
      double x2 = 2, y2 = 3.0;

      util::point p1(x1, y1, 0.0);
      util::point p2(x2, y2, 0.0);

      double x3 = 0.0, y3 = x3+1;

      double x4 = 2.2, y4 = x4+1;

      double x5 = -3.2, y5 = x5+1;

      util::point p_left(x3, y3 + 0.2, 0.0);

      util::point p_right(x3, y3 - 0.2, 0.0);

      util::point p_on(x3, y3, 0.0);

      util::point p_on_left_tol(x3, y3 + 0.1*mesh_size/10.0, 0.0);

      util::point p_on_right_tol(x3, y3 - 0.1*mesh_size/10.0, 0.0);

      util::point p_up(x4, y4, 0.0);

      util::point p_down(x5, y5, 0.0);

      // util::point p_left_up(-0.1, 1.2, 0.0);

      // util::point p_left_down(-0.1, -0.2, 0.0);

      // util::point p_right_up(0.1, 1.2, 0.0);
      // util::point p_right_down(-0.1, -0.2, 0.0);

      // check
      std::cout<<"*********** p_left\n";
      int loc = util::methods::locationComparedToCrackLine(p_left, p1, p2, mesh_size, true);
      if (loc != 0) {
         std::cout<<"Error: debug failed for p_left locationComparedToCrackLine().\n";
      }
      std::cout<<"***********\n";

      std::cout<<"\n";
      std::cout<<"*********** p_right\n";
      loc = util::methods::locationComparedToCrackLine(p_right, p1, p2, mesh_size, true);
      if (loc != 0) {
         std::cout<<"Error: debug failed for p_right locationComparedToCrackLine().\n";
      }
      std::cout<<"***********\n";

      std::cout<<"\n";
      std::cout<<"*********** p_on\n";      
      loc = util::methods::locationComparedToCrackLine(p_on, p1, p2, mesh_size, true);
      if (loc != 1) {
         std::cout<<"Error: debug failed for p_on locationComparedToCrackLine().\n";
      }
      std::cout<<"***********\n";

      std::cout<<"\n";
      std::cout<<"*********** p_on_left_tol\n";
      loc = util::methods::locationComparedToCrackLine(p_on_left_tol, p1, p2, mesh_size, true);
      if (loc != -1) {
         std::cout<<"Error: debug failed for p_on_left_tol locationComparedToCrackLine().\n";
      }
      std::cout<<"***********\n";

      std::cout<<"\n";
      std::cout<<"*********** p_on_right_tol\n";
      loc = util::methods::locationComparedToCrackLine(p_on_right_tol, p1, p2, mesh_size, true);
      if (loc != 1) {
         std::cout<<"Error: debug failed for p_on_right_tol locationComparedToCrackLine().\n";
      }
      std::cout<<"***********\n";

      std::cout<<"\n";
      std::cout<<"*********** p_up\n";
      loc = util::methods::locationComparedToCrackLine(p_up, p1, p2, mesh_size, true);
      if (loc != 0) {
         std::cout<<"Error: debug failed for p_up locationComparedToCrackLine().\n";
      }
      std::cout<<"***********\n";

      std::cout<<"\n";
      std::cout<<"*********** p_down\n";
      loc = util::methods::locationComparedToCrackLine(p_down, p1, p2, mesh_size, true);
      if (loc != 0) {
         std::cout<<"Error: debug failed for p_down locationComparedToCrackLine().\n";
      }
      std::cout<<"***********\n";

      // std::cout<<"\n";
      // std::cout<<"*********** p_left_up\n";
      // loc = util::methods::locationComparedToCrackLine(p_left_up, p1, p2, mesh_size, true);
      // if (loc != 0) {
      //    std::cout<<"Error: debug failed for p_left_up locationComparedToCrackLine().\n";
      // }
      // std::cout<<"***********\n";

      // std::cout<<"\n";
      // std::cout<<"*********** p_left_down\n";
      // loc = util::methods::locationComparedToCrackLine(p_left_down, p1, p2, mesh_size, true);
      // if (loc != 0) {
      //    std::cout<<"Error: debug failed for p_left_down locationComparedToCrackLine().\n";
      // }
      // std::cout<<"***********\n";

      // std::cout<<"\n";
      // std::cout<<"*********** p_right_down\n";
      // loc = util::methods::locationComparedToCrackLine(p_right_down, p1, p2, mesh_size, true);
      // if (loc != 0) {
      //    std::cout<<"Error: debug failed for p_right_down locationComparedToCrackLine().\n";
      // }
      // std::cout<<"***********\n";

      // std::cout<<"\n";
      // std::cout<<"*********** p_right_up\n";
      // loc = util::methods::locationComparedToCrackLine(p_right_up, p1, p2, mesh_size, true);
      // if (loc != 0) {
      //    std::cout<<"Error: debug failed for p_right_up locationComparedToCrackLine().\n";
      // }
   } 


   // 
   // first get data corresponding to Dirichlet bc
   //
   double tol = 1.0E-12;
   if (create_dirichlet_files == true) {

      //
      // open necessary files
      //
      std::ofstream F_dir_geom;
      std::ofstream F_dir_u_l;
      std::ofstream F_dir_u_r;
      std::ofstream F_dir_u_t;
      std::ofstream F_dir_u_b;
      std::ofstream F_dir_u_ic;
      std::ofstream F_dir_v_ic;
      std::ofstream F_dir_frac_flag;
      std::ofstream F_dir_u_l_wc;
      std::ofstream F_dir_u_r_wc;

      // std::string file_test;
      // file_test.append("./a/");
      // file_test.append(geom_diri_filename);
      // F_dir_geom.open(file_test);

      F_dir_geom.open(geom_diri_filename);
      // header
      F_dir_geom<<"id,x,y,volume\n";

      
      F_dir_u_l.open(displ_bc_left_diri_filename);
      // header
      F_dir_u_l<<"id\n";
      
      
      F_dir_u_r.open(displ_bc_right_diri_filename);
      // header
      F_dir_u_r<<"id\n";


      F_dir_u_t.open(displ_bc_top_diri_filename);
      // header
      F_dir_u_t<<"id\n";
      
      
      F_dir_u_b.open(displ_bc_bottom_diri_filename);
      // header
      F_dir_u_b<<"id\n";      

      
      if (displ_ic_flag != 0) {
         F_dir_u_ic.open(displ_ic_diri_filename);
         // header
         F_dir_u_ic<<"id,x,y\n";
      }

      if (vel_ic_flag != 0) {
         F_dir_v_ic.open(vel_ic_diri_filename);
         // header
         F_dir_v_ic<<"id,x,y\n";   
      }   

      if (fracture_data.fracture_enabled == true) {
         F_dir_frac_flag.open(fracture_flag_diri_filename);
         // header
         F_dir_frac_flag<<"id,fracture_flag\n";
      }

      F_dir_u_l_wc.open(displ_bc_left_diri_wc_filename);
      // header
      F_dir_u_l_wc<<"id\n";
      
      
      F_dir_u_r_wc.open(displ_bc_right_diri_wc_filename);
      // header
      F_dir_u_r_wc<<"id\n";

      // // old calc
      // // discretize D union D_nonlocal
      // double x_left = domain.first[0] - horizon - 2.0 * mesh_size;
      // double x_right = domain.second[0] + horizon + 2.0 * mesh_size;
      // double y_bottom = domain.first[1] - horizon - 2.0 * mesh_size;
      // double y_top = domain.second[1] + horizon + 2.0 * mesh_size;

      // int n_nodes = 0; // node counter
      // int nx = int((x_right - x_left)/mesh_size);
      // int ny = int((y_top - y_bottom)/mesh_size);

      // discretize D union D_nonlocal
      // modify domain parameters 
      std::pair<std::vector<double>, std::vector<double>> nonlocal_domain;
      nonlocal_domain.first.push_back(0.0);
      nonlocal_domain.first.push_back(0.0);
      nonlocal_domain.second.push_back(0.0);
      nonlocal_domain.second.push_back(0.0);

      int eps_h = ratio;

      //
      // modify x right boundary
      //
      int ba_h = int( (domain.second[0] - domain.first[0])/mesh_size );

      // difference in length of interval
      double diff = (domain.second[0] - domain.first[0]) - (ba_h * mesh_size);

      // check if this difference is positive 
      // by design diff will always be nonnegative
      if (compare::definitelyGreaterThan(diff, 0.0)) {

         // difference is positive meaning original interval is greater
         // than ba_h * mesh_size

         // check if the difference is above half og mesh_size
         // if yes then redefine the end of interval 
         if (compare::definitelyGreaterThan(diff, 0.5 * mesh_size)) 
            domain.second[0] = domain.first[0] + (ba_h+1) * mesh_size;
         else 
            domain.second[0] = domain.first[0] + ba_h * mesh_size;
      }  
      
      //
      // modify y top boundary
      //
      ba_h = int( (domain.second[1] - domain.first[1])/mesh_size );

      // difference in length of interval
      diff = (domain.second[1] - domain.first[1]) - (ba_h * mesh_size);

      // check if this difference is positive 
      // by design diff will always be nonnegative
      if (compare::definitelyGreaterThan(diff, 0.0)) {

         // difference is positive meaning original interval is greater
         // than ba_h * mesh_size

         // check if the difference is above half og mesh_size
         // if yes then redefine the end of interval 
         if (compare::definitelyGreaterThan(diff, 0.5 * mesh_size)) 
            domain.second[1] = domain.first[1] + (ba_h+1) * mesh_size;
         else 
            domain.second[1] = domain.first[1] + ba_h * mesh_size;
      }

      //
      // nonlocal boundary
      //

      // nonlocal x left
      nonlocal_domain.first[0] = domain.first[0] - double(eps_h)*mesh_size;

      // nonlocal x right
      nonlocal_domain.second[0] = domain.second[0] + double(eps_h)*mesh_size;

      // nonlocal y bottom 
      nonlocal_domain.first[1] = domain.first[1] - double(eps_h)*mesh_size;

      // nonlocal y top 
      nonlocal_domain.second[1] = domain.second[1] + double(eps_h)*mesh_size;

      // number of divisions in x direction 
      int nx = int((nonlocal_domain.second[0] - nonlocal_domain.first[0])/mesh_size);

      // number of divisions in y direction 
      int ny = int((nonlocal_domain.second[1] - nonlocal_domain.first[1])/mesh_size);

      // total number of nodes
      int total_num_nodes = (nx+1)*(ny+1);

      int n_nodes = 0;

      // For dirichlet, we have changed the boundary so we need
      // to adjust the crack location as well.
      std::cout<<"***            ***\n";
      std::cout<<"Warning: For Dirichlet data, modify the location of crack accordingly.\n"; 
      std::cout<<"***            ***\n";
      
      // we have two methods to create a crack. We choose one as follows.
      int method = 0;
      if (fracture_data.fracture_enabled == true and fracture_data.method == 1)
         method = 1;

      if (method == 1) {

      	// in this method, we do duplicate the nodes on crack line. We mark
      	// one set of duplicate nodes as -1 and other as +1.

	      for (int i=0;i<=nx;i++) {
	      	for (int j=0;j<=ny;j++) {

	      		util::point x = util::point();
	      		x.x = nonlocal_domain.first[0] + double(i) * mesh_size;
	      		x.y = nonlocal_domain.first[1] + double(j) * mesh_size;

	      		double vol_i = mesh_size * mesh_size;

	            // flag = 0 node is not on fracture edge
	            // flag = -1 node is on fracture edge and is on left side
	            // flag = 1 node is on fracture edge and is on right side
	            int fracture_flag = 0;
	            bool duplicate = false;

	            // crack is a line between two points
               // find the location of x wrt to crack line
               bool crack_along_x_axis = false;
               fracture_flag = util::methods::locationComparedToCrackLineAlongAxis(x, fracture_data.p1, fracture_data.p2, mesh_size, crack_along_x_axis);

               if (fracture_flag == -1 or fracture_flag == 1)
                  duplicate = true;

	            if (duplicate == false) {

	            	// node is not on crack edge

	            	// compute the volume of node correctly
		      		if (i == 0 or i == nx) {

		      			if (j == 0 or j == ny) {
		      				// corner node
		      				vol_i = mesh_size * mesh_size / 4.0;
		      			}
		      			else {
		      				// not a corner but a boundary node
		      				vol_i = mesh_size * mesh_size / 2.0;
		      			}
		      		}

		      		if (j == 0 or j == ny) {

		      			if (i == 0 or i == nx) {
		      				// corner node
		      				vol_i = mesh_size * mesh_size / 4.0;
		      			}
		      			else {
		      				// not a corner but a boundary node
		      				vol_i = mesh_size * mesh_size / 2.0;
		      			}
		      		}

		      		// write to geometry file
		      		F_dir_geom<<n_nodes<<","<<x.x<<","<<x.y<<","<<vol_i<<"\n";

                  F_dir_frac_flag<<n_nodes<<","<<fracture_flag<<"\n";
	            }
	            else {

	            	// create two set of nodes with same coordinate, mark one 
	            	// as left side of crack and other as right side of crack.

	            	// compute volume correctly
	            	vol_i = mesh_size * mesh_size / 2.0;

	            	// if node is on the bottom edge, modify the volume further
	            	if (compare::definitelyLessThan(x.y, 1.0E-10))
	            		vol_i = vol_i / 2.0;

	            	// first write left side node
	            	fracture_flag = -1;

		      		F_dir_geom<<n_nodes<<","<<x.x<<","<<x.y<<","<<vol_i<<"\n";

                  F_dir_frac_flag<<n_nodes<<","<<fracture_flag<<"\n";

		      		// write right side node
		      		fracture_flag = 1;

		      		F_dir_geom<<n_nodes + 1<<","<<x.x<<","<<x.y<<","<<vol_i<<"\n";

                  F_dir_frac_flag<<n_nodes + 1<<","<<fracture_flag<<"\n";
	            }

	      		// get initial condition on nodes
	      		util::point u = util::point();
	      		util::point v = util::point();
	   			
	   			fd_fracture_2d::computeIC(x, u, v, domain, test_flag,displ_ic_flag, displ_ic_params,	vel_ic_flag, vel_ic_params);

	      		// check whether node is inside, in left nonlocal boundary,
	      		// right nonlocal boundary, top nonlocal boundary
	      		// and bottom nonlocal boundary
	      		if (x.y <= domain.first[1]) {

	      			if (duplicate == false)
                     F_dir_u_b<<n_nodes<<"\n";
	      			else {
	      				F_dir_u_b<<n_nodes<<"\n";
	      				F_dir_u_b<<n_nodes + 1<<"\n";
	      			}
	      		}
	      		
               if (x.y >= domain.second[1]) 
	      			F_dir_u_t<<n_nodes<<"\n";
	      		
               if (x.x <= domain.first[0] and (x.y > domain.first[1] and x.y < domain.second[1])) 
	      			F_dir_u_l<<n_nodes<<"\n";
	      		
               if (x.x >= domain.second[0] and (x.y > domain.first[1] and x.y < domain.second[1])) 
	      			F_dir_u_r<<n_nodes<<"\n";

               if (x.x <= domain.first[0]) 
                  F_dir_u_l_wc<<n_nodes<<"\n";
               
               if (x.x >= domain.second[0]) 
                  F_dir_u_r_wc<<n_nodes<<"\n";

	      		if (duplicate == false) {
	            	if (displ_ic_flag != 0)
                     F_dir_u_ic<<n_nodes<<","<<u.x<<","<<u.y<<"\n";
                  if (vel_ic_flag != 0)
                     F_dir_v_ic<<n_nodes<<","<<v.x<<","<<v.y<<"\n";
	            }
	            else {

                  if (displ_ic_flag != 0)
                     F_dir_u_ic<<n_nodes<<","<<u.x<<","<<u.y<<"\n";
                  if (vel_ic_flag != 0)
                     F_dir_v_ic<<n_nodes<<","<<v.x<<","<<v.y<<"\n";

                  if (displ_ic_flag != 0)
                     F_dir_u_ic<<n_nodes+1<<","<<u.x<<","<<u.y<<"\n";
                  if (vel_ic_flag != 0)
                     F_dir_v_ic<<n_nodes+1<<","<<v.x<<","<<v.y<<"\n";
	            }
	      	
	         	// increment the mesh node counter
	         	if (duplicate == false)
	         		n_nodes = n_nodes + 1;
	         	else
	         		n_nodes = n_nodes + 2;
	      	} // loop over j
	      } // loop over i
	   } // if method 1
	   else if(method == 0) {

	   	// in this method, we do not remove or add node. We simply shift the
	   	// the crack line to position such that it does not intersect any
	   	// node. We then indicate the node which lies in left and right
	   	// side of crack line as -1 and +1.

	   	if (fracture_data.fracture_enabled == true) {
            // check loc_x
            int r_loc_x_h = fracture_data.p1.x/mesh_size;

            if (std::abs(fracture_data.p1.x - double(r_loc_x_h)*mesh_size) < 1.0E-12) {

               // loc_x is multiple of mesh size thus we shift it a bit
               fracture_data.p1.x += 5.0E-9;

               fracture_data.p2.x += 5.0E-9;
            }
         }

	      for (int i=0;i<=nx;i++) {
	      	for (int j=0;j<=ny;j++) {

	      		util::point x = util::point();
	      		x.x = nonlocal_domain.first[0] + double(i) * mesh_size;
	      		x.y = nonlocal_domain.first[1] + double(j) * mesh_size;

	      		double vol_i = mesh_size * mesh_size;


	            // flag = 0 node is not on fracture edge
               // flag = -1 node is on fracture edge and is on left side
               // flag = 1 node is on fracture edge and is on right side
               int fracture_flag = 0;
               bool duplicate = false;
               
               if (fracture_data.fracture_enabled == true) {
               
                  // crack is a line between two points
                  // find the location of x wrt to crack line
                  bool crack_along_x_axis = false;
                  fracture_flag = util::methods::locationComparedToCrackLineAlongAxis(x, fracture_data.p1, fracture_data.p2, mesh_size, crack_along_x_axis);
               }

            	// compute the volume of node correctly
	      		if (i == 0 or i == nx) {

	      			if (j == 0 or j == ny) {
	      				// corner node
	      				vol_i = mesh_size * mesh_size / 4.0;
	      			}
	      			else {
	      				// not a corner but a boundary node
	      				vol_i = mesh_size * mesh_size / 2.0;
	      			}
	      		}

	      		if (j == 0 or j == ny) {

	      			if (i == 0 or i == nx) {
	      				// corner node
	      				vol_i = mesh_size * mesh_size / 4.0;
	      			}
	      			else {
	      				// not a corner but a boundary node
	      				vol_i = mesh_size * mesh_size / 2.0;
	      			}
	      		}

		      	// write to geometry file
		      	F_dir_geom<<n_nodes<<","<<x.x<<","<<x.y<<","<<vol_i<<"\n";

               if (fracture_data.fracture_enabled == true)
                  F_dir_frac_flag<<n_nodes<<","<<fracture_flag<<"\n";

	      		// get initial condition on nodes
	      		util::point u = util::point();
	      		util::point v = util::point();
	   			
	   			fd_fracture_2d::computeIC(x, u, v, domain, test_flag,displ_ic_flag, displ_ic_params,	vel_ic_flag, vel_ic_params);

	      		// check whether node is inside, in left nonlocal boundary,
	      		// right nonlocal boundary, top nonlocal boundary
	      		// and bottom nonlocal boundary
	      		if (x.y <= domain.first[1])
      				F_dir_u_b<<n_nodes<<"\n";

               if (x.y >= domain.second[1]) 
	      			F_dir_u_t<<n_nodes<<"\n";
	      		
               if (x.x <= domain.first[0] 
                  and (compare::definitelyGreaterThan(x.y, domain.first[1]) 
                  and compare::definitelyLessThan(x.y, domain.second[1]) ) )
                  F_dir_u_l<<n_nodes<<"\n";

               if (x.x >= domain.second[0] 
                  and ( compare::definitelyGreaterThan(x.y, domain.first[1])
                  and   compare::definitelyLessThan(x.y, domain.second[1]) ) ) 
                  F_dir_u_r<<n_nodes<<"\n";

               if (x.x <= domain.first[0])
                  F_dir_u_l_wc<<n_nodes<<"\n";
               
               if (x.x >= domain.second[0])
                  F_dir_u_r_wc<<n_nodes<<"\n";

               if (displ_ic_flag != 0)
                  F_dir_u_ic<<n_nodes<<","<<u.x<<","<<u.y<<"\n";
               if (vel_ic_flag != 0)
                  F_dir_v_ic<<n_nodes<<","<<v.x<<","<<v.y<<"\n";
	      	
         		n_nodes = n_nodes + 1;
	      	} // loop over j
	      } // loop over i
	   } // if method 0

      std::cout<<"Total number of nodes (Dirichlet Data) = "<<n_nodes<<"\n";

      if (fracture_data.fracture_enabled == true) {

         std::cout<<"Fracture p1 = ("<<fracture_data.p1.x<<","<<fracture_data.p1.y<<"), p2 = ("<<fracture_data.p2.x<<","<<fracture_data.p2.y<<") \n";
      }

      std::cout<<"Domain = ["<<domain.first[0]<<","<<domain.first[1]<<"]x["<<domain.second[0]<<","<<domain.second[1]<<"]\n";
      std::cout<<"Nonlocal domain = ["<<nonlocal_domain.first[0]<<","<<nonlocal_domain.first[1]<<"]x["<<nonlocal_domain.second[0]<<","<<nonlocal_domain.second[1]<<"]\n";

      // close all files
      F_dir_geom.close();
      F_dir_u_l.close();
      F_dir_u_r.close();
      F_dir_u_t.close();
      F_dir_u_b.close();
      if (displ_ic_flag != 0) 
         F_dir_u_ic.close();
      if (vel_ic_flag != 0) 
         F_dir_v_ic.close();
      if (fracture_data.fracture_enabled == true)
         F_dir_frac_flag.close();
      F_dir_u_l_wc.close();
      F_dir_u_r_wc.close();
   }

   if (create_neuman_files == true) {

      // open files for corresponding data Neuman
      std::ofstream F_neu_geom;
      std::ofstream F_neu_u_l;
      std::ofstream F_neu_u_r;
      std::ofstream F_neu_u_t;
      std::ofstream F_neu_u_b;
      std::ofstream F_neu_u_ic;
      std::ofstream F_neu_v_ic;
      std::ofstream F_neu_frac_flag;
      std::ofstream F_neu_u_l_wc;
      std::ofstream F_neu_u_r_wc;

      F_neu_geom.open(geom_neum_filename);
      // header
      F_neu_geom<<"id,x,y,volume\n";

      
      F_neu_u_l.open(displ_bc_left_neum_filename);
      // header
      F_neu_u_l<<"id\n";
      
      
      F_neu_u_r.open(displ_bc_right_neum_filename);
      // header
      F_neu_u_r<<"id\n";


      F_neu_u_t.open(displ_bc_top_neum_filename);
      // header
      F_neu_u_t<<"id\n";
      
      
      F_neu_u_b.open(displ_bc_bottom_neum_filename);
      // header
      F_neu_u_b<<"id\n"; 
      
      if (displ_ic_flag != 0) {
         F_neu_u_ic.open(displ_ic_neum_filename);
         // header
         F_neu_u_ic<<"id,x,y\n";
      }

      if (vel_ic_flag != 0) {
         F_neu_v_ic.open(vel_ic_neum_filename);
         // header
         F_neu_v_ic<<"id,x,y\n";
      }

      if (fracture_data.fracture_enabled == true) {

         F_neu_frac_flag.open(fracture_flag_neum_filename);
         // header
         F_neu_frac_flag<<"id,fracture_flag\n";
      }

      F_neu_u_l_wc.open(displ_bc_left_neum_wc_filename);
      // header
      F_neu_u_l_wc<<"id\n";
      
      
      F_neu_u_r_wc.open(displ_bc_right_neum_wc_filename);
      // header
      F_neu_u_r_wc<<"id\n";

      //
      // modify boundary so that discretization is good
      //

      //
      // modify x right
      //
      int ba_h = int( (domain.second[0] - domain.first[0])/mesh_size );
      
      // difference in length of interval
      double diff = (domain.second[0] - domain.first[0]) - (ba_h * mesh_size);

      // check if this difference is positive 
      // by design diff will always be nonnegative
      if (compare::definitelyGreaterThan(diff, 0.0)) {

         // difference is positive meaning original interval is greater
         // than ba_h * mesh_size

         // check if the difference is above half og mesh_size
         // if yes then redefine the end of interval 
         if (compare::definitelyGreaterThan(diff, 0.5 * mesh_size)) 
            domain.second[0] = domain.first[0] + (ba_h+1) * mesh_size;
         else 
            domain.second[0] = domain.first[0] + ba_h * mesh_size;
      }

      //
      // modify y top
      //
      ba_h = int( (domain.second[1] - domain.first[1])/mesh_size );

      // difference in length of interval
      diff = (domain.second[1] - domain.first[1]) - (ba_h * mesh_size);

      // check if this difference is positive 
      // by design diff will always be nonnegative
      if (compare::definitelyGreaterThan(diff, 0.0)) {

         // difference is positive meaning original interval is greater
         // than ba_h * mesh_size

         // check if the difference is above half og mesh_size
         // if yes then redefine the end of interval 
         if (compare::definitelyGreaterThan(diff, 0.5 * mesh_size)) 
            domain.second[1] = domain.first[1] + (ba_h+1) * mesh_size;
         else 
            domain.second[1] = domain.first[1] + ba_h * mesh_size;
      }

      // number of divisions in x direction 
      int nx = int((domain.second[0] - domain.first[0])/mesh_size);

      // number of divisions in y direction 
      int ny = int((domain.second[1] - domain.first[1])/mesh_size);

      // total number of nodes
      int total_num_nodes = (nx+1)*(ny+1);

      int n_nodes = 0;

      // we have two method to create a crack. We choose one as follows.
      int method = 0;
      if (fracture_data.fracture_enabled == true and fracture_data.method == 1)
         method = 1;

      if (method == 1) {

      	// in this method, we do duplicate the nodes on crack line. We mark
      	// one set of duplicate nodes as -1 and other as +1.

	      for (int i=0;i<=nx;i++) {
	         for (int j=0;j<=ny;j++) {

	            util::point x = util::point();
	            x.x = domain.first[0] + double(i) * mesh_size;
	            x.y = domain.first[1] + double(j) * mesh_size;

	            double vol_i = mesh_size * mesh_size;

	            // flag = 0 node is not on fracture edge
	            // flag = -1 node is on fracture edge and is on left side
	            // flag = 1 node is on fracture edge and is on right side
	            int fracture_flag = 0;
	            bool duplicate = false;

               // crack is a line between two points
               // find the location of x wrt to crack line
               bool crack_along_x_axis = false;
               fracture_flag = util::methods::locationComparedToCrackLineAlongAxis(x, fracture_data.p1, fracture_data.p2, mesh_size, crack_along_x_axis);

               if (fracture_flag == -1 or fracture_flag == 1)
                  duplicate = true;

	            if (duplicate == false) {

	            	// node is not on crack edge

	            	// compute the volume of node correctly
		      		if (i == 0 or i == nx) {

		      			if (j == 0 or j == ny) {
		      				// corner node
		      				vol_i = mesh_size * mesh_size / 4.0;
		      			}
		      			else {
		      				// not a corner but a boundary node
		      				vol_i = mesh_size * mesh_size / 2.0;
		      			}
		      		}

		      		if (j == 0 or j == ny) {

		      			if (i == 0 or i == nx) {
		      				// corner node
		      				vol_i = mesh_size * mesh_size / 4.0;
		      			}
		      			else {
		      				// not a corner but a boundary node
		      				vol_i = mesh_size * mesh_size / 2.0;
		      			}
		      		}

		      		// write to geometry file
		      		F_neu_geom<<n_nodes<<","<<x.x<<","<<x.y<<","<<vol_i<<"\n";
                  F_neu_frac_flag<<n_nodes<<","<<fracture_flag<<"\n";   
	            }
	            else {

	            	// create two set of nodes with same coordinate, mark one 
	            	// as left side of crack and other as right side of crack.

	            	// compute volume correctly
	            	vol_i = mesh_size * mesh_size / 2.0;

	            	// if node is on the bottom edge, modify the volume further
	            	if (x.y < 1.0E-10)
	            		vol_i = vol_i / 2.0;

	            	// first write left side node
	            	fracture_flag = -1;

		      		F_neu_geom<<n_nodes<<","<<x.x<<","<<x.y<<","<<vol_i<<"\n";
                  F_neu_frac_flag<<n_nodes<<","<<fracture_flag<<"\n";

		      		// write right side node
		      		fracture_flag = 1;

		      		F_neu_geom<<n_nodes + 1<<","<<x.x<<","<<x.y<<","<<vol_i<<"\n";
                  F_neu_frac_flag<<n_nodes + 1<<","<<fracture_flag<<"\n";
	            }

	            // get initial condition on nodes
	            util::point u = util::point();
	            util::point v = util::point();
	            
	            fd_fracture_2d::computeIC(x, u, v, domain, test_flag,
	               displ_ic_flag, displ_ic_params,
	               vel_ic_flag, vel_ic_params);

	            // check whether node is inside, in left nonlocal boundary,
	            // right nonlocal boundary, top nonlocal boundary
	            // and bottom nonlocal boundary
	            if (j == 0) {

	               if (duplicate == false )
	               	F_neu_u_b<<n_nodes<<"\n";
	               else {
	               	F_neu_u_b<<n_nodes<<"\n";
	               	F_neu_u_b<<n_nodes+1<<"\n";
	               }
	            }

	            if (j == ny) 
	               F_neu_u_t<<n_nodes<<"\n";

	            if (i == 0) 
	               if (j !=0 and j != ny)
	                  F_neu_u_l<<n_nodes<<"\n";

	            if (i == nx) 
	               if (j !=0 and j != ny)
	                  F_neu_u_r<<n_nodes<<"\n";

               if (i == 0) 
                  F_neu_u_l_wc<<n_nodes<<"\n";

               if (i == nx) 
                  F_neu_u_r_wc<<n_nodes<<"\n";

	            // x is inside material domain
	            if (duplicate == false) {
	            	
                  if (displ_ic_flag != 0)
                     F_neu_u_ic<<n_nodes<<","<<u.x<<","<<u.y<<"\n";
                  if (vel_ic_flag != 0)
	            	   F_neu_v_ic<<n_nodes<<","<<v.x<<","<<v.y<<"\n";
	            }
	            else {

                  if (displ_ic_flag != 0)
                     F_neu_u_ic<<n_nodes<<","<<u.x<<","<<u.y<<"\n";
                  if (vel_ic_flag != 0)
                     F_neu_v_ic<<n_nodes<<","<<v.x<<","<<v.y<<"\n";

                  if (displ_ic_flag != 0)
                     F_neu_u_ic<<n_nodes+1<<","<<u.x<<","<<u.y<<"\n";
                  if (vel_ic_flag != 0)
                     F_neu_v_ic<<n_nodes+1<<","<<v.x<<","<<v.y<<"\n";
	            }
	         
	         	// increment the mesh node counter
	         	if (duplicate == false)
	         		n_nodes = n_nodes + 1;
	         	else
	         		n_nodes = n_nodes + 2;
	         } // loop over j
	      } // loop over i
	   }// if method 1
	   else if (method == 0) { 

	   	// in this method, we do not remove or add node. We simply shift the
	   	// the crack line to position such that it does not intersect any
	   	// node. We then indicate the node which lies in left and right
	   	// side of crack line as -1 and +1.

         if (fracture_data.fracture_enabled == true) {

            if (fracture_data.along_x_axis == true) {

               // check loc_y
               int r_loc_y_h = fracture_data.p1.y/mesh_size;

               if (std::abs(fracture_data.p1.y - double(r_loc_y_h)*mesh_size) < 1.0E-12) {

                  // loc_y is multiple of mesh size thus we shift it a bit
                  fracture_data.p1.y += 5.0E-9;

                  fracture_data.p2.y += 5.0E-9;
               }
            }
            else {
               // check loc_x
               int r_loc_x_h = fracture_data.p1.x/mesh_size;

               if (std::abs(fracture_data.p1.x - double(r_loc_x_h)*mesh_size) < 1.0E-12) {

                  // loc_x is multiple of mesh size thus we shift it a bit
                  fracture_data.p1.x += 5.0E-9;

                  fracture_data.p2.x += 5.0E-9;
               }
            }
         }

	      for (int i=0;i<=nx;i++) {
	         for (int j=0;j<=ny;j++) {

	            util::point x = util::point();
	            x.x = domain.first[0] + double(i) * mesh_size;
	            x.y = domain.first[1] + double(j) * mesh_size;

	            double vol_i = mesh_size * mesh_size;

               // flag = 0 node is not on fracture edge
               // flag = -1 node is on fracture edge and is on left side
               // flag = 1 node is on fracture edge and is on right side
               int fracture_flag = 0;
               bool duplicate = false;
               
               if (fracture_data.fracture_enabled == true) {
	            
                  // crack is a line between two points
                  // find the location of x wrt to crack line
                  fracture_flag = util::methods::locationComparedToCrackLineAlongAxis(x, fracture_data.p1, fracture_data.p2, mesh_size, fracture_data.along_x_axis);
               }

            	// compute the volume of node correctly
	      		if (i == 0 or i == nx) {

	      			if (j == 0 or j == ny) {
	      				// corner node
	      				vol_i = mesh_size * mesh_size / 4.0;
	      			}
	      			else {
	      				// not a corner but a boundary node
	      				vol_i = mesh_size * mesh_size / 2.0;
	      			}
	      		}

	      		if (j == 0 or j == ny) {

	      			if (i == 0 or i == nx) {
	      				// corner node
	      				vol_i = mesh_size * mesh_size / 4.0;
	      			}
	      			else {
	      				// not a corner but a boundary node
	      				vol_i = mesh_size * mesh_size / 2.0;
	      			}
	      		}

		      	// write to geometry file
		      	F_neu_geom<<n_nodes<<","<<x.x<<","<<x.y<<","<<vol_i<<"\n";

               if (fracture_data.fracture_enabled == true) 
                  F_neu_frac_flag<<n_nodes<<","<<fracture_flag<<"\n";

	            // get initial condition on nodes
	            util::point u = util::point();
	            util::point v = util::point();
	            
	            fd_fracture_2d::computeIC(x, u, v, domain, test_flag,
	               displ_ic_flag, displ_ic_params,
	               vel_ic_flag, vel_ic_params);

	            // check whether node is inside, in left nonlocal boundary,
	            // right nonlocal boundary, top nonlocal boundary
	            // and bottom nonlocal boundary
	            if (j == 0) 
               	F_neu_u_b<<n_nodes<<"\n";

	            if (j == ny) 
	               F_neu_u_t<<n_nodes<<"\n";

	            if (i == 0) 
	               if (j !=0 and j != ny)
	                  F_neu_u_l<<n_nodes<<"\n";

	            if (i == nx) 
	               if (j !=0 and j != ny)
	                  F_neu_u_r<<n_nodes<<"\n";

               if (i == 0) 
                  F_neu_u_l_wc<<n_nodes<<"\n";

               if (i == nx) 
                  F_neu_u_r_wc<<n_nodes<<"\n";

	            // x is inside material domain
               if (displ_ic_flag != 0)
                  F_neu_u_ic<<n_nodes<<","<<u.x<<","<<u.y<<"\n";
               if (vel_ic_flag != 0)
                  F_neu_v_ic<<n_nodes<<","<<v.x<<","<<v.y<<"\n";
	         
	         	n_nodes = n_nodes + 1;
	         } // loop over j
	      } // loop over i
	   } // if method 0


      if (fracture_data.fracture_enabled == true) {
         std::cout<<"Fracture p1 = ("<<fracture_data.p1.x<<","<<fracture_data.p1.y<<"), p2 = ("<<fracture_data.p2.x<<","<<fracture_data.p2.y<<") \n";
      }

      std::cout<<"Total number of nodes (Neumann Data) = "<<n_nodes<<"\n";
      std::cout<<"Domain = ["<<domain.first[0]<<","<<domain.first[1]<<"]x["<<domain.second[0]<<","<<domain.second[1]<<"]\n";

      // close all files
      F_neu_geom.close();
      F_neu_u_l.close();
      F_neu_u_r.close();
      F_neu_u_t.close();
      F_neu_u_b.close();
      if (displ_ic_flag != 0) 
         F_neu_u_ic.close();
      if (vel_ic_flag != 0) 
         F_neu_v_ic.close();
      if (fracture_data.fracture_enabled == true) 
         F_neu_frac_flag.close();
      F_neu_u_l_wc.close();
      F_neu_u_r_wc.close();
   }   
}
