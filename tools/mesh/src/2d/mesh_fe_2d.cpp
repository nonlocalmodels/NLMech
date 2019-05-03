// Copyright (c)		2017 Prashant K. Jha
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txts

#include "mesh_twod.hpp"

namespace fe_2d {

//
// read data from mesh_fd_2d.yaml to generate initial configuration data
// 
void
readDataFile(std::string &path_dir,
      std::string &mesh_dir_file,
      std::string &bc_l_dir_file,
      std::string &bc_r_dir_file,
      std::string &bc_t_dir_file,
      std::string &bc_b_dir_file,
      std::string &bc_l_dir_wc_file,
      std::string &bc_r_dir_wc_file,
      std::string &path_neu,
      std::string &mesh_neu_file,
      std::string &bc_l_neu_file,
      std::string &bc_r_neu_file,
      std::string &bc_t_neu_file,
      std::string &bc_b_neu_file,
      std::string &bc_l_neu_wc_file,
      std::string &bc_r_neu_wc_file,
      bool & create_dirichlet_files,
      bool & create_neuman_files,
      std::pair<std::vector<double>, std::vector<double>>   &Domain,
      double                            &Horizon,
      int                               &Ratio,
      double                            &meshSize,
      std::string 							 &mesh_type,
      size_t &u_ic_flag,
      std::vector<double> &u_ic_params,
      size_t &v_ic_flag,
      std::vector<double> &v_ic_params,
      bool & testFlag,
      bool & bending_out,
      YAML::Node config) {

   // // YAML file
   // config = YAML::LoadFile(config_file);

   //
   // local variable
   //
   std::vector<double> d;
   std::string dummy;

   // read
   if (config["OutputFile"]["Dirichlet_Filename"]) {
      create_dirichlet_files = true;

      util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename", "Path");
      path_dir = config["OutputFile"]["Dirichlet_Filename"]["Path"].as<std::string>();

      util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename", "Fe_Mesh_File");
      dummy = config["OutputFile"]["Dirichlet_Filename"]["Fe_Mesh_File"].as<std::string>();
      // mesh_dir_file.append(path);
      // mesh_dir_file.append("/");

      // check if file_name has .vtu at end if yes then remove it
      std::string ext = ".vtu";
      size_t found = dummy.find(ext);
      if (found != std::string::npos)
         dummy.erase(found, ext.length());

      mesh_dir_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename", "BC_Displ_ID_Left_File");
      dummy = config["OutputFile"]["Dirichlet_Filename"]["BC_Displ_ID_Left_File"].as<std::string>();
      bc_l_dir_file.append(path_dir);
      bc_l_dir_file.append("/");
      bc_l_dir_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename", "BC_Displ_ID_Right_File");
      dummy = config["OutputFile"]["Dirichlet_Filename"]["BC_Displ_ID_Right_File"].as<std::string>();
      bc_r_dir_file.append(path_dir);
      bc_r_dir_file.append("/");
      bc_r_dir_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename", "BC_Displ_ID_Top_File");
      dummy = config["OutputFile"]["Dirichlet_Filename"]["BC_Displ_ID_Top_File"].as<std::string>();
      bc_t_dir_file.append(path_dir);
      bc_t_dir_file.append("/");
      bc_t_dir_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename", "BC_Displ_ID_Bottom_File");
      dummy = config["OutputFile"]["Dirichlet_Filename"]["BC_Displ_ID_Bottom_File"].as<std::string>();
      bc_b_dir_file.append(path_dir);
      bc_b_dir_file.append("/");
      bc_b_dir_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename", "BC_Displ_ID_Left_File");
      dummy = config["OutputFile"]["Dirichlet_Filename"]["BC_Displ_ID_Left_File"].as<std::string>();
      bc_l_dir_wc_file.append(path_dir);
      bc_l_dir_wc_file.append("/");

      // remove csv from dummy
      ext = ".csv";
      found = dummy.find(ext);
      if (found != std::string::npos)
         dummy.erase(found, ext.length());

      // add wc to dummy
      dummy.append("_wc");

      // add csv back
      dummy.append(ext);

      // add dummy to filename
      bc_l_dir_wc_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename", "BC_Displ_ID_Right_File");
      dummy = config["OutputFile"]["Dirichlet_Filename"]["BC_Displ_ID_Right_File"].as<std::string>();
      bc_r_dir_wc_file.append(path_dir);
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

      // add dummy to filename
      bc_r_dir_wc_file.append(dummy);
   }

   if (config["OutputFile"]["Neuman_Filename"]) {
      create_neuman_files = true;

      util::methods::checkForData(config, "OutputFile", "Neuman_Filename", "Path");
      path_neu = config["OutputFile"]["Neuman_Filename"]["Path"].as<std::string>();

      util::methods::checkForData(config, "OutputFile", "Neuman_Filename", "Fe_Mesh_File");
      dummy = config["OutputFile"]["Neuman_Filename"]["Fe_Mesh_File"].as<std::string>();
      // mesh_neu_file.append(path);
      // mesh_neu_file.append("/");

      // check if file_name has .vtu at end. If yes then remove it
      std::string ext = ".vtu";
      size_t found = dummy.find(ext);
      if (found != std::string::npos)
         dummy.erase(found, ext.length());

      mesh_neu_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Neuman_Filename", "BC_Displ_ID_Left_File");
      dummy = config["OutputFile"]["Neuman_Filename"]["BC_Displ_ID_Left_File"].as<std::string>();
      bc_l_neu_file.append(path_neu);
      bc_l_neu_file.append("/");
      bc_l_neu_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Neuman_Filename", "BC_Displ_ID_Right_File");
      dummy = config["OutputFile"]["Neuman_Filename"]["BC_Displ_ID_Right_File"].as<std::string>();
      bc_r_neu_file.append(path_neu);
      bc_r_neu_file.append("/");
      bc_r_neu_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Neuman_Filename", "BC_Displ_ID_Top_File");
      dummy = config["OutputFile"]["Neuman_Filename"]["BC_Displ_ID_Top_File"].as<std::string>();
      bc_t_neu_file.append(path_neu);
      bc_t_neu_file.append("/");
      bc_t_neu_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Neuman_Filename", "BC_Displ_ID_Bottom_File");
      dummy = config["OutputFile"]["Neuman_Filename"]["BC_Displ_ID_Bottom_File"].as<std::string>();
      bc_b_neu_file.append(path_neu);
      bc_b_neu_file.append("/");
      bc_b_neu_file.append(dummy);

      util::methods::checkForData(config, "OutputFile", "Neuman_Filename", "BC_Displ_ID_Left_File");
      dummy = config["OutputFile"]["Neuman_Filename"]["BC_Displ_ID_Left_File"].as<std::string>();
      bc_l_neu_wc_file.append(path_neu);
      bc_l_neu_wc_file.append("/");

      // remove csv from dummy
      ext = ".csv";
      found = dummy.find(ext);
      if (found != std::string::npos)
         dummy.erase(found, ext.length());

      // const std::string ext(".csv");
      // if ( dummy != ext and 
      //      dummy.size() > ext.size() and
      //      dummy.substr(dummy.size() - ext.size()) == ".gz" ) {

      //    // if so then strip them off
      //    dummy = dummy.substr(0, dummy.size() - ext.size());
      // }

      // add wc to dummy
      dummy.append("_wc");

      // add csv back
      dummy.append(ext);

      // add dummy to filename
      bc_l_neu_wc_file.append(dummy);

      // std::cout<<bc_l_neu_wc_file<<"\n";

      util::methods::checkForData(config, "OutputFile", "Neuman_Filename", "BC_Displ_ID_Right_File");
      dummy = config["OutputFile"]["Neuman_Filename"]["BC_Displ_ID_Right_File"].as<std::string>();
      bc_r_neu_wc_file.append(path_neu);
      bc_r_neu_wc_file.append("/");

      // remove csv from dummy
      ext = ".csv";
      found = dummy.find(ext);
      if (found != std::string::npos)
         dummy.erase(found, ext.length());

      // if ( dummy != ext and 
      //      dummy.size() > ext.size() and
      //      dummy.substr(dummy.size() - ext.size()) == ".gz" ) {

      //    // if so then strip them off
      //    dummy = dummy.substr(0, dummy.size() - ext.size());
      // }

      // add wc to dummy
      dummy.append("_wc");

      // add csv back
      dummy.append(ext);

      // add dummy to filename
      bc_r_neu_wc_file.append(dummy);

      // std::cout<<bc_r_neu_wc_file<<"\n";
   } 

   util::methods::checkForData(config, "MeshData");
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

      // mesh type: either Uniform_Square or Uniform_Triangle
      if (!config["MeshData"]["Mesh_Type"]) {
      	std::cerr<<"Mesh_Type not specified.\n";
      	exit(1);
      }
      mesh_type = config["MeshData"]["Mesh_Type"].as<std::string>();
   }

   // 
   // check if bending data is enabled
   //
   if (config["BendingData"]) 
      bending_out =  true;

   // for fe, no test is implemented yet
   testFlag = false;

   // read ic data only if it is not test
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

   if (test_flag == false) {
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
   }   
   else {
      return;
   }

   return;
}


//
// for bending loading
//
struct bendingData{
   bool th_pt_out;
   bool fr_pt_out;
   util::point sp_pt_l;
   util::point sp_pt_r;
   util::point ld_pt_l;
   util::point ld_pt_r;
};

struct bendingDataFE{
   bool th_pt_out;
   bool fr_pt_out;
   util::point sp_pt_l;
   util::point sp_pt_r;
   util::point sp_pt_ll;
   util::point sp_pt_rr;
   util::point ld_pt_l;
   util::point ld_pt_r;
   util::point ld_pt_ll;
   util::point ld_pt_rr;
};
#define BENDINGDATA_INIT {false, false, util::point(), util::point(), util::point(), util::point()}

#define BENDINGDATAFE_INIT {false, false, util::point(), util::point(), util::point(), util::point(), util::point(), util::point(), util::point(), util::point()}

struct findBendingPointData{
   size_t id_sp_pt_l;
   double dist_sp_pt_l;

   size_t id_sp_pt_r;
   double dist_sp_pt_r;

   size_t id_ld_pt_l;
   double dist_ld_pt_l;

   size_t id_ld_pt_r;
   double dist_ld_pt_r;
};

struct findBendingPointDataFE{
   size_t id_sp_pt_l;
   double dist_sp_pt_l;

   size_t id_sp_pt_ll;
   double dist_sp_pt_ll;

   size_t id_sp_pt_r;
   double dist_sp_pt_r;

   size_t id_sp_pt_rr;
   double dist_sp_pt_rr;

   size_t id_ld_pt_l;
   double dist_ld_pt_l;

   size_t id_ld_pt_ll;
   double dist_ld_pt_ll;

   size_t id_ld_pt_r;
   double dist_ld_pt_r;

   size_t id_ld_pt_rr;
   double dist_ld_pt_rr;
};

#define FINDBENDINGPOINTDATA_INIT {0, 1000.0, 0, 1000.0, 0, 1000.0, 0, 1000.0}

#define FINDBENDINGPOINTDATAFE_INIT {0, 1000.0, 0, 1000.0, 0, 1000.0, 0, 1000.0, 0, 1000.0, 0, 1000.0, 0, 1000.0, 0, 1000.0}


void checkForBendingData(YAML::Node config,
   std::pair<std::vector<double>, std::vector<double>> domain,
   size_t id,
   util::point x,
   bool output = false) {

   static struct bendingData bdData = BENDINGDATA_INIT;
   // bdData.th_pt_out =  false;
   // bdData.fr_pt_out =  false;
   // bdData.sp_pt_l = util::point();
   // bdData.sp_pt_r = util::point();
   // bdData.th_ld_pt = util::point();
   // bdData.fr_ld_pt_l = util::point();
   // bdData.fr_ld_pt_r = util::point();

   static std::ofstream F_sp_l_bd;
   static std::ofstream F_sp_r_bd;
   static std::ofstream F_ld_bd;

   // read data
   static int init = -1;
   if (init == -1) {
      if (config["BendingData"]["Three_Point_Bending"]) {

         bdData.th_pt_out = true;

         // read distance of support point from left edge
         double dist = config["BendingData"]["Three_Point_Bending"]["Distance_Support_Point"].as<double>();

         bdData.sp_pt_l.x = domain.first[0] + dist;
         bdData.sp_pt_l.y = domain.first[1];

         bdData.sp_pt_r.x = domain.second[0] - dist;
         bdData.sp_pt_r.y = domain.first[1]; // it should be first

         // read distance of load point from left edge
         dist = config["BendingData"]["Three_Point_Bending"]["Distance_Load_Point"].as<double>();

         bdData.ld_pt_l.x = domain.first[0] + dist;
         bdData.ld_pt_l.y = domain.second[1];

         //
         // first read the path
         //
         std::string path = config["OutputFile"]["Neuman_Filename"]["Path"].as<std::string>();

         // get filename to output node id corresponding to left support point
         std::string filename;
         filename.append(path);
         filename.append("/");
         std::string dummy = config["BendingData"]["Three_Point_Bending"]["Filename"]["Support_Point_Left_File"].as<std::string>();
         filename.append(dummy);

         // create a file stream and write header
         F_sp_l_bd.open(filename);
         // header
         F_sp_l_bd<<"id\n";     

         // get filename to output node id corresponding to right support point
         dummy = config["BendingData"]["Three_Point_Bending"]["Filename"]["Support_Point_Right_File"].as<std::string>();
         std::string filename2;
         filename2.append(path);
         filename2.append("/");
         filename2.append(dummy);

         // create a file stream and write header
         F_sp_r_bd.open(filename2);
         // header
         F_sp_r_bd<<"id\n";     

         // get filename to output node id corresponding to load point
         dummy = config["BendingData"]["Three_Point_Bending"]["Filename"]["Load_Point_File"].as<std::string>();
         std::string filename3;
         filename3.append(path);
         filename3.append("/");
         filename3.append(dummy);

         // create a file stream and write header
         F_ld_bd.open(filename3);
         // header
         F_ld_bd<<"id\n";
      }

      if (config["BendingData"]["Four_Point_Bending"]) {

         bdData.fr_pt_out = true;

         // read distance of support point from left edge
         double dist = config["BendingData"]["Four_Point_Bending"]["Distance_Support_Point"].as<double>();

         bdData.sp_pt_l.x = domain.first[0] + dist;
         bdData.sp_pt_l.y = domain.first[1];

         bdData.sp_pt_r.x = domain.second[0] - dist;
         bdData.sp_pt_r.y = domain.first[1]; // it should be first

         // read distance of load point from left edge
         dist = config["BendingData"]["Four_Point_Bending"]["Distance_Load_Point"].as<double>();

         bdData.ld_pt_l.x = domain.first[0] + dist;
         bdData.ld_pt_l.y = domain.second[1];

         bdData.ld_pt_r.x = domain.second[0] - dist;
         bdData.ld_pt_r.y = domain.second[1];

         //
         // first read the path
         //
         std::string path = config["OutputFile"]["Neuman_Filename"]["Path"].as<std::string>();

         // get filename to output node id corresponding to left support point
         std::string dummy = config["BendingData"]["Four_Point_Bending"]["Filename"]["Support_Point_Left_File"].as<std::string>();
         std::string filename;
         filename.append(path);
         filename.append("/");
         filename.append(dummy);

         // create a file stream and write header
         F_sp_l_bd.open(filename);
         // header
         F_sp_l_bd<<"id\n";     

         // get filename to output node id corresponding to right support point
         dummy = config["BendingData"]["Four_Point_Bending"]["Filename"]["Support_Point_Right_File"].as<std::string>();
         std::string filename2;
         filename2.append(path);
         filename2.append("/");
         filename2.append(dummy);

         // create a file stream and write header
         F_sp_r_bd.open(filename2);
         // header
         F_sp_r_bd<<"id\n";     

         // get filename to output node id corresponding to load point
         dummy = config["BendingData"]["Four_Point_Bending"]["Filename"]["Load_Point_File"].as<std::string>();
         std::string filename3;
         filename3.append(path);
         filename3.append("/");
         filename3.append(dummy);

         // create a file stream and write header
         F_ld_bd.open(filename3);
         // header
         F_ld_bd<<"id\n";
      }

      // check if three point and four point both data are provided
      if (bdData.th_pt_out ==  true and bdData.fr_pt_out ==  true) {
         std::cerr<<"Check input file. Either three point data or four data point data can be processed not both at the same time.\n";
         exit(1);
      }

      init = 0;
   }

   // to find the closest node to the support points and load points,
   // we need to measure the distance 
   static struct findBendingPointData ptData = FINDBENDINGPOINTDATA_INIT;

   if (output == false) {

      //
      if (bdData.th_pt_out == true) {

         util::point d = x - bdData.sp_pt_l;
         if (d.length() < ptData.dist_sp_pt_l and d.y == 0.0) {
            ptData.id_sp_pt_l = id;
            ptData.dist_sp_pt_l = d.length();
         }

         d = x - bdData.sp_pt_r;
         if (d.length() < ptData.dist_sp_pt_r and d.y == 0.0) {
            ptData.id_sp_pt_r = id;
            ptData.dist_sp_pt_r = d.length();
         }

         d = x - bdData.ld_pt_l;
         if (d.length() < ptData.dist_ld_pt_l and d.y == 0.0) {
            ptData.id_ld_pt_l = id;
            ptData.dist_ld_pt_l = d.length();
         }
      }
      else if (bdData.fr_pt_out == true) {

         util::point d = x - bdData.sp_pt_l;
         if (d.length() < ptData.dist_sp_pt_l and d.y == 0.0) {
            ptData.id_sp_pt_l = id;
            ptData.dist_sp_pt_l = d.length();
         }

         d = x - bdData.sp_pt_r;
         if (d.length() < ptData.dist_sp_pt_r and d.y == 0.0) {
            ptData.id_sp_pt_r = id;
            ptData.dist_sp_pt_r = d.length();
         }

         d = x - bdData.ld_pt_l;
         if (d.length() < ptData.dist_ld_pt_l and d.y == 0.0) {
            ptData.id_ld_pt_l = id;
            ptData.dist_ld_pt_l = d.length();
         }

         d = x - bdData.ld_pt_r;
         if (d.length() < ptData.dist_ld_pt_r and d.y == 0.0) {
            ptData.id_ld_pt_r = id;
            ptData.dist_ld_pt_r = d.length();
         }
      }
   }
   else {

      // we are at the end of node generation, therefore, write the
      // data stored in ptData to corresponding files
      F_sp_l_bd<<ptData.id_sp_pt_l<<"\n";
      F_sp_r_bd<<ptData.id_sp_pt_r<<"\n";

      if (bdData.th_pt_out ==  true) {
         F_ld_bd<<ptData.id_ld_pt_l<<"\n";
      }
      else if (bdData.fr_pt_out ==  true) {
         F_ld_bd<<ptData.id_ld_pt_l<<"\n";
         F_ld_bd<<ptData.id_ld_pt_r<<"\n";
      }
   }
}

void checkForBendingDataFE(YAML::Node config,
   std::pair<std::vector<double>, std::vector<double>> domain,
   double mesh_size,
   size_t id,
   util::point x,
   bool output = false) {

   static struct bendingDataFE bdData = BENDINGDATAFE_INIT;
   // bdData.th_pt_out =  false;
   // bdData.fr_pt_out =  false;
   // bdData.sp_pt_l = util::point();
   // bdData.sp_pt_r = util::point();
   // bdData.th_ld_pt = util::point();
   // bdData.fr_ld_pt_l = util::point();
   // bdData.fr_ld_pt_r = util::point();

   static std::ofstream F_sp_l_bd;
   static std::ofstream F_sp_r_bd;
   static std::ofstream F_ld_bd;

   // read data
   static int init = -1;
   if (init == -1) {
      if (config["BendingData"]["Three_Point_Bending"]) {

         bdData.th_pt_out = true;

         // read distance of support point from left edge
         double dist = config["BendingData"]["Three_Point_Bending"]["Distance_Support_Point"].as<double>();

         bdData.sp_pt_l.x = domain.first[0] + dist;
         bdData.sp_pt_l.y = domain.first[1];

         bdData.sp_pt_r.x = domain.second[0] - dist;
         bdData.sp_pt_r.y = domain.first[1]; // it should be first

         bdData.sp_pt_ll.x = domain.first[0] + dist - mesh_size;
         bdData.sp_pt_ll.y = domain.first[1];

         bdData.sp_pt_rr.x = domain.second[0] - dist + mesh_size;
         bdData.sp_pt_rr.y = domain.first[1]; // it should be first

         // read distance of load point from left edge
         dist = config["BendingData"]["Three_Point_Bending"]["Distance_Load_Point"].as<double>();

         bdData.ld_pt_l.x = domain.first[0] + dist;
         bdData.ld_pt_l.y = domain.second[1];

         bdData.ld_pt_ll.x = domain.first[0] + dist - mesh_size;
         bdData.ld_pt_ll.y = domain.second[1];

         //
         // first read the path
         //
         std::string path = config["OutputFile"]["Neuman_Filename"]["Path"].as<std::string>();

         // get filename to output node id corresponding to left support point
         std::string filename;
         filename.append(path);
         filename.append("/");
         std::string dummy = config["BendingData"]["Three_Point_Bending"]["Filename"]["Support_Point_Left_File"].as<std::string>();
         filename.append(dummy);

         // create a file stream and write header
         F_sp_l_bd.open(filename);
         // header
         F_sp_l_bd<<"id\n";     

         // get filename to output node id corresponding to right support point
         dummy = config["BendingData"]["Three_Point_Bending"]["Filename"]["Support_Point_Right_File"].as<std::string>();
         std::string filename2;
         filename2.append(path);
         filename2.append("/");
         filename2.append(dummy);

         // create a file stream and write header
         F_sp_r_bd.open(filename2);
         // header
         F_sp_r_bd<<"id\n";     

         // get filename to output node id corresponding to load point
         dummy = config["BendingData"]["Three_Point_Bending"]["Filename"]["Load_Point_File"].as<std::string>();
         std::string filename3;
         filename3.append(path);
         filename3.append("/");
         filename3.append(dummy);

         // create a file stream and write header
         F_ld_bd.open(filename3);
         // header
         F_ld_bd<<"id\n";
      }

      if (config["BendingData"]["Four_Point_Bending"]) {

         bdData.fr_pt_out = true;

         // read distance of support point from left edge
         double dist = config["BendingData"]["Four_Point_Bending"]["Distance_Support_Point"].as<double>();

         bdData.sp_pt_l.x = domain.first[0] + dist;
         bdData.sp_pt_l.y = domain.first[1];

         bdData.sp_pt_r.x = domain.second[0] - dist;
         bdData.sp_pt_r.y = domain.first[1]; // it should be first

         bdData.sp_pt_ll.x = domain.first[0] + dist - mesh_size;
         bdData.sp_pt_ll.y = domain.first[1];

         bdData.sp_pt_rr.x = domain.second[0] - dist + mesh_size;
         bdData.sp_pt_rr.y = domain.first[1]; // it should be first

         // read distance of load point from left edge
         dist = config["BendingData"]["Four_Point_Bending"]["Distance_Load_Point"].as<double>();

         bdData.ld_pt_l.x = domain.first[0] + dist;
         bdData.ld_pt_l.y = domain.second[1];

         bdData.ld_pt_r.x = domain.second[0] - dist;
         bdData.ld_pt_r.y = domain.second[1];

         bdData.ld_pt_ll.x = domain.first[0] + dist - mesh_size;
         bdData.ld_pt_ll.y = domain.second[1];

         bdData.ld_pt_rr.x = domain.second[0] - dist + mesh_size;
         bdData.ld_pt_rr.y = domain.second[1];

         //
         // first read the path
         //
         std::string path = config["OutputFile"]["Neuman_Filename"]["Path"].as<std::string>();

         // get filename to output node id corresponding to left support point
         std::string dummy = config["BendingData"]["Four_Point_Bending"]["Filename"]["Support_Point_Left_File"].as<std::string>();
         std::string filename;
         filename.append(path);
         filename.append("/");
         filename.append(dummy);

         // create a file stream and write header
         F_sp_l_bd.open(filename);
         // header
         F_sp_l_bd<<"id\n";     

         // get filename to output node id corresponding to right support point
         dummy = config["BendingData"]["Four_Point_Bending"]["Filename"]["Support_Point_Right_File"].as<std::string>();
         std::string filename2;
         filename2.append(path);
         filename2.append("/");
         filename2.append(dummy);

         // create a file stream and write header
         F_sp_r_bd.open(filename2);
         // header
         F_sp_r_bd<<"id\n";     

         // get filename to output node id corresponding to load point
         dummy = config["BendingData"]["Four_Point_Bending"]["Filename"]["Load_Point_File"].as<std::string>();
         std::string filename3;
         filename3.append(path);
         filename3.append("/");
         filename3.append(dummy);

         // create a file stream and write header
         F_ld_bd.open(filename3);
         // header
         F_ld_bd<<"id\n";
      }

      // check if three point and four point both data are provided
      if (bdData.th_pt_out ==  true and bdData.fr_pt_out ==  true) {
         std::cerr<<"Check input file. Either three point data or four data point data can be processed not both at the same time.\n";
         exit(1);
      }

      init = 0;
   }

   // to find the closest node to the support points and load points,
   // we need to measure the distance 
   static struct findBendingPointDataFE ptData = FINDBENDINGPOINTDATAFE_INIT;

   if (output == false) {

      //
      if (bdData.th_pt_out == true) {

         util::point d = x - bdData.sp_pt_l;
         if (d.length() < ptData.dist_sp_pt_l and d.y == 0.0) {
            ptData.id_sp_pt_l = id;
            ptData.dist_sp_pt_l = d.length();
         }

         d = x - bdData.sp_pt_ll;
         if (d.length() < ptData.dist_sp_pt_ll and d.y == 0.0) {
            ptData.id_sp_pt_ll = id;
            ptData.dist_sp_pt_ll = d.length();
         }

         d = x - bdData.sp_pt_r;
         if (d.length() < ptData.dist_sp_pt_r and d.y == 0.0) {
            ptData.id_sp_pt_r = id;
            ptData.dist_sp_pt_r = d.length();
         }

         d = x - bdData.sp_pt_rr;
         if (d.length() < ptData.dist_sp_pt_rr and d.y == 0.0) {
            ptData.id_sp_pt_rr = id;
            ptData.dist_sp_pt_rr = d.length();
         }

         d = x - bdData.ld_pt_l;
         if (d.length() < ptData.dist_ld_pt_l and d.y == 0.0) {
            ptData.id_ld_pt_l = id;
            ptData.dist_ld_pt_l = d.length();
         }

         d = x - bdData.ld_pt_ll;
         if (d.length() < ptData.dist_ld_pt_ll and d.y == 0.0) {
            ptData.id_ld_pt_ll = id;
            ptData.dist_ld_pt_ll = d.length();
         }
      }
      else if (bdData.fr_pt_out == true) {

         util::point d = x - bdData.sp_pt_l;
         if (d.length() < ptData.dist_sp_pt_l and d.y == 0.0) {
            ptData.id_sp_pt_l = id;
            ptData.dist_sp_pt_l = d.length();
         }

         d = x - bdData.sp_pt_ll;
         if (d.length() < ptData.dist_sp_pt_ll and d.y == 0.0) {
            ptData.id_sp_pt_ll = id;
            ptData.dist_sp_pt_ll = d.length();
         }

         d = x - bdData.sp_pt_r;
         if (d.length() < ptData.dist_sp_pt_r and d.y == 0.0) {
            ptData.id_sp_pt_r = id;
            ptData.dist_sp_pt_r = d.length();
         }

         d = x - bdData.sp_pt_rr;
         if (d.length() < ptData.dist_sp_pt_rr and d.y == 0.0) {
            ptData.id_sp_pt_rr = id;
            ptData.dist_sp_pt_rr = d.length();
         }

         d = x - bdData.ld_pt_l;
         if (d.length() < ptData.dist_ld_pt_l and d.y == 0.0) {
            ptData.id_ld_pt_l = id;
            ptData.dist_ld_pt_l = d.length();
         }

         d = x - bdData.ld_pt_ll;
         if (d.length() < ptData.dist_ld_pt_ll and d.y == 0.0) {
            ptData.id_ld_pt_ll = id;
            ptData.dist_ld_pt_ll = d.length();
         }

         d = x - bdData.ld_pt_r;
         if (d.length() < ptData.dist_ld_pt_r and d.y == 0.0) {
            ptData.id_ld_pt_r = id;
            ptData.dist_ld_pt_r = d.length();
         }

         d = x - bdData.ld_pt_rr;
         if (d.length() < ptData.dist_ld_pt_rr and d.y == 0.0) {
            ptData.id_ld_pt_rr = id;
            ptData.dist_ld_pt_rr = d.length();
         }
      }
   }
   else {

      // we are at the end of node generation, therefore, write the
      // data stored in ptData to corresponding files
      F_sp_l_bd<<ptData.id_sp_pt_l<<"\n";
      F_sp_l_bd<<ptData.id_sp_pt_ll<<"\n";

      F_sp_r_bd<<ptData.id_sp_pt_r<<"\n";
      F_sp_r_bd<<ptData.id_sp_pt_rr<<"\n";

      if (bdData.th_pt_out ==  true) {
         F_ld_bd<<ptData.id_ld_pt_l<<"\n";
         F_ld_bd<<ptData.id_ld_pt_ll<<"\n";
      }
      else if (bdData.fr_pt_out ==  true) {
         F_ld_bd<<ptData.id_ld_pt_l<<"\n";
         F_ld_bd<<ptData.id_ld_pt_ll<<"\n";

         F_ld_bd<<ptData.id_ld_pt_r<<"\n";
         F_ld_bd<<ptData.id_ld_pt_rr<<"\n";
      }
   }
}


//
// create uniform square mesh
//
void UniformSquareMesh(YAML::Node config) {

   // //  local variables
   // YAML::Node config;
   // std::string config_filename = "mesh_fe_2d.yaml";

   std::pair<std::vector<double>, std::vector<double>> domain;
   double horizon;
   int ratio;
   double mesh_size;
   size_t displ_ic_flag;
   std::vector<double> displ_ic_params;
   size_t vel_ic_flag;
   std::vector<double> vel_ic_params;
   int center[3];

   bool test_flag;  // flag if input file is for testing the code
	
	std::string mesh_type;
   
   bool create_dirichlet_files = false;
   std::string path_diri;
   std::string mesh_diri_filename;
   std::string displ_bc_left_diri_filename;
   std::string displ_bc_right_diri_filename;
   std::string displ_bc_top_diri_filename;
   std::string displ_bc_bottom_diri_filename;
   std::string displ_bc_left_diri_wc_filename; // with corner points
   std::string displ_bc_right_diri_wc_filename; // with corner points

   bool create_neuman_files = false;
   std::string path_neum;
   std::string mesh_neum_filename;
   std::string displ_bc_left_neum_filename;
   std::string displ_bc_right_neum_filename;
   std::string displ_bc_top_neum_filename;
   std::string displ_bc_bottom_neum_filename;   
   std::string displ_bc_left_neum_wc_filename; // with corner points
   std::string displ_bc_right_neum_wc_filename; // with corner points  

   // bending loading related data (Implemented only for neuman type)
   bool bd_dout = false;

   // read data from config file
   readDataFile(path_diri,
      mesh_diri_filename,
      displ_bc_left_diri_filename,
      displ_bc_right_diri_filename,
      displ_bc_top_diri_filename,
      displ_bc_bottom_diri_filename,
      displ_bc_left_diri_wc_filename,
      displ_bc_right_diri_wc_filename,
      path_neum,
      mesh_neum_filename,
      displ_bc_left_neum_filename,
      displ_bc_right_neum_filename,
      displ_bc_top_neum_filename,
      displ_bc_bottom_neum_filename,
      displ_bc_left_neum_wc_filename,
      displ_bc_right_neum_wc_filename,
      create_dirichlet_files,
      create_neuman_files,
      domain,
      horizon,
      ratio,
      mesh_size,
      mesh_type,
      displ_ic_flag,
      displ_ic_params,
      vel_ic_flag,
      vel_ic_params,
      test_flag,
      bd_dout,
      config);

   // 
   // first get data corresponding to Dirichlet bc
   //
   double tol = 1.0E-12;
   if (create_dirichlet_files == true) {

      //
      //	open necessary files
      //
      std::ofstream File_disp_diri_left;
      File_disp_diri_left.open(displ_bc_left_diri_filename);
      // header
      File_disp_diri_left<<"id\n";
      
      std::ofstream File_disp_diri_right;
      File_disp_diri_right.open(displ_bc_right_diri_filename);
      // header
      File_disp_diri_right<<"id\n";

      std::ofstream File_disp_diri_top;
      File_disp_diri_top.open(displ_bc_top_diri_filename);
      // header
      File_disp_diri_top<<"id\n";
      
      std::ofstream File_disp_diri_bottom;
      File_disp_diri_bottom.open(displ_bc_bottom_diri_filename);
      // header
      File_disp_diri_bottom<<"id\n";

      std::ofstream File_disp_diri_left_wc;
      File_disp_diri_left_wc.open(displ_bc_left_diri_wc_filename);
      // header
      File_disp_diri_left_wc<<"id\n";
      
      std::ofstream File_disp_diri_right_wc;
      File_disp_diri_right_wc.open(displ_bc_right_diri_wc_filename);
      // header
      File_disp_diri_right_wc<<"id\n";


      // discretize D union D_nonlocal
      // modify domain parameters 
      std::pair<std::vector<double>, std::vector<double>> nonlocal_domain;
      nonlocal_domain.first.push_back(0.0);
      nonlocal_domain.first.push_back(0.0);
      nonlocal_domain.second.push_back(0.0);
      nonlocal_domain.second.push_back(0.0);

      int eps_h = int(horizon/mesh_size);

      // modify x_left boundary
      if (double(eps_h)*mesh_size < horizon)
      	eps_h++;
      nonlocal_domain.first[0] = domain.first[0] - double(eps_h)*mesh_size;

      // modify domain's x end
      int ba_h = int( (domain.second[0] - domain.first[0])/mesh_size );
      double B = domain.first[0] + ba_h * mesh_size;
      if (B < domain.second[0])
      	domain.second[0] = domain.first[0] + ba_h * mesh_size;

      // modify domain's x_right boundary
      nonlocal_domain.second[0] = domain.second[0] + double(eps_h)*mesh_size;

      // modify y_bottom boundary
      nonlocal_domain.first[1] = domain.first[1] - double(eps_h)*mesh_size;

      // modify domain's y end
      ba_h = int( (domain.second[1] - domain.first[1])/mesh_size );
      B = domain.first[1] + ba_h * mesh_size;
      if (B < domain.second[1])
      	domain.second[1] = domain.first[1] + ba_h * mesh_size;

      // modify domain's x_right boundary
      nonlocal_domain.second[1] = domain.second[1] + double(eps_h)*mesh_size;

      // number of divisions in x direction 
      int nx = int((nonlocal_domain.second[0] - nonlocal_domain.first[0])/mesh_size);

      // number of divisions in y direction 
      int ny = int((nonlocal_domain.second[1] - nonlocal_domain.first[1])/mesh_size);

   	// total number of nodes
   	int total_num_nodes = (nx+1)*(ny+1);

   	// total number of elements
   	int total_num_elems = nx * ny;

   	// resize the nodes and elements data
   	std::vector<util::node> nodes;
   	nodes.resize(total_num_nodes);

   	std::vector<util::element> elements;
   	elements.resize(total_num_elems);

   	int node_counter = 0;
   	for (int j=0; j<=ny; j++) {
   		// j is along y direction
   		for(int i=0; i<=nx; i++) {
   			// i is along x direction

   			// node number
   			// int n = j * Nx + i;
   			int n = node_counter;
   			util::point x = util::point();
   			x.x = nonlocal_domain.first[0] + double(i) * mesh_size;
   			x.y = nonlocal_domain.first[1] + double(j) * mesh_size;

   			// create node data
   			util::node n_node = util::node(x.x,x.y);

   			n_node.n = n;

   			// for velocity and displacement
   			util::point v = util::point();
   			util::point u = util::point();

   			// set fixity mask of n_node and also add the id to
   			// list of ids on nonlocal boundaries
   			if (x.y <= domain.first[1]) {

   				// within bottom nonlocal boundary
   				n_node.fixity |= NONLOCAL_MASK;
   				n_node.fixity |= NONLOCAL_Y_BEGIN_MASK;

               if (x.y == domain.first[1]) {
                  n_node.fixity |= SURFACE_MASK;
                  n_node.fixity |= SURFACE_Y_BEGIN_MASK;                  
               }

      			// write to displ_bc_bottom_filename
      			File_disp_diri_bottom<<n<<"\n";

      			// compute initial condition
      			// if (shape_bottom != "Fixed") {
   				// 	computeIC(x, u, v, domain, test_flag, test_string,
   				// 		displ_ic_flag, displ_ic_params,
   				// 		vel_ic_flag, vel_ic_params);
   				// }
   			}
   			
            if (x.y >= domain.second[1]) {

   				// within top nonlocal boundary
   				n_node.fixity |= NONLOCAL_MASK;
   				n_node.fixity |= NONLOCAL_Y_END_MASK;

               if (x.y == domain.second[1]) {
                  n_node.fixity |= SURFACE_MASK;
                  n_node.fixity |= SURFACE_Y_END_MASK;
               }               

      			// write to displ_bc_bottom_filename
      			File_disp_diri_top<<n<<"\n";

      			// compute initial condition
      			// if (shape_top != "Fixed") {
   				// 	computeIC(x, u, v, domain, test_flag, test_string,
   				// 		displ_ic_flag, displ_ic_params,
   				// 		vel_ic_flag, vel_ic_params);
   				// }
   			}
   			
            if (x.x <= domain.first[0] and (x.y > domain.first[1] and x.y < domain.second[1])) {

   				// within left nonlocal boundary
   				n_node.fixity |= NONLOCAL_MASK;
   				n_node.fixity |= NONLOCAL_X_BEGIN_MASK;
               
               if (x.x == domain.first[0]) {
                  n_node.fixity |= SURFACE_MASK;
                  n_node.fixity |= SURFACE_X_BEGIN_MASK;
               }

      			// write to displ_bc_bottom_filename
      			File_disp_diri_left<<n<<"\n";

      			// compute initial condition
      			// if (shape_left != "Fixed") {
   				// 	computeIC(x, u, v, domain, test_flag, test_string,
   				// 		displ_ic_flag, displ_ic_params,
   				// 		vel_ic_flag, vel_ic_params);
   				// }
   			}
   			
            if (x.x >= domain.second[0] and (x.y > domain.first[1] and x.y < domain.second[1])) {
   				
   				// within right nonlocal boundary
   				n_node.fixity |= NONLOCAL_MASK;
   				n_node.fixity |= NONLOCAL_X_END_MASK;

               if (x.x == domain.second[0]) {
                  n_node.fixity |= SURFACE_MASK;
                  n_node.fixity |= SURFACE_X_END_MASK;                  
               }

      			// write to displ_bc_bottom_filename
      			File_disp_diri_right<<n<<"\n";

      			// compute initial condition
      			// if (shape_right != "Fixed") {
   				// 	computeIC(x, u, v, domain, test_flag, test_string,
   				// 		displ_ic_flag, displ_ic_params,
   				// 		vel_ic_flag, vel_ic_params);
   				// }
   			}

            if (x.x <= domain.first[0]) {

               // within left nonlocal boundary
               n_node.fixity |= NONLOCAL_MASK;
               n_node.fixity |= NONLOCAL_X_BEGIN_MASK;
               
               if (x.x == domain.first[0]) {
                  n_node.fixity |= SURFACE_MASK;
                  n_node.fixity |= SURFACE_X_BEGIN_MASK;
               }

               // write to displ_bc_bottom_filename
               File_disp_diri_left_wc<<n<<"\n";

               // compute initial condition
               // if (shape_left != "Fixed") {
               //    computeIC(x, u, v, domain, test_flag, test_string,
               //       displ_ic_flag, displ_ic_params,
               //       vel_ic_flag, vel_ic_params);
               // }
            }
            
            if (x.x >= domain.second[0]) {
               
               // within right nonlocal boundary
               n_node.fixity |= NONLOCAL_MASK;
               n_node.fixity |= NONLOCAL_X_END_MASK;

               if (x.x == domain.second[0]) {
                  n_node.fixity |= SURFACE_MASK;
                  n_node.fixity |= SURFACE_X_END_MASK;                  
               }

               // write to displ_bc_bottom_filename
               File_disp_diri_right_wc<<n<<"\n";

               // compute initial condition
               // if (shape_right != "Fixed") {
               //    computeIC(x, u, v, domain, test_flag, test_string,
               //       displ_ic_flag, displ_ic_params,
               //       vel_ic_flag, vel_ic_params);
               // }
            }

            // compute initial condition
  				computeIC(x, u, v, domain, test_flag,
   						displ_ic_flag, displ_ic_params,
   						vel_ic_flag, vel_ic_params);

   			// apply displacement to node
   			n_node.x = n_node.X + u;

   			// velocity
            n_node.v = v;

   			// add node to nodes vector
   			nodes[n] = n_node;

   			// increment the counter
   			node_counter++;
   		}
   	}

   	if (node_counter != total_num_nodes) {
   		std::cerr<<"Error: Final node counter = "<<node_counter<<" and total number of nodes = "<<total_num_nodes<<" not matching. \n";
   		std::cerr<<"Error: Check UniformTwoDSquareMesh().\n";
   		exit(1);
   	}	

   	int element_counter = 0;
   	// another loop for elements
   	for (int j=0; j<ny; j++) {
   		// j is along y direction
   		for(int i=0; i<nx; i++) {
   			// i is along x direction

   			// element number
   			// int n = j * Nx + i;
   			int n = element_counter;

   			// create element data
   			util::element n_elem = util::element(size_t(n));

   			// set type
   			n_elem.type = VTK_TYPE_QUAD;

   			// element node connectivity (put it in anti clockwise order)
   			// correct formula for node number
   			int n1 = j * (nx+1) + i;
   			int n2 = j * (nx+1) + i + 1;
   			int n3 = (j+1) * (nx+1) + i;
   			int n4 = (j+1) * (nx+1) + i + 1;

   			// update in element

   			// // anti clockwise ordering of nodes
   			n_elem.nodes.push_back(size_t(n1));
   			n_elem.nodes.push_back(size_t(n2));
   			n_elem.nodes.push_back(size_t(n4));
   			n_elem.nodes.push_back(size_t(n3));

   			// update in node data too
   			nodes[n1].elems.push_back(n);
   			nodes[n2].elems.push_back(n);
   			nodes[n3].elems.push_back(n);
   			nodes[n4].elems.push_back(n);

   			// write area to element data
   			n_elem.area = mesh_size * mesh_size;

   			// set fixity data of element
   			if (util::bit_match(nodes[n3].fixity, NONLOCAL_Y_BEGIN_MASK)) {

   				// top layer of element is inside bottom nonlocal boundary
   				// thus mark the element as nonlocal and also as
   				// nonlocal Y begin
   				n_elem.fixity |= NONLOCAL_MASK;
   				n_elem.fixity |= NONLOCAL_Y_BEGIN_MASK;
   			}
   			
   			if (util::bit_match(nodes[n2].fixity, NONLOCAL_Y_END_MASK)) {

   				// bottom layer of element is inside top nonlocal boundary
   				// thus mark the element as nonlocal and also as
   				// nonlocal Y end
   				n_elem.fixity |= NONLOCAL_MASK;
   				n_elem.fixity |= NONLOCAL_Y_END_MASK;
   			}
   			
   			double n2_x = nodes[n2].X[0];
   			if (n2_x <= domain.first[0]) {

   				// so the element is within left nonlocal boundary
   				double n2_y = nodes[n2].X[1];
   				double n3_y = nodes[n3].X[1];

   				if (n2_y >= domain.first[1] and n3_y <= domain.second[1]) {
   					// this element is inside left nonlocal boundary
   					n_elem.fixity |= NONLOCAL_MASK;
   					n_elem.fixity |= NONLOCAL_X_BEGIN_MASK;
   				}
   			}
   			
   			double n1_x = nodes[n1].X[0];
   			if (n1_x >= domain.second[0]) {

   				// so the element is within left nonlocal boundary
   				double n2_y = nodes[n2].X[1];
   				double n3_y = nodes[n3].X[1];

   				if (n2_y >= domain.first[1] and n3_y <= domain.second[1]) {
   					// this element is inside left nonlocal boundary
   					n_elem.fixity |= NONLOCAL_MASK;
   					n_elem.fixity |= NONLOCAL_X_END_MASK;
   				}
   			}

   			//
   			// n_elem data is complete so add it to elements vector
   			//
   			elements[n] = n_elem;

   			// increment the element counter
   			element_counter++;
   		}
   	}

   	// check if total number of elements and element_counter is same
   	if (element_counter != total_num_elems) {
   		std::cerr<<"Error: Final element counter = "<<element_counter<<" and total number of elementes = "<<total_num_elems<<" not matching. \n";
   		std::cerr<<"Error: Check UniformTwoDSquareMesh().\n";
   		exit(1);
   	}	

   	// check if total num of nodes calculated and size of nodes is same
   	if (nodes.size() != total_num_nodes) {
   		std::cerr<<"Error: Nodes data size = "<<nodes.size()<<" and total number of nodes = "<<total_num_nodes<<" not matching.\n";
   		std::cerr<<"Error: Check UniformTwoDSquareMesh().\n";
   		exit(1);
   	}

   	if (elements.size() != total_num_elems) {
   		std::cerr<<"Error: Elements data size = "<<elements.size()<<" and total number of elements = "<<total_num_elems<<" not matching.\n";
   		std::cerr<<"Error: Check UniformTwoDSquareMesh().\n";
   		exit(1);
   	}	

   	std::cout<<"Total number of nodes (Dirichlet Data) = "<<total_num_nodes<<"\n";

      // call VtkWriter to output the data
      IO::VtkWriter  writer = IO::VtkWriter();
      writer.open(path_diri,mesh_diri_filename,false);
      writer.appendMesh(&nodes, &elements);
      // writer.appendData("Displacement", displacement);
      // writer.appendData("Velocity", velocity);
      writer.addTimeStep(0.0);
      writer.close();

      // close csv files
      File_disp_diri_left.close();
      File_disp_diri_right.close();
      File_disp_diri_bottom.close();
      File_disp_diri_top.close();
      File_disp_diri_left_wc.close();
      File_disp_diri_right_wc.close();
   } // end of dirichlet data

   // 
   // get data corresponding to Neuman bc
   //
   // double tol = 1.0E-12;
   if (create_neuman_files == true) {

      //
      // open necessary files
      //
      std::ofstream File_disp_neum_left;
      File_disp_neum_left.open(displ_bc_left_neum_filename);
      // header
      File_disp_neum_left<<"id\n";
      
      std::ofstream File_disp_neum_right;
      File_disp_neum_right.open(displ_bc_right_neum_filename);
      // header
      File_disp_neum_right<<"id\n";

      std::ofstream File_disp_neum_top;
      File_disp_neum_top.open(displ_bc_top_neum_filename);
      // header
      File_disp_neum_top<<"id\n";
      
      std::ofstream File_disp_neum_bottom;
      File_disp_neum_bottom.open(displ_bc_bottom_neum_filename);
      // header
      File_disp_neum_bottom<<"id\n";

      std::ofstream File_disp_neum_left_wc;
      File_disp_neum_left_wc.open(displ_bc_left_neum_wc_filename);
      // header
      File_disp_neum_left_wc<<"id\n";
      
      std::ofstream File_disp_neum_right_wc;
      File_disp_neum_right_wc.open(displ_bc_right_neum_wc_filename);
      // header
      File_disp_neum_right_wc<<"id\n";

      //
      // modify boundary so that discretization is good
      //
      // modify domain's x end
      int ba_h = int( (domain.second[0] - domain.first[0])/mesh_size );
      double B = domain.first[0] + ba_h * mesh_size;
      if (B < domain.second[0])
         domain.second[0] = domain.first[0] + ba_h * mesh_size;

      // modify domain's y end
      ba_h = int( (domain.second[1] - domain.first[1])/mesh_size );
      B = domain.first[1] + ba_h * mesh_size;
      if (B < domain.second[1])
         domain.second[1] = domain.first[1] + ba_h * mesh_size;

      // number of divisions in x direction 
      int nx = int((domain.second[0] - domain.first[0])/mesh_size);

      // number of divisions in y direction 
      int ny = int((domain.second[1] - domain.first[1])/mesh_size);

      // total number of nodes
      int total_num_nodes = (nx+1)*(ny+1);

      // total number of elements
      int total_num_elems = nx * ny;

      // resize the nodes and elements data
      std::vector<util::node> nodes;
      nodes.resize(total_num_nodes);

      std::vector<util::element> elements;
      elements.resize(total_num_elems);

      int node_counter = 0;
      for (int j=0; j<=ny; j++) {
         // j is along y direction
         for(int i=0; i<=nx; i++) {
            // i is along x direction

            // node number
            // int n = j * Nx + i;
            int n = node_counter;
            util::point x = util::point();
            x.x = domain.first[0] + double(i) * mesh_size;
            x.y = domain.first[1] + double(j) * mesh_size;

            // create node data
            util::node n_node = util::node(x.x,x.y);

            n_node.n = n;

            // for velocity and displacement
            util::point v = util::point();
            util::point u = util::point();

            // set fixity mask of n_node and also add the id to
            // list of ids on nonlocal boundaries
            if (j == 0) {

               // within bottom nonlocal boundary
               n_node.fixity |= SURFACE_MASK;
               n_node.fixity |= SURFACE_Y_BEGIN_MASK;

               // write to displ_bc_bottom_filename
               File_disp_neum_bottom<<n<<"\n";

               // compute initial condition
               // if (shape_bottom != "Fixed") {
               //    computeIC(x, u, v, domain, test_flag, test_string,
               //       displ_ic_flag, displ_ic_params,
               //       vel_ic_flag, vel_ic_params);
               // }
            }

            if (j == ny) {

               // within top nonlocal boundary
               n_node.fixity |= SURFACE_MASK;
               n_node.fixity |= SURFACE_Y_END_MASK;

               // write to displ_bc_bottom_filename
               File_disp_neum_top<<n<<"\n";

               // compute initial condition
               // if (shape_top != "Fixed") {
               //    computeIC(x, u, v, domain, test_flag, test_string,
               //       displ_ic_flag, displ_ic_params,
               //       vel_ic_flag, vel_ic_params);
               // }
            }

            if (i == 0) {

               n_node.fixity |= SURFACE_MASK;
               n_node.fixity |= SURFACE_X_BEGIN_MASK;

               // we exclude corner points of domain in the list
               if (j !=0 and j != ny) {
                  // write to displ_bc_bottom_filename
                  File_disp_neum_left<<n<<"\n";
               }
            }

            if (i == nx) {

               n_node.fixity |= SURFACE_MASK;
               n_node.fixity |= SURFACE_X_END_MASK;

               // we exclude corner points of domain in the list
               if (j !=0 and j != ny) {               
                  // write to displ_bc_bottom_filename
                  File_disp_neum_right<<n<<"\n";
               }
            }

            if (i == 0) {

               n_node.fixity |= SURFACE_MASK;
               n_node.fixity |= SURFACE_X_BEGIN_MASK;

               // // we exclude corner points of domain in the list
               // if (j !=0 and j != ny) {
               //    // write to displ_bc_bottom_filename
               File_disp_neum_left_wc<<n<<"\n";
               // }
            }

            if (i == nx) {

               n_node.fixity |= SURFACE_MASK;
               n_node.fixity |= SURFACE_X_END_MASK;

               // // we exclude corner points of domain in the list
               // if (j !=0 and j != ny) {               
               //    // write to displ_bc_bottom_filename
               File_disp_neum_right_wc<<n<<"\n";
               // }
            }

            // compute initial condition
            computeIC(x, u, v, domain, test_flag,
                     displ_ic_flag, displ_ic_params,
                     vel_ic_flag, vel_ic_params);

            //
            // we now check for bending data
            //
            if (bd_dout == true) {
               checkForBendingDataFE(config, domain, mesh_size, node_counter, n_node.X);
            }    

            // apply displacement to node
            n_node.x = n_node.X + u;

            // velocity
            n_node.v = v;

            // add node to nodes vector
            nodes[n] = n_node;

            // increment the counter
            node_counter++;
         }
      }

      if (node_counter != total_num_nodes) {
         std::cerr<<"Error: Final node counter = "<<node_counter<<" and total number of nodes = "<<total_num_nodes<<" not matching. \n";
         std::cerr<<"Error: Check UniformTwoDSquareMesh().\n";
         exit(1);
      }  

      int element_counter = 0;
      // another loop for elements
      for (int j=0; j<ny; j++) {
         // j is along y direction
         for(int i=0; i<nx; i++) {
            // i is along x direction

            // element number
            // int n = j * Nx + i;
            int n = element_counter;

            // create element data
            util::element n_elem = util::element(size_t(n));

            // set type
            n_elem.type = VTK_TYPE_QUAD;

            // element node connectivity (put it in anti clockwise order)
            // correct formula for node number
            int n1 = j * (nx+1) + i;
            int n2 = j * (nx+1) + i + 1;
            int n3 = (j+1) * (nx+1) + i;
            int n4 = (j+1) * (nx+1) + i + 1;

            // update in element

            // // anti clockwise ordering of nodes
            n_elem.nodes.push_back(size_t(n1));
            n_elem.nodes.push_back(size_t(n2));
            n_elem.nodes.push_back(size_t(n4));
            n_elem.nodes.push_back(size_t(n3));

            // update in node data too
            nodes[n1].elems.push_back(n);
            nodes[n2].elems.push_back(n);
            nodes[n3].elems.push_back(n);
            nodes[n4].elems.push_back(n);

            // write area to element data
            n_elem.area = mesh_size * mesh_size;

            //
            // n_elem data is complete so add it to elements vector
            //
            elements[n] = n_elem;

            // increment the element counter
            element_counter++;
         }
      }

      // check if total number of elements and element_counter is same
      if (element_counter != total_num_elems) {
         std::cerr<<"Error: Final element counter = "<<element_counter<<" and total number of elementes = "<<total_num_elems<<" not matching. \n";
         std::cerr<<"Error: Check UniformTwoDSquareMesh().\n";
         exit(1);
      }  

      // check if total num of nodes calculated and size of nodes is same
      if (nodes.size() != total_num_nodes) {
         std::cerr<<"Error: Nodes data size = "<<nodes.size()<<" and total number of nodes = "<<total_num_nodes<<" not matching.\n";
         std::cerr<<"Error: Check UniformTwoDSquareMesh().\n";
         exit(1);
      }

      if (elements.size() != total_num_elems) {
         std::cerr<<"Error: Elements data size = "<<elements.size()<<" and total number of elements = "<<total_num_elems<<" not matching.\n";
         std::cerr<<"Error: Check UniformTwoDSquareMesh().\n";
         exit(1);
      }  

      //
      // write the bending data to file
      //
      if (bd_dout == true) {
         checkForBendingDataFE(config, domain, mesh_size, total_num_nodes, util::point(), true);
      } 

      std::cout<<"Total number of nodes (Neumann Data) = "<<total_num_nodes<<"\n";

      // call VtkWriter to output the data
      IO::VtkWriter  writer = IO::VtkWriter();
      writer.open(path_neum, mesh_neum_filename, false);
      writer.appendMesh(&nodes, &elements);
      // writer.appendData("Displacement", displacement);
      // writer.appendData("Velocity", velocity);
      writer.addTimeStep(0.0);
      writer.close();

      // close csv files
      File_disp_neum_left.close();
      File_disp_neum_right.close();
      File_disp_neum_bottom.close();
      File_disp_neum_top.close();
      File_disp_neum_left_wc.close();
      File_disp_neum_right_wc.close();
   } // end of neuman data

   return;
}

//
//	create uniform triangle mesh
//
void UniformTriangleMesh(YAML::Node config) {

   // //  local variables
   // YAML::Node config;
   // std::string config_filename = "mesh_fe_2d.yaml";

   std::pair<std::vector<double>, std::vector<double>> domain;
   double horizon;
   int ratio;
   double mesh_size;
   size_t displ_ic_flag;
   std::vector<double> displ_ic_params;
   size_t vel_ic_flag;
   std::vector<double> vel_ic_params;
   int center[3];

   bool test_flag;  // flag if input file is for testing the code
	
	std::string mesh_type;
   
   bool create_dirichlet_files = false;
   std::string path_diri;
   std::string mesh_diri_filename;
   std::string displ_bc_left_diri_filename;
   std::string displ_bc_right_diri_filename;
   std::string displ_bc_top_diri_filename;
   std::string displ_bc_bottom_diri_filename;
   std::string displ_bc_left_diri_wc_filename;
   std::string displ_bc_right_diri_wc_filename;

   bool create_neuman_files = false;
   std::string path_neum;
   std::string mesh_neum_filename;
   std::string displ_bc_left_neum_filename;
   std::string displ_bc_right_neum_filename;
   std::string displ_bc_top_neum_filename;
   std::string displ_bc_bottom_neum_filename;  
   std::string displ_bc_left_neum_wc_filename;
   std::string displ_bc_right_neum_wc_filename;  

   // bending loading related data (Implemented only for neuman type)
   bool bd_dout = false; 

   // read data from config file
   readDataFile(path_diri,
      mesh_diri_filename,
      displ_bc_left_diri_filename,
      displ_bc_right_diri_filename,
      displ_bc_top_diri_filename,
      displ_bc_bottom_diri_filename,
      displ_bc_left_diri_wc_filename,
      displ_bc_right_diri_wc_filename,
      path_neum,
      mesh_neum_filename,
      displ_bc_left_neum_filename,
      displ_bc_right_neum_filename,
      displ_bc_top_neum_filename,
      displ_bc_bottom_neum_filename,
      displ_bc_left_neum_wc_filename,
      displ_bc_right_neum_wc_filename,
      create_dirichlet_files,
      create_neuman_files,
      domain,
      horizon,
      ratio,
      mesh_size,
      mesh_type,
      displ_ic_flag,
      displ_ic_params,
      vel_ic_flag,
      vel_ic_params,
      test_flag,
      bd_dout,
      config);

   // 
   // first get data corresponding to Dirichlet bc
   //
   double tol = 1.0E-12;
   if (create_dirichlet_files == true) {

      //
      // open necessary files
      //
      std::ofstream File_disp_diri_left;
      File_disp_diri_left.open(displ_bc_left_diri_filename);
      // header
      File_disp_diri_left<<"id\n";
      
      std::ofstream File_disp_diri_right;
      File_disp_diri_right.open(displ_bc_right_diri_filename);
      // header
      File_disp_diri_right<<"id\n";

      std::ofstream File_disp_diri_top;
      File_disp_diri_top.open(displ_bc_top_diri_filename);
      // header
      File_disp_diri_top<<"id\n";
      
      std::ofstream File_disp_diri_bottom;
      File_disp_diri_bottom.open(displ_bc_bottom_diri_filename);
      // header
      File_disp_diri_bottom<<"id\n";

      std::ofstream File_disp_diri_left_wc;
      File_disp_diri_left_wc.open(displ_bc_left_diri_wc_filename);
      // header
      File_disp_diri_left_wc<<"id\n";
      
      std::ofstream File_disp_diri_right_wc;
      File_disp_diri_right_wc.open(displ_bc_right_diri_wc_filename);
      // header
      File_disp_diri_right_wc<<"id\n";


      // discretize D union D_nonlocal
      // modify domain parameters 
      std::pair<std::vector<double>, std::vector<double>> nonlocal_domain;
      nonlocal_domain.first.push_back(0.0);
      nonlocal_domain.first.push_back(0.0);
      nonlocal_domain.second.push_back(0.0);
      nonlocal_domain.second.push_back(0.0);

      int eps_h = int(horizon/mesh_size);

      // modify x_left boundary
      if (double(eps_h)*mesh_size < horizon)
      	eps_h++;
      nonlocal_domain.first[0] = domain.first[0] - double(eps_h)*mesh_size;

      // modify domain's x end
      int ba_h = int( (domain.second[0] - domain.first[0])/mesh_size );
      double B = domain.first[0] + ba_h * mesh_size;
      if (B < domain.second[0])
      	domain.second[0] = domain.first[0] + ba_h * mesh_size;

      // modify domain's x_right boundary
      nonlocal_domain.second[0] = domain.second[0] + double(eps_h)*mesh_size;

      // modify y_bottom boundary
      nonlocal_domain.first[1] = domain.first[1] - double(eps_h)*mesh_size;

      // modify domain's y end
      ba_h = int( (domain.second[1] - domain.first[1])/mesh_size );
      B = domain.first[1] + ba_h * mesh_size;
      if (B < domain.second[1])
      	domain.second[1] = domain.first[1] + ba_h * mesh_size;

      // modify domain's x_right boundary
      nonlocal_domain.second[1] = domain.second[1] + double(eps_h)*mesh_size;

      // number of divisions in x direction 
      int nx = int((nonlocal_domain.second[0] - nonlocal_domain.first[0])/mesh_size);

      // number of divisions in y direction 
      int ny = int((nonlocal_domain.second[1] - nonlocal_domain.first[1])/mesh_size);

   	// total number of nodes
   	int total_num_nodes = (nx+1)*(ny+1);

   	// total number of elements
   	int total_num_elems = 2 * nx * ny;

   	// resize the nodes and elements data
   	std::vector<util::node> nodes;
   	nodes.resize(total_num_nodes);

   	std::vector<util::element> elements;
   	elements.resize(total_num_elems);

   	int node_counter = 0;
   	for (int j=0; j<=ny; j++) {
   		// j is along y direction
   		for(int i=0; i<=nx; i++) {
   			// i is along x direction

   			// node number
   			// int n = j * Nx + i;
   			int n = node_counter;
   			util::point x = util::point();
   			x.x = nonlocal_domain.first[0] + double(i) * mesh_size;
   			x.y = nonlocal_domain.first[1] + double(j) * mesh_size;

   			// create node data
   			util::node n_node = util::node(x.x,x.y);

   			n_node.n = n;

   			// for velocity and displacement
   			util::point v = util::point();
   			util::point u = util::point();

   			// set fixity mask of n_node and also add the id to
   			// list of ids on nonlocal boundaries
   			if (x.y <= domain.first[1]) {

   				// within bottom nonlocal boundary
   				n_node.fixity |= NONLOCAL_MASK;
   				n_node.fixity |= NONLOCAL_Y_BEGIN_MASK;

               if (x.y == domain.first[1]) {
                  n_node.fixity |= SURFACE_MASK;
                  n_node.fixity |= SURFACE_Y_BEGIN_MASK;                  
               }

               // write to displ_bc_bottom_filename
               File_disp_diri_bottom<<n<<"\n";

      			// compute initial condition
      			// if (shape_bottom != "Fixed") {
   				// 	computeIC(x, u, v, domain, test_flag, test_string,
   				// 		displ_ic_flag, displ_ic_params,
   				// 		vel_ic_flag, vel_ic_params);
   				// }
   			}
   			
            if (x.y >= domain.second[1]) {

   				// within top nonlocal boundary
   				n_node.fixity |= NONLOCAL_MASK;
   				n_node.fixity |= NONLOCAL_Y_END_MASK;

               if (x.y == domain.second[1]) {
                  n_node.fixity |= SURFACE_MASK;
                  n_node.fixity |= SURFACE_Y_END_MASK;
               }               

               // write to displ_bc_bottom_filename
               File_disp_diri_top<<n<<"\n";

      			// compute initial condition
      			// if (shape_top != "Fixed") {
   				// 	computeIC(x, u, v, domain, test_flag, test_string,
   				// 		displ_ic_flag, displ_ic_params,
   				// 		vel_ic_flag, vel_ic_params);
   				// }
   			}
   			
            if (x.x <= domain.first[0] and (x.y > domain.first[1] and x.y < domain.second[1])) {

   				// within left nonlocal boundary
   				n_node.fixity |= NONLOCAL_MASK;
   				n_node.fixity |= NONLOCAL_X_BEGIN_MASK;

               if (x.x == domain.first[0]) {
                  n_node.fixity |= SURFACE_MASK;
                  n_node.fixity |= SURFACE_X_BEGIN_MASK;
               }

               // write to displ_bc_bottom_filename
               File_disp_diri_left<<n<<"\n";

      			// compute initial condition
      			// if (shape_left != "Fixed") {
   				// 	computeIC(x, u, v, domain, test_flag, test_string,
   				// 		displ_ic_flag, displ_ic_params,
   				// 		vel_ic_flag, vel_ic_params);
   				// }
   			}
   			
            if (x.x >= domain.second[0] and (x.y > domain.first[1] and x.y < domain.second[1])) {
   				
   				// within right nonlocal boundary
   				n_node.fixity |= NONLOCAL_MASK;
   				n_node.fixity |= NONLOCAL_X_END_MASK;

               if (x.x == domain.second[0]) {
                  n_node.fixity |= SURFACE_MASK;
                  n_node.fixity |= SURFACE_X_END_MASK;                  
               }

               // write to displ_bc_bottom_filename
               File_disp_diri_right<<n<<"\n";

      			// compute initial condition
      			// if (shape_right != "Fixed") {
   				// 	computeIC(x, u, v, domain, test_flag, test_string,
   				// 		displ_ic_flag, displ_ic_params,
   				// 		vel_ic_flag, vel_ic_params);
   				// }
   			}

            if (x.x <= domain.first[0]) {

               // within left nonlocal boundary
               n_node.fixity |= NONLOCAL_MASK;
               n_node.fixity |= NONLOCAL_X_BEGIN_MASK;

               if (x.x == domain.first[0]) {
                  n_node.fixity |= SURFACE_MASK;
                  n_node.fixity |= SURFACE_X_BEGIN_MASK;
               }

               // write to displ_bc_bottom_filename
               File_disp_diri_left_wc<<n<<"\n";

               // compute initial condition
               // if (shape_left != "Fixed") {
               //    computeIC(x, u, v, domain, test_flag, test_string,
               //       displ_ic_flag, displ_ic_params,
               //       vel_ic_flag, vel_ic_params);
               // }
            }
            
            if (x.x >= domain.second[0]) {
               
               // within right nonlocal boundary
               n_node.fixity |= NONLOCAL_MASK;
               n_node.fixity |= NONLOCAL_X_END_MASK;

               if (x.x == domain.second[0]) {
                  n_node.fixity |= SURFACE_MASK;
                  n_node.fixity |= SURFACE_X_END_MASK;                  
               }

               // write to displ_bc_bottom_filename
               File_disp_diri_right_wc<<n<<"\n";

               // compute initial condition
               // if (shape_right != "Fixed") {
               //    computeIC(x, u, v, domain, test_flag, test_string,
               //       displ_ic_flag, displ_ic_params,
               //       vel_ic_flag, vel_ic_params);
               // }
            }

            computeIC(x, u, v, domain, test_flag,
   						displ_ic_flag, displ_ic_params,
   						vel_ic_flag, vel_ic_params);

   			// apply displacement to node
   			n_node.x = n_node.X + u;

   			// velocity
            n_node.v = v;
            
   			// add node to nodes vector
   			nodes[n] = n_node;

   			// increment the counter
   			node_counter++;
   		}
   	}

   	if (node_counter != total_num_nodes) {
   		std::cerr<<"Error: Final node counter = "<<node_counter<<" and total number of nodes = "<<total_num_nodes<<" not matching. \n";
   		std::cerr<<"Error: Check UniformTwoDTriangleMesh().\n";
   		exit(1);
   	}	

   	int element_counter = 0; // see below for element type
   	// another loop for elements
   	for (int j=0; j<ny; j++) {
   		// j is along y direction
   		for(int i=0; i<nx; i++) {
   			// i is along x direction

   			// // in MeshTwoD.hpp the description of various types
   			// // is given. 
   			// // Currently we have implemented Type 1
   			// switch (element_type) {
   			// 	case 1: 
   			// 	{
   			// if i is even then we create elements with nodes n1,n2,n4,n3
   			// if i is odd then we create elements with nodes n2,n3,n5,n4

   			bool is_even = false;
   			if (i%2 == 0)
   				is_even = true;

   			// get node numbers
   			int n1 = j * (nx+1) +  i;
   			int n2 = j * (nx+1) + i + 1;
   			int n3 = (j+1) * (nx+1) + i;
   			int n4 = (j+1) * (nx+1) + i + 1;

   			// get element counter 
   			// T1 = j * 2 * Nx + 2*i
   			// T2 = j * 2 * Nx + 2*i + 1
   			int T1 = element_counter;
   			element_counter++;
   			int T2 = element_counter;
   			element_counter++;

   			// create triangle 1 and 2
   			util::element elem_1 = util::element(size_t(T1));
   			util::element elem_2 = util::element(size_t(T2));

   			// set type
   			elem_1.type = VTK_TYPE_TRIANGLE;
   			elem_2.type = VTK_TYPE_TRIANGLE;

   			// element node connectivity
   			if (is_even == true) {
   				// T1 (anticlockwise order)
   				elem_1.nodes.push_back(size_t(n1));
   				elem_1.nodes.push_back(size_t(n2));
   				elem_1.nodes.push_back(size_t(n4));

   				// T2 (anticlockwise order)
   				elem_2.nodes.push_back(size_t(n1));
   				elem_2.nodes.push_back(size_t(n4));
   				elem_2.nodes.push_back(size_t(n3));

   				// also put the element number in node data
   				nodes[n1].elems.push_back(T1);
   				nodes[n2].elems.push_back(T1);
   				nodes[n4].elems.push_back(T1);

   				nodes[n1].elems.push_back(T2);
   				nodes[n4].elems.push_back(T2);
   				nodes[n3].elems.push_back(T2);

   				if ( 	util::bit_match(nodes[n1].fixity, NONLOCAL_MASK) and 
   						util::bit_match(nodes[n2].fixity, NONLOCAL_MASK) and 
   						util::bit_match(nodes[n4].fixity, NONLOCAL_MASK) ) {

   					elem_1.fixity |= NONLOCAL_MASK;

   					double n1_x = nodes[n1].X[0];
   					double n2_x = nodes[n2].X[0];
   					double n2_y = nodes[n2].X[1];
   					double n4_y = nodes[n4].X[1];

   					if (n4_y <= domain.first[1]) 
   						elem_1.fixity |= NONLOCAL_Y_BEGIN_MASK;

   					if (n2_y >= domain.second[1])
   						elem_1.fixity |= NONLOCAL_Y_END_MASK;

   					if (n2_x <= domain.first[0])
   						if (n2_y >= domain.first[1] and n4_y <= domain.second[1])
   							elem_1.fixity |= NONLOCAL_X_BEGIN_MASK;

   					if (n1_x >= domain.second[0])
   						if (n2_y >= domain.first[1] and n4_y <= domain.second[1])
   							elem_1.fixity |= NONLOCAL_X_END_MASK;
   				}

   				if ( 	util::bit_match(nodes[n1].fixity, NONLOCAL_MASK) and 
   						util::bit_match(nodes[n4].fixity, NONLOCAL_MASK) and 
   						util::bit_match(nodes[n3].fixity, NONLOCAL_MASK) ) {

   					elem_2.fixity |= NONLOCAL_MASK;

   					double n1_x = nodes[n1].X[0];
   					double n1_y = nodes[n1].X[1];
   					double n4_x = nodes[n4].X[0];
   					double n4_y = nodes[n4].X[1];

   					if (n4_y <= domain.first[1]) 
   						elem_2.fixity |= NONLOCAL_Y_BEGIN_MASK;

   					if (n1_y >= domain.second[1])
   						elem_2.fixity |= NONLOCAL_Y_END_MASK;

   					if (n4_x <= domain.first[0])
   						if (n1_y >= domain.first[1] and n4_y <= domain.second[1])
   							elem_2.fixity |= NONLOCAL_X_BEGIN_MASK;

   					if (n1_x >= domain.second[0])
   						if (n1_y >= domain.first[1] and n4_y <= domain.second[1])
   							elem_2.fixity |= NONLOCAL_X_END_MASK;
   				}
   			}
   			else {
   				// T1 (anticlockwise order)
   				elem_1.nodes.push_back(size_t(n1));
   				elem_1.nodes.push_back(size_t(n2));
   				elem_1.nodes.push_back(size_t(n3));

   				// T2 (anticlockwise order)
   				elem_2.nodes.push_back(size_t(n2));
   				elem_2.nodes.push_back(size_t(n4));
   				elem_2.nodes.push_back(size_t(n3));

   				// also put the element number in node data
   				nodes[n1].elems.push_back(T1);
   				nodes[n2].elems.push_back(T1);
   				nodes[n3].elems.push_back(T1);

   				nodes[n2].elems.push_back(T2);
   				nodes[n4].elems.push_back(T2);
   				nodes[n3].elems.push_back(T2);

   				if ( 	util::bit_match(nodes[n1].fixity, NONLOCAL_MASK) and 
   						util::bit_match(nodes[n2].fixity, NONLOCAL_MASK) and 
   						util::bit_match(nodes[n3].fixity, NONLOCAL_MASK) ) {

   					elem_1.fixity |= NONLOCAL_MASK;

   					double n1_x = nodes[n1].X[0];
   					double n2_x = nodes[n2].X[0];
   					double n2_y = nodes[n2].X[1];
   					double n3_y = nodes[n3].X[1];

   					if (n3_y <= domain.first[1]) 
   						elem_1.fixity |= NONLOCAL_Y_BEGIN_MASK;

   					if (n2_y >= domain.second[1])
   						elem_1.fixity |= NONLOCAL_Y_END_MASK;

   					if (n2_x <= domain.first[0])
   						if (n2_y >= domain.first[1] and n3_y <= domain.second[1])
   							elem_1.fixity |= NONLOCAL_X_BEGIN_MASK;

   					if (n1_x >= domain.second[0])
   						if (n2_y >= domain.first[1] and n3_y <= domain.second[1])
   							elem_1.fixity |= NONLOCAL_X_END_MASK;
   				}

   				if ( 	util::bit_match(nodes[n2].fixity, NONLOCAL_MASK) and 
   						util::bit_match(nodes[n4].fixity, NONLOCAL_MASK) and 
   						util::bit_match(nodes[n3].fixity, NONLOCAL_MASK) ) {

   					elem_2.fixity |= NONLOCAL_MASK;

   					double n2_x = nodes[n2].X[0];
   					double n2_y = nodes[n2].X[1];
   					double n3_x = nodes[n3].X[0];
   					double n3_y = nodes[n3].X[1];

   					if (n3_y <= domain.first[1]) 
   						elem_2.fixity |= NONLOCAL_Y_BEGIN_MASK;

   					if (n2_y >= domain.second[1])
   						elem_2.fixity |= NONLOCAL_Y_END_MASK;

   					if (n2_x <= domain.first[0])
   						if (n2_y >= domain.first[1] and n3_y <= domain.second[1])
   							elem_2.fixity |= NONLOCAL_X_BEGIN_MASK;

   					if (n3_x >= domain.second[0])
   						if (n2_y >= domain.first[1] and n3_y <= domain.second[1])
   							elem_2.fixity |= NONLOCAL_X_END_MASK;
   				}
   			}

   			// write area to element data
   			// can also call util::methods::triangleArea()
   			elem_1.area = 0.5*mesh_size*mesh_size; 
   			elem_2.area = elem_1.area;

   			// n_elem data is complete so add it to elements vector
   			elements[T1] = elem_1;
   			elements[T2] = elem_2;

   			// Note element counter has already been incremented above

   			// 	} // case 1 bracket
   			// 	break;
   			// } // switch bracket
   		} // loop over i
   	} // loop over j

   	// check if total number of elements and element_counter is same
   	if (element_counter != total_num_elems) {
   		std::cerr<<"Error: Final element counter = "<<element_counter<<" and total number of elementes = "<<total_num_elems<<" not matching. \n";
   		std::cerr<<"Error: Check UniformTwoDTriangleMesh().\n";
   		exit(1);
   	}	

   	// check if total num of nodes calculated and size of nodes is same
   	if (nodes.size() != total_num_nodes) {
   		std::cerr<<"Error: Nodes data size = "<<nodes.size()<<" and total number of nodes = "<<total_num_nodes<<" not matching.\n";
   		std::cerr<<"Error: Check UniformTwoDTriangleMesh().\n";
   		exit(1);
   	}

   	if (elements.size() != total_num_elems) {
   		std::cerr<<"Error: Elements data size = "<<elements.size()<<" and total number of elements = "<<total_num_elems<<" not matching.\n";
   		std::cerr<<"Error: Check UniformTwoDTriangleMesh().\n";
   		exit(1);
   	}	

   	std::cout<<"Total number of nodes (Dirichlet Data) = "<<total_num_nodes<<"\n";

      // call VtkWriter to output the data
      IO::VtkWriter  writer = IO::VtkWriter();
      writer.open(path_diri,mesh_diri_filename, false);
      writer.appendMesh(&nodes, &elements);
      // writer.appendData("Displacement", displacement);
      // writer.appendData("Velocity", velocity);
      writer.addTimeStep(0.0);
      writer.close();

      // close csv files
      File_disp_diri_left.close();
      File_disp_diri_right.close();
      File_disp_diri_bottom.close();
      File_disp_diri_top.close();
      File_disp_diri_left_wc.close();
      File_disp_diri_right_wc.close();
   } // end of dirichlet data

   // 
   // get data corresponding to Neuman bc
   //
   if (create_neuman_files == true) {

      //
      // open necessary files
      //
      std::ofstream File_disp_neum_left;
      File_disp_neum_left.open(displ_bc_left_neum_filename);
      // header
      File_disp_neum_left<<"id\n";
      
      std::ofstream File_disp_neum_right;
      File_disp_neum_right.open(displ_bc_right_neum_filename);
      // header
      File_disp_neum_right<<"id\n";

      std::ofstream File_disp_neum_top;
      File_disp_neum_top.open(displ_bc_top_neum_filename);
      // header
      File_disp_neum_top<<"id\n";
      
      std::ofstream File_disp_neum_bottom;
      File_disp_neum_bottom.open(displ_bc_bottom_neum_filename);
      // header
      File_disp_neum_bottom<<"id\n";

      std::ofstream File_disp_neum_left_wc;
      File_disp_neum_left_wc.open(displ_bc_left_neum_wc_filename);
      // header
      File_disp_neum_left_wc<<"id\n";
      
      std::ofstream File_disp_neum_right_wc;
      File_disp_neum_right_wc.open(displ_bc_right_neum_wc_filename);
      // header
      File_disp_neum_right_wc<<"id\n";


      //
      // modify boundary so that discretization is good
      //
      // modify domain's x end
      int ba_h = int( (domain.second[0] - domain.first[0])/mesh_size );
      double B = domain.first[0] + ba_h * mesh_size;
      if (B < domain.second[0])
         domain.second[0] = domain.first[0] + ba_h * mesh_size;

      // modify domain's y end
      ba_h = int( (domain.second[1] - domain.first[1])/mesh_size );
      B = domain.first[1] + ba_h * mesh_size;
      if (B < domain.second[1])
         domain.second[1] = domain.first[1] + ba_h * mesh_size;

      // number of divisions in x direction 
      int nx = int((domain.second[0] - domain.first[0])/mesh_size);

      // number of divisions in y direction 
      int ny = int((domain.second[1] - domain.first[1])/mesh_size);

      // total number of nodes
      int total_num_nodes = (nx+1)*(ny+1);

      // total number of elements
      int total_num_elems = 2 * nx * ny;

      // resize the nodes and elements data
      std::vector<util::node> nodes;
      nodes.resize(total_num_nodes);

      std::vector<util::element> elements;
      elements.resize(total_num_elems);

      int node_counter = 0;
      for (int j=0; j<=ny; j++) {
         // j is along y direction
         for(int i=0; i<=nx; i++) {
            // i is along x direction

            // node number
            // int n = j * Nx + i;
            int n = node_counter;
            util::point x = util::point();
            x.x = domain.first[0] + double(i) * mesh_size;
            x.y = domain.first[1] + double(j) * mesh_size;

            // create node data
            util::node n_node = util::node(x.x,x.y);

            n_node.n = n;

            // for velocity and displacement
            util::point v = util::point();
            util::point u = util::point();

            // set fixity mask of n_node and also add the id to
            // list of ids on nonlocal boundaries
            if (j == 0) {

               // within bottom nonlocal boundary
               n_node.fixity |= SURFACE_MASK;
               n_node.fixity |= SURFACE_Y_BEGIN_MASK;

               // write to displ_bc_bottom_filename
               File_disp_neum_bottom<<n<<"\n";

               // compute initial condition
               // if (shape_bottom != "Fixed") {
               //    computeIC(x, u, v, domain, test_flag, test_string,
               //       displ_ic_flag, displ_ic_params,
               //       vel_ic_flag, vel_ic_params);
               // }
            }

            if (j == ny) {

               // within top nonlocal boundary
               n_node.fixity |= SURFACE_MASK;
               n_node.fixity |= SURFACE_Y_END_MASK;

               // write to displ_bc_bottom_filename
               File_disp_neum_top<<n<<"\n";

               // compute initial condition
               // if (shape_top != "Fixed") {
               //    computeIC(x, u, v, domain, test_flag, test_string,
               //       displ_ic_flag, displ_ic_params,
               //       vel_ic_flag, vel_ic_params);
               // }
            }

            if (i == 0) {

               n_node.fixity |= SURFACE_MASK;
               n_node.fixity |= SURFACE_X_BEGIN_MASK;

               // we exclude corner points of domain in the list
               if (j !=0 and j != ny) {
                  // write to displ_bc_bottom_filename
                  File_disp_neum_left<<n<<"\n";
               }
            }

            if (i == nx) {

               n_node.fixity |= SURFACE_MASK;
               n_node.fixity |= SURFACE_X_END_MASK;

               // we exclude corner points of domain in the list
               if (j !=0 and j != ny) {               
                  // write to displ_bc_bottom_filename
                  File_disp_neum_right<<n<<"\n";
               }
            }

            if (i == 0) {

               n_node.fixity |= SURFACE_MASK;
               n_node.fixity |= SURFACE_X_BEGIN_MASK;

               // // we exclude corner points of domain in the list
               // if (j !=0 and j != ny) {
               //    // write to displ_bc_bottom_filename
               File_disp_neum_left_wc<<n<<"\n";
               // }
            }

            if (i == nx) {

               n_node.fixity |= SURFACE_MASK;
               n_node.fixity |= SURFACE_X_END_MASK;

               // // we exclude corner points of domain in the list
               // if (j !=0 and j != ny) {               
               //    // write to displ_bc_bottom_filename
               File_disp_neum_right_wc<<n<<"\n";
               // }
            }

            computeIC(x, u, v, domain, test_flag,
                     displ_ic_flag, displ_ic_params,
                     vel_ic_flag, vel_ic_params);

            //
            // we now check for bending data
            //
            if (bd_dout == true) {
               checkForBendingDataFE(config, domain, mesh_size, node_counter, n_node.X);
            }

            // apply displacement to node
            n_node.x = n_node.X + u;

            // velocity
            n_node.v = v;
            
            // add node to nodes vector
            nodes[n] = n_node;

            // increment the counter
            node_counter++;
         }
      }

      if (node_counter != total_num_nodes) {
         std::cerr<<"Error: Final node counter = "<<node_counter<<" and total number of nodes = "<<total_num_nodes<<" not matching. \n";
         std::cerr<<"Error: Check UniformTwoDTriangleMesh().\n";
         exit(1);
      }  

      int element_counter = 0; // see below for element type
      // another loop for elements
      for (int j=0; j<ny; j++) {
         // j is along y direction
         for(int i=0; i<nx; i++) {
            // i is along x direction

            // // in MeshTwoD.hpp the description of various types
            // // is given. 
            // // Currently we have implemented Type 1
            // switch (element_type) {
            //    case 1: 
            //    {
            // if i is even then we create elements with nodes n1,n2,n4,n3
            // if i is odd then we create elements with nodes n2,n3,n5,n4

            bool is_even = false;
            if (i%2 == 0)
               is_even = true;

            // get node numbers
            int n1 = j * (nx+1) +  i;
            int n2 = j * (nx+1) + i + 1;
            int n3 = (j+1) * (nx+1) + i;
            int n4 = (j+1) * (nx+1) + i + 1;

            // get element counter 
            // T1 = j * 2 * Nx + 2*i
            // T2 = j * 2 * Nx + 2*i + 1
            int T1 = element_counter;
            element_counter++;
            int T2 = element_counter;
            element_counter++;

            // create triangle 1 and 2
            util::element elem_1 = util::element(size_t(T1));
            util::element elem_2 = util::element(size_t(T2));

            // set type
            elem_1.type = VTK_TYPE_TRIANGLE;
            elem_2.type = VTK_TYPE_TRIANGLE;

            // element node connectivity
            if (is_even == true) {
               // T1 (anticlockwise order)
               elem_1.nodes.push_back(size_t(n1));
               elem_1.nodes.push_back(size_t(n2));
               elem_1.nodes.push_back(size_t(n4));

               // T2 (anticlockwise order)
               elem_2.nodes.push_back(size_t(n1));
               elem_2.nodes.push_back(size_t(n4));
               elem_2.nodes.push_back(size_t(n3));

               // also put the element number in node data
               nodes[n1].elems.push_back(T1);
               nodes[n2].elems.push_back(T1);
               nodes[n4].elems.push_back(T1);

               nodes[n1].elems.push_back(T2);
               nodes[n4].elems.push_back(T2);
               nodes[n3].elems.push_back(T2);
            }
            else {
               // T1 (anticlockwise order)
               elem_1.nodes.push_back(size_t(n1));
               elem_1.nodes.push_back(size_t(n2));
               elem_1.nodes.push_back(size_t(n3));

               // T2 (anticlockwise order)
               elem_2.nodes.push_back(size_t(n2));
               elem_2.nodes.push_back(size_t(n4));
               elem_2.nodes.push_back(size_t(n3));

               // also put the element number in node data
               nodes[n1].elems.push_back(T1);
               nodes[n2].elems.push_back(T1);
               nodes[n3].elems.push_back(T1);

               nodes[n2].elems.push_back(T2);
               nodes[n4].elems.push_back(T2);
               nodes[n3].elems.push_back(T2);
            }

            // write area to element data
            // can also call util::methods::triangleArea()
            elem_1.area = 0.5*mesh_size*mesh_size; 
            elem_2.area = elem_1.area;

            // n_elem data is complete so add it to elements vector
            elements[T1] = elem_1;
            elements[T2] = elem_2;

            // Note element counter has already been incremented above

            //    } // case 1 bracket
            //    break;
            // } // switch bracket
         } // loop over i
      } // loop over j

      // check if total number of elements and element_counter is same
      if (element_counter != total_num_elems) {
         std::cerr<<"Error: Final element counter = "<<element_counter<<" and total number of elementes = "<<total_num_elems<<" not matching. \n";
         std::cerr<<"Error: Check UniformTwoDTriangleMesh().\n";
         exit(1);
      }  

      // check if total num of nodes calculated and size of nodes is same
      if (nodes.size() != total_num_nodes) {
         std::cerr<<"Error: Nodes data size = "<<nodes.size()<<" and total number of nodes = "<<total_num_nodes<<" not matching.\n";
         std::cerr<<"Error: Check UniformTwoDTriangleMesh().\n";
         exit(1);
      }

      if (elements.size() != total_num_elems) {
         std::cerr<<"Error: Elements data size = "<<elements.size()<<" and total number of elements = "<<total_num_elems<<" not matching.\n";
         std::cerr<<"Error: Check UniformTwoDTriangleMesh().\n";
         exit(1);
      }  

      //
      // write the bending data to file
      //
      if (bd_dout == true) {
         checkForBendingDataFE(config, domain, mesh_size, total_num_nodes, util::point(), true);
      } 

      std::cout<<"Total number of nodes (Neumann Data) = "<<total_num_nodes<<"\n";

      // call VtkWriter to output the data
      IO::VtkWriter  writer = IO::VtkWriter();
      writer.open(path_neum,mesh_neum_filename, false);
      writer.appendMesh(&nodes, &elements);
      // writer.appendData("Displacement", displacement);
      // writer.appendData("Velocity", velocity);
      writer.addTimeStep(0.0);
      writer.close();

      // close csv files
      File_disp_neum_left.close();
      File_disp_neum_right.close();
      File_disp_neum_bottom.close();
      File_disp_neum_top.close();
      File_disp_neum_left_wc.close();
      File_disp_neum_right_wc.close();
   } // end of neuman data   

   return;	
}

//
// create uniform triangle mesh
//
// Mesh:
//
//         o---------o---------o---------o---------o
//         | \     / | \     / | \     / | \     / |
//         |  \  /   |  \  /   |  \  /   |  \  /   |
//         |    o    |    o    |    o    |    o    |
//         |  /  \   |  /  \   |  /  \   |  /  \   |
//         | /     \ | /     \ | /     \ | /     \ |
//         o---------o---------o---------o---------o
//         | \     / | \     / | \     / | \     / |
//         |  \  /   |  \  /   |  \  /   |  \  /   |
//         |    o    |    o    |    o    |    o    |
//         |  /  \   |  /  \   |  /  \   |  /  \   |
//         | /     \ | /     \ | /     \ | /     \ |
//         o---------o---------o---------o---------o
//         | \     / | \     / | \     / | \     / |
//         |  \  /   |  \  /   |  \  /   |  \  /   |
//         |    o    |    o    |    o    |    o    |
//         |  /12 \  |  /  \   |  /  \   |  /  \   |
//       11| /     \ | /     \ | /     \ | /     \ |
//     --- o---------o---------o---------o---------o
//      |  | \     / | \     / | \     / | \     / |
//      |  |  \  /   |  \  /   |  \  /   |  \  /   |
//     2h  |    o    |    o    |    o    |    o    |
//      |  |  / 2\   |  / 5\   |  / 7\   |  / 9\   |
//      |  | /     \ | /     \ | /     \ | /     \ |
//     --- o---------o---------o---------o---------o
//         1         3         6         8         10
//
//         |--- 2h --|
//
//
void UniformTriangleSymMesh(YAML::Node config) {

   // //  local variables
   // YAML::Node config;
   // std::string config_filename = "mesh_fe_2d.yaml";

   std::pair<std::vector<double>, std::vector<double>> domain;
   double horizon;
   int ratio;
   double mesh_size;
   size_t displ_ic_flag;
   std::vector<double> displ_ic_params;
   size_t vel_ic_flag;
   std::vector<double> vel_ic_params;
   int center[3];

   bool test_flag;  // flag if input file is for testing the code
   
   std::string mesh_type;
   
   bool create_dirichlet_files = false;
   std::string path_diri;
   std::string mesh_diri_filename;
   std::string displ_bc_left_diri_filename;
   std::string displ_bc_right_diri_filename;
   std::string displ_bc_top_diri_filename;
   std::string displ_bc_bottom_diri_filename;
   std::string displ_bc_left_diri_wc_filename;
   std::string displ_bc_right_diri_wc_filename;

   bool create_neuman_files = false;
   std::string path_neum;
   std::string mesh_neum_filename;
   std::string displ_bc_left_neum_filename;
   std::string displ_bc_right_neum_filename;
   std::string displ_bc_top_neum_filename;
   std::string displ_bc_bottom_neum_filename;  
   std::string displ_bc_left_neum_wc_filename;
   std::string displ_bc_right_neum_wc_filename;  

   // bending loading related data (Implemented only for neuman type)
   bool bd_dout = false; 

   // read data from config file
   readDataFile(path_diri,
      mesh_diri_filename,
      displ_bc_left_diri_filename,
      displ_bc_right_diri_filename,
      displ_bc_top_diri_filename,
      displ_bc_bottom_diri_filename,
      displ_bc_left_diri_wc_filename,
      displ_bc_right_diri_wc_filename,
      path_neum,
      mesh_neum_filename,
      displ_bc_left_neum_filename,
      displ_bc_right_neum_filename,
      displ_bc_top_neum_filename,
      displ_bc_bottom_neum_filename,
      displ_bc_left_neum_wc_filename,
      displ_bc_right_neum_wc_filename,
      create_dirichlet_files,
      create_neuman_files,
      domain,
      horizon,
      ratio,
      mesh_size,
      mesh_type,
      displ_ic_flag,
      displ_ic_params,
      vel_ic_flag,
      vel_ic_params,
      test_flag,
      bd_dout,
      config);


   double tol = 1.0E-12;
   // 
   // get data corresponding to Neuman bc
   //
   if (create_neuman_files == true) {

      //
      // open necessary files
      //
      std::ofstream File_disp_neum_left;
      File_disp_neum_left.open(displ_bc_left_neum_filename);
      // header
      File_disp_neum_left<<"id\n";
      
      std::ofstream File_disp_neum_right;
      File_disp_neum_right.open(displ_bc_right_neum_filename);
      // header
      File_disp_neum_right<<"id\n";

      std::ofstream File_disp_neum_top;
      File_disp_neum_top.open(displ_bc_top_neum_filename);
      // header
      File_disp_neum_top<<"id\n";
      
      std::ofstream File_disp_neum_bottom;
      File_disp_neum_bottom.open(displ_bc_bottom_neum_filename);
      // header
      File_disp_neum_bottom<<"id\n";

      std::ofstream File_disp_neum_left_wc;
      File_disp_neum_left_wc.open(displ_bc_left_neum_wc_filename);
      // header
      File_disp_neum_left_wc<<"id\n";
      
      std::ofstream File_disp_neum_right_wc;
      File_disp_neum_right_wc.open(displ_bc_right_neum_wc_filename);
      // header
      File_disp_neum_right_wc<<"id\n";


      //
      // modify boundary so that discretization is good
      //
      //
      // Here we want rectangle to be multiple of 2h
      //
      double mesh_size_grid = 2.0 * mesh_size;

      // modify domain's x end
      int ba_h = int( (domain.second[0] - domain.first[0])/(mesh_size_grid) );
      double B = domain.first[0] + ba_h * mesh_size_grid;
      if (B < domain.second[0])
         domain.second[0] = domain.first[0] + ba_h * mesh_size_grid;

      // modify domain's y end
      ba_h = int( (domain.second[1] - domain.first[1])/mesh_size_grid );
      B = domain.first[1] + ba_h * mesh_size_grid;
      if (B < domain.second[1])
         domain.second[1] = domain.first[1] + ba_h * mesh_size_grid;

      // number of divisions in x direction 
      int nx_grid = int((domain.second[0] - domain.first[0])/mesh_size_grid);

      // number of divisions in y direction 
      int ny_grid = int((domain.second[1] - domain.first[1])/mesh_size_grid);

      // total number of nodes
      int total_num_nodes = (nx_grid+1)*(ny_grid+1) + nx_grid*ny_grid;

      // total number of elements
      int total_num_elems = 4 * nx_grid * ny_grid;

      // resize the nodes and elements data
      std::vector<util::node> nodes;
      nodes.resize(total_num_nodes);

      std::vector<util::element> elements;
      elements.resize(total_num_elems);

      int node_counter = 0;
      for (int j=0; j<=ny_grid; j++) {
         // j is along y direction
         for(int i=0; i<=nx_grid; i++) {
            // i is along x direction

            //
            // first node (grid point)
            //
            int n = node_counter;
            util::point x = util::point();
            x.x = domain.first[0] + double(i) * mesh_size_grid;
            x.y = domain.first[1] + double(j) * mesh_size_grid;

            // create node data
            util::node n_node = util::node(x.x,x.y);

            n_node.n = n;

            // for velocity and displacement
            util::point v = util::point();
            util::point u = util::point();

            // set fixity mask of n_node and also add the id to
            // list of ids on nonlocal boundaries
            if (j == 0) {

               // within bottom nonlocal boundary
               n_node.fixity |= SURFACE_MASK;
               n_node.fixity |= SURFACE_Y_BEGIN_MASK;

               // write to displ_bc_bottom_filename
               File_disp_neum_bottom<<n<<"\n";

               // compute initial condition
               // if (shape_bottom != "Fixed") {
               //    computeIC(x, u, v, domain, test_flag, test_string,
               //       displ_ic_flag, displ_ic_params,
               //       vel_ic_flag, vel_ic_params);
               // }
            }

            if (j == ny_grid) {

               // within top nonlocal boundary
               n_node.fixity |= SURFACE_MASK;
               n_node.fixity |= SURFACE_Y_END_MASK;

               // write to displ_bc_bottom_filename
               File_disp_neum_top<<n<<"\n";

               // compute initial condition
               // if (shape_top != "Fixed") {
               //    computeIC(x, u, v, domain, test_flag, test_string,
               //       displ_ic_flag, displ_ic_params,
               //       vel_ic_flag, vel_ic_params);
               // }
            }

            if (i == 0) {

               n_node.fixity |= SURFACE_MASK;
               n_node.fixity |= SURFACE_X_BEGIN_MASK;

               // we exclude corner points of domain in the list
               if (j !=0 and j != ny_grid) {
                  // write to displ_bc_bottom_filename
                  File_disp_neum_left<<n<<"\n";
               }
            }

            if (i == nx_grid) {

               n_node.fixity |= SURFACE_MASK;
               n_node.fixity |= SURFACE_X_END_MASK;

               // we exclude corner points of domain in the list
               if (j !=0 and j != ny_grid) {               
                  // write to displ_bc_bottom_filename
                  File_disp_neum_right<<n<<"\n";
               }
            }

            if (i == 0) {

               n_node.fixity |= SURFACE_MASK;
               n_node.fixity |= SURFACE_X_BEGIN_MASK;

               // // we exclude corner points of domain in the list
               // if (j !=0 and j != ny) {
               //    // write to displ_bc_bottom_filename
               File_disp_neum_left_wc<<n<<"\n";
               // }
            }

            if (i == nx_grid) {

               n_node.fixity |= SURFACE_MASK;
               n_node.fixity |= SURFACE_X_END_MASK;

               // // we exclude corner points of domain in the list
               // if (j !=0 and j != ny) {               
               //    // write to displ_bc_bottom_filename
               File_disp_neum_right_wc<<n<<"\n";
               // }
            }

            computeIC(x, u, v, domain, test_flag,
                     displ_ic_flag, displ_ic_params,
                     vel_ic_flag, vel_ic_params);

            //
            // we now check for bending data
            //
            if (bd_dout == true) {
               checkForBendingDataFE(config, domain, mesh_size, node_counter, n_node.X);
            }

            // apply displacement to node
            n_node.x = n_node.X + u;

            // velocity
            n_node.v = v;
            
            // add node to nodes vector
            nodes[n] = n_node;

            // increment the counter
            node_counter++;

            //
            // second node (center of grid)
            //
            bool create_second_node = true;

            if (i == nx_grid)
               create_second_node = false;

            if (j == ny_grid)
               create_second_node = false; 

            if (create_second_node == false)
               continue;

            n = node_counter;
            x = util::point();
            x.x = domain.first[0] + double(i) * mesh_size_grid + 0.5 * mesh_size_grid;
            x.y = domain.first[1] + double(j) * mesh_size_grid + 0.5 * mesh_size_grid;

            // create node data
            util::node n_node_2 = util::node(x.x,x.y);

            n_node_2.n = n;

            // for velocity and displacement
            v = util::point();
            u = util::point();

            computeIC(x, u, v, domain, test_flag,
                     displ_ic_flag, displ_ic_params,
                     vel_ic_flag, vel_ic_params);

            //
            // we now check for bending data
            //
            if (bd_dout == true) {
               checkForBendingDataFE(config, domain, mesh_size, node_counter, n_node_2.X);
            }

            // apply displacement to node
            n_node_2.x = n_node_2.X + u;

            // velocity
            n_node_2.v = v;
            
            // add node to nodes vector
            nodes[n] = n_node_2;

            // increment the counter
            node_counter++;
         }
      }

      if (node_counter != total_num_nodes) {
         std::cerr<<"Error: Final node counter = "<<node_counter<<" and total number of nodes = "<<total_num_nodes<<" not matching. \n";
         std::cerr<<"Error: Check UniformTriangleSymMesh().\n";
         exit(1);
      }  

      int element_counter = 0;
      // another loop for elements
      for (int j=0; j<ny_grid; j++) {
         // j is along y direction
         for(int i=0; i<nx_grid; i++) {
            // i is along x direction

            // 
            // Each grid consists of 4 triangles and 5 nodes
            //       n4          n3
            //         o---------o
            //         | \ T3  / |
            //         |  \  /   |
            //         |T0  o T2 |
            //         |  /  \   |
            //         | / T1  \ |
            //         o---------o
            //        n0        n2
            //
            // center node: n1
            //
            // size of grid is 2h
            //
            //
            // get node numbers
            std::vector<size_t> ns(5, 0);

            // first node
            ns[0] = j * (nx_grid + 1 + nx_grid) +  2*i;

            // second node
            ns[1] = j * (nx_grid + 1 + nx_grid) +  2*i + 1;

            // third node
            ns[2] = j * (nx_grid + 1 + nx_grid) +  2*i + 2;

            if (j+1 < ny_grid) {

               // fourth node
               ns[3] = (j+1) * (nx_grid + 1 + nx_grid) +  2*i + 2;

               // fifth node
               ns[4] = (j+1) * (nx_grid + 1 + nx_grid) +  2*i;
            }
            else {

               // fourth node
               ns[3] = (j+1) * (nx_grid + 1 + nx_grid) +  i + 1;

               // fifth node
               ns[4] = (j+1) * (nx_grid + 1 + nx_grid) +  i;
            }
            
            // verify 
            for (size_t k=0; k<5; k++)
               if (ns[k] != nodes[ns[k]].n) {
                  std::cerr << "Error: Node number incorrect, k = "<<k<<" ns[k] = "<<ns[k]<<" n = "<<nodes[ns[k]].n<<", i = "<<i<<" j = "<<j<<" nx = "<<nx_grid<<" ny = "<<ny_grid<<"\n";
                  exit(1);
               }

            //
            // first element
            //
            for (size_t k =0; k<4; k++) {

               int T = j * 4 * nx_grid + 4 * i + k;
               util::element elem = util::element(T);

               // set type
               elem.type = VTK_TYPE_TRIANGLE;

               // element-node connectivity
               if (k==0) {
                  elem.nodes.push_back(ns[0]);
                  elem.nodes.push_back(ns[1]);
                  elem.nodes.push_back(ns[4]);

                  // node-element connectivity
                  nodes[ns[0]].elems.push_back(T);
                  nodes[ns[1]].elems.push_back(T);
                  nodes[ns[4]].elems.push_back(T);
               }
               else if (k==1) {
                  elem.nodes.push_back(ns[0]);
                  elem.nodes.push_back(ns[2]);
                  elem.nodes.push_back(ns[1]);

                  // node-element connectivity
                  nodes[ns[0]].elems.push_back(T);
                  nodes[ns[2]].elems.push_back(T);
                  nodes[ns[1]].elems.push_back(T);
               }
               else if (k==2) {
                  elem.nodes.push_back(ns[1]);
                  elem.nodes.push_back(ns[2]);
                  elem.nodes.push_back(ns[3]);

                  // node-element connectivity
                  nodes[ns[1]].elems.push_back(T);
                  nodes[ns[2]].elems.push_back(T);
                  nodes[ns[3]].elems.push_back(T);
               }
               else if (k==3) {
                  elem.nodes.push_back(ns[1]);
                  elem.nodes.push_back(ns[3]);
                  elem.nodes.push_back(ns[4]);

                  // node-element connectivity
                  nodes[ns[1]].elems.push_back(T);
                  nodes[ns[3]].elems.push_back(T);
                  nodes[ns[4]].elems.push_back(T);
               }


               // area
               elem.area = 0.5 * mesh_size * mesh_size; 

               // store the element
               elements[T] = elem;
            }

            element_counter += 4;
         } // loop over i
      } // loop over j

      // check if total number of elements and element_counter is same
      if (element_counter != total_num_elems) {
         std::cerr<<"Error: Final element counter = "<<element_counter<<" and total number of elementes = "<<total_num_elems<<" not matching. \n";
         std::cerr<<"Error: Check UniformTriangleSymMesh().\n";
         exit(1);
      }  

      // check if total num of nodes calculated and size of nodes is same
      if (nodes.size() != total_num_nodes) {
         std::cerr<<"Error: Nodes data size = "<<nodes.size()<<" and total number of nodes = "<<total_num_nodes<<" not matching.\n";
         std::cerr<<"Error: Check UniformTriangleSymMesh().\n";
         exit(1);
      }

      if (elements.size() != total_num_elems) {
         std::cerr<<"Error: Elements data size = "<<elements.size()<<" and total number of elements = "<<total_num_elems<<" not matching.\n";
         std::cerr<<"Error: Check UniformTriangleSymMesh().\n";
         exit(1);
      }  

      //
      // write the bending data to file
      //
      if (bd_dout == true) {
         checkForBendingDataFE(config, domain, mesh_size, total_num_nodes, util::point(), true);
      } 

      std::cout<<"Total number of nodes (Neumann Data) = "<<total_num_nodes<<"\n";

      // call VtkWriter to output the data
      IO::VtkWriter  writer = IO::VtkWriter();
      writer.open(path_neum,mesh_neum_filename, false);
      writer.appendMesh(&nodes, &elements);
      // writer.appendData("Displacement", displacement);
      // writer.appendData("Velocity", velocity);
      writer.addTimeStep(0.0);
      writer.close();

      // close csv files
      File_disp_neum_left.close();
      File_disp_neum_right.close();
      File_disp_neum_bottom.close();
      File_disp_neum_top.close();
      File_disp_neum_left_wc.close();
      File_disp_neum_right_wc.close();
   } // end of neuman data   

   return;  
}

} // namespace fe_2d

void mesh::fe_twod(YAML::Node config) {

	// YAML::Node config = YAML::LoadFile("mesh_fe_2d.yaml");

	std::string mesh_type;
	if (!config["MeshData"]["Mesh_Type"]) {
      	std::cerr<<"Mesh_Type not specified.\n";
      	exit(1);
   }
   mesh_type = config["MeshData"]["Mesh_Type"].as<std::string>();

   if (mesh_type == "Uniform_Square")
   	fe_2d::UniformSquareMesh(config);
   else if (mesh_type == "Uniform_Triangle")
   	fe_2d::UniformTriangleMesh(config);
   else if (mesh_type == "Uniform_Triangle_Sym")
      fe_2d::UniformTriangleSymMesh(config);
   else {
   	std::cerr<<"Error: Check Mesh_Type data. Currently only Uniform_Square and Uniform_Triangle is implemented.\n";
   	exit(1);
   }

}