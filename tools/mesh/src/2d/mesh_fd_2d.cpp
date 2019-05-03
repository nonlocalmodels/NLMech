// Copyright (c)     2017
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

#include "mesh_twod.hpp"

static int counter=0;

namespace fd_2d {

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
      std::string &geom_neu_file,
      std::string &bc_l_neu_file,
      std::string &bc_r_neu_file,
      std::string &bc_t_neu_file,
      std::string &bc_b_neu_file,
      std::string &bc_l_neu_wc_file,
      std::string &bc_r_neu_wc_file,
      std::string &u_ic_neu_file,
      std::string &v_ic_neu_file,
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
      size_t & testFlag,
      std::string &testString,
      std::vector<double> &testParam,
      bool & void_enable,
      std::string &void_shape,
      bool & bending_out,
      YAML::Node config) {

   // YAML file
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

   if (config["Is_It_Test"]) {
      testFlag = 1;
      testString = config["Is_It_Test"].as<std::string>();
      auto a = config["Test_Param"];
      for (auto j : a) 
         testParam.push_back(j.as<double>());

      if (testParam.size() != 3) {
         std::cerr<<"Check number of parameters in Test_Param. It requires three parameters.\n";
         exit(1);
      }
   }
   else {
      testFlag = 0;
      testString.append("No_Test");
      testParam.push_back(0.);
   }

   // read ic data only if it is not test
   if (testFlag == 0) {
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

   //
   // check if void is enables
   //
   if (config["VoidData"]) {
      void_enable = true;

      if (config["VoidData"]["Rectangle"]) 
         void_shape = "Rectangle";
      else if (config["VoidData"]["Circle"]) 
         void_shape = "Circle";
      else if (config["VoidData"]["Ellipse"]) 
         void_shape = "Ellipse";

      std::cout<<void_enable<<","<<void_shape<<"\n";
   }

   // 
   // check if bending data is enabled
   //
   if (config["BendingData"]) 
      bending_out =  true;

   // perform checks
   if (testString != "No_Test") 
      if ( testString != "two_d_sine" and testString != "two_d_tquadratic") {
         std::cerr<<"Is_It_Test = "<<testString<<" in config file is not implemented."<<"\n";
         exit(1);
      }

   // check if material domain is [0,1] \times [0,1] for test simulation
   if (testFlag == 1) {
      if (Domain.first[0] <= -0.00001 or Domain.first[0] >= 0.00001 
         or Domain.first[1] <= -0.00001 or Domain.first[1] >= 0.00001) {
         std::cerr<<"Test_Flag = "<<testFlag<<" of config file is implemented only for domain [0,1]x[0,1]"<<"\n";
         exit(1);
      }
      if (Domain.second[0] <= 0.99999 or Domain.second[0] >= 1.00001 
         or Domain.second[1] <= 0.99999 or Domain.second[1] >= 1.00001) {
         std::cerr<<"Test_Flag = "<<testFlag<<" of config file is implemented only for domain [0,1]x[0,1]"<<"\n";
         exit(1);
      } 
   }   

   return;
} 

//
// getGuassianField()
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
   size_t test_flag,
   std::string test_string,
   std::vector<double> test_param,
   size_t u_ic_flag,
   std::vector<double> u_ic_params,
   size_t v_ic_flag,
   std::vector<double> v_ic_params) {

   if (test_flag == 0 ) {
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
      if (test_string == "two_d_sine") {
         // initial condition for testing the peridynamics code
         // see Wiki for more information
         const double tol = 1.0E-12;
         double dx = test_param[0];
         double dy = test_param[1];
         double a = test_param[2];

         double alpha_xi = 0.;
         if (compare::definitelyGreaterThan(xi.x, 0.0) 
            and compare::definitelyGreaterThan(xi.y, 0.0) 
            and compare::definitelyLessThan(xi.x, 1.0) 
            and compare::definitelyLessThan(xi.y, 1.0) 
            )
            alpha_xi = a * xi.x * xi.y * (1. -xi.x) * (1. -xi.y);

         double xi_dot_d = xi.x * dx + xi.y * dy;
         ui.x =  alpha_xi * std::sin(M_PI * xi_dot_d) * dx;
         ui.y =  alpha_xi * std::sin(M_PI * xi_dot_d) * dy;

         vi.x = M_PI * alpha_xi * std::cos(M_PI * xi_dot_d) * dx;
         vi.y = M_PI * alpha_xi * std::cos(M_PI * xi_dot_d) * dy;
      }

      if (test_string == "two_d_tquadratic") {
         // initial condition for testing the peridynamics code
         // see Wiki for more information
         const double tol = 1.0E-12;
         double dx = test_param[0];
         double dy = test_param[1];
         double a = test_param[2];

         double alpha_xi = 0.;
         if (compare::definitelyGreaterThan(xi.x, 0.0) 
            and compare::definitelyGreaterThan(xi.y, 0.0) 
            and compare::definitelyLessThan(xi.x, 1.0) 
            and compare::definitelyLessThan(xi.y, 1.0) 
            )
            alpha_xi = a * xi.x * xi.y * (1. -xi.x) * (1. -xi.y);

         ui.x =  0.;
         ui.y =  0.;

         vi.x = alpha_xi * dx;
         vi.y = alpha_xi * dy;
      }      
   }

   return;
}


//
// checkVoid
//
bool checkForVoidRectangle(YAML::Node config,
   util::point x) {

   static util::point p1 = util::point();
   static util::point p2 = util::point();

   static int init = -1;
   if(init == -1) {

      // read data
      auto e = config["VoidData"]["Rectangle"]["LB_Corner"];
      size_t k = 0;
      for (auto j : e) {
         if (k==0)
            p1.x = j.as<double>();
         if (k==1)
            p1.y = j.as<double>();

         k++;
      }

      e = config["VoidData"]["Rectangle"]["RT_Corner"];
      k = 0;
      for (auto j : e) {
         if (k==0)
            p2.x = j.as<double>();
         if (k==1)
            p2.y = j.as<double>();

         k++;
      }

      init = 0;
   }

   // check if point is inside the void
   return util::methods::isPointInsideRectangle(x, p1.x, p2.x, p1.y, p2.y);
}

bool checkForVoidCircle(YAML::Node config,
   util::point x) {

   static util::point center = util::point();
   static double radius = 0.0;

   static int init = -1;
   if(init == -1) {

      // read data
      auto e = config["VoidData"]["Circle"]["Center"];
      size_t k = 0;
      for (auto j : e) {
         center[k] = j.as<double>();
         k++;
      }

      radius = config["VoidData"]["Circle"]["Radius"].as<double>();

      init = 0;
   }

   // check if point is inside the void
   const double tol = 1.0E-12;
   util::point d = x - center;

   if (compare::definitelyGreaterThan(d.length(), radius - tol))
      return true;
   else 
      return false;
}

bool checkForVoidEllipse(YAML::Node config,
   util::point x) {

   // not implemented yet
   return false;
}

bool checkForVoid(YAML::Node config,
   std::string void_shape,
   util::point x) {

   if (void_shape == "Rectangle") 
      return checkForVoidRectangle(config, x);
   else if (void_shape == "Circle") 
      return checkForVoidCircle(config, x);
   else if (void_shape == "Ellipse") 
      return checkForVoidEllipse(config, x);
   else {
      std::cerr<<"Error: Check input file for void shape\n";
      exit(1);
   }
}

struct bendingData{
   bool th_pt_out;
   bool fr_pt_out;
   util::point sp_pt_l;
   util::point sp_pt_r;
   util::point ld_pt_l;
   util::point ld_pt_r;
};
#define BENDINGDATA_INIT {false, false, util::point(), util::point(), util::point(), util::point()}

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
#define FINDBENDINGPOINTDATA_INIT {0, 1000.0, 0, 1000.0, 0, 1000.0, 0, 1000.0}


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

} // end of namespace fd_2d

void mesh::fd_twod(YAML::Node config) {

   //  local variables
   // YAML::Node config;
   // std::string config_filename = "mesh_fd_2d.yaml";

   std::pair<std::vector<double>, std::vector<double>> domain;
   double horizon;
   int ratio;
   double mesh_size;
   size_t displ_ic_flag;
   std::vector<double> displ_ic_params;
   size_t vel_ic_flag;
   std::vector<double> vel_ic_params;
   int center[3];

   size_t test_flag;  // flag if input file is for testing the code
   std::string test_string;
   std::vector<double> test_param;

   bool create_dirichlet_files = false;
   std::string geom_diri_filename;
   std::string displ_bc_left_diri_filename;
   std::string displ_bc_right_diri_filename;
   std::string displ_bc_top_diri_filename;
   std::string displ_bc_bottom_diri_filename;
   std::string displ_ic_diri_filename;
   std::string vel_ic_diri_filename;
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
   std::string displ_bc_left_neum_wc_filename;
   std::string displ_bc_right_neum_wc_filename;

   // bending loading related data
   bool bd_dout = false;

   // void related data
   bool void_enable = false;
   std::string void_shape;

   // flag for volume correction
   bool vol_cor_flag = true;

   // read data from config file
   fd_2d::readDataFile(geom_diri_filename,
      displ_bc_left_diri_filename,
      displ_bc_right_diri_filename,
      displ_bc_top_diri_filename,
      displ_bc_bottom_diri_filename,
      displ_bc_left_diri_wc_filename,
      displ_bc_right_diri_wc_filename,
      displ_ic_diri_filename,
      vel_ic_diri_filename,
      geom_neum_filename,
      displ_bc_left_neum_filename,
      displ_bc_right_neum_filename,
      displ_bc_top_neum_filename,
      displ_bc_bottom_neum_filename,
      displ_bc_left_neum_wc_filename,
      displ_bc_right_neum_wc_filename,
      displ_ic_neum_filename,
      vel_ic_neum_filename,
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
      test_flag,
      test_string,
      test_param,
      void_enable,
      void_shape,
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
      std::ofstream F_dir_geom;
      std::ofstream F_dir_u_l;
      std::ofstream F_dir_u_r;
      std::ofstream F_dir_u_t;
      std::ofstream F_dir_u_b;
      std::ofstream F_dir_u_ic;
      std::ofstream F_dir_v_ic;
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

      if (test_flag != 0 or displ_ic_flag != 0) {
         F_dir_u_ic.open(displ_ic_diri_filename);
         // header
         F_dir_u_ic<<"id,x,y\n";
      }

      if (test_flag != 0 or vel_ic_flag != 0) {
         F_dir_v_ic.open(vel_ic_diri_filename);
         // header
         F_dir_v_ic<<"id,x,y\n";   
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
      for (int i=0;i<=nx;i++) {
         for (int j=0;j<=ny;j++) {

            util::point x = util::point();
            x.x = nonlocal_domain.first[0] + double(i) * mesh_size;
            x.y = nonlocal_domain.first[1] + double(j) * mesh_size;

            double vol_i = mesh_size * mesh_size;

            if (vol_cor_flag == true) {

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
            }

            // write to geometry file
            F_dir_geom<<n_nodes<<","<<x.x<<","<<x.y<<","<<vol_i<<"\n";

            // get initial condition on nodes
            util::point u = util::point();
            util::point v = util::point();
            
            fd_2d::computeIC(x, u, v, domain, test_flag, test_string, test_param,
               displ_ic_flag, displ_ic_params,
               vel_ic_flag, vel_ic_params);

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

            // x is inside material domain
            if (test_flag != 0 or displ_ic_flag != 0)
               F_dir_u_ic<<n_nodes<<","<<u.x<<","<<u.y<<"\n";
            if (test_flag != 0 or vel_ic_flag != 0) 
               F_dir_v_ic<<n_nodes<<","<<v.x<<","<<v.y<<"\n";           
         
            // increment the mesh node counter
            n_nodes++;
         } // loop over j
      } // loop over i

      if (total_num_nodes != n_nodes) {
         std::cerr<<"Error: Check total number of nodes.\n";
         exit(1);
      }

      std::cout<<"Number of nodes (Dirichlet Data) = "<<n_nodes<<"\n";

      std::cout<<"Domain = ["<<domain.first[0]<<","<<domain.first[1]<<"]x["<<domain.second[0]<<","<<domain.second[1]<<"]\n";
      std::cout<<"Nonlocal domain = ["<<nonlocal_domain.first[0]<<","<<nonlocal_domain.first[1]<<"]x["<<nonlocal_domain.second[0]<<","<<nonlocal_domain.second[1]<<"]\n";

      // close all files
      F_dir_geom.close();
      F_dir_u_l.close();
      F_dir_u_r.close();
      F_dir_u_t.close();
      F_dir_u_b.close();
      if (test_flag != 0 or displ_ic_flag != 0)
         F_dir_u_ic.close();
      if (test_flag != 0 or vel_ic_flag != 0) 
         F_dir_v_ic.close();
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
      
      if (test_flag != 0 or displ_ic_flag != 0) {
         F_neu_u_ic.open(displ_ic_neum_filename);
         // header
         F_neu_u_ic<<"id,x,y\n";
      }

      if (test_flag != 0 or vel_ic_flag != 0) {
         F_neu_v_ic.open(vel_ic_neum_filename);
         // header
         F_neu_v_ic<<"id,x,y\n";
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
      for (int i=0;i<=nx;i++) {
         for (int j=0;j<=ny;j++) {

            util::point x = util::point();
            x.x = domain.first[0] + double(i) * mesh_size;
            x.y = domain.first[1] + double(j) * mesh_size;

            double vol_i = mesh_size * mesh_size;

            if (vol_cor_flag == true) {
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
            }

            //
            // check for void, if inside void then skip this point
            //
            bool inside_void = false;
            if (void_enable == true)
               inside_void = fd_2d::checkForVoid(config, void_shape, x);

            if (inside_void == true)
               continue;

            // write to geometry file
            F_neu_geom<<n_nodes<<","<<x.x<<","<<x.y<<","<<vol_i<<"\n";

            // get initial condition on nodes
            util::point u = util::point();
            util::point v = util::point();
            
            fd_2d::computeIC(x, u, v, domain, test_flag, test_string, test_param,
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
            if (test_flag != 0 or displ_ic_flag != 0) 
               F_neu_u_ic<<n_nodes<<","<<u.x<<","<<u.y<<"\n";
            if (test_flag != 0 or vel_ic_flag != 0) 
               F_neu_v_ic<<n_nodes<<","<<v.x<<","<<v.y<<"\n";       

            //
            // we now check for bending data
            //
            if (bd_dout == true) {
               fd_2d::checkForBendingData(config, domain, n_nodes, x);
            }    
         
            // increment the mesh node counter
            n_nodes++;
         } // loop over j
      } // loop over i

      //
      // write the bending data to file
      //
      if (bd_dout == true) {
         fd_2d::checkForBendingData(config, domain, n_nodes, util::point(), true);
      } 

      if (total_num_nodes != n_nodes) {
         std::cerr<<"Error: Check total number of nodes.\n";
         exit(1);
      }
      std::cout<<"Number of nodes (Neumann Data) = "<<n_nodes<<"\n";
      std::cout<<"Domain = ["<<domain.first[0]<<","<<domain.first[1]<<"]x["<<domain.second[0]<<","<<domain.second[1]<<"]\n";

      // close all files
      F_neu_geom.close();
      F_neu_u_l.close();
      F_neu_u_r.close();
      F_neu_u_t.close();
      F_neu_u_b.close();
      if (test_flag != 0 or displ_ic_flag != 0) 
         F_neu_u_ic.close();
      if (test_flag != 0 or vel_ic_flag != 0) 
         F_neu_v_ic.close();
      F_neu_u_l_wc.close();
      F_neu_u_r_wc.close();
   }   
}