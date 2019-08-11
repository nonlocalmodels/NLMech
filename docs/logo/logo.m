clear all; clc;

% size
sim_Ly = 1.5;
sim_Lx = 1.0;
sim_Lz = 1.0;

% width of letter in x-axis
sim_Wx = 0.25;

% gap in vertical axis
sim_gapy = 0.5;

% width of letter M
sim_WM = 0.75 * sim_Lx;

% height of mid point in letter 'M'
sim_Mh = 0.25;

% width of mid horizontal line in letter M
sim_WMx = sim_gapy;

% mesh size
sim_h = 0.025;

% start line of letter M
sim_startMx = 2 * sim_Lx - 2 * sim_Wx;

% +ve slope of slanted line in N, M is
sim_slope = (sim_Ly - sim_gapy) / (sim_Lx - 2 * sim_Wx);

% choose mid point in letter M so that it has slope equal to sim_slope
sim_midPtMy = sim_Ly - sim_slope * (0.5 * sim_WM);

% choose consistent point for lower slanted line in letter M
sim_ptMliney = (sim_Ly - sim_gapy) - sim_slope * (0.5 * sim_WM + 0.5 * sim_WMx);


%
% create .geo file for gmsh
%
printf("Creating .geo file for gmsh\n");
geof = fopen('logo.geo','w');
fprintf(geof,"cl__1 = 1;\n");

%
% points
%
fprintf(geof,"Point(1) = {%4.6e, %4.6e, 0, %4.6e};\n", 0.0, 0.0, sim_h);
fprintf(geof,"Point(2) = {%4.6e, %4.6e, 0, %4.6e};\n", 0.0, sim_Ly, sim_h);
fprintf(geof,"Point(3) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_Wx, 0.0, sim_h);
fprintf(geof,"Point(4) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_Wx, sim_Ly - sim_gapy, sim_h);
fprintf(geof,"Point(5) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_Wx, sim_Ly, sim_h);
fprintf(geof,"Point(6) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_Lx - sim_Wx, 0.0, sim_h);
fprintf(geof,"Point(7) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_Lx - sim_Wx, sim_gapy, sim_h);
fprintf(geof,"Point(8) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_Lx - sim_Wx, sim_Ly, sim_h);
fprintf(geof,"Point(9) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_Lx, sim_Ly, sim_h);
fprintf(geof,"Point(10) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_Lx, 0.5 * sim_gapy, sim_h);
fprintf(geof,"Point(11) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_startMx, 0.5 * sim_gapy, sim_h);
fprintf(geof,"Point(12) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_startMx, sim_Ly, sim_h);
fprintf(geof,"Point(13) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_startMx + sim_Wx, 0.0, sim_h);
fprintf(geof,"Point(14) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_startMx + sim_Wx, sim_Ly - sim_gapy, sim_h);
fprintf(geof,"Point(15) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_startMx + sim_Wx, sim_Ly, sim_h);
fprintf(geof,"Point(16) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_startMx + sim_Wx + 0.5 * sim_WM, sim_midPtMy - sim_gapy, sim_h);
fprintf(geof,"Point(17) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_startMx + sim_Wx + 0.5 * sim_WM, sim_midPtMy, sim_h);
fprintf(geof,"Point(18) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_startMx + sim_Wx + sim_WM, 0.0, sim_h);
fprintf(geof,"Point(19) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_startMx + sim_Wx + sim_WM + sim_Wx, 0.0, sim_h);
fprintf(geof,"Point(20) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_startMx + sim_Wx + sim_WM, sim_Ly - sim_gapy, sim_h);
fprintf(geof,"Point(21) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_startMx + sim_Wx + sim_WM, sim_Ly, sim_h);
fprintf(geof,"Point(22) = {%4.6e, %4.6e, 0, %4.6e};\n", sim_startMx + sim_Wx + sim_WM + sim_Wx, sim_Ly, sim_h);

%
% lines
%
fprintf(geof,"Line(1) = {1, 2};\n");
fprintf(geof,"Line(2) = {2, 5};\n");
fprintf(geof,"Line(3) = {5, 7};\n");
fprintf(geof,"Line(4) = {7, 8};\n");
fprintf(geof,"Line(5) = {8, 9};\n");
fprintf(geof,"Line(6) = {9, 10};\n");
fprintf(geof,"Line(7) = {10, 11};\n");
fprintf(geof,"Line(8) = {11, 12};\n");
fprintf(geof,"Line(9) = {12, 15};\n");
fprintf(geof,"Line(10) = {15, 17};\n");
fprintf(geof,"Line(11) = {17, 21};\n");
fprintf(geof,"Line(12) = {21, 22};\n");
fprintf(geof,"Line(13) = {22, 19};\n");
fprintf(geof,"Line(14) = {19, 18};\n");
fprintf(geof,"Line(15) = {18, 20};\n");
fprintf(geof,"Line(16) = {20, 16};\n");
fprintf(geof,"Line(17) = {16, 14};\n");
fprintf(geof,"Line(18) = {14, 13};\n");
fprintf(geof,"Line(19) = {13, 6};\n");
fprintf(geof,"Line(20) = {6, 4};\n");
fprintf(geof,"Line(21) = {4, 3};\n");
fprintf(geof,"Line(22) = {3, 1};\n");

%
% surfaces
%
fprintf(geof,"Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22};\n");


%
% plane surface
%
fprintf(geof,"Plane Surface(1) = {1};\n");

% close file
fclose(geof);
