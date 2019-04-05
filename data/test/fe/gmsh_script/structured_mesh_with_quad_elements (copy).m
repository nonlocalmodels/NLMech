clear all; clc;

% Domain
Lx = 1.0;
Ly = 1.0;
h = 0.2;

% print to file
printf("cl__1 = 1;\n");
% printf("Mesh.Algorithm = 8;\n");

% four corners of domain
p = 1;
printf("Point(%d) = {0, 0, 0, %4.6e};\n", p, h);
p++;
printf("Point(%d) = {%4.6e, 0, 0, %4.6e};\n", p, Lx, h);
p++;
printf("Point(%d) = {%4.6e, %4.6e, 0, %4.6e};\n", p, Lx, Ly, h);
p++;
printf("Point(%d) = {0, %4.6e, 0, %4.6e};\n", p, Ly, h);
p++;

% create lines on lower boundary 
printf("Line(1) = {1, 2};\n");
printf("Line(2) = {2, 3};\n");
printf("Line(3) = {3, 4};\n");
printf("Line(4) = {4, 1};\n");

% line loop to form surface
printf("Line Loop(5) = {1, 2, 3, 4};\n");
printf("Plane Surface(6) = {5};\n");

% for structured mesh
printf("Transfinite Surface {6} = {1, 2, 3, 4};\n");
printf("Transfinite Line {4, 2} = 10 Using Progression 1;\n");
printf("Transfinite Line {1, 3} = 15 Using Progression 1;\n");
printf("Recombine Surface {6};\n");