function new_point = mapPointBetweenTriangles3D(point,T1,T2)

% This script does a barycentric coordinate mapping
% T1 and T2 are faces, 'point' is a location inside T1.
% The output is 'new_point', i.e. the location of 'point' inside T2.

new_point=[];

a = T1(1,:);
b = T1(2,:);
c = T1(3,:);

point = point';

pac = triangle_area([point;a;c]);
pab = triangle_area([point;a;b]);
pcb = triangle_area([point;c;b]);

total_a = pac+pab+pcb;

pac = pac/total_a;
pab = pab/total_a;
pcb = pcb/total_a;

new_point(1,1) = ( pcb * T2(1,1) ) + (pac * T2(2,1)) + (pab * T2(3,1)) ;
new_point(2,1) = ( pcb * T2(1,2) ) + (pac * T2(2,2)) + (pab * T2(3,2));
new_point(3,1) = ( pcb * T2(1,3) ) + (pac * T2(2,3)) + (pab * T2(3,3));

