function new_point = mapPointBetweenTriangles2D(point,T1,T2)

% calculate areas

new_point=[];

%keyboard;

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
%new_point(3,1) = ( pcb * T2(1,3) ) + (pac * T2(2,3)) + (pab * T2(3,3));



%{
if 1

    % New
    new_point=[];

    %keyboard;
    %xsi = triangle_barycentric_2d ( [these_coords(1:2:end,1) these_coords(2:2:end,1)]', [new_mesh_coords_(ii(kk),1);new_mesh_coords_(ii(kk),2)] )
    %p = triangle_xsi_to_xy_2d ( [these_verts(1:2:end,1) these_verts(2:2:end,1)]', xsi )

    xsi = triangle_barycentric_2d ( T1', point );
    
    
    new_point=zeros(3,1);
    
    b2 = xsi(1,2);
    b3 = xsi(1,3);
    
    new_point(1,1) = ( (1 - b2 - b3) * T2(1,1) ) + (b2 * T2(2,1)) + (b3 * T2(3,1)) ;
    new_point(2,1) = ( (1 - b2 - b3) * T2(1,2) ) + (b2 * T2(2,2)) + (b3 * T2(3,2));
    new_point(3,1) = ( (1 - b2 - b3) * T2(1,3) ) + (b2 * T2(2,3)) + (b3 * T2(3,3));
    
    % Px = ( (1 - b2 - b3) * V1x ) + (b2 * V2x) + (b3 * V3x) ;
    %Py = ( (1 - b2 - b3) * V1y ) + (b2 * V2y) + (b3 * V3y);
    %Pz = ( (1 - b2 - b3) * V1z ) + (b2 * V2z) + (b3 * V3z);
    
    %new_point = triangle_xsi_to_xy_2d ( T2', xsi );

    
else

    % Old

    x = point(1,1);
    y = point(2,1);

    % Get point from T2 triangle and calculate its position in T1 triangle.

    xDash1 = T2(1,1);
    xDash2 = T2(2,1);
    xDash3 = T2(3,1);
    yDash1 = T2(1,2);
    yDash2 = T2(2,2);
    yDash3 = T2(3,2);

    denominator = -xDash2*yDash3 + xDash2*yDash1 + xDash1*yDash3 + xDash3*yDash2 - xDash3*yDash1 - xDash1*yDash2;

    beta = (y*xDash3 - xDash1*y - xDash3*yDash1 - yDash3*x +xDash1*yDash3 + x*yDash1)/denominator;

    gamma = (x*yDash2 - x*yDash1 - xDash1*yDash2 - xDash2*y + xDash2*yDash1 + xDash1*y)/denominator;

    alpha = 1 - (beta + gamma);

    newX = round(alpha * T1(1,1) + beta * T1(2,1) + gamma * T1(3,1));
    newY = round(alpha * T1(1,2) + beta * T1(2,2) + gamma * T1(3,2));

    new_point = [newX;newY];

end
%}