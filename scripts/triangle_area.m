function [area]=triangle_area(P)
% This function gives the area of a triangle
%
% [area]=triangle_area(Points)
%
% Points: The Points should be a numeric array, of size 3xn,
%         thus the points can be 2D, 3D... nD

%

[k,m]=size(P); if(k~=3), error('Points are not a 3xn array'); end

% Length of edges
L=[sqrt(sum((P(1,:)-P(2,:)).^2)) sqrt(sum((P(2,:)-P(3,:)).^2)) sqrt(sum((P(3,:)-P(1,:)).^2))];

% Area calculation with Heron's formula
s = ((L(1)+L(2)+L(3))/2);
area = sqrt(s*(s-L(1))*(s-L(2))*(s-L(3)));


