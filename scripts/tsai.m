function p1 = tsai(Cal_1C,v1)

%
% function p1 = tsai(Cal_1C,v1)
%
% Given a 3D point v1, find its 2D projection p1.
% Requires camera paraneter matrix Cal_1C .
%
% (c) Darren Cosker, 2008.
%
% Based on an earlier C++ implementation by 3DMD Ltd
%

v1.x = v1.x - Cal_1C.X;
v1.y = v1.y - Cal_1C.Y;
v1.z = v1.z - Cal_1C.Z;

p1.x = (Cal_1C.M(1,1)*v1.x) + (Cal_1C.M(2,1)*v1.y) + (Cal_1C.M(3,1)*v1.z);
p1.y = (Cal_1C.M(1,2)*v1.x) + (Cal_1C.M(2,2)*v1.y) + (Cal_1C.M(3,2)*v1.z);
p1.z = (Cal_1C.M(1,3)*v1.x) + (Cal_1C.M(2,3)*v1.y) + (Cal_1C.M(3,3)*v1.z);

p1.x = (Cal_1C.f * p1.x)/p1.z;
p1.y = (Cal_1C.f * p1.y)/p1.z;

R = sqrt((p1.x.*p1.x)+(p1.y.*p1.y));

r=R;
r2 = r*r;
epsilon = 0.01/1000.0;

fprime = 1+r2*(3*Cal_1C.K + r2*5*Cal_1C.K2);

if(fprime<1.0e-4)
    disp('fprime is low');
end
error=10000000;

while((abs(R)>0.001) & (abs(error/R)>epsilon))
    distortion = r*(1 + r2*(Cal_1C.K+r2*Cal_1C.K2));
    r = r - (distortion - R)/fprime;
    error = (R - distortion)/fprime;
end

p1.x = p1.x * (r/R);
p1.y = p1.y * (r/R);

p1.x = p1.x/Cal_1C.x_;
p1.y = p1.y/Cal_1C.y_;

p1.x = p1.x + Cal_1C.a;
p1.y = p1.y + Cal_1C.b;
p1.y = Cal_1C.is(1,2) - p1.y;