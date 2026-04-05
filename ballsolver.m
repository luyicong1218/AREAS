% function [xa,yb,zc,ra] = ballsolver(x)
% syms a b c r
% [a,b,c,r]= solve((x(1,1)-a)^2+(x(2,1)-b)^2+(x(3,1)-c)^2 == r^2,(x(1,2)-a)^2+(x(2,2)-b)^2+(x(3,2)-c)^2 == r^2,(x(1,3)-a)^2+(x(2,3)-b)^2+(x(3,3)-c)^2 == r^2,(x(1,4)-a)^2+(x(2,4)-b)^2+(x(3,4)-c)^2 == r^2);
% xa=a(1,1);
% yb=b(1,1);
% zc=c(1,1);
% ra=max(r(1,1),r(2,1));

A=[1
    0
    0];
B=[0
    1
    0];
C=[0
    0
    1];
D=[0
    0
    -1];
mat=[transpose([A B C D]),transpose([1,1,1,1])];
mat1=[transpose([-sum(A.^2),-sum(B.^2),-sum(C.^2),-sum(D.^2)]),mat(:,2:4)];
mat2=[mat(:,1),transpose([-sum(A.^2),-sum(B.^2),-sum(C.^2),-sum(D.^2)]),mat(:,3:4)];
mat3=[mat(:,1:2),transpose([-sum(A.^2),-sum(B.^2),-sum(C.^2),-sum(D.^2)]),mat(:,4)];
mat4=[mat(:,1:3),transpose([-sum(A.^2),-sum(B.^2),-sum(C.^2),-sum(D.^2)])];

temp_a=-det(mat1)/det(mat);
temp_b=-det(mat2)/det(mat);
temp_c=-det(mat3)/det(mat);
temp_r=det(mat4)/det(mat)

a=temp_a/-2;
b=temp_b/-2;
c=temp_c/-2;
r=sqrt(a^2+b^2+c^2-temp_r);


