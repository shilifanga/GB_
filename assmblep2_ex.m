function M1 = assmblep2_ex(N,h,quad)
%ASSMBLEP2_EX 此处显示有关此函数的摘要
M = zeros(N+3,N+3);

for n1 = 4:N

    % M_{\alpha,\beta): alpha = beta-2;
    % test function
    index = 1 ;
    phi = GD_basis_p2(index);
    % trial function
    index = 3 ;
    [psi] =  GD_basis_p2leftshift4(index);
    % assmble position
    Nphi0phim2 = h./2.* sum(quad.weights.*phi(quad.points-2).*psi(quad.points-2));
    M(n1,n1-2) = Nphi0phim2;

    % M_{\alpha,\beta): alpha = beta-1;
    % test function
    index= 1;
    phi_1 = GD_basis_p2(index);
    % trial function
    index= 2;
    psi_1 = GD_basis_p2leftshift2(index);
    % test function
    index= 2;
    phi_2 = GD_basis_p2(index);
    % trial function
    index= 3;
    psi_2 = GD_basis_p2leftshift2(index);
    % assmble position
    M(n1,n1-1) = h./2.* sum(quad.weights.*phi_1(quad.points-2).*psi_1(quad.points-2)) + ...
        h./2.* sum(quad.weights.*phi_2(quad.points).*psi_2(quad.points)) ;
    Nphi0phim1 = M(n1,n1-1) ;
    % M_{\alpha,\beta): alpha = beta;
    % test function
    index= 1;
    phi_1 = GD_basis_p2(index);
    % trial function
    index= 1;
    psi_1 = GD_basis_p2(index);
    % test function
    index= 2;
    phi_2 = GD_basis_p2(index);
    % trial function
    index= 2;
    psi_2 = GD_basis_p2(index);
    % test function
    index= 3;
    phi_3 = GD_basis_p2(index);
    % trial function
    index= 3;
    psi_3 = GD_basis_p2(index);
    % assmble position
    M(n1,n1) =  h./2.* sum(quad.weights.*phi_1(quad.points-2).*psi_1(quad.points-2)) + ...
        h./2.* sum(quad.weights.*phi_2(quad.points).*psi_2(quad.points)) + ...
        h./2.* sum(quad.weights.*phi_3(quad.points+2).*psi_3(quad.points+2));
    Nphi0phi0 = M(n1,n1);
    % M_{\alpha,\beta): alpha = beta+1;
    % test function
    index= 2;
    phi_1 = GD_basis_p2(index);
    % trial function
    index= 1;
    psi_1 = GD_basis_p2rightshift2(index);
    % test function
    index= 3;
    phi_2 = GD_basis_p2(index);
    % trial function
    index= 2;
    psi_2 = GD_basis_p2rightshift2(index);
    % assmble position
    M(n1,n1+1) = h./2.* sum(quad.weights.*phi_1(quad.points).*psi_1(quad.points)) + ...
        h./2.* sum(quad.weights.*phi_2(quad.points+2).*psi_2(quad.points+2)) ;
    Nphi0phip1 =  M(n1,n1+1);

    % M_{\alpha,\beta): alpha = beta+2;
    % test function
    index = 3 ;
    phi = GD_basis_p2(index);
    % trial function
    index = 1 ;
    [psi] =  GD_basis_p2rightshift4(index);
    % assmble position
    M(n1,n1+2) = h./2.* sum(quad.weights.*phi(quad.points+2).*psi(quad.points+2));
    Nphi0phip2 =  M(n1,n1+2);
end

M1 = zeros(N+1,N+1);% \phi_{0,1,2,3,N}----->N+1行,N+1列
M1(2:end-1,:) = M(3:end-2,2:end-1);

% 利用外插法处理边界条件
% 左边界
% \phi_3 需要修改的项
% \int_{0}^{1} \phi_0\phi_{-1}
index= 2;
phi_1 = GD_basis_p2(index);
index= 3;
psi_1 = GD_basis_p2leftshift2(index);
phi0phim1 =  h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points+1)).*psi_1(1/2.*(quad.points+1)))./2;

index= 3;
phi_1 = GD_basis_p2(index);
phim1phim1 =  h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points+5)).*phi_1(1/2.*(quad.points+5)))./2;
M1(3,1) =  M1(3,1) - phi0phim1 - 4*phim1phim1;

index= 1;
phi_1 = GD_basis_p2(index);
index= 3;
psi_1 = GD_basis_p2leftshift4(index);
phi1phim1 =  h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points-3)).*psi_1(1/2.*(quad.points-3)))./2;
M1(3,2) =  M1(3,2) - phi1phim1 + 6*phim1phim1;
M1(3,3) =  M1(3,3) - 4*phim1phim1;
M1(3,4) =  M1(3,4) + 1*phim1phim1;

% \phi_2
M1(2,1) = Nphi0phim2 + 4*phi0phim1 + 16*phim1phim1 ;
M1(2,2) = Nphi0phim1 + 4*phim1phim1 - 24*phim1phim1 ;
M1(2,3) = Nphi0phi0 + 16*phim1phim1 ;
M1(2,4) = Nphi0phip1 - 4*phim1phim1 ;
M1(2,5) = Nphi0phip2;

% \phi_0
index= 2;
phi_1 = GD_basis_p2(index);
index= 2;
psi_1 = GD_basis_p2(index);
phi0phi0_1 =  h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points+1)).*psi_1(1/2.*(quad.points+1)))./2;
index= 3;
phi_1 = GD_basis_p2(index);
index= 3;
psi_1 = GD_basis_p2(index);
phi0phi0_2 =  h./2.* sum(quad.weights.*phi_1(quad.points+2).*psi_1(quad.points+2));

M1(1,1) = phi0phi0_1 + phi0phi0_2 + 8*phi0phim1 + 16*phim1phim1 ;

index= 1;
phi_1 = GD_basis_p2(index);
index= 2;
psi_1 = GD_basis_p2leftshift2(index);
phi1phi0_1 =  h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points+1)).*psi_1(1/2.*(quad.points+1)))./2;
index= 2;
phi_1 = GD_basis_p2(index);
index= 3;
psi_1 = GD_basis_p2leftshift2(index);
phi1phi0_2 =  h./2.* sum(quad.weights.*phi_1(quad.points+2).*psi_1(quad.points+2));

index= 1;
phi_1 = GD_basis_p2(index);
index= 3;
psi_1 = GD_basis_p2leftshift4(index);
phi1phim1 =  h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points+1)).*psi_1(1/2.*(quad.points+1)));

M1(1,2) = phi1phi0_1 + phi1phi0_2 + phi1phim1 -6*phi0phim1 - 24*phim1phim1;
M1(1,3) = Nphi0phim2 + 4*phi0phim1 + 16*phim1phim1;
M1(1,4) = -phi0phim1 - 4*phim1phim1;


%%% 利用对称性修改最后三行
len = flip(M1(3,1:5));
M1(N-1,N+1-4:N+1) = len;

len = flip(M1(2,1:5));
M1(N,N+1-4:N+1) = len;

len = flip(M1(1,1:5));
M1(N+1,N+1-4:N+1) = len;







end




function [phi] = GD_basis_p2(index)%,indexL,indexR

if index == 1
    phi = @(x) 1/8.*(x+4).*(x+2);
elseif index == 2
    phi = @(x) -1/4.*(x+2).*(x-2);
elseif index == 3
    phi = @(x) 1/8.*(x-2).*(x-4);
end


end

function [phi] = GD_basis_p2leftshift2(index)%,indexL,indexR
f = GD_basis_p2(index);
phi =  @(x) f(x + 2);
end

function [phi] = GD_basis_p2leftshift4(index)%,indexL,indexR
f = GD_basis_p2(index);
phi =  @(x) f(x + 4);
end

function [phi] = GD_basis_p2rightshift2(index)%,indexL,indexR
f = GD_basis_p2(index);
phi =  @(x) f(x - 2);
end

function [phi] = GD_basis_p2rightshift4(index)%,indexL,indexR
f = GD_basis_p2(index);
phi =  @(x) f(x - 4);
end


