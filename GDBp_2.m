% This function is used to assmble the mass matrices
% We assume that the computational domain [a, b] is uniformly! divided into N subintervals.
% For boundary, we simply consider the ghost basis methods.
% The boundary conditions we use are the Drichlet boundary conditions,
% which imply u(a)= ua; u(b)= ub with ua,ub are given.
% for p=2, the Extrapolation method is also concluded.

a = 0; b =1;
N = 1000;
h = (b-a)/N;
x = a:h:b;
xc = a+h/2:h:b-h/2;
%Drichlet boundary
ua = 1; ub= 1;
GDbasistype = 2;

boundarycase = 'ghost';
  boundarycase = 'extrapolation';
% 积分点
% quad = GaussQuadratureRule_line(5, 101);% [-1,1]
% 积分点
quad.points = [-0.9061798459386640; ...
    -0.5384693101056830; ...
    0; ...
    0.5384693101056830; ...
    0.9061798459386640];

quad.weights = [0.2369268850561890; ...
    0.4786286704993660; ...
    0.5688888888888889; ...
    0.4786286704993660; ...
    0.2369268850561890];

if GDbasistype==2
    switch boundarycase
        case 'ghost'
            A = assmblep2(N,h,quad);

            A = modifyboundaryp2(A,h,quad,N);

        case 'extrapolation'
            A = assmblep2_ex(N,h,quad);
    end
    % elseif GDbasistype==4
    %
    %     A = assmblep4(N,h,quad);
    %
    %     A = modifyboundaryp4(A,h,quad,N);

end


function M =  assmblep2(N,h,quad)

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
    M(n1,n1-2) = h./2.* sum(quad.weights.*phi(quad.points-2).*psi(quad.points-2));

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


    % M_{\alpha,\beta): alpha = beta+2;
    % test function
    index = 3 ;
    phi = GD_basis_p2(index);
    % trial function
    index = 1 ;
    [psi] =  GD_basis_p2rightshift4(index);
    % assmble position
    M(n1,n1+2) = h./2.* sum(quad.weights.*phi(quad.points+2).*psi(quad.points+2));

end





end



function M = modifyboundaryp2(M,h,quad,N)
% M(:,:) = 0;
% ghostmethod
% M_{\alpha,\beta): alpha = 2 (\phi_1); beta = alpha -2(\phi_{-1});
% test function
index = 1 ;
phi = GD_basis_p2(index);
% trial function
index = 3 ;
[psi] =  GD_basis_p2leftshift4(index);
% assmble position
M(3,1) = h./2.* sum(quad.weights.*phi(1/2.*(quad.points-3)).*psi(1/2.*(quad.points-3)))./2;

% % 检查
% f = @(r)  1/8.*(r).*(r+2) .* 1/8.*(r+2).*(r+4) ;
% result = integral(@(x)f(x), -2, -1, 'AbsTol', 1e-12)/2.*h;
% er = abs(result-M(3,1))
% M_{\alpha,\beta): alpha 3 (\phi_1); beta = alpha -1(\phi_{0});
% test function
index= 1;
phi_1 = GD_basis_p2(index);
% trial function
index= 2;
psi_1 = GD_basis_p2leftshift2(index);
% test function(正常区间)
index= 2;
phi_2 = GD_basis_p2(index);
% trial function
index= 3;
psi_2 = GD_basis_p2leftshift2(index);
% assmble position
M(3,2) = h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points-3)).*psi_1(1/2.*(quad.points-3)))./2 + ...
    h./2.* sum(quad.weights.*phi_2(quad.points).*psi_2(quad.points)) ;
% 检查
% f = @(r)  -1/4.*(r).*(r+4) .* (1/8).*(r+2).*(r+4) ;
% result1 = integral(@(x)f(x), -2, -1, 'AbsTol', 1e-12)/2.*h;
% g = @(r)  -1/4.*(r-2).*(r+2) .* (1/8).*(r-2+2).*(r-4+2) ;
% result2 = integral(@(x)g(x), -1, 1, 'AbsTol', 1e-12)/2.*h;
% result = result1 + result2;
% er = abs(result-M(3,2))

% M_{\alpha,\beta): alpha = beta;
% test function
index= 1;
phi_1 = GD_basis_p2(index);
% trial function (最左半个区间)
index= 1;
psi_1 = GD_basis_p2(index);
% test function (正常)
index= 2;
phi_2 = GD_basis_p2(index);
% trial function
index= 2;
psi_2 = GD_basis_p2(index);
% test function (正常)
index= 3;
phi_3 = GD_basis_p2(index);
% trial function
index= 3;
psi_3 = GD_basis_p2(index);
% assmble position
M(3,3) =  h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points-3)).*psi_1(1/2.*(quad.points-3)))./2 + ...
    h./2.* sum(quad.weights.*phi_2(quad.points).*psi_2(quad.points)) + ...
    h./2.* sum(quad.weights.*phi_3(quad.points+2).*psi_3(quad.points+2));
% % 检查
% f = @(r)  1/8.*(r+2).*(r+4) .* (1/8).*(r+2).*(r+4) ;
% result1 = integral(@(x)f(x), -2, -1, 'AbsTol', 1e-12)/2.*h;
% g = @(r)  -1/4.*(r-2).*(r+2) .* (-1/4).*(r-2).*(r+2) ;
% result2 = integral(@(x)g(x), -1, 1, 'AbsTol', 1e-12)/2.*h;
% hf = @(r)  1/8.*(r-2).*(r-4) .* (1/8).*(r-2).*(r-4) ;
% result3 = integral(@(x)hf(x), 1, 3, 'AbsTol', 1e-12)/2.*h;
% result = result1 + result2 + result3;
% er = abs(result-M(3,3))

%%% 对于test function: phi_1,与实验函数\phi_2,\phi_3是正常可以计算；从而在此先不计算
% M_{\alpha,\beta): alpha = beta+1;
M(3,4) = M(4,5) ;
M(3,5) = M(4,6) ;


%%%
% ghostmethod
% M_{\alpha,\beta): alpha = 1 (\phi_0); beta = alpha -1(\phi_{-1});
% test function
index = 2 ;
phi = GD_basis_p2(index);
% trial function
index = 3 ;
[psi] =  GD_basis_p2leftshift2(index);
% assmble position
M(2,1) = h./2.* sum(quad.weights.*phi(1/2.*(quad.points+1)).*psi(1/2.*(quad.points+1)))./2;

% 检查
% f = @(r)  1/8.*(r).*(r-2) .* (-1/4).*(r+2).*(r-2) ;
% result = integral(@(x)f(x), 0, 1, 'AbsTol', 1e-12)/2.*h;
% er = abs(result-M(2,1))
% M_{\alpha,\beta): alpha 1 (\phi_0); beta = alpha -1(\phi_{0});
% test function
index= 2;
phi_1 = GD_basis_p2(index);
% trial function
index= 2;
psi_1 = GD_basis_p2(index);
% test function(正常区间)
index= 3;
phi_2 = GD_basis_p2(index);
% trial function
index= 3;
psi_2 = GD_basis_p2(index);
% assmble position
M(2,2) = h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points+1)).*psi_1(1/2.*(quad.points+1)))./2 + ...
    h./2.* sum(quad.weights.*phi_2(quad.points+2).*psi_2(quad.points+2)) ;
% 检查
% f = @(r)  -1/4.*(r-2).*(r+2) .* (-1/4).*(r-2).*(r+2) ;
% result1 = integral(@(x)f(x), 0, 1, 'AbsTol', 1e-12)/2.*h;
% g = @(r)  1/8.*(r-2).*(r-4) .* (1/8).*(r-2).*(r-4) ;
% result2 = integral(@(x)g(x), 1, 3, 'AbsTol', 1e-12)/2.*h;
% result = result1 + result2;
% er = abs(result-M(2,2))

% M_{\alpha,\beta): alpha = beta;
% test function
index= 2;
phi_1 = GD_basis_p2(index);
% trial function
index= 1;
psi_1 = GD_basis_p2rightshift2(index);
% test function (正常)
index= 3;
phi_2 = GD_basis_p2(index);
% trial function
index= 2;
psi_2 = GD_basis_p2rightshift2(index);

% assmble position
M(2,3) =  h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points+1)).*psi_1(1/2.*(quad.points+1)))./2 + ...
    h./2.* sum(quad.weights.*phi_2(quad.points+2).*psi_2(quad.points+2)) ;

% 检查
% f = @(r)  -1/4.*(r-2).*(r+2) .* (1/8).*(r+2).*(r) ;
% result1 = integral(@(x)f(x), 0, 1, 'AbsTol', 1e-12)/2.*h;
% g = @(r)  1/8.*(r-2).*(r-4) .* (-1/4).*(r).*(r+4) ;
% result2 = integral(@(x)g(x), 1, 3, 'AbsTol', 1e-12)/2.*h;
% result = result1 + result2  ;
% er = abs(result-M(2,3))

% M_{\alpha,\beta): alpha = beta;
% test function
index= 3;
phi_1 = GD_basis_p2(index);
% trial function
index= 1;
psi_1 = GD_basis_p2rightshift4(index);

% assmble position
M(2,4) =  h./2.* sum(quad.weights.*phi_1(quad.points+2).*psi_1(quad.points+2)) ;

% 检查
% g = @(r)  1/8.*(r-2).*(r-4) .* (1/8).*(r).*(r-2) ;
% result2 = integral(@(x)g(x), 1, 3, 'AbsTol', 1e-12)/2.*h;
% result =   result2  ;
% er = abs(result-M(2,4))

%%% 修改右边界部分
% test function： \PHI_9
M(N+1,N+1-2) = M(N,N+1-3);% 正常
M(N+1,N+1-1) = M(N,N+1-2);% 正常
% M_{\alpha,\beta): alpha = beta;
% phi_9\phi_9
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
M(N+1,N+1) =  h./2.* sum(quad.weights.*phi_1(quad.points-2).*psi_1(quad.points-2)) + ...
    h./2.* sum(quad.weights.*phi_2(quad.points).*psi_2(quad.points)) + ...
    h./2.* sum(quad.weights.*phi_3(1/2.*(quad.points+3)).*psi_3(1/2.*(quad.points+3)))./2;
% 检查
% f = @(r)  1/8.*(r+2).*(r+4) .* (1/8).*(r+2).*(r+4) ;
% result1 = integral(@(x)f(x), -3, -1, 'AbsTol', 1e-12)/2.*h;
% g = @(r)  -1/4.*(r+2).*(r-2) .* (-1/4).*(r+2).*(r-2) ;
% result2 = integral(@(x)g(x), -1, 1, 'AbsTol', 1e-12)/2.*h;
% fr = @(r)  1/8.*(r-4).*(r-2) .* (1/8).*(r-4).*(r-2) ;
% result3 = integral(@(x)fr(x), 1, 2, 'AbsTol', 1e-12)/2.*h;
% result = result1 + result2 + result3;
% er = abs(result-M(N+1,N+1))

% M_{\alpha,\beta): alpha = beta + 1;
% phi_9\phi_10
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
M(N+1,N+2) =  h./2.* sum(quad.weights.*phi_1(quad.points).*psi_1(quad.points)) + ...
    h./2.* sum(quad.weights.*phi_2(1/2.*(quad.points+3)).*psi_2(1/2.*(quad.points+3)))./2;

% 检查
% f = @(r)  -1/4.*(r+2).*(r-2) .* (1/8).*(r+2).*(r) ;
% result1 = integral(@(x)f(x), -1, 1, 'AbsTol', 1e-12)/2.*h;
% g = @(r)  (1/8).*(r-2).*(r-4).* (-1/4).*(r-4).*(r) ;
% result2 = integral(@(x)g(x), 1, 2, 'AbsTol', 1e-12)/2.*h;
% result = result1 + result2  ;
% er = abs(result-M(N+1,N+2))

% M_{\alpha,\beta): alpha = beta + 2;
% phi_9\phi_11
% test function
index= 3;
phi_1 = GD_basis_p2(index);
% trial function
index= 1;
psi_1 = GD_basis_p2rightshift4(index);

% assmble position
M(N+1,N+3) =  h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points+3)).*psi_1(1/2.*(quad.points+3)))./2;

% % 检查
% f = @(r)  1/8.*(r).*(r-2) .* (1/8).*(r-2).*(r-4) ;
% result1 = integral(@(x)f(x), 1, 2, 'AbsTol', 1e-12)/2.*h;
%
% result = result1  ;
% er = abs(result-M(N+1,N+3))


%%% 修改右边界部分
% test function： \PHI_10
M(N+2,N+2-2) = M(N+1,N+2-3);% 正常\PHI_8\PHI_10
% M_{\alpha,\beta): alpha = beta;
% phi_9\phi_10
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
M(N+2,N+2-1) =  h./2.* sum(quad.weights.*phi_1(quad.points-2).*psi_1(quad.points-2)) + ...
    h./2.* sum(quad.weights.*phi_2(1/2.*(quad.points-1)).*psi_2(1/2.*(quad.points-1)))./2;
% 检查
% f = @(r)  1/8.*(r+2).*(r+4) .* (-1/4).*(r).*(r+4) ;
% result1 = integral(@(x)f(x), -3, -1, 'AbsTol', 1e-12)/2.*h;
% g = @(r)  -1/4.*(r+2).*(r-2) .* (1/8).*(r).*(r-2) ;
% result2 = integral(@(x)g(x), -1, 0, 'AbsTol', 1e-12)/2.*h;
% result = result1 + result2;
% er = abs(result-M(N+2,N+1))

% M_{\alpha,\beta): alpha = beta + 1;
% phi_10\phi_10
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

% assmble position
M(N+2,N+2) =  h./2.* sum(quad.weights.*phi_1(quad.points-2).*psi_1(quad.points-2)) + ...
    h./2.* sum(quad.weights.*phi_2(1/2.*(quad.points-1)).*psi_2(1/2.*(quad.points-1)))./2;

% 检查
% f = @(r)  1/8.*(r+2).*(r+4) .* (1/8).*(r+2).*(r+4) ;
% result1 = integral(@(x)f(x), -3, -1, 'AbsTol', 1e-12)/2.*h;
% g = @(r)  (-1/4).*(r-2).*(r+2).* (-1/4).*(r-2).*(r+2) ;
% result2 = integral(@(x)g(x), -1, 0, 'AbsTol', 1e-12)/2.*h;
% result = result1 + result2  ;
% er = abs(result-M(N+2,N+2))

% M_{\alpha,\beta): alpha = beta + 2;
% phi_11\phi_10
% test function
index= 2;
phi_1 = GD_basis_p2(index);
% trial function
index= 1;
psi_1 = GD_basis_p2rightshift2(index);
% assmble position
M(N+2,N+3) =  h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points-1)).*psi_1(1/2.*(quad.points-1)))./2;

% 检查
% f = @(r)  1/8.*(r).*(r+2) .* (-1/4).*(r-2).*(r+2) ;
% result1 = integral(@(x)f(x), -1, 0, 'AbsTol', 1e-12)/2.*h;
%
% result = result1  ;
% er = abs(result-M(N+2,N+3))

M(1,1) = 1;
M(end,end) = 1;




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