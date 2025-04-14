% This function is used to assmble the mass matrices
% We assume that the computational domain [a, b] is uniformly! divided into N subintervals.
% For boundary, we simply consider the ghost basis methods.
% The boundary conditions we use are the Drichlet boundary conditions,
% which imply u(a)= ua; u(b)= ub with ua,ub are given.
a = 0; b =1;
N = 10;
h = (b-a)/N;
x = a:h:b;
xc = a+h/2:h:b-h/2;
%Drichlet boundary
ua = 1; ub= 1;
GDbasistype = 4;

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


A = assmblep4(N,h,quad);

A(1,1) = 1;
A(2,2) = 1;
A(end-1,end-1) = 1;
A(end,end) = 1;

function M =  assmblep4(N,h,quad)

M = zeros(N+5,N+5);

%% STEP1 : 内部点
for n1 = 6:N

    % M_{\alpha,\beta): alpha = beta-4;
    % test function
    index = 1 ;
    phi = GD_basis_p4(index);
    % trial function
    index = 5 ;
    [psi] =  GD_basis_p4leftshift8(index);
    % assmble position
    M(n1,n1-4) = h./2.* sum(quad.weights.*phi(quad.points-4).*psi(quad.points-4));

    % M_{\alpha,\beta): alpha = beta-3;
    % test function
    index= 1;
    phi_1 = GD_basis_p4(index);
    % trial function
    index= 4;
    psi_1 = GD_basis_p4leftshift6(index);
    % test function
    index= 2;
    phi_2 = GD_basis_p4(index);
    % trial function
    index= 5;
    psi_2 = GD_basis_p4leftshift6(index);
    % assmble position
    N1_N1M3_1 = h./2.* sum(quad.weights.*phi_1(quad.points-4).*psi_1(quad.points-4)) ;
    N1_N1M3_2 = h./2.* sum(quad.weights.*phi_2(quad.points-2).*psi_2(quad.points-2)) ;

    M(n1,n1-3) =  N1_N1M3_1+ N1_N1M3_2;
    
    % M_{\alpha,\beta): alpha = beta-2;
    % test function
    index= 1;
    phi_1 = GD_basis_p4(index);
    % trial function
    index= 3;
    psi_1 = GD_basis_p4leftshift4(index);
    % test function
    index= 2;
    phi_2 = GD_basis_p4(index);
    % trial function
    index= 4;
    psi_2 = GD_basis_p4leftshift4(index);
    % test function
    index= 3;
    phi_3 = GD_basis_p4(index);
    % trial function
    index= 5;
    psi_3 = GD_basis_p4leftshift4(index);
    % assmble position
    N1_N1M2_1 =  h./2.* sum(quad.weights.*phi_1(quad.points-4).*psi_1(quad.points-4));
    N1_N1M2_2 =  h./2.* sum(quad.weights.*phi_2(quad.points-2).*psi_2(quad.points-2));
    N1_N1M2_3 =  h./2.* sum(quad.weights.*phi_3(quad.points).*psi_3(quad.points));
    M(n1,n1-2) = N1_N1M2_1 + N1_N1M2_2 + N1_N1M2_3;
   

    % M_{\alpha,\beta): alpha = beta-1;
    % test function
    index= 1;
    phi_1 = GD_basis_p4(index);
    % trial function
    index= 2;
    psi_1 = GD_basis_p4leftshift2(index);
    % test function
    index= 2;
    phi_2 = GD_basis_p4(index);
    % trial function
    index= 3;
    psi_2 = GD_basis_p4leftshift2(index);
    % test function
    index= 3;
    phi_3 = GD_basis_p4(index);
    % trial function
    index= 4;
    psi_3 = GD_basis_p4leftshift2(index);
    % test function
    index= 4;
    phi_4 = GD_basis_p4(index);
    % trial function
    index= 5;
    psi_4 = GD_basis_p4leftshift2(index);
    % assmble position
    N1_N1M1_1 =  h./2.* sum(quad.weights.*phi_1(quad.points-4).*psi_1(quad.points-4))  ;
    N1_N1M1_2 =    h./2.* sum(quad.weights.*phi_2(quad.points-2).*psi_2(quad.points-2)) ;
    N1_N1M1_3 =   h./2.* sum(quad.weights.*phi_3(quad.points).*psi_3(quad.points)) ;
    N1_N1M1_4 =    h./2.* sum(quad.weights.*phi_4(quad.points+2).*psi_4(quad.points+2));
    M(n1,n1-1) =  N1_N1M1_1 +  N1_N1M1_2 +  N1_N1M1_3 +  N1_N1M1_4;
   

    % M_{\alpha,\beta): alpha = beta;
    % test function
    index= 1;
    phi_1 = GD_basis_p4(index);
    % trial function
    index= 1;
    psi_1 = GD_basis_p4(index);
    % test function
    index= 2;
    phi_2 = GD_basis_p4(index);
    % trial function
    index= 2;
    psi_2 = GD_basis_p4(index);
    % test function
    index= 3;
    phi_3 = GD_basis_p4(index);
    % trial function
    index= 3;
    psi_3 = GD_basis_p4(index);
    % test function
    index= 4;
    phi_4 = GD_basis_p4(index);
    % trial function
    index= 4;
    psi_4 = GD_basis_p4(index);
    % test function
    index= 5;
    phi_5 = GD_basis_p4(index);
    % trial function
    index= 5;
    psi_5 = GD_basis_p4(index);
    % assmble position
    N1_N1_1 =   h./2.* sum(quad.weights.*phi_1(quad.points-4).*psi_1(quad.points-4)) ;
    N1_N1_2 =   h./2.* sum(quad.weights.*phi_2(quad.points-2).*psi_2(quad.points-2)) ;
    N1_N1_3 =   h./2.* sum(quad.weights.*phi_3(quad.points).*psi_3(quad.points)) ;
    N1_N1_4 =   h./2.* sum(quad.weights.*phi_4(quad.points+2).*psi_4(quad.points+2));
    N1_N1_5 =  h./2.* sum(quad.weights.*phi_5(quad.points+4).*psi_5(quad.points+4));
    M(n1,n1) =   N1_N1_1 +  N1_N1_2 +  N1_N1_3 +  N1_N1_4 +  N1_N1_5;
    

    % M_{\alpha,\beta): alpha = beta+1;
    % test function
    index= 2;
    phi_1 = GD_basis_p4(index);
    % trial function
    index= 1;
    psi_1 = GD_basis_p4rightshift2(index);
    % test function
    index= 3;
    phi_2 = GD_basis_p4(index);
    % trial function
    index= 2;
    psi_2 = GD_basis_p4rightshift2(index);
    % test function
    index= 4;
    phi_3 = GD_basis_p4(index);
    % trial function
    index= 3;
    psi_3 = GD_basis_p4rightshift2(index);
    % test function
    index= 5;
    phi_4 = GD_basis_p4(index);
    % trial function
    index= 4;
    psi_4 = GD_basis_p4rightshift2(index);

    % assmble position
    N1_N1P1_1 =   h./2.* sum(quad.weights.*phi_1(quad.points-2).*psi_1(quad.points-2));
    N1_N1P1_2 =      h./2.* sum(quad.weights.*phi_2(quad.points).*psi_2(quad.points));
    N1_N1P1_3 =      h./2.* sum(quad.weights.*phi_3(quad.points+2).*psi_3(quad.points+2));
    N1_N1P1_4 =      h./2.* sum(quad.weights.*phi_4(quad.points+4).*psi_4(quad.points+4));
    M(n1,n1+1) =  N1_N1P1_1 +  N1_N1P1_2 +  N1_N1P1_3 +  N1_N1P1_4;
    
    % M_{\alpha,\beta): alpha = beta+2;
    % test function
    index= 3;
    phi_1 = GD_basis_p4(index);
    % trial function
    index= 1;
    psi_1 = GD_basis_p4rightshift4(index);
    % test function
    index= 4;
    phi_2 = GD_basis_p4(index);
    % trial function
    index= 2;
    psi_2 = GD_basis_p4rightshift4(index);
    % test function
    index= 5;
    phi_3 = GD_basis_p4(index);
    % trial function
    index= 3;
    psi_3 = GD_basis_p4rightshift4(index);


    % assmble position
    N1_N1P2_1 =  h./2.* sum(quad.weights.*phi_1(quad.points).*psi_1(quad.points));
    N1_N1P2_2 =     h./2.* sum(quad.weights.*phi_2(quad.points+2).*psi_2(quad.points+2));
    N1_N1P2_3 =      h./2.* sum(quad.weights.*phi_3(quad.points+4).*psi_3(quad.points+4));
    M(n1,n1+2) = N1_N1P2_1  + N1_N1P2_2 + N1_N1P2_3;
    

    % M_{\alpha,\beta): alpha = beta+3;
    % test function
    index= 4;
    phi_1 = GD_basis_p4(index);
    % trial function
    index= 1;
    psi_1 = GD_basis_p4rightshift6(index);
    % test function
    index= 5;
    phi_2 = GD_basis_p4(index);
    % trial function
    index= 2;
    psi_2 = GD_basis_p4rightshift6(index);

    % assmble position
    N1_N1P3_1 = h./2.* sum(quad.weights.*phi_1(quad.points+2).*psi_1(quad.points+2));
    N1_N1P3_2 = h./2.* sum(quad.weights.*phi_2(quad.points+4).*psi_2(quad.points+4));
    M(n1,n1+3) =  N1_N1P3_1 +  N1_N1P3_2 ;
    

    % M_{\alpha,\beta): alpha = beta+4;
    % test function
    index= 5;
    phi_1 = GD_basis_p4(index);
    % trial function
    index= 1;
    psi_1 = GD_basis_p4rightshift8(index);

    % assmble position
    M(n1,n1+4) =  h./2.* sum(quad.weights.*phi_1(quad.points+4).*psi_1(quad.points+4));
    

end


%% STEP2 : 边界点
index = 1 ;
phi = GD_basis_p4(index);
% trial function
index = 5 ;
[psi] =  GD_basis_p4leftshift8(index);
% assmble position
M(5,1) = h./2.* sum(quad.weights.*phi(1/2.*(quad.points-7)).*psi(1/2.*(quad.points-7)))./2;
 
% test function
index= 1;
phi_1 = GD_basis_p4(index);
% trial function
index= 4;
psi_1 = GD_basis_p4leftshift6(index);

% assmble position
M(5,2) = h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points-7)).*psi_1(1/2.*(quad.points-7)))./2 + ...
    N1_N1M3_2 ;
 

% M_{\alpha,\beta): alpha = beta-2;
% test function
index= 1;
phi_1 = GD_basis_p4(index);
% trial function
index= 3;
psi_1 = GD_basis_p4leftshift4(index);
% assmble position
M(5,3) =  h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points-7)).*psi_1(1/2.*(quad.points-7)))./2+ ...
    N1_N1M2_2 + N1_N1M2_3;

 
% test function
index= 1;
phi_1 = GD_basis_p4(index);
% trial function
index= 2;
psi_1 = GD_basis_p4leftshift2(index);
% assmble position
M(5,4) =   h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points-7)).*psi_1(1/2.*(quad.points-7)))./2 + ...
    N1_N1M1_2 + N1_N1M1_3 + N1_N1M1_4;
 

index= 1;
phi_1 = GD_basis_p4(index);
% trial function
index= 1;
psi_1 = GD_basis_p4(index);
% assmble position
M(5,5) =   h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points-7)).*psi_1(1/2.*(quad.points-7)))./2 + ...
    N1_N1_2 + N1_N1_3 + N1_N1_4  + N1_N1_5;
 
%%% 对于test function: phi_1,与实验函数\phi_2,\phi_3,\phi_4是正常可以计算；从而在此先不计算
% M_{\alpha,\beta): alpha = beta+1;
M(5,6) = M(6,7) ;
M(5,7) = M(6,8) ;
M(5,8) = M(6,9) ;
M(5,9) = M(6,10) ;

%%% \phi_1
index = 2 ;
phi = GD_basis_p4(index);
% trial function
index = 5 ;
[psi] =  GD_basis_p4leftshift6(index);
% assmble position
M(4,1) = h./2.* sum(quad.weights.*phi(1/2.*(quad.points-3)).*psi(1/2.*(quad.points-3)))./2;
 

% test function
index= 2;
phi_1 = GD_basis_p4(index);
% trial function
index= 4;
psi_1 = GD_basis_p4leftshift4(index);

% assmble position
M(4,2) = h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points-3)).*psi_1(1/2.*(quad.points-3)))./2 + ...
    N1_N1M2_3 ;
 

% test function
index= 2;
phi_1 = GD_basis_p4(index);
% trial function
index= 3;
psi_1 = GD_basis_p4leftshift2(index);
% assmble position
M(4,3) =  h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points-3)).*psi_1(1/2.*(quad.points-3)))./2+ ...
    N1_N1M1_3 + N1_N1M1_4;

 

% test function
index= 2;
phi_1 = GD_basis_p4(index);
% trial function
index= 2;
psi_1 = GD_basis_p4(index);
% assmble position
M(4,4) =   h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points-3)).*psi_1(1/2.*(quad.points-3)))./2 + ...
    N1_N1_3 + N1_N1_4 + N1_N1_5;
 

index= 2;
phi_1 = GD_basis_p4(index);
% trial function
index= 1;
psi_1 = GD_basis_p4rightshift2(index);
% assmble position
M(4,5) =   h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points-3)).*psi_1(1/2.*(quad.points-3)))./2 + ...
    N1_N1P1_2 + N1_N1P1_3 + N1_N1P1_4  ;

 
%%% 对于test function: phi_1,与实验函数\phi_2,\phi_3,\phi_4是正常可以计算；从而在此先不计算
% M_{\alpha,\beta): alpha = beta+1;
M(4,6) = M(5,7) ;
M(4,7) = M(5,8) ;
M(4,8) = M(5,9) ;

%%% \phi_0
index = 3 ;
phi = GD_basis_p4(index);
% trial function
index = 5 ;
[psi] =  GD_basis_p4leftshift4(index);
% assmble position
M(3,1) = h./2.* sum(quad.weights.*phi(1/2.*(quad.points+1)).*psi(1/2.*(quad.points+1)))./2;
 

% test function
index= 3;
phi_1 = GD_basis_p4(index);
% trial function
index= 4;
psi_1 = GD_basis_p4leftshift2(index);

% assmble position
M(3,2) = h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points+1)).*psi_1(1/2.*(quad.points+1)))./2 + ...
    N1_N1M1_4 ;
 
% test function
index= 3;
phi_1 = GD_basis_p4(index);
% trial function
index= 3;
psi_1 = GD_basis_p4(index);
% assmble position
M(3,3) =  h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points+1)).*psi_1(1/2.*(quad.points+1)))./2+ ...
    N1_N1_4 + N1_N1_5;

 
% test function
index= 3;
phi_1 = GD_basis_p4(index);
% trial function
index= 2;
psi_1 = GD_basis_p4rightshift2(index);
% assmble position
M(3,4) =   h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points+1)).*psi_1(1/2.*(quad.points+1)))./2 + ...
    N1_N1P1_3 + N1_N1P1_4 ;
 

index= 3;
phi_1 = GD_basis_p4(index);
% trial function
index= 1;
psi_1 = GD_basis_p4rightshift4(index);
% assmble position
M(3,5) =   h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points+1)).*psi_1(1/2.*(quad.points+1)))./2 + ...
    N1_N1P2_2 +  N1_N1P2_3   ;
 
%%% 对于test function: phi_1,与实验函数\phi_2,\phi_3,\phi_4是正常可以计算；从而在此先不计算
% M_{\alpha,\beta): alpha = beta+1;
M(3,6) = M(4,7) ;
M(3,7) = M(4,8) ;


%%% 右边界的处理
% \phi_8
M(N+1,N-3) = M(N,N-4) ;
M(N+1,N-2) = M(N,N-3) ;
M(N+1,N-1) = M(N,N-2) ;
M(N+1,N) = M(N,N-1) ;

index= 5;
phi_1 = GD_basis_p4(index);
% trial function
index= 5;
psi_1 = GD_basis_p4(index);
% assmble position
M(N+1,N+1) =   h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points+7)).*psi_1(1/2.*(quad.points+7)))./2 + ...
    N1_N1_1 + N1_N1_2 + N1_N1_3  + N1_N1_4;
 

index = 5 ;
phi = GD_basis_p4(index);
% trial function
index = 4 ;
[psi] =  GD_basis_p4rightshift2(index);
% assmble position
M(N+1,N+2) = h./2.* sum(quad.weights.*phi(1/2.*(quad.points+7)).*psi(1/2.*(quad.points+7)))./2 + ...
    N1_N1P1_1 + N1_N1P1_2 + N1_N1P1_3;

  
index = 5 ;
phi = GD_basis_p4(index);
% trial function
index = 3 ;
[psi] =  GD_basis_p4rightshift4(index);
% assmble position
M(N+1,N+3) = h./2.* sum(quad.weights.*phi(1/2.*(quad.points+7)).*psi(1/2.*(quad.points+7)))./2 + ...
    N1_N1P2_1 + N1_N1P2_2;

  
index = 5 ;
phi = GD_basis_p4(index);
% trial function
index = 2 ;
[psi] =  GD_basis_p4rightshift6(index);
% assmble position
M(N+1,N+4) = h./2.* sum(quad.weights.*phi(1/2.*(quad.points+7)).*psi(1/2.*(quad.points+7)))./2 + ...
    N1_N1P3_1;
 

index = 5 ;
phi = GD_basis_p4(index);
% trial function
index = 1 ;
[psi] =  GD_basis_p4rightshift8(index);
% assmble position
M(N+1,N+5) = h./2.* sum(quad.weights.*phi(1/2.*(quad.points+7)).*psi(1/2.*(quad.points+7)))./2;

  
% \phi_9
M(N+2,N-2) = M(N,N-4) ;
M(N+2,N-1) = M(N,N-3) ;
M(N+2,N) = M(N,N-2) ;

index= 4;
phi_1 = GD_basis_p4(index);
% trial function
index= 5;
psi_1 = GD_basis_p4leftshift2(index);
% assmble position
M(N+2,N+1) =   h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points+3)).*psi_1(1/2.*(quad.points+3)))./2 + ...
    N1_N1M1_1 + N1_N1M1_2 + N1_N1M1_3;

  
index= 4;
phi_1 = GD_basis_p4(index);
% trial function
index= 4;
psi_1 = GD_basis_p4(index);
% assmble position
M(N+2,N+2) =   h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points+3)).*psi_1(1/2.*(quad.points+3)))./2 + ...
    N1_N1_1 + N1_N1_2 + N1_N1_3;

 

index= 4;
phi_1 = GD_basis_p4(index);
% trial function
index= 3;
psi_1 = GD_basis_p4rightshift2(index);
% assmble position
M(N+2,N+3) =   h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points+3)).*psi_1(1/2.*(quad.points+3)))./2 + ...
    N1_N1P1_1 + N1_N1P1_2;

 

index= 4;
phi_1 = GD_basis_p4(index);
% trial function
index= 2;
psi_1 = GD_basis_p4rightshift4(index);
% assmble position
M(N+2,N+4) =   h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points+3)).*psi_1(1/2.*(quad.points+3)))./2 + ...
    N1_N1P2_1;

 

index= 4;
phi_1 = GD_basis_p4(index);
% trial function
index= 1;
psi_1 = GD_basis_p4rightshift6(index);
% assmble position
M(N+2,N+5) =   h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points+3)).*psi_1(1/2.*(quad.points+3)))./2;

 

% \phi_10
M(N+3,N-1) = M(N,N-4) ;
M(N+3,N) = M(N,N-3) ;


index= 3;
phi_1 = GD_basis_p4(index);
% trial function
index= 5;
psi_1 = GD_basis_p4leftshift4(index);
% assmble position
M(N+3,N+1) =   h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points-1)).*psi_1(1/2.*(quad.points-1)))./2 + ...
    N1_N1M2_1 + N1_N1M2_2;

 

index= 3;
phi_1 = GD_basis_p4(index);
% trial function
index= 4;
psi_1 = GD_basis_p4leftshift2(index);
% assmble position
M(N+3,N+2) =   h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points-1)).*psi_1(1/2.*(quad.points-1)))./2 + ...
    N1_N1M1_1 + N1_N1M1_2  ;

 
index= 3;
phi_1 = GD_basis_p4(index);
% trial function
index= 3;
psi_1 = GD_basis_p4(index);
% assmble position
M(N+3,N+3) =   h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points-1)).*psi_1(1/2.*(quad.points-1)))./2 + ...
    N1_N1_1 +  N1_N1_2 ;

 

index= 3;
phi_1 = GD_basis_p4(index);
% trial function
index= 2;
psi_1 = GD_basis_p4rightshift2(index);
% assmble position
M(N+3,N+4) =   h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points-1)).*psi_1(1/2.*(quad.points-1)))./2 + ...
    N1_N1P1_1   ;
 

index= 3;
phi_1 = GD_basis_p4(index);
% trial function
index= 1;
psi_1 = GD_basis_p4rightshift4(index);
% assmble position
M(N+3,N+5) =   h./2.* sum(quad.weights.*phi_1(1/2.*(quad.points-1)).*psi_1(1/2.*(quad.points-1)))./2;

%T = readtable('xiadate.txt');
%xia = table2array(T);

%for k=3:13
%   err =  max(abs(h.*xia(k,:) - M(k,:)));
%   if err>1e-15
%       disp('mass matrix is not corrcet!');
%   end
%end

end

%%% GD basis
function [phi] = GD_basis_p4(index)%,indexL,indexR

if index == 1
    phi = @(x) 1/384.*(x+8).*(x+6).*(x+4).*(x+2);
elseif index == 2
    phi = @(x) -1/96.*(x+6).*(x+4).*(x+2).*(x-2);
elseif index == 3
    phi = @(x) 1/64.*(x+4).*(x+2).*(x-2).*(x-4);
elseif index == 4
    phi = @(x)  -1/96.*(x+2).*(x-2).*(x-4).*(x-6);
elseif index == 5
    phi = @(x) 1/384.*(x-2).*(x-4).*(x-6).*(x-8);
end


end

function [phi] = GD_basis_p4leftshift2(index)%,indexL,indexR
f = GD_basis_p4(index);
phi =  @(x) f(x + 2);
end

function [phi] = GD_basis_p4leftshift4(index)%,indexL,indexR
f = GD_basis_p4(index);
phi =  @(x) f(x + 4);
end

function [phi] = GD_basis_p4leftshift6(index)%,indexL,indexR
f = GD_basis_p4(index);
phi =  @(x) f(x + 6);
end

function [phi] = GD_basis_p4leftshift8(index)%,indexL,indexR
f = GD_basis_p4(index);
phi =  @(x) f(x + 8);
end


function [phi] = GD_basis_p4rightshift2(index)%,indexL,indexR
f = GD_basis_p4(index);
phi =  @(x) f(x - 2);
end

function [phi] = GD_basis_p4rightshift4(index)%,indexL,indexR
f = GD_basis_p4(index);
phi =  @(x) f(x - 4);
end

function [phi] = GD_basis_p4rightshift6(index)%,indexL,indexR
f = GD_basis_p4(index);
phi =  @(x) f(x - 6);
end

function [phi] = GD_basis_p4rightshift8(index)%,indexL,indexR
f = GD_basis_p4(index);
phi =  @(x) f(x - 8);
end



