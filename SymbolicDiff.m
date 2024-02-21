clc
clear all

X = sym('X', [5 1]) ;
P = sym('P', [9 1]) ;
 syms Ess Gss Iss HGP QG QI QE HGP QG QI QE  

syms t SmoothStep
%Glucagon Molar mass: 3482.747 g/mol
 
kon  = 6e-5*60; %1/pmol/hour
koff = 0.24*60;
 

 
 Vii = P(2);
 mu  =  P(3);
 Kid = P(4);
 DegI = P(5);
 DegE = P(6);
  kin = P(7);
 K1 = P(8);
 P0 = P(9);
 
 krec = 0.003*60;
 kpin = 0;
 Kii = 0;
 Rss = krec/( (krec+kpin)+(krec+kin)/(koff+kin)*mu*kon*Ess );
 REss = kon*mu*Ess*Rss/(koff+kin);
 Eh = mu*X(3) ;
 Vid = (HGP - Vii*Gss/(Kii+Gss))*((Kid+Gss)/(Gss*Iss));
 
 HC=2;
 
r(1,1) =  QG +  P0 + P(1)*( X(5).^HC./(K1^HC+ X(5).^HC) ) ...
              - Vii*X(1)/(Kii+X(1))  - Vid*X(1).*X(2)/(Kid+X(1))  ;
r(2,1) =   QI  - DegI*X(2);                                          %I
r(3,1) =   QE  - DegE* X(3)  ;                  %E PicoMolar/L/min
r(4,1) =   -kon*Eh.*X(4) +koff*X(5) -kpin *X(4) + krec*(1-X(4)-X(5)); % picomolar
r(5,1) =    kon*Eh.*X(4) -koff*X(5) - kin*X(5); 
 

 
 %%
% ODE function for solver
FID=fopen('gigrfunode.m','w');
    fprintf(FID, '%s\n\n',' function r = gigrfunode(t,X,P)');
        fprintf(FID, '%s\n',' global Ess Gss Iss HGP QG QI QE ');

        for i=1:length(r)
        fprintf(FID, ' r(%d,1) = %s; \n\n ',i, (r(i)));
        end
        
    fprintf(FID, '%s\n','end');
        
        
 FID=fopen('gigrfun.m','w');
    fprintf(FID, '%s\n\n',' function r = gigrfun(t,fd_cell,P)');
     fprintf(FID, '%s\n',' global Ess Gss Iss HGP QG QI QE ');
     fprintf(FID, '%s\n\n','X = eval_fdcell(t,fd_cell,0);');

        for i=1:length(r)
        fprintf(FID, ' r{%d} = %s; \n\n ',i, (r(i)));
        end
        
    fprintf(FID, '%s\n','end');
 
 

% Jacobian dfdx

dfdx = cell(length(X),length(X));
FID=fopen('gigrdfdx.m','w');
    fprintf(FID, '%s\n\n',' function dfdx = gigrdfdx(t,fd_cell,P)');
            fprintf(FID, '%s\n',' global Ess Gss Iss HGP QG QI QE ');
    fprintf(FID, ' dfdx = cell(%d,%d); \n\n ',length(X),length(X));
    fprintf(FID, ' X =  (eval_fdcell(t,fd_cell,0)); \n ');
 
 
    

for i=1:length(r)
    for j=1:length(X)
        dfdx{i,j} = diff(r(i),X(j));
        fprintf(FID, 'dfdx{%d,%d}=%s ; \n',i,j,  (dfdx{i,j}));
    end
end
    fprintf(FID, '%s\n','end');
    
  
% dfdp
dfdp = cell(length(X),length(P));
FID=fopen('gigrdfdp.m','w');
    fprintf(FID, '%s\n',' function dfdp = gigrdfdp(t,fd_cell,P) ');
    fprintf(FID, '%s\n',' global Ess Gss Iss HGP QG QI QE ');
    fprintf(FID, ' dfdp = cell(%d,%d); \n\n ',length(X),length(P));
    fprintf(FID, ' X =  (eval_fdcell(t,fd_cell,0)); \n ');


for i=1:length(r)
    for j=1:length(P)
        dfdp{i,j} = diff(r(i),P(j));
        fprintf(FID, 'dfdp{%d,%d}=%s ; \n',i,j,  (dfdp{i,j}));
    end
end
    fprintf(FID, '%s\n','end');

%d2fdx2
d2fdx2 = cell(length(X),length(X),length(X));
FID=fopen('gigrd2fdx2.m','w');
    fprintf(FID, '%s\n',' function d2fdx2 = gigrd2fdx2(t,fd_cell,P) ');
    fprintf(FID, '%s\n',' global Ess Gss Iss HGP QG QI QE ');
    fprintf(FID, ' d2fdx2 = cell(%d,%d,%d); \n\n ',length(X),length(X),length(X));
    fprintf(FID, ' X =  (eval_fdcell(t,fd_cell,0)); \n ');

for i=1:length(r)
    for j=1:length(X)
        for k=1:length(X)
        d2fdx2{i,j,k} = diff(r(i),X(j),X(k));
        fprintf(FID, 'd2fdx2{%d,%d,%d}=%s ; \n',i,j,k,  (d2fdx2{i,j,k}));
        end
    end
end
    fprintf(FID, '%s\n','end');
 

%d2fdxdp
d2fdxdp = cell(length(X),length(X),length(P));
FID=fopen('gigrd2fdxdp.m','w');
    fprintf(FID, '%s\n',' function d2fdxdp = gigrd2fdxdp(t,fd_cell,P) ');
    fprintf(FID, '%s\n',' global Ess Gss Iss HGP QG QI QE ');
    fprintf(FID, ' d2fdxdp = cell(%d,%d,%d); \n\n ',length(X),length(X),length(P));
    fprintf(FID, ' X =  (eval_fdcell(t,fd_cell,0)); \n ');


for i=1:length(r)
    for j=1:length(X)
        for k=1:length(P)
        d2fdxdp{i,j,k} = diff(r(i),X(j),P(k));
        fprintf(FID, 'd2fdxdp{%d,%d,%d}=%s ; \n',i,j,k,  (d2fdxdp{i,j,k}));
        end
    end
end
    fprintf(FID, '%s\n','end');

 
%d2fdp2
d2fdp2 = cell(length(X),length(P),length(P));
FID=fopen('gigrd2fdp2.m','w');
    fprintf(FID, '%s\n',' function d2fdp2 = gigrd2fdp2(t,fd_cell,P) ');
    fprintf(FID, '%s\n',' global Ess Gss Iss HGP QG QI QE ');
    fprintf(FID, ' d2fdp2 = cell(%d,%d,%d); \n\n ',length(X),length(P),length(P));
    fprintf(FID, ' X = (eval_fdcell(t,fd_cell,0)); \n ');


for i=1:length(r)
    for j=1:length(P)
        for k=1:length(P)
        d2fdp2{i,j,k} = diff(r(i),P(j),P(k));
        fprintf(FID, 'd2fdp2{%d,%d,%d}=%s ; \n',i,j,k,  (d2fdp2{i,j,k}));
        end
    end
end
    fprintf(FID, '%s\n','end');

    
    
        for i=1:length(X)
        s_old = sprintf('X%d',i); 
        s_new = sprintf('X{%d}',i);
        find_and_replace('gigrdfdx.m',s_old,s_new)
        find_and_replace('gigrdfdp.m',s_old,s_new)
        find_and_replace('gigrd2fdx2.m',s_old,s_new)
        find_and_replace('gigrd2fdxdp.m',s_old,s_new)
        find_and_replace('gigrd2fdp2.m',s_old,s_new)
        find_and_replace('gigrfun.m',s_old,s_new)       
        
        s1_old = sprintf('X%d',i); 
        s1_new = sprintf('X(%d)',i);
                find_and_replace('gigrfunode.m',s1_old,s1_new)

        
        
        end
        %%
    for index=1:(length(P))
        
        s_old = sprintf('P%d',(index)) ;
        s_new = sprintf('P(%d)',(index));
        find_and_replace('gigrdfdx.m',s_old,s_new)
        find_and_replace('gigrdfdp.m',s_old,s_new)
        find_and_replace('gigrd2fdx2.m',s_old,s_new)
        find_and_replace('gigrd2fdxdp.m',s_old,s_new)
        find_and_replace('gigrd2fdp2.m',s_old,s_new)
        find_and_replace('gigrfunode.m',s_old,s_new)
        find_and_replace('gigrfun.m',s_old,s_new)
        
    end
    
    %%
%         find_and_replace('gigrdfdx.m','SmoothStep','SmoothStep(t)')
%         find_and_replace('gigrdfdp.m','SmoothStep','SmoothStep(t)')
%         find_and_replace('gigrd2fdx2.m','SmoothStep','SmoothStep(t)')
%         find_and_replace('gigrd2fdxdp.m','SmoothStep','SmoothStep(t)')
%         find_and_replace('gigrd2fdp2.m','SmoothStep','SmoothStep(t)')
%         find_and_replace('gigrfunode.m','SmoothStep','SmoothStep(t)')
%         find_and_replace('gigrfun.m','SmoothStep','SmoothStep(t)')
% 
% str='(1-exp(-heaviside(t-180)*(t-180)))';
% 
%         find_and_replace('gigrdfdx.m','SmoothStep',str)
%         find_and_replace('gigrdfdp.m','SmoothStep',str)
%         find_and_replace('gigrd2fdx2.m','SmoothStep',str)
%         find_and_replace('gigrd2fdxdp.m','SmoothStep',str)
%         find_and_replace('gigrd2fdp2.m','SmoothStep',str)
%         find_and_replace('gigrfunode.m','SmoothStep',str)
%         find_and_replace('gigrfun.m','SmoothStep',str)
%  
 
 
        find_and_replace('gigrdfdx.m',1,1,1 )
        find_and_replace('gigrdfdp.m',1,1,1 )
        find_and_replace('gigrd2fdx2.m',1,1,1 )
        find_and_replace('gigrd2fdxdp.m',1,1,1 )
        find_and_replace('gigrd2fdp2.m',1,1,1 )
        find_and_replace('gigrfunode.m',1,1,1 )
        find_and_replace('gigrfun.m',1,1,1 )

