 function dfdx = gigrdfdx(t,fd_cell,P)

 global Ess Gss Iss HGP QG QI QE 
 dfdx = cell(5,5); 

  X =  (eval_fdcell(t,fd_cell,0)); 
 dfdx{1,1}=(X{1}.*X{2}.*(HGP - P(2)).*(Gss + P(4)))./(Gss.*Iss.*(P(4) + X{1}).^2) - (X{2}.*(HGP - P(2)).*(Gss + P(4)))./(Gss.*Iss.*(P(4) + X{1})) ; 
dfdx{1,2}=-(X{1}.*(HGP - P(2)).*(Gss + P(4)))./(Gss.*Iss.*(P(4) + X{1})) ; 
dfdx{1,3}=0 ; 
dfdx{1,4}=0 ; 
dfdx{1,5}=(2.*P(1).*X{5})./(P(8).^2 + X{5}.^2) - (2.*P(1).*X{5}.^3)./(P(8).^2 + X{5}.^2).^2 ; 
dfdx{2,1}=0 ; 
dfdx{2,2}=-P(5) ; 
dfdx{2,3}=0 ; 
dfdx{2,4}=0 ; 
dfdx{2,5}=0 ; 
dfdx{3,1}=0 ; 
dfdx{3,2}=0 ; 
dfdx{3,3}=-P(6) ; 
dfdx{3,4}=0 ; 
dfdx{3,5}=0 ; 
dfdx{4,1}=0 ; 
dfdx{4,2}=0 ; 
dfdx{4,3}=-(9.*P(3).*X{4})./2500 ; 
dfdx{4,4}=- (9.*P(3).*X{3})./2500 - 9./50 ; 
dfdx{4,5}=711./50 ; 
dfdx{5,1}=0 ; 
dfdx{5,2}=0 ; 
dfdx{5,3}=(9.*P(3).*X{4})./2500 ; 
dfdx{5,4}=(9.*P(3).*X{3})./2500 ; 
dfdx{5,5}=- P(7) - 72./5 ; 
end
