 function d2fdxdp = gigrd2fdxdp(t,fd_cell,P) 
 global Ess Gss Iss HGP QG QI QE 
 d2fdxdp = cell(5,5,9); 

  X =  (eval_fdcell(t,fd_cell,0)); 
 d2fdxdp{1,1,1}=0 ; 
d2fdxdp{1,1,2}=(X{2}.*(Gss + P(4)))./(Gss.*Iss.*(P(4) + X{1})) - (X{1}.*X{2}.*(Gss + P(4)))./(Gss.*Iss.*(P(4) + X{1}).^2) ; 
d2fdxdp{1,1,3}=0 ; 
d2fdxdp{1,1,4}=(X{1}.*X{2}.*(HGP - P(2)))./(Gss.*Iss.*(P(4) + X{1}).^2) - (X{2}.*(HGP - P(2)))./(Gss.*Iss.*(P(4) + X{1})) + (X{2}.*(HGP - P(2)).*(Gss + P(4)))./(Gss.*Iss.*(P(4) + X{1}).^2) - (2.*X{1}.*X{2}.*(HGP - P(2)).*(Gss + P(4)))./(Gss.*Iss.*(P(4) + X{1}).^3) ; 
d2fdxdp{1,1,5}=0 ; 
d2fdxdp{1,1,6}=0 ; 
d2fdxdp{1,1,7}=0 ; 
d2fdxdp{1,1,8}=0 ; 
d2fdxdp{1,1,9}=0 ; 
d2fdxdp{1,2,1}=0 ; 
d2fdxdp{1,2,2}=(X{1}.*(Gss + P(4)))./(Gss.*Iss.*(P(4) + X{1})) ; 
d2fdxdp{1,2,3}=0 ; 
d2fdxdp{1,2,4}=(X{1}.*(HGP - P(2)).*(Gss + P(4)))./(Gss.*Iss.*(P(4) + X{1}).^2) - (X{1}.*(HGP - P(2)))./(Gss.*Iss.*(P(4) + X{1})) ; 
d2fdxdp{1,2,5}=0 ; 
d2fdxdp{1,2,6}=0 ; 
d2fdxdp{1,2,7}=0 ; 
d2fdxdp{1,2,8}=0 ; 
d2fdxdp{1,2,9}=0 ; 
d2fdxdp{1,3,1}=0 ; 
d2fdxdp{1,3,2}=0 ; 
d2fdxdp{1,3,3}=0 ; 
d2fdxdp{1,3,4}=0 ; 
d2fdxdp{1,3,5}=0 ; 
d2fdxdp{1,3,6}=0 ; 
d2fdxdp{1,3,7}=0 ; 
d2fdxdp{1,3,8}=0 ; 
d2fdxdp{1,3,9}=0 ; 
d2fdxdp{1,4,1}=0 ; 
d2fdxdp{1,4,2}=0 ; 
d2fdxdp{1,4,3}=0 ; 
d2fdxdp{1,4,4}=0 ; 
d2fdxdp{1,4,5}=0 ; 
d2fdxdp{1,4,6}=0 ; 
d2fdxdp{1,4,7}=0 ; 
d2fdxdp{1,4,8}=0 ; 
d2fdxdp{1,4,9}=0 ; 
d2fdxdp{1,5,1}=(2.*X{5})./(P(8).^2 + X{5}.^2) - (2.*X{5}.^3)./(P(8).^2 + X{5}.^2).^2 ; 
d2fdxdp{1,5,2}=0 ; 
d2fdxdp{1,5,3}=0 ; 
d2fdxdp{1,5,4}=0 ; 
d2fdxdp{1,5,5}=0 ; 
d2fdxdp{1,5,6}=0 ; 
d2fdxdp{1,5,7}=0 ; 
d2fdxdp{1,5,8}=(8.*P(1).*P(8).*X{5}.^3)./(P(8).^2 + X{5}.^2).^3 - (4.*P(1).*P(8).*X{5})./(P(8).^2 + X{5}.^2).^2 ; 
d2fdxdp{1,5,9}=0 ; 
d2fdxdp{2,1,1}=0 ; 
d2fdxdp{2,1,2}=0 ; 
d2fdxdp{2,1,3}=0 ; 
d2fdxdp{2,1,4}=0 ; 
d2fdxdp{2,1,5}=0 ; 
d2fdxdp{2,1,6}=0 ; 
d2fdxdp{2,1,7}=0 ; 
d2fdxdp{2,1,8}=0 ; 
d2fdxdp{2,1,9}=0 ; 
d2fdxdp{2,2,1}=0 ; 
d2fdxdp{2,2,2}=0 ; 
d2fdxdp{2,2,3}=0 ; 
d2fdxdp{2,2,4}=0 ; 
d2fdxdp{2,2,5}=-1 ; 
d2fdxdp{2,2,6}=0 ; 
d2fdxdp{2,2,7}=0 ; 
d2fdxdp{2,2,8}=0 ; 
d2fdxdp{2,2,9}=0 ; 
d2fdxdp{2,3,1}=0 ; 
d2fdxdp{2,3,2}=0 ; 
d2fdxdp{2,3,3}=0 ; 
d2fdxdp{2,3,4}=0 ; 
d2fdxdp{2,3,5}=0 ; 
d2fdxdp{2,3,6}=0 ; 
d2fdxdp{2,3,7}=0 ; 
d2fdxdp{2,3,8}=0 ; 
d2fdxdp{2,3,9}=0 ; 
d2fdxdp{2,4,1}=0 ; 
d2fdxdp{2,4,2}=0 ; 
d2fdxdp{2,4,3}=0 ; 
d2fdxdp{2,4,4}=0 ; 
d2fdxdp{2,4,5}=0 ; 
d2fdxdp{2,4,6}=0 ; 
d2fdxdp{2,4,7}=0 ; 
d2fdxdp{2,4,8}=0 ; 
d2fdxdp{2,4,9}=0 ; 
d2fdxdp{2,5,1}=0 ; 
d2fdxdp{2,5,2}=0 ; 
d2fdxdp{2,5,3}=0 ; 
d2fdxdp{2,5,4}=0 ; 
d2fdxdp{2,5,5}=0 ; 
d2fdxdp{2,5,6}=0 ; 
d2fdxdp{2,5,7}=0 ; 
d2fdxdp{2,5,8}=0 ; 
d2fdxdp{2,5,9}=0 ; 
d2fdxdp{3,1,1}=0 ; 
d2fdxdp{3,1,2}=0 ; 
d2fdxdp{3,1,3}=0 ; 
d2fdxdp{3,1,4}=0 ; 
d2fdxdp{3,1,5}=0 ; 
d2fdxdp{3,1,6}=0 ; 
d2fdxdp{3,1,7}=0 ; 
d2fdxdp{3,1,8}=0 ; 
d2fdxdp{3,1,9}=0 ; 
d2fdxdp{3,2,1}=0 ; 
d2fdxdp{3,2,2}=0 ; 
d2fdxdp{3,2,3}=0 ; 
d2fdxdp{3,2,4}=0 ; 
d2fdxdp{3,2,5}=0 ; 
d2fdxdp{3,2,6}=0 ; 
d2fdxdp{3,2,7}=0 ; 
d2fdxdp{3,2,8}=0 ; 
d2fdxdp{3,2,9}=0 ; 
d2fdxdp{3,3,1}=0 ; 
d2fdxdp{3,3,2}=0 ; 
d2fdxdp{3,3,3}=0 ; 
d2fdxdp{3,3,4}=0 ; 
d2fdxdp{3,3,5}=0 ; 
d2fdxdp{3,3,6}=-1 ; 
d2fdxdp{3,3,7}=0 ; 
d2fdxdp{3,3,8}=0 ; 
d2fdxdp{3,3,9}=0 ; 
d2fdxdp{3,4,1}=0 ; 
d2fdxdp{3,4,2}=0 ; 
d2fdxdp{3,4,3}=0 ; 
d2fdxdp{3,4,4}=0 ; 
d2fdxdp{3,4,5}=0 ; 
d2fdxdp{3,4,6}=0 ; 
d2fdxdp{3,4,7}=0 ; 
d2fdxdp{3,4,8}=0 ; 
d2fdxdp{3,4,9}=0 ; 
d2fdxdp{3,5,1}=0 ; 
d2fdxdp{3,5,2}=0 ; 
d2fdxdp{3,5,3}=0 ; 
d2fdxdp{3,5,4}=0 ; 
d2fdxdp{3,5,5}=0 ; 
d2fdxdp{3,5,6}=0 ; 
d2fdxdp{3,5,7}=0 ; 
d2fdxdp{3,5,8}=0 ; 
d2fdxdp{3,5,9}=0 ; 
d2fdxdp{4,1,1}=0 ; 
d2fdxdp{4,1,2}=0 ; 
d2fdxdp{4,1,3}=0 ; 
d2fdxdp{4,1,4}=0 ; 
d2fdxdp{4,1,5}=0 ; 
d2fdxdp{4,1,6}=0 ; 
d2fdxdp{4,1,7}=0 ; 
d2fdxdp{4,1,8}=0 ; 
d2fdxdp{4,1,9}=0 ; 
d2fdxdp{4,2,1}=0 ; 
d2fdxdp{4,2,2}=0 ; 
d2fdxdp{4,2,3}=0 ; 
d2fdxdp{4,2,4}=0 ; 
d2fdxdp{4,2,5}=0 ; 
d2fdxdp{4,2,6}=0 ; 
d2fdxdp{4,2,7}=0 ; 
d2fdxdp{4,2,8}=0 ; 
d2fdxdp{4,2,9}=0 ; 
d2fdxdp{4,3,1}=0 ; 
d2fdxdp{4,3,2}=0 ; 
d2fdxdp{4,3,3}=-(9.*X{4})./2500 ; 
d2fdxdp{4,3,4}=0 ; 
d2fdxdp{4,3,5}=0 ; 
d2fdxdp{4,3,6}=0 ; 
d2fdxdp{4,3,7}=0 ; 
d2fdxdp{4,3,8}=0 ; 
d2fdxdp{4,3,9}=0 ; 
d2fdxdp{4,4,1}=0 ; 
d2fdxdp{4,4,2}=0 ; 
d2fdxdp{4,4,3}=-(9.*X{3})./2500 ; 
d2fdxdp{4,4,4}=0 ; 
d2fdxdp{4,4,5}=0 ; 
d2fdxdp{4,4,6}=0 ; 
d2fdxdp{4,4,7}=0 ; 
d2fdxdp{4,4,8}=0 ; 
d2fdxdp{4,4,9}=0 ; 
d2fdxdp{4,5,1}=0 ; 
d2fdxdp{4,5,2}=0 ; 
d2fdxdp{4,5,3}=0 ; 
d2fdxdp{4,5,4}=0 ; 
d2fdxdp{4,5,5}=0 ; 
d2fdxdp{4,5,6}=0 ; 
d2fdxdp{4,5,7}=0 ; 
d2fdxdp{4,5,8}=0 ; 
d2fdxdp{4,5,9}=0 ; 
d2fdxdp{5,1,1}=0 ; 
d2fdxdp{5,1,2}=0 ; 
d2fdxdp{5,1,3}=0 ; 
d2fdxdp{5,1,4}=0 ; 
d2fdxdp{5,1,5}=0 ; 
d2fdxdp{5,1,6}=0 ; 
d2fdxdp{5,1,7}=0 ; 
d2fdxdp{5,1,8}=0 ; 
d2fdxdp{5,1,9}=0 ; 
d2fdxdp{5,2,1}=0 ; 
d2fdxdp{5,2,2}=0 ; 
d2fdxdp{5,2,3}=0 ; 
d2fdxdp{5,2,4}=0 ; 
d2fdxdp{5,2,5}=0 ; 
d2fdxdp{5,2,6}=0 ; 
d2fdxdp{5,2,7}=0 ; 
d2fdxdp{5,2,8}=0 ; 
d2fdxdp{5,2,9}=0 ; 
d2fdxdp{5,3,1}=0 ; 
d2fdxdp{5,3,2}=0 ; 
d2fdxdp{5,3,3}=(9.*X{4})./2500 ; 
d2fdxdp{5,3,4}=0 ; 
d2fdxdp{5,3,5}=0 ; 
d2fdxdp{5,3,6}=0 ; 
d2fdxdp{5,3,7}=0 ; 
d2fdxdp{5,3,8}=0 ; 
d2fdxdp{5,3,9}=0 ; 
d2fdxdp{5,4,1}=0 ; 
d2fdxdp{5,4,2}=0 ; 
d2fdxdp{5,4,3}=(9.*X{3})./2500 ; 
d2fdxdp{5,4,4}=0 ; 
d2fdxdp{5,4,5}=0 ; 
d2fdxdp{5,4,6}=0 ; 
d2fdxdp{5,4,7}=0 ; 
d2fdxdp{5,4,8}=0 ; 
d2fdxdp{5,4,9}=0 ; 
d2fdxdp{5,5,1}=0 ; 
d2fdxdp{5,5,2}=0 ; 
d2fdxdp{5,5,3}=0 ; 
d2fdxdp{5,5,4}=0 ; 
d2fdxdp{5,5,5}=0 ; 
d2fdxdp{5,5,6}=0 ; 
d2fdxdp{5,5,7}=-1 ; 
d2fdxdp{5,5,8}=0 ; 
d2fdxdp{5,5,9}=0 ; 
end