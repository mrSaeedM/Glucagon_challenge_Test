 function r = gigrfun(t,fd_cell,P)

 global Ess Gss Iss HGP QG QI QE 
X = eval_fdcell(t,fd_cell,0);

 r{1} = P(9) - P(2) + QG + (P(1).*X{5}.^2)./(P(8).^2 + X{5}.^2) - (X{1}.*X{2}.*(HGP - P(2)).*(Gss + P(4)))./(Gss.*Iss.*(P(4) + X{1})); 

  r{2} = QI - P(5).*X{2}; 

  r{3} = QE - P(6).*X{3}; 

  r{4} = (711.*X{5})./50 - (9.*X{4})./50 - (9.*P(3).*X{3}.*X{4})./2500 + 9./50; 

  r{5} = (9.*P(3).*X{3}.*X{4})./2500 - P(7).*X{5} - (72.*X{5})./5; 

 end
