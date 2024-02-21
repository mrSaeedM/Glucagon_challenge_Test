function Q=PtoQ(P)
global alpha
 

for i=1:length(P)
    Q(i)=1/alpha(i)*log(1+exp(alpha(i)*P(i)));
end