function P=QtoP(Q)
global alpha
for i=1:length(Q)
        P(i)= log(exp(alpha(i)*Q(i))-1)/alpha(i);
end
