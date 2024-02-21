function out=SmoothStep(t)
alpha = .3;
if (t<180)
    out=0;
elseif t>=18 && t<360
    out = 1- exp(-alpha*(t-180));
else
    out= (1- exp(-alpha*(360-180)))*(exp(-alpha*(t-360)));
end

end