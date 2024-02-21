clc
clear all
addpath Auxilary
savepath
format shortG
% Demonstration of Profiled Estimation of Differential Equations
odefn    =@(t,X) gigrfunode(t,X,pars) ;      % Function for ODE solver (exact)
% RHS Functions
fn.fn       = @gigrfun  ;       % RHS function
% Now derivatives of the function with respect to system components and  parameters 
fn.dfdx     = @gigrdfdx;      % Derivative wrt inputs (Jacobian)
fn.dfdp     = @gigrdfdp;      % Derviative wrt parameters
% Now we need functions to compute all three sets of second derivatives:
fn.d2fdx2   = @gigrd2fdx2;    % Hessian wrt inputs
fn.d2fdxdp  = @gigrd2fdxdp;   % Hessian wrt inputs and parameters
fn.d2fdp2   = @gigrd2fdp2;    % Hessian wrt parameters.    

global Ess Gss Iss HGP QG QI QE  

SubjectNumber = '31b400';

filename1 = sprintf('Subject%s.xlsx',SubjectNumber);
[num,txt,raw]  = xlsread(filename1);
filename2 = sprintf('Subject%s_HGP.xlsx',SubjectNumber); 
[num1,txt1,raw1]  = xlsread(filename2);
 HepaticGP = num1(:,6)*180/4.44*1e-6*60; %(micromol/min)*(180 microg/micromol)*(1e-6g/microg)/4.44L = g/hour/L
 HepaticGD = num1(:,7)*180/4.44*1e-6*60; %(micromol/min)*(180 microg/micromol)*(1e-6g/microg)/4.44L = g/hour/L
  
  HGP = nanmean(HepaticGP(1:3))

%% 
L=size(num,1);
time = num([1:L],4)/60;  % hours
E    = num([1:L],6)/3482.747*1000 ;  %PicoMol/L  V_E=19.6
I    = num([1:L],7);    %mU/L  V_I=1.52
G    = num([1:L],8)/4.44;%*1000;    %becomes g/L
R    = zeros(L,1);
RE   = zeros(L,1);

full_time=[180/60;time([6:L])];

 
  % and parameters
  L1=length(time(time<3))-1;
 Ess = nanmean(E([2:L1])); Iss = nanmean(I([2:L1])); Gss = nanmean(G([2:L1]));
 
QG= 0.00005*80/4.44*60;%*1000; %g/L/hour
QE = 3*80/19.6/3482.747*1000*60 ; %E PicoMolar/L/hour
QI = 4*2/1.52*60 ;   % mU/L/hour

DegI = 4*2/1.52/(nanmean(I([7:end])))*60 ;%  r(2,1) =   4*2/1.52   - P(3)*X(2);                                          %I
DegE = 3*80/19.6/3482.747*1000/(nanmean(E([7:end])))*60;% r(3,1) =   QE  -P(4)* X(3)  ;   

startpars=[   3.8
        .4*HGP
        4
        20
       DegI
       DegE
       .3*60
    0.004
   1e-4]

           

     wts=2*[1  .1 .1 .1 .1]; 
lambda =  1*[2800  ,  10, 10, 1200, 900]    ;% [1e7 ,  10000, 1000, 1e9]    ;
lambda0 = 0.1;
  
full_path=[[Gss;G([6:L])] [Iss;I([6:L])] [Ess;E([6:L])] R([5:L]) RE([5:L])];
 P = startpars;
 pars = startpars;
 Vii = P(2);
 mu  = P(3);
 Kid = P(4);
 krec = 0.003*60;
 kin = P(7);
 K1=P(8);
 P0=P(9);
 
kon  = 6e-5*60; %1/pmol/hour
koff = 0.24*60;
Kii  = 0; kpin=0;
HC=2;
 Rss = krec/((krec+kpin)+(krec+kin)/(koff+kin)*mu*kon*Ess);
 REss = kon*mu*Ess*Rss/(koff+kin)
 
 if  1.2*HGP - (P0+ P(1)* REss^HC/(K1^2+REss^HC))<0
      1.2*HGP - (P0+ P(1)* REss^HC/(K1^2+REss^HC))
     msg = 'Error occurred.';
error(msg)
 end
  Vid = (HGP - Vii*Gss/(Kii+Gss))*((Kid+Gss)/(Gss*Iss))
 

%%
% Now we need to specify the times at which we observe the system.

tspan = full_time;
  
obs_pts = cell(1,5);

tG=[];
tI=[];
tE=[];

     for i=1:length(full_time)
        if ~isnan(full_path(i,1))
            tG = [tG ,i];
        end
     end
                 
     for i=1:length(full_time)
        if ~isnan(full_path(i,2))
            tI = [tI ,i];
        end
     end    
            
     for i=1:length(full_time)
        if ~isnan(full_path(i,3))
            tE = [tE ,i];
        end
     end      
    
 

obs_pts{1} = tG;       
obs_pts{2} = tI;  
obs_pts{3} = tE;

% Finally, we will want to be able to plot what a true solution looks like
% on a fairly fine grid. This specifies that grid. 

tfine = full_time(1):0.1:full_time(end);     

% Calculate trajectories
%
% We use the MATLAB routine |ode45| to solve the differential equation.

% First we need to set up a convergence tolerance for the numerical
% solution:

odeopts = odeset('RelTol',1e-13);



 
% Set up observations
%
% Finally, we set up MATLAB objects for the observations. These will be the
% objects |Tcell| and |Ycell| which will be cell arrays
% containing the observation times and observation values respectively.
% Each element of the cell array corresponds to one component of the
% GIGR system. 

% We start by defining cell arrays:

Tcell = cell(1,size(full_path,2));
path_cell = Tcell;

% We take the data from the solution of the differential equation and put
% it into the appropriate component. 

for i = 1:length(obs_pts)
    Tcell{i} = full_time(obs_pts{i});
    path_cell{i} = full_path(obs_pts{i},i);
end

% Finally, we add random observational noise to the 'path' variable. 

Ycell = path_cell;                
for i = 1:length(path_cell)
    Ycell{i} = path_cell{i} ;
end

 
       

% Now, we need to define some meta-parameters of a B-spline basis. First of
% these is the number of knots:

nknots =51;%1*length(tspan);    

% Then the order of B-spline that we will employ:

norder =5;

% Finally, since we will need to evaluate a non-linear penaly, we need the
% number of quadrature points between knots that will be used to perform
% numerical integration:

nquad =6;     

% Profiling optimisation control
%
% Here, we set up some control parameters for performing non-linear
% optimization. All optimization is doen by a Gauss-Newton method. However,
% there are two distinct levels of optimization. These control values take
% the form of the options in the MATLAB nonlinear optimization toolbox. 
%
% In all cases, the options should be set to use a Jacobian. I have also
% specified some tolerances and the 'Display' variable. 

% Firstly, there is the outer-optimization of the structural parameters.
% We need to be less concerned about very fine convergence, but it is
% appropriate to display progress every iteration.  

lsopts_out = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','iter','MaxIter',1000,'TolFun',1e-14,'TolX',1e-14);

% Then, there is the inner optimization loop to perform a smooth with a
% non-linear penalty. Here it is important to have tight convergence
% tolerances, but displaying the progress of the optimization will make the
% output messy. The Gauss-Newton optimization can be speeded up using a
% set of sparse matrix multiplication routines which is specified by
% setting 'JacobMult' to '@SparseJMfun'. 

lsopts_in = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','off','MaxIter',1000,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun);

% Finally, sometimes we want to just do the smoothing. The following
% options are the same as for the inner optimization, but the 'Display'
% option is set to output a summary when the routine finishes.  

lsopts_other = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','off','MaxIter',1000,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun);


% Setting up Functional Data Objects

% Firstly, we need to produce a basis function for each component of the
% system. To begin with, they all need to have the same range specified:

range = [min(full_time),max(full_time)];

% Now we create a cell-array containing the knots for B-spline bases. In
% this case, each basis will contain ?? equally spaced knots. Note
% however, that we could have different bases. 

knots_cell = cell(1,size(Ycell,2));
knots_cell(:) = {linspace(range(1),range(2),nknots)};
% myKNOTS=[linspace(range(1),(range(1)+range(2))/2,floor(nknots/4)),...
%  linspace((range(1)+range(2))/2,range(2),nknots-floor(nknots/4))];
% knots_cell(:) = {myKNOTS};

% The bases are also contained in a cell array:

basis_cell = cell(1,length(Ycell)); 

% At the same time, we will create a cell array of Lfd objects. This will
% be used to perform an initial smooth of the data:

Lfd_cell = cell(1,length(Ycell));

% We will also need to calculate the number of basis functions for each
% component of the system:

 
% Finally, quadrature points will have to be common to all bases in the
% system. The following code creates 'bigknots', a large vector
% containing all the knots for all the bases, it also calculates the number
% of basis functions for each component:  

bigknots = knots_cell{1};               
nbasis(1) = length(knots_cell{1}) + norder - 2;          
 
for i = 2:length(Ycell)
    bigknots = [bigknots knots_cell{i}];
    nbasis(i) =  length(knots_cell{i}) + norder -2;
end

% We can now use 'bigknots' to create a set of quadrature points.
% 'quadvals' contains Simpson's Rule quadrature points and quadrature
% weights based on 'nquad' quadrature points between each successive
% knot value. 

quadvals = MakeQuadPoints(bigknots,nquad);   

% Now we can create the basis and Lfd objects to popluate 'basis_cell'
% and 'Lfd_cell'. In this case 'MakeBasis' simply creates a B-spline basis
% and attaches 'quadvals' as quadrature points. We have chosen 'Lfd_cell'
% to contain objects penalizing the first derivative of a smooth.    

for i = 1:length(Ycell)
    basis_cell{i} = MakeBasis(range,nbasis(i),norder,...  
        knots_cell{i},quadvals,1);                        
    Lfd_cell{i} = fdPar(basis_cell{i},1,lambda0);         
end

% Smooth the data

% As a first step, we create a smooth of the data without reference to the
% differential equation. 'smoothfd_cell' does this for each component
% individually using the penalty specified in 'Lfd_cell'. 

 
% First we start with a non-parametric smooth of the data. 

DEfd = smoothfd_cell(Ycell,Tcell,Lfd_cell);

coefs = getcellcoefs(DEfd);
 
 %%
ind =[4 5];
% pcoefs = kron(ones(nbasis(ind),1),1);
coefs = [.6*ones(nbasis(4),1); .005*ones(nbasis(5),1)];
coefs1 = lsqnonlin(@SplineCoefErr_DEfit,coefs,[],[],lsopts_other,DEfd,ind,fn,pars,[],[]);

DEfd = update_fdcell(coefs1,ind,DEfd);



coefs = getcellcoefs(DEfd); % This is what we'll need to feed into model-based
                            % smooths. 
 

%% Re-smoothing with model-based penalty
%
% Now we get into some of the grunt work of the method. First, we will
% smooth the data using the differential equation as a penalty, but with
% the jittered parameters. 

% This is done by a call to the MATLAB optimizer 'lsqnonlin' which
% gives out the optimized coefficients of the basis. 

% options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','FunctionTolerance',1e-12);

lambda1= [1e8 ,  1000, 1000, 1e5, 1e5] ;
[newcoefs,resnorm2] = lsqnonlin(@SplineCoefErr,coefs,[],[],lsopts_in,...
        basis_cell,Ycell,Tcell,wts,lambda1,fn,[],pars,[]);
    
    

tDEfd = Make_fdcell(newcoefs,basis_cell);
%  size(newcoefs)

%% plot results along with exact solution
figure(3)
clf
devals = eval_fdcell(tfine,tDEfd,0);
% temp=cell2mat(devals(4));
% devals{4}=devals{4}-(temp(1,1)-Rss)



for i = 1:length(Ycell)
    subplot(length(Ycell),1,i);
    plot(tfine,devals{i},'r','LineWidth',2);
    hold on;
    plot(Tcell{i},Ycell{i},'b.');
    hold off
    if i==1
        ylabel('\fontsize{13} G')
        title(['\fontsize{13} Raw data (.), ', ...
               'exact solution (r-) and true path (g-)'])
    elseif i==2
        xlabel('\fontsize{13} t')
        ylabel('\fontsize{13} I')
    elseif i==3
        xlabel('\fontsize{13} t')
        ylabel('\fontsize{13} E')        
    else
        ylabel('\fontsize{13} R')
    end
end

drawnow
 
%% Perform the Profiled Estimation
%
% |Profile_GausNewt| runs the Guass-Newton iteration for the outer
% optimization in profiling. It outputs the new parameter estimates along
% with a cell-array of functional data objects that give the model-based
% smooth to the data 
%

%wts=[10 .10 1 1]
%  [newpars,newDEfd_cell] = Profile_GausNewt(pars,lsopts_out,tDEfd,fn,...
%     lambda,Ycell,Tcell,wts,[],lsopts_in,[],[],[],[]);
%   [newpars,newDEfd_cell] = My_Profile_GausNewt(startpars,lsopts_out,tDEfd,fn,...
%      lambda,Ycell,Tcell,wts,[],lsopts_in,[],@pen,@dpen,[]);
 
%   lambda=1e6*[1 1 1 1000]
%  lambda =wts
startpars';
  [newpars,newDEfd_cell] = My_Profile_GausNewt(startpars,lsopts_out,tDEfd,fn,...
     lambda,Ycell,Tcell,wts,[],lsopts_in,[],[],[],[]);
                             
  finalcoefs = getcellcoefs(newDEfd_cell)  ;
  size(finalcoefs)
disp(['Initial parameter values: '])
 (startpars)

disp(['New parameter estimates: ']);
 (newpars)




%% plot smooth with profile-estimated parameters

% newpars=[   3.8
%         .4*HGP
%         4
%         20
%        DegI
%        DegE
%        .3*60
%     0.004
%    1e-4]


   P = newpars;
 Vii = P(2);
 mu  = P(3);
 Kid = P(4);
 kin = P(7);
 K1=P(8);
 P0=P(9);
 
 HC = 2;
 Rss = krec/((krec+kpin)+(krec+kin)/(koff+kin)*mu*kon*Ess);
 REss = kon*mu*Ess*Rss/(koff+kin)
  1.2*HGP - (P0+ P(1)* REss^HC/(K1^2+REss^HC))
 Vid = (HGP - Vii*Gss/(Kii+Gss))*((Kid+Gss)/(Gss*Iss))
 y0=[Gss Iss Ess Rss REss];
 [t,X] = ode15s(@gigrfunode,tspan(1):.01:tspan(end),y0,odeopts,newpars);
 

figure(4)
clf
% devals = eval_fdcell(tfine,newDEfd_cell,0);
Data = [G,I,E];

for i = 1:4
    subplot(2,2,i)
    if i<4
%     plot(tfine,devals{i},'r-','LineWidth',2);
    end
    hold on;
    plot(t,X(:,i),'g-','LineWidth',2);
    set(gca, 'FontSize', 14);

     
    if i==1
         plot(time,Data(:,i),'b*');
         xlabel('  t (hour)')
         ylabel('  G (g/L)')
         str=sprintf('Subject %s',SubjectNumber);
         title(str) 
         plot([-5/60, 180/60 ],[Gss Gss],'g-','LineWidth',2)    
              axis([-10/60  360/60  0 4 ])
              plot([180/60 ,180/60 ],[ 0 4 ], 'k-.')
%                   legend('Profiled','ODE','Data','location','northwest')

    elseif i==2
            plot(time,Data(:,i),'b*');
            xlabel('  t (hour)')
            ylabel('  I (mU/L)')
            axis([-10/60  360/60  0 30])
            plot([180/60 ,180/60 ],[ 0 30], 'k-.')
            plot([-5/60 , 180/60 ],[Iss Iss],'g-','LineWidth',2)

    elseif i==3
            plot(time,Data(:,i),'b*');
            xlabel('  t (hour)')
            ylabel('  E (PicoMolar/L)')   
            axis([-10/60  360/60  0 100])
            plot([180/60 ,180/60 ],[ 0 100], 'k-.')
            plot([-5/60 , 180/60 ],[Ess Ess],'g-','LineWidth',2)          

    else
        ylabel('R/Rtot')
        axis([-10/60  360/60  0 1])
    end
end

figure(4)
subplot(2,2,4)
RE =X(:,5); % Eh*R/Kd
Ri = 1 - RE-X(:,4);
plot(t, 10*  RE,'b','LineWidth',2);
plot(t,   Ri,'r','LineWidth',2);
% legend('R','10* RE',' Ri','Location','southwest')
        plot([180/60 ,180/60 ],[ 0 100], 'k-.')
                plot([-5/60 , 180/60 ],[Rss Rss],'g-','LineWidth',2)          


drawnow
 %


figure(5)
clf
plot(t,P0+P(1)*X(:,5).^HC./(K1^2+X(:,5).^HC),'b','LineWidth',2)
hold on
plot(t,  Vii*X(:,1)./(Kii+X(:,1)) + Vid*X(:,1).*X(:,2)./(Kid+X(:,1)),'r','LineWidth',2)
plot(time(end-12:end),HepaticGP(end-12:end),'bo-','linewidth',1)
plot(time(end-12:end),HepaticGD(end-12:end),'rs-','linewidth',1)

% legend('1','2','3')
 axis([0 6 0 12])
%%
str = sprintf('Full_model_Result_Subject%s.mat',SubjectNumber);
save(str)

%% Squared Error Performance

% Squared error for estimated parameters

newpreds = eval_fdcell(Tcell,newDEfd_cell,0);
new_err = cell(length(newpreds));
for i = 1:4
    new_err{i} = wts(i)*(newpreds{i} - Ycell{i}).^2;
end

new_err = mean(cell2mat(new_err))

% Squared error for true parameters
