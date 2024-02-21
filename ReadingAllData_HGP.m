clc
clear all
  %  close all

filename = 'CHDRdata2.xlsx';
[num,txt,raw]  = xlsread(filename);
j=1;
ind=1;
index=[];
txt(1,:)=[]; % title of the table
SubjectNumber=[]
    
    for i=1:size(num,1)
         if num(i,4)<0 && isequal(char(txt(i,end)),'400 mg') ...
                 && isequal(char(txt(i,3)),'Before treatment')      
        %if num(i,4)<0 && isequal(char(txt(i,end)),'400 mg')

        index(j,1) = i;
        SubjectNumber(ind)=num(i,1);ind=ind+1;
        
         elseif  num(i,4)==360 && isequal(char(txt(i,end)),'400 mg') ...
                  && isequal(char(txt(i,3)),'Before treatment') 
         %elseif  num(i,4)==360 && isequal(char(txt(i,end)),'400 mg')
            
            index(j,2)=i;
            j=j+1;
        end

    end
    
 raw(1,:)=[]; % clearing the title row
    
    for k=  1: length(SubjectNumber)
        
%      time = num([index(2*k-1):index(2*k)],4);
%      E    = num([index(2*k-1):index(2*k)],6);
%      I    = num([index(2*k-1):index(2*k)],7);
%      G    = num([index(2*k-1):index(2*k)],8);
% 
filename=sprintf('Subject%db400.xlsx',SubjectNumber(k))
     
xlswrite(filename,raw([index(k,1):index(k,2)],:))
% writetable(T,filename,'Sheet',1,'Range','D1')
     
%      i=1;
%       while time(i)<=180
%          StableG1(i)=G(i);
%          i=i+1;
%      end
%      StableG1(isnan(StableG1),:)=[] ;
%      Gss1 = mean(StableG1);
%      stableG = G([end-5:end]);
%      stableG(isnan(stableG),:)=[];
%      Gss2 = mean(stableG);
%      ratioG(k)=Gss2/Gss1;
%      
%      
%           i=1;
%      while time(i)<=180
%          StableI1(i)=I(i);
%          i=i+1;
%      end
%      StableI1(isnan(StableI1))=[];
%      Iss1 = mean(StableI1)
%      stableI = I([end-8:end]);
%      stableI(isnan(stableI),:)=[];
%      Iss2 = mean(stableI)
%      ratioI(k)=Iss2/Iss1;
%      
%       i=1;
%      while time(i)<=180
%          StableE1(i)=E(i);
%          i=i+1;
%      end
%      StableE1(isnan(StableE1))=[];
%      Ess1 = mean(StableE1);
%      stableE = E([end-8:end]);
%      stableE(isnan(stableE),:)=[];
%      Ess2 = mean(stableE);
%      ratioE(k)=Ess2/Ess1;    
%      
%      
%  
% 
% figure(k)
% clf
%  subplot(1,3,1)
% scatter(time,E,'filled')
% title('Glucagon')
% set(gca, 'FontSize', 18);
% axis([-20 400 0 250])
% hold on
% plot([180,180],[ 0 250], 'k-.')
% plot([360,360],[ 0 250], 'k-.')
% plot([time(1), 180],[Ess1 Ess1],'r-')
% plot([time(end-8), 360],[Ess2 Ess2],'r-')
% 
% 
% 
% subplot(1,3,2)
% scatter(time,I,'filled')
% str= sprintf('Subject %d \n insulin',SubjectNumber(k))
% title(str)
% set(gca, 'FontSize', 18);
% axis([-20 400 0 50])
% hold on
% plot([180,180],[ 0 250], 'k-.')
% plot([360,360],[ 0 250], 'k-.')
% plot([time(1), 180],[Iss1 Iss1],'r-')
% plot([time(end-8), 360],[Iss2 Iss2],'r-')
% 
% 
% subplot(1,3,3)
% scatter(time,G,'filled')
% title('Glucose')
% set(gca, 'FontSize', 18);
% axis([-20 400 0 20])
% hold on
% plot([180,180],[ 0 250], 'k-.')
% plot([360,360],[ 0 250], 'k-.')
% plot([0, 180],[Gss1 Gss1],'r-')
% plot([time(end-5), 360],[Gss2 Gss2],'r-')
 
end  
  SubjectNumber
  
  
  
      %%
%     figure(25)
%     clf
%     hist(ratioG./ratioE)
%     figure(1)
%     clf 
%     sz=2*ones(k,1);
%     scatter([1:k],ratioE,50,'b','filled')  
%      hold on
%      scatter([1:k],ratioG,50,'r','filled')
%         legend('Es2/Es1','Gs2/Gs1')
%     grid on
    %grid minor
%     xticks([ 1:k])
%   xticklabels({  ' 6','',' 8','',' 11','',' 15','',' 18','',' 22','',' 23','',' 27','',' 29','',' 34','',' 39','',' 137',''})
% xticks([ 1:1:k])
 %    xticklabels({   '      S6','','     S8','','      S11','','      S15','','       S18','','       S22','','       S23','','       S27','','       S29','','       S34','','       S39','','       S137'})
  
%    xticks([ 1:1:k])
%    xticklabels({SubjectNumber})
%  
% axis([0 k+1 1 4])
%    xlabel('Subject Number')
%    ylabel('Ratio of steady states')
% ax = gca; % current axes
% ax.FontSize = 16;%     