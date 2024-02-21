clc
clear all
format shortG
DataBefore=...
    [7.2	1.1417	3.9993	20.002	36.565	4.4162	18.001	0.0057993	0	;...
3.7999	0.78835	4.9919	15.001	38.172	4.9611	12.007	0.0069925	0	;...
4.1094	0.83083	4.497	11.352	43.412	5.0656	12.183	0.0058401	0	;...
4.1001	0.69452	5.9997	5.0001	30.471	4.0489	9.0007	0.0099991	0	;...
7.218	1.1414	6.2155	2.4197	26.928	5.6045	15.441	0.0070713	0	;...
8.7975	0.35076	2.8855	16.24	20.084	9.7802	27.139	0.0031891	0	;...
4.9009	1	9.9996	20.031	26.721	6.8956	17.998	0.0045001	0	;...
6.3778	0.97329	1.6301	1.8176	21.491	5.0594	56.409	0.0011275	4.00E-08	];
Vid= [ 3.3495    1.3577    1.4150    0.7582    0.5578    2.0569    2.7105    0.4462]';
DataBefore = [DataBefore Vid]

DataAfter =[...
6.4852	0.52993	3.9619	4.0678	38.53	4.5422	9.5907	0.0099973	0.54396	;...
4.4797	1.0518	3.0747	11.076	37.701	4.7787	18.181	0.0050748	0.41361	;...
3.5491	0.60377	3.5315	10.455	43.94	4.2161	9.8511	0.0041815	0.10265	;...
7.8496	0.90667	5.0001	5.0011	28.241	5.0072	12.001	0.014998	0	;...
6.5004	1.5363	6.001	2.9994	26.593	5.6757	12.005	0.0089956	0	;...
9.9989	0.98766	5.0999	17.988	22.267	7.7188	12.008	0.01199	2.60E-07	;...
4.9	0.89079	6.3	18	25.355	8.3402	13.8	0.0072	0	;...
8.1001	0.93301	7.9989	9.9795	24.812	5.7891	8.963	0.018996	1.60E-07];	

Vid=[ 1.5136    1.6047    2.2563    0.7691    0.4817    2.4782    2.0534    1.2301]';
DataAfter = [DataAfter Vid]

MEANBefore = mean(DataBefore)
STDBefore = std(DataBefore)

MEANAfter = mean(DataAfter)
STDAfter  = std(DataAfter)

P=[];
for i=1:size(DataAfter,2)
p = signrank(DataBefore(:,i),DataAfter(:,i));
 
P=[P; p];
end
P'
 