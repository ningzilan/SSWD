clear;clc;
tab=readtable('..\code\Biase_densty.csv');
genename=tab{:,1}; % choose the top 1000 genes
data=tab(:,2:end);
data=table2array(data);
Y_optimal = K_means_3(data);
tab=table(genename,Y_optimal);
writetable(tab,'..\code\Biase_subcluster.txt');
 
