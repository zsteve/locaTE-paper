clear all
close all
addpath(genpath('./'))%Add the necessary path

path_name=['./']; 
data_matrix = dlmread([path_name 'datamatrix.csv'],',')';
pseudotime_array = dlmread([path_name 'dpt.csv'],',');
A = dlmread([path_name 'J.csv'],',');%The real binary graph from the literature

t_array=pseudotime_array;%Choose the time label (real or pseudo) 
X=[t_array,data_matrix'];

[~,I]=sort(X(:,1));
X=X(I,:);%X is the spacetime matrix of the data (of size Cx(1+G)), with the
%chosen (real or pseudo) time in the first column and the gene expression in the others 
%% TESTING GRISLI
Alpha=@(Kx,Dt,sigx,sigt)exp(-(Kx.^2)/(2*sigx^2)).*exp(-(Dt.^2)/(2*sigt^2)).*(Dt.^2);%The kernel we use
R=100;
L_array=[1 5 10 25 50 75 100] % 70;
alpha_min=.3;
saveResultsBool=true;
nbtry=1;

saveFileName='output_GRISLI_Gvelo.txt';
[AUROC_GRISLI, AUPR_GRISLI, elapsedTime, TPR_array_area, FPR_array_area,...
    PPV_array_area,A_app_array_Rnk0]=Test_GRISLI_realdata(A,X,L_array,Alpha,R,alpha_min,saveResultsBool, saveFileName,false);%,V0,nbtry
% ROCPlot(FPR_array_area,TPR_array_area,"GRISLI");
for i=1:length(L_array)
    csvwrite(sprintf('A_grisli_L_%d.csv', L_array(i)), A_app_array_Rnk0(:, :, i))
end
