
edit_components_GUI(F_dff,A_keep,C_keep,b,f,Cn,S_dec,options);
Coor = plot_contours(A_keep,Cn,options,1);


pause;
load tmp.mat


A_keep = A_keep(:,~reject);
C_keep = C_keep(~reject,:);
C_dec = C_dec(~reject,:);
F0 = F0(~reject,:);
F_dff = F_dff(~reject,:);
R_keep = R_keep(~reject,:);
S_dec = S_dec(~reject,:); 
S_keep = S_keep(~reject,:);

svfile = matfile(strcat(scan,'_cnmf_results.mat'),'Writable',true);
svfile.A_keep = A_keep;
svfile.C_keep = C_keep;
svfile.C_dec = C_dec;
svfile.F0 = F0;
svfile.F_dff = F_dff;
svfile.R_keep = R_keep;
svfile.S_dec = S_dec; 
svfile.S_keep = S_keep;
svfile.f = f;
svfile.b = b;
svfile.P = P;
% svfile.RESULTS = RESULTS;
svfile.bl = bl;

% save(strcat(scan,'_cnmf_results.mat'),...
%     'options','A_keep','C_keep','R_keep','S_dec','C_dec','b','f','P','RESULTS',...
%     'Cn','F_dff','F0','C_dec','bl')