% complete pipeline for calcium imaging data pre-processing
% function convert_sbx_motion_correct_run_cnmf(foldername,scan)
gcp; % start a parallel engine
% raw_data_folder
% output_folder
%% convert sbx to hdf5

cd(output_folder); 
try
    mkdir(output_folder,'motionCorrect');
    cd('motionCorrect');
catch
    cd('motionCorrect');
end
h5fname = fullfile(raw_data_folder,'motionCorrect',strcat(scan,'_green.h5'));
if ~exist(h5fname,'file')
    info = sbx2hdf5(scan);
else
    info = load([scan '.mat']);
    info = info.info;
end
FOV = size(read_file(h5fname,1,1));

%% motion correct (and save registered h5 files as 2d matrices (to be used in the end)..)
% register files one by one. use template obtained from file n to
% initialize template of file n + 1; 

options_mc = NoRMCorreSetParms('d1',FOV(1),'d2',FOV(2),'grid_size',[128,128],'init_batch',200,...
                'overlap_pre',64,'mot_uf',4,'bin_width',200,'max_shift',24,'max_dev',8,'us_fac',50,...
                'output_type','h5');
%%

append = '_mc';
options_mc.h5_filename = fullfile(foldername,[h5fname(1:end-3),append,'.h5']);
if ~exist(options_mc.h5_filename,'file')
    [M,shifts,template,options_mc,col_shift] = normcorre_batch(h5fname,options_mc,[]);
    save(fullfile(foldername,[scan,'_green_shifts',append,'.mat']),'shifts','-v7.3');           
    % save shifts of each file at the respective folder  
end


%% downsample h5 files and save into a single memory mapped matlab file
    

fr = info.resfreq/info.recordsPerBuffer; % frame rate
batch_size = 2000;                    % read chunks of that size

%% now run CNMF on patches on the downsampled file, set parameters first

data = options_mc.h5_filename;
h= h5info(data);
sizY= h.Datasets(2).Dataspace.Size;
patch_size = [40,40];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [8,8];                        % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = 10;                                            % number of components to be found
tau = 10;                                          % std of gaussian kernel (half size of neuron) 
p = 2;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.7;                                  % merging threshold


options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'deconv_method','constrained_foopsi',...    % neural activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'ssub',12,...                                % spatial downsampling when processing
    'tsub',1,...                                % further temporal downsampling when processing
    'merge_thr',merge_thr,...                   % merging threshold
    'max_size_thr',300,'min_size_thr',20,...    % max/min acceptable size for each component
    'spatial_method','regularized',...          % method for updating spatial components
    'df_prctile',50,...                         % take the median of background fluorescence to compute baseline fluorescence 
    'fr',fr,...                            % downsamples
    'space_thresh',0.7,...                      % space correlation acceptance threshold
    'min_SNR',5.0,...                           % trace SNR acceptance threshold
    'nb',1,...                                  % number of background components per patch
    'gnb',3,...                                 % number of global background components
    'cnn_thr',0.4,...
    'decay_time',1.0);                       % length of typical transient for the indicator used



%% Run on patches (the main work is done here)

[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,0,options);  % do not perform deconvolution here since
                                                                             % we are operating on downsampled data

%% classify components
rval_space = classify_comp_corr(data,A,C,b,f,options);
ind_corr = rval_space > options.space_thresh;           % components that pass the correlation test
                                        % this test will keep processes
                                        
%% further classification with cnn_classifier
try  % matlab 2017b or later is needed
    [ind_cnn,value] = cnn_classifier(A,sizY(1:2),'cnn_model',options.cnn_thr);
catch
    disp("cnn failed");
    ind_cnn = true(size(A,2),1);                         % components that pass the CNN classifier
end     
                            
%% event exceptionality

fitness = compute_event_exceptionality(C+YrA,options.N_samples_exc,options.robust_std);
ind_exc = (fitness < options.min_fitness);

%% select components

keep = (ind_corr | ind_cnn) & ind_exc;


%% keep only the active components    
A_tmp = A(:,keep);
C_tmp = C(keep,:);
R_tmp = YrA(keep,:);
%clear A C YrA
A= []; C=[]; YrA=[];
% % % kill deadband ROIs
cm = com(A_tmp,sizY(1),sizY(2));
cm_keep = cm(:,2)>75 & cm(:,2)<720;

A_keep = A_tmp(:,cm_keep);
C_keep = C_tmp(cm_keep,:);
R_keep = R_tmp(cm_keep,:);

%clear A_tmp C_tmp R_tmp
A_tmp = []; C_tmp = []; R_tmp =[];


%% save results
svfile = matfile(strcat(scan,'_cnmf_results_pre.mat'),'Writable',true);
svfile.options = options; svfile.A=A_keep; svfile.C=C_keep;
svfile.R=R_keep; svfile.b = b; svfile.f = f;

template_sc = template/max(max(template));
template_sc(template_sc>.25) = .25;
svfile.template = 4*template_sc; 
