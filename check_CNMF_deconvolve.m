load(fullfile(outputfolder,strcat(scan,'_cnmf_results_pre.mat')));

template_sc = template/max(max(template));
template_sc(template_sc>.2) = .2;

%
cm = com(A_keep,options.d1,options.d2,options.d3);
if exist('edge1','var') && exist('edge2','var')
    cm_keep = cm(:,2)>edge1 & cm(:,2)<edge2;
    A_keep = A_keep(:,cm_keep);
    C_keep = C_keep(cm_keep,:);
    C_dec = C_dec(cm_keep,:);
    R_keep = R_keep(cm_keep,:);
    cm = cm(cm_keep,:);
end

    

[A_final, C_final] = viewNeurons_patches(A_keep,C_keep,b,f,5*template_sc,cm,options) ; 

%% extract fluorescence on native temporal resolution

N = size(C_final,1);                             % total number of components
T = size(C_final,2);                                    % total number of timesteps
% F = C_final + R;                       % full fluorescence
S = zeros(N,T);


%% extract DF/F and deconvolve DF/F traces

[F_dff,F0] = detrend_df_f(A_final,b,C_final,f,[],options);
options.p=2;
p = 2;
%%
C_dec = zeros(N,T);         % deconvolved DF/F traces
S_dec = zeros(N,T);         % deconvolved neural activity
bl = zeros(N,1);            % baseline for each trace (should be close to zero since traces are DF/F)
neuron_sn = zeros(N,1);     % noise level at each trace
g = cell(N,1);              % discrete time constants for each trace
if p == 1; model_ar = 'ar1'; elseif p == 2; model_ar = 'ar2'; else; error('This order of dynamics is not supported'); end
 
for i = 1:N
    spkmin = options.spk_SNR*GetSn(F_dff(i,:));
    lam = choose_lambda(exp(-1/(options.fr*options.decay_time)),GetSn(F_dff(i,:)),options.lam_pr);
    [cc,spk,opts_oasis] = deconvolveCa(F_dff(i,:),model_ar,'method','thresholded','optimize_pars',true,'maxIter',20,...
                                'window',150,'lambda',lam,'smin',spkmin);
    bl(i) = opts_oasis.b;
    C_dec(i,:) = cc(:)' + bl(i);
    S_dec(i,:) = spk(:);
    neuron_sn(i) = opts_oasis.sn;
    g{i} = opts_oasis.pars(:)';
    disp(['Performing deconvolution. Trace ',num2str(i),' out of ',num2str(N),' finished processing.'])
end


%%

svfile = matfile(strcat(scan,'_cnmf_results.mat'),'Writable',true);
svfile.options = options; svfile.A=A_final; svfile.C=C_final;
svfile.C_dec = C_dec; svfile.S_dec=S_dec;
svfile.F_dff = F_dff; svfile.F0 = F0;
svfile.b = b; svfile.f = f;
svfile.template = template; 
