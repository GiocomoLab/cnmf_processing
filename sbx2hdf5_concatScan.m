function [infos, max_idx, fr] = sbx2hdf5_concatScan(fnames,output_name)


gcp;

% load each mat file and make sure the nchan is the same and fr is the same

% save nsamples for each scan and max idx
max_idx = 0; n_samples = 0; infos = cell(2,1);
for i  = 1:numel(fnames)
    fnames{i}
    fname = load([fnames{i} '.mat']);
    info = fname.info;
    % find number of channels
    if info.channels == 1
        info.nchan = 2; factor = 1;
    else
        info.nchan = 1; factor = 2;
    end
    
    % number of planes
    info.nplanes = max([1,length(info.otparam)-1]);
    
    if ~isempty(info.otparam)
        fr = info.resfreq/info.recordsPerBuffer/(length(info.otparam)-1);
    else
        fr = info.resfreq/info.recordsPerBuffer;
    end
    
    
    % number of frames
    d = dir([fnames{i} '.sbx']);
    
    info.nsamples = (info.sz(2) * info.recordsPerBuffer * 2 * info.nchan);   % bytes per record
    info.max_idx =  d.bytes/info.recordsPerBuffer/info.sz(2)*factor/4 -1;
    max_idx = max_idx + info.max_idx;
    n_samples = n_samples + info.nsamples;
    infos{i} = info;
end



% create h5 files to write


if infos{1}.nplanes == 1
    shp = [infos{1}.sz(1), infos{1}.sz(2), max_idx+1];
elseif infos{1}.nplanes > 1
    shp = [info{1}.sz(1)*info{1}.nplanes, info{1}.sz(2), ceil((max_idx+1)/info.nplanes)];
end

shp
h5create([output_name '_green.h5'],'/mov',shp);
if info.nchan == 2
    h5create([output_name '_red.h5'], '/mov',shp);
end


% write files in manageable chunks
step_size = 1000;

tmp_max_idx = 0;
for f = 1:numel(fnames)
    fname = fnames{f}; info = infos{f};
    tmp_min_idx = tmp_max_idx; tmp_max_idx = tmp_max_idx+info.max_idx;
    
    h5_first_frames = tmp_min_idx:step_size:tmp_max_idx;
    sbx_first_frames = 0:step_size:info.max_idx;
    for i = 1:length(sbx_first_frames)
        disp(i)
        sbx_first_frame = sbx_first_frames(i); h5_first_frame = h5_first_frames;
        last_frame = min([sbx_first_frame+step_size info.max_idx+1]);
        k = last_frame-sbx_first_frame;
        tmp_dat = single(sbxread(fname,sbx_first_frame,k));
        
        
        if info.nplanes>1
            for plane_ind = 1:nplanes
                disp('in loop');
                plane = tmp_dat(:,:,:,plane_ind:nplanes:k);
                plane = squeeze(plane(1,:,:,:));
                plnsz = size(plane);
                
                h5write([output_name '_green.h5'],'/mov',plane,[(plane_ind-1)*plnsz(1)+1,1, h5_first_frame/nplanes+1],[plnsz(1) plnsz(2) plnsz(3)]);
                
            end
        else
            
            if ~isempty(tmp_dat)
%                 h5_first_frame+1
%                 shp(1)
%                 shp(2)
%                 size(tmp_dat,4)
                h5write([output_name '_green.h5'],'/mov',squeeze(tmp_dat(1,:,:,:)),[1,1,h5_first_frame(i)+1],[shp(1) shp(2) size(tmp_dat,4)]);
            end
        end
        
    end
    
    
end



