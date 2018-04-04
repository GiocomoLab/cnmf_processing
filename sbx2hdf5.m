function info = sbx2hdf5(fname)


gcp;
tmp = load([fname '.mat']);
info = tmp.info;
% find number of channels
if info.channels == 1
    info.nchan = 2; factor = 1;
else
    info.nchan = 1; factor = 2;
end

% number of planes
info.nplanes = max([1,length(info.otwave_um)]);

% number of frames
d = dir([fname '.sbx']);

info.nsamples = (info.sz(2) * info.recordsPerBuffer * 2 * info.nchan);   % bytes per record     
info.max_idx =  d.bytes/info.recordsPerBuffer/info.sz(2)*factor/4 -1;


% create h5 files to write 

if info.nplanes == 1
    shp = [info.sz(1), info.sz(2), info.max_idx+1];
elseif info.nplanes > 1
    shp = [info.sz(1), info.sz(2), info.nplanes, ceil((info.max_idx+1)/info.nplanes)];
end

h5create([fname '_green.h5'],'/mov',shp);
if info.nchan == 2
    h5create([fname '_red.h5'], '/mov',shp);
end


% tmp_dat = single(sbxread(fname,0,info.max_idx+1));
% if info.nplanes>1
%     for plane_ind = 1:info.nplanes
%         plane = tmp_dat(:,:,:,plane_ind:info.nplanes:k);
%         h5write([fname '_green.h5'],'/mov',squeeze(plane(1,:,:,:)),[1,1,plane_ind, 1],[shp(1) shp(2) 1 size(plane,4)]);
%         if info.nchan >1
%             h5write([fname '_red.h5'],'/mov',squeeze(plane(2,:,:,:)),[1,1,plane_ind, 1],[shp(1) shp(2) 1 size(plane,4)]);
%         end
%     end
% else
%     h5write([fname '_green.h5'],'/mov', squeeze(tmp_dat(1,:,:,:)), [1 1 1],[shp(1) shp(2) size(tmp_dat,4)]);
%     if info.nchan>1
%         h5write([fname '_red.h5'], '/mov',squeeze(tmp_dat(2,:,:,:)),[1 1 1], [shp(1) shp(2) size(tmp_dat,4)]);
%     end
% end


% write files in manageable chunks
step_size = 1000;
first_frames = 0:step_size:info.max_idx;
max_idx = info.max_idx; nplanes = info.nplanes;
nchan = info.nchan;
disp([fname '_green.h5'])
pwd
parfor i = 1:length(first_frames)
    disp(i)
    first_frame = first_frames(i);
    
    last_frame = min([first_frame+step_size max_idx+1]);
    k = last_frame-first_frame;
    tmp_dat = single(sbxread(fname,first_frame,k));
    
    
    if nplanes>1
        for plane_ind = 1:nplanes
            plane = tmp_dat(:,:,:,plane_ind:nplanes:k);
            h5write([fname '_green.h5'],'/mov',squeeze(plane(1,:,:,:)),[1,1,plane_ind, first_frame/nplanes+1],[shp(1) shp(2) 1 size(plane,4)]);
%             if nchan >1
%                 h5write([fname '_red.h5'],'/mov',squeeze(plane(2,:,:,:)),[1,1,plane_ind, first_frame/nplanes+1],[shp(1) shp(2) 1 size(plane,4)]);
%             end
        end
    else
        
        if ~isempty(tmp_dat)
            h5write([fname '_green.h5'],'/mov',squeeze(tmp_dat(1,:,:,:)),[1,1,first_frame+1],[shp(1) shp(2) size(tmp_dat,4)]);
%         if nchan >1
%             h5write([fname '_red.h5'],'/mov',squeeze(tmp_dat(2,:,:,:)),[1,1,first_frame/nplanes+1],[shp(1) shp(2) size(tmp_dat,4)]);
%         end
        end
    end
        
end
                
            
            
        
    

