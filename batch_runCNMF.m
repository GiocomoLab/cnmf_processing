gcp; % start a parallel engine

serverDir = 'Z:\2P_data\2AFC_V3';

scanStruct = struct();
% scanStruct(1).mouse = '4054010.5';
% % scanStruct(1).scans = {'TT_to1/TT_to1_000_005','TT_to2/TT_to2_000_002','TT_to4/TT_to4_000_005',...
% %     'TT_train3/TT_train3_000_001'};
% scanStruct(1).scans = {'TT_to4/TT_to4_000_005',...
%     'TT_train3/TT_train3_000_001'};
% 
scanStruct(1).mouse = '4139190.1';
scanStruct(1).scans = {'TT_to12/TT_to12_000_001'};

% scanStruct(1).mouse = '4139190.3';
% scanStruct(3).scans = {'TT_train1/TT_train1_000_006','TT_train1/TT_train1_001_001',...
%     'TT_train3/TT_train3_000_002','TT_to2/TT_to2_000_002'};
% scanStruct(1).scans = {'TT_to10/TT_to10_000_003'};



% scanStruct(4).mouse = '4054011.1';
% scanStruct(4).scans = {'TT_to4/TT_to4_000_003'};


for i = 1:numel(scanStruct)
    mouse = scanStruct(i).mouse;
    cd(fullfile('E:\',mouse));
    for scanName = scanStruct(i).scans
        scan = scanName{:}
        
        %if ~exist(fullfile('E:\',mouse,strcat(scan,'.sbx')),'file')
        [pathstr,name,ext]= fileparts(fullfile('E:\',mouse,strcat(scan,'.sbx')));
        if ~exist(pathstr,'file')
            mkdir(pathstr)
        end
        if ~exist(fullfile('E:\',mouse,strcat(scan,'.sbx')),'file')
            copyfile(fullfile(serverDir,mouse,strcat(scan,'.sbx')),fullfile('E:\',mouse,strcat(scan,'.sbx')))
        end
        if ~exist(fullfile('E:\',mouse,strcat(scan,'.mat')),'file')
            copyfile(fullfile(serverDir,mouse,strcat(scan,'.mat')),fullfile('E:\',mouse,strcat(scan,'.mat')))
        end
        %end
        
        foldername = fullfile('E:\',mouse);
        edge1 = 25;
        edge2 = 780;
            
        runCNMF;
        
        delete(strcat(scan,'_green_mc.h5'));
        delete(strcat(scan,'_green.h5'));
        delete(fullfile('E:\',mouse,strcat(scan,'.sbx')));
        
        
            
        
        
        
        
    end
end

  
  
