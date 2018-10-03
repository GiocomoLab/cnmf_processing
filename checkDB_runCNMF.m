

% open database
dbname = fullfile( 'G:\My Drive\VR_Data\TwoTower\','behavior.sqlite');
conn = sqlite(dbname,'readonly');

% 1-MouseName 2-DateFolder 3-SessionNumber 4-Track 5-RewardCount 6-Imaging
data = fetch(conn,'SELECT MouseName, DateFolder, SessionNumber, Track, RewardCount, Imaging FROM sessions');


% find indices where imaging == 1
imaging = logical(cell2mat(data(:,6)));
data_imaging= data(imaging,:);

edge1 = 25;
edge2 = 780;


for s =  1:size(data_imaging,1)
    % check if [session]_cnmf_results_pre.mat exists
    scan_folder = fullfile('G:\My Drive\2P_Data\TwoTower',data_imaging{s,1},data_imaging{s,2},...
        data_imaging{s,4});
    cnmf_file_name = fullfile(scan_folder,strcat(data_imaging{s,4},'_*',num2str(data_imaging{s,3}),'_*','_cnmf_results_pre.mat'));
    match_files = dir(cnmf_file_name);
    % if it doesn't exist, runCNMF
    if isempty(match_files)
        % find scan file
        sbx_file_dir = dir(fullfile(scan_folder,strcat(data_imaging{s,4},'_*',num2str(data_imaging{s,3}),'_*','.sbx')));
%         if ~isempty(sbx_file_dir)
            if length(sbx_file_dir) == 1 % check how many sbx files there are, assuming one scan file per session
                sbx_file_name = sbx_file_dir.name;
            else
                disp('more than one sbx file for this session. taking the last modified one');
                file_datenums = [sbx_file_dir(:).datenum];
                [~,I]=max(file_datenums);
                sbx_file_name = sbx_file_dir(I).name;
            end
            
            % prep folders for running on local SSD
            [~,scan,~]= fileparts(sbx_file_name);
            foldername = fullfile('E:\',data_imaging{s,1},data_imaging{s,2},data_imaging{s,4});
            
            if ~exist(foldername,'dir'); mkdir(foldername); end
            
            copyfile(fullfile(scan_folder,sbx_file_name),fullfile(foldername,[scan '.sbx']));
            copyfile(fullfile(scan_folder,[scan '.mat']),fullfile(foldername,[scan '.mat']));
            
            runCNMF; % workhorse - need to work on functionalizing
            
            cd('Z:\scripts\cnmf_processing');
            
            % copy results
            copyfile(fullfile(foldername,[scan '_cnmf_results_pre.mat']),fullfile(scan_folder,[scan '_cnmf_results_pre.mat']));
            
            % clear space on SSD
            [status,message,message_id]=rmdir(foldername,'s');
%         end
    end
    
    
    
    
    
    
end

