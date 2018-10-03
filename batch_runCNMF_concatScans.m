

% scanStruct(1).mouse = '4139190.1';
% scanStruct(1).date_folder = '15_06_2018';
% scanStruct(1).scan = {'TwoTower_noTimeout/TwoTower_noTimeout_1_000';
%                       'TwoTower_noTimeout/TwoTower_noTimeout_3_000'};
% scanStruct(1).output = 'TwoTower_noTimeout_1__TwoTower_noTimeout_3';
                  
scanStruct(2).mouse = '4139190.1';
scanStruct(2).date_folder = '16_06_2018';
scanStruct(2).scan = {'TwoTower_noTimeout/TwoTower_noTimeout_2_000';
                      'TwoTower_Timeout/TwoTower_Timeout_2_000'};
scanStruct(2).output = 'TwoTower_noTimeout_2__TwoTower_Timeout_2';
                  

scanStruct(3).mouse = '4139190.1';
scanStruct(3).date_folder = '11_07_2018';
scanStruct(3).scan = {'TwoTower_Timeout/TwoTower_Timeout_017_023';
                    'TwoTower_BlackWhite_noTimeout/TwoTower_BlackWhite_noTimeout_001_026'};
scanStruct(3).output = 'TwoTower_Timeout_17__TwoTower_BlackWhite_noTimeout_1';

scanStruct(4).mouse = '4139190.1';
scanStruct(4).date_folder = '21_06_2018';
scanStruct(4).scan = {'TwoTower_Timeout/TwoTower_Timeout_2_000';
                      'TwoTower_Timeout/TwoTower_Timeout_3_000'};
scanStruct(4).output = 'TwoTower_Timeout_2__TwoTower_Timeout_3';
                  
                  
scanStruct(5).mouse = '4139190.3';
scanStruct(5).date_folder = '12_07_2018';
scanStruct(5).scan = {'TwoTower_Timeout/TwoTower_Timeout_001_004';
                      'TwoTower_BlackWhite_noTimeout/TwoTower_BlackWhite_noTimeout_001_006'};
scanStruct(5).output = 'TwoTower_Timeout_1__TwoTower_BlackWhite_noTimeout_1';
                  
scanStruct(6).mouse = '4139190.3';
scanStruct(6).date_folder = '14_07_2018';
scanStruct(6).scan = {'TwoTower_BlackWhite_noTimeout/TwoTower_BlackWhite_noTimeout_002_003';
                      'TwoTower_Timeout/TwoTower_Timeout_016_014'};
scanStruct(6).output = 'TwoTower_BlackWhite_noTimeout_2__TwoTower_Timeout_16';
                  

scanStruct(7).mouse = '4139190.3';
scanStruct(7).date_folder = '15_06_2018';
scanStruct(7).scan = {'TwoTower_noTimeout/TwoTower_noTimeout_1_000';
                      'TwoTower_noTimeout/TwoTower_noTimeout_3_000'};
scanStruct(7).output = 'TwoTower_noTimeout_1__TwoTower_noTimeout_3';
                  
scanStruct(8).mouse = '4139202.2';
scanStruct(8).date_folder = '20_07_2018';
scanStruct(8).scan = {'TwoTower_noTimeout/TwoTower_noTimeout_001_002';
                      'TwoTower_noTimeout/TwoTower_noTimeout_002_004'};
scanStruct(8).output = 'TwoTower_noTimeout_1__TwoTower_noTimeout_2';
                 

edge1 = 25;
edge2 = 780;


for s =  2:numel(scanStruct)
    ss = scanStruct(s);
    % check if [session]_cnmf_results_pre.mat exists
    scan_folder = fullfile('G:\My Drive\2P_Data\TwoTower',ss.mouse,ss.date_folder);
    
    %copy date folder to SpeedProcess
    ssd_folder = fullfile('X:\',ss.mouse,ss.date_folder);
    mkdir(ssd_folder);
    copyfile(fullfile(scan_folder,'*'),ssd_folder);
    
    % run CNMF
    cd(ssd_folder)
    foldername = ssd_folder;
    scan = ss.scan;
    output_name = ss.output;
    runCNMF_concatScan; % workhorse - need to work on functionalizing
    
    % copy results to G: ... date folder
    copyfile(fullfile(foldername,[output_name '_cnmf_results_pre.mat']),fullfile(scan_folder,[output_name '_cnmf_results_pre.mat']));
    cd('G:\My Drive\2P_Data\TwoTower');
    % clear space on SSD
    [status,message,message_id]=rmdir(foldername,'s');
    
    
    
    
end

