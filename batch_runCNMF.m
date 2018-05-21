scans = {'1Port_000_002','1Port_000_003','1Port_000_014','1Port_000_029','1Port_002_002'};
foldername = 'E:\\4054010.5\1Port';
outputfolder = foldername;

cd(foldername);
for scancell = scans
    scan = scancell{:}
    edge1 = 1; edge2 = 796;
%     copyfile(fullfile('X:\4058541.3\rawDat',strcat(scan,'.sbx')),strcat(scan,'.sbx'));
%     copyfile(fullfile('X:\4058541.3\rawDat',strcat(scan,'.mat')),strcat(scan,'.mat'));
    
    runCNMF;
    
    delete *.h5
end

