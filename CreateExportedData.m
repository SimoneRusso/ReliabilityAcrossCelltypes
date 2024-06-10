
clear 

% Last Noise 1 run is block number 3 - rename files in folders instead of running again

%MainFolder = 'C:\Users\simon\PycharmProjects\pythonProject1\allensdk\CellType\Data\cell_types\';
MainFolder = 'D:\AllenInstitute\CellType\cell_types\';
Table = readtable('C:\Users\simon\PycharmProjects\pythonProject1\allensdk\CellType\Data\cell_types\cell_types_specimen_details.csv');

SelectedSweepName = 'Noise 1';

% Rerun Noise 1 _ 1

%WindowOfInterest_Noise1_LowInt = [403965 1004010];
%WindowOfInterest_Noise1_MidInt = [2004010 2604010];
%WindowOfInterest_Noise1_HighInt = [3604010 4204010];
%WindowOfInterest_4604001 = WindowOfInterest_Noise1_MidInt;


%WindowOfInterest_Noise1_LowInt = [101003 251002];
%WindowOfInterest_Noise1_MidInt = [501005 651005];
%WindowOfInterest_Noise1_HighInt = [901005 1051010];
%WindowOfInterest_1401000_1401050 = WindowOfInterest_Noise1_MidInt;


%WindowOfInterest_Noise1_LowInt = [403979 1004010];
%WindowOfInterest_Noise1_MidInt = [2004010 2604010];
%WindowOfInterest_Noise1_HighInt = [3604010 4204010];
%WindowOfInterest_5004001_5604001_5604000_6004001 = WindowOfInterest_Noise1_MidInt;

BlockOfInterest = 1;

DirMainFolder = dir(MainFolder);
DirMainFolder(1:2) = [];
DirMainFolder = DirMainFolder([DirMainFolder.isdir]);

DirNames = regexprep({DirMainFolder.name},'specimen_','');

DirNames = cellfun(@(x) str2num(x), DirNames, 'UniformOutput', true);
%%.
AllOutputsTable = [];


%
%Noise 2 - file 170 wrong??
for counterCell = 1%:size(DirNames,2)%1:size(Table,1) --- 252
    try
        currSpecimenID = DirNames(counterCell);%Table.specimen__id(1);
        currFolderName = ['specimen_' num2str(currSpecimenID)];
        currNwbFile = fullfile(MainFolder,currFolderName,'ephys.nwb');
        currIndexTable = find(Table.specimen__id==currSpecimenID);
        
        %if Table(currIndexTable). %check if already loaded and in case exclude
        
        %CurrInfo = h5info(currNwbFile);
        
        %
        CurrInfo = h5info(currNwbFile,'/acquisition/timeseries');
        SweepsPaths = {CurrInfo.Groups.Name};
        
        SelectedSweeps = [];
        SelectedSweepData = [];
        for counterSweep = 1:size(SweepsPaths,2)
            currSweep = SweepsPaths{counterSweep};
            
            currSweepName = h5read(currNwbFile,[currSweep,'/aibs_stimulus_name']);
            
            if strcmp(currSweepName, SelectedSweepName)
                currSweepData = h5read(currNwbFile,[currSweep,'/data']);
                currSweepStim = h5read(currNwbFile,[regexprep(currSweep,'acquisition/timeseries','stimulus/presentation'),'/data']);
                
                tmpStartEndStimuli = find(diff(currSweepStim~=0));
                WindowOfInterest = [tmpStartEndStimuli(2*BlockOfInterest+1) tmpStartEndStimuli(2*BlockOfInterest+2)];
                
                
                Date = h5read(currNwbFile,'/file_create_date');
                if str2double(Date{1}(end-3:end))>=2016
                    sRate = 50000;
                else
                    sRate = 200000;
                end
                if ~exist(fullfile(MainFolder,currFolderName,'sRate.mat'))
                    save(fullfile(MainFolder,currFolderName,'sRate.mat'),'sRate')
                end
                
                %tmpStartStimuli = find(diff(currSweepStim>0)>0);
                %tmpEndStimuli = find(diff(currSweepStim>0)<0);
                %WindowOfInterest = [tmpStartStimuli(BlockOfInterest+1) tmpEndStimuli(BlockOfInterest+1)];
                
                
                %if size(currSweepData,1)==4604001
                %    WindowOfInterest = WindowOfInterest_4604001;
                %elseif any([size(currSweepData,1)==1401000,size(currSweepData,1)==1401050])
                %    WindowOfInterest = WindowOfInterest_1401000_1401050;
                %elseif any([size(currSweepData,1)==5604001,size(currSweepData,1)==6004001, size(currSweepData,1) == 5604000, size(currSweepData,1) == 5004001])
                %    WindowOfInterest = WindowOfInterest_5004001_5604001_5604000_6004001;
                %else
                %    'Error? Unexpected window of interest?'
                %    Error
                %end
                
                SelectedSweepData(end+1,:) = currSweepData(WindowOfInterest(1):WindowOfInterest(2));
                SelectedSweeps(end+1) = counterSweep;
            end
            
        end
        
        if size(SelectedSweeps,2)>1
            [Output,SpikesTs] = CTAllen_AnalysisCycle(SelectedSweepData,sRate);
            
            %save(fullfile(MainFolder,currFolderName,'SpikesTs.mat'),'SpikesTs')
            save(fullfile(MainFolder,currFolderName,['SpikesTs_' SelectedSweepName '_' num2str(BlockOfInterest) '.mat']),'SpikesTs')
            
            %Output.SessID = currSpecimenID;
            Output = struct2table(Output);
            
            if isempty(AllOutputsTable)
                AllOutputsTable = Output;
                Table(:,end+1:end+size(Output,2))=array2table(nan([size(Table,1) size(Output,2)]));
                
                tmpFields = fields(Output);
                Table.Properties.VariableNames(end+1-size(Output,2):end) = tmpFields(1:end-3);
                %Table(currIndexTable,end+1:end+size(Output,2)) = Output;
                %Table.Properties.VariableNames(end+1)
            else
                AllOutputsTable = [AllOutputsTable;Output];
                Table(currIndexTable,end+1-size(Output,2):end) = Output;
            end
            %plot(SelectedSweepData')
            %title(num2str(Output.PercMatchingSpikes))
            %pause()
            %cla
        end
        
    catch
        'Error'
    end
    counterCell
end
