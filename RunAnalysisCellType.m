
clear 

%MainFolder = 'C:\Users\simon\PycharmProjects\pythonProject1\allensdk\CellType\Data\cell_types\';
MainFolder = 'D:\AllenInstitute\CellType\cell_types\';
Table = readtable('C:\Users\simon\PycharmProjects\pythonProject1\allensdk\CellType\Data\cell_types\cell_types_specimen_details.csv');

SelectedSweepName = 'Noise 1';


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

BlockOfInterest = 2;

DirMainFolder = dir(MainFolder);
DirMainFolder(1:2) = [];
DirMainFolder = DirMainFolder([DirMainFolder.isdir]);

DirNames = regexprep({DirMainFolder.name},'specimen_','');

DirNames = cellfun(@(x) str2num(x), DirNames, 'UniformOutput', true);
%%.
AllOutputsTable = [];


%
%Noise 2 - file 170 wrong??
All_Spikes = cell(1,size(Table,1));
for counterCell = 1%:size(DirNames,2)%1:size(Table,1) --- 252
    %try=
        currSpecimenID = DirNames(counterCell);%Table.specimen__id(1);
        currFolderName = ['specimen_' num2str(currSpecimenID)];
        currNwbFile = fullfile(MainFolder,currFolderName,'ephys.nwb');
        currIndexTable = find(Table.specimen__id==currSpecimenID);
        
        
        CurrInfo = h5info(currNwbFile,'/acquisition/timeseries');
        SweepsPaths = {CurrInfo.Groups.Name};
        SweepsPaths_Spikes = regexprep(SweepsPaths,'/acquisition/timeseries','/analysis/spike_times');
        
        counterSweep = 1;
        while counterSweep>0
            try
            currSweep_StartTime = h5read(currNwbFile,[regexprep(SweepsPaths{counterSweep},'acquisition/timeseries','epochs'),'/start_time']);
            currSweep_StopTime = h5read(currNwbFile,[regexprep(SweepsPaths{counterSweep},'acquisition/timeseries','epochs'),'/stop_time']);
            currSweep_nSamples = h5read(currNwbFile,[regexprep(SweepsPaths{counterSweep},'acquisition/timeseries','epochs'),'/stimulus/timeseries/num_samples']);
            counterSweep = 0;
            catch
                counterSweep = counterSweep+1;
            end
        end
        
        
        ApproxSamplingFreq = round(1/((currSweep_StopTime-currSweep_StartTime)/double(currSweep_nSamples)));
        
        if ApproxSamplingFreq>100000
            sRate = 200000;
        else
            sRate = 50000;
        end
        
        SelectedSweeps = [];
        SelectedSweepData = [];
        SelectedSweepData_Spikes = {};
        SelectedSweepData_AllSpikes = {};
        for counterSweep = 1:size(SweepsPaths,2)
            currSweep = SweepsPaths{counterSweep};
            currSweep_Spikes = SweepsPaths_Spikes{counterSweep};
            
            currSweepName = h5read(currNwbFile,[currSweep,'/aibs_stimulus_name']);
            
            if strcmp(currSweepName, SelectedSweepName)
                
                
                %currSweepData = h5read(currNwbFile,[currSweep,'/data']);
                currSweepStim = h5read(currNwbFile,[regexprep(currSweep,'acquisition/timeseries','stimulus/presentation'),'/data']);
                currSweepData_Spikes = h5read(currNwbFile,currSweep_Spikes)*sRate;
                
                
                
                tmpStartEndStimuli = find(diff(currSweepStim~=0));
                WindowOfInterest = [tmpStartEndStimuli(2*BlockOfInterest+1) tmpStartEndStimuli(2*BlockOfInterest+2)];
                
                
                %figure, hold on
                %plot(currSweepStim)
                %scatter(round(currSweepData_Spikes),ones(size(currSweepData_Spikes))/100000000000)
                %scatter(WindowOfInterest,ones(size(WindowOfInterest))/100000000000)
                
                
                
                currSweepData_SelectedSpikes = currSweepData_Spikes(find((currSweepData_Spikes>WindowOfInterest(1)).*(currSweepData_Spikes<WindowOfInterest(2))));
                
                SelectedSweepData_Spikes{end+1} = currSweepData_SelectedSpikes/sRate; 
                SelectedSweepData_AllSpikes{end+1} = currSweepData_Spikes/sRate; 
                
                SelectedSweeps(end+1) = counterSweep;
            end
            
        end
        
        All_Spikes{currIndexTable}= SelectedSweepData_AllSpikes;
        
        
        if size(SelectedSweeps,2)>1
            
            %[Output] = CTAllen_AnalysisCycle_Fast(SelectedSweepData_Spikes);
            
            
            %save(fullfile(MainFolder,currFolderName,'SpikesTs.mat'),'SpikesTs')
            %save(fullfile(MainFolder,currFolderName,['SelectedSweepDataSpikes_' SelectedSweepName '_' num2str(BlockOfInterest) '.mat']),'SelectedSweepData_Spikes')
            %save(fullfile(MainFolder,currFolderName,['SelectedSweeps_' SelectedSweepName '_' num2str(BlockOfInterest) '.mat']),'SelectedSweeps')
            
            Output.PercMatchingSpikes = nanmean(nanmean(Output.PercMatchingSpikes));
            Output.nSpikes = nanmean(nanmean(Output.nSpikes));
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
        
    %catch
    %    'Error'
    %end
    counterCell
end




%save(['Table_' SelectedSweepName '_' num2str(BlockOfInterest)],'Table')
%save(['AllSpikes_' SelectedSweepName],'All_Spikes')

%%


IndexesNonNan = find((not(isnan(Table.nSpikes))).*(not((Table.nSpikes==0))).*(not(isnan(Table.PercMatchingSpikes))).*(not((Table.PercMatchingSpikes==0))));

figure
scatter(((Table.nSpikes(IndexesNonNan))),(Table.PercMatchingSpikes(IndexesNonNan)))
%Model = fitlm(log10(Table.nSpikes(IndexesNonNan)),log10(Table.PercMatchingSpikes(IndexesNonNan)))

%figure
%plot(Model)
%%

SpeciesFilter = {'Mus musculus_Shuff' 'Mus musculus'};
%SpeciesFilter = {'Mus musculus_Shuff' 'Mus musculus' 'Homo Sapiens_Shuff' 'Homo Sapiens'};
%SpeciesFilter = {'Homo Sapiens_Shuff' 'Homo Sapiens'};

TablesPerArea = {};

%SelectedVariable = 'structure__acronym';
%SelectedVariable = 'structure__layer';
%SelectedVariable = 'tag__dendrite_type';
%SelectedVariable = 'donor__sex';
%SelectedVariable = 'donor__disease_state';
%SelectedVariable = 'donor__species';
SelectedVariable = 'line_name';
%SelectedVariable = 'specimen__hemisphere';
%SelectedVariable = 'structure__name';
%SelectedVariable = 'cell_reporter_status';
%SelectedVariable = 'tag__apical';
%SelectedVariable = 'donor__age';



UniqueStructAcronyms = unique(Table.(SelectedVariable));
UniqueStructAcronyms(find(strcmp(UniqueStructAcronyms,'_Shuff'))) = [];
if ~iscell(UniqueStructAcronyms)
    UniqueStructAcronyms(isnan(UniqueStructAcronyms)) = [];
end
for counterStructAcronym = 1:size(UniqueStructAcronyms,1)
    if iscell(UniqueStructAcronyms)
        currStructAcronym = UniqueStructAcronyms{counterStructAcronym};
        currStructAcronymCorrected = regexprep(currStructAcronym,'-','_');
        currStructAcronymCorrected = regexprep(currStructAcronymCorrected,'/','');
        currStructAcronymCorrected = regexprep(currStructAcronymCorrected,',','');
        currStructAcronymCorrected = regexprep(currStructAcronymCorrected,'"','');
        currStructAcronymCorrected = regexprep(currStructAcronymCorrected,' ','_');
        currStructAcronymCorrected = erase(currStructAcronymCorrected,'|');
        if ~isempty(currStructAcronymCorrected)
            
            if ~isnan(str2double(currStructAcronymCorrected(1)))
                currStructAcronymCorrected = ['a' currStructAcronymCorrected];
            end
            
            Filters = ones(size(Table.(SelectedVariable)));
            
            %Filters = Filters .* strcmp(Table.('cell_reporter_status'), 'positive');
            if strcmp(SelectedVariable,'line_name') %%% CHECK IF THIS MAKES SENSE -  What is cell_reporter_status ????
                Filters = Filters .* strcmp(Table.('cell_reporter_status'), 'positive');
            end
            if ~isempty(SpeciesFilter)
                Filters = Filters .* ismember(Table.('donor__species'), SpeciesFilter);
            end
            
            
            TablesPerArea.(currStructAcronymCorrected) = Table(find(strcmp(Table.(SelectedVariable), currStructAcronym).*Filters),:);
            
        end
    else
        currStructAcronym = UniqueStructAcronyms(counterStructAcronym);
        currStructAcronymCorrected = ['L' num2str(currStructAcronym)];
        TablesPerArea.(regexprep(currStructAcronymCorrected,'/','')) = Table(find((Table.(SelectedVariable)==currStructAcronym)),:);

    end
    
    
end

%


SelectedAreasToPlot = fields(TablesPerArea)';%{'VISp6a','VISp5','VISp4','VISp23','VISp1'};
ToBoxplot = {};
ToBoxplotAreas = {};
ToKruskal = [];
ToKruskalAreas = [];

%figure, hold on

for counterAreas = 1:size(SelectedAreasToPlot,2)
    currAreaToPlot = SelectedAreasToPlot{counterAreas};
    tmpValue = TablesPerArea.(currAreaToPlot).PercMatchingSpikes;
    %tmpValue = TablesPerArea.(currAreaToPlot).DistanceBetweenSpikes_Matched_MeanAbs;
    %tmpValue = TablesPerArea.(currAreaToPlot).DistanceBetweenSpikes_Matched_Std;
    %tmpValue = TablesPerArea.(currAreaToPlot).DistanceBetweenSpikes_All_Std;
    %tmpValue = TablesPerArea.(currAreaToPlot).DistanceBetweenSpikes_All_MeanAbs;
    %tmpValue = TablesPerArea.(currAreaToPlot).CorrFull_Spearman_R;
    %tmpValue = TablesPerArea.(currAreaToPlot).CorrFull_Pearson_R;
    %tmpValue = TablesPerArea.(currAreaToPlot).CorrTrim_Spearman_R;
    %tmpValue = TablesPerArea.(currAreaToPlot).CorrTrim_Pearson_R;
    %tmpValue = TablesPerArea.(currAreaToPlot).ISI_Std;
    %tmpValue = TablesPerArea.(currAreaToPlot).ISI_Mean;
    %tmpValue = TablesPerArea.(currAreaToPlot).nSpikes;
    
    
    tmpValueNSpikes = TablesPerArea.(currAreaToPlot).nSpikes;
    
    
    % Filter max n spikes
    tmpValue = tmpValue(find((TablesPerArea.(currAreaToPlot).nSpikes>2).*(TablesPerArea.(currAreaToPlot).nSpikes<500)));
    %tmpValue = tmpValue(find((TablesPerArea.(currAreaToPlot).nSpikes>=00).*(TablesPerArea.(currAreaToPlot).nSpikes<50)));
    
    tmpValueNSpikes = tmpValueNSpikes(find((TablesPerArea.(currAreaToPlot).nSpikes>2).*(TablesPerArea.(currAreaToPlot).nSpikes<500)));
    
    tmpValueNSpikes(isnan(tmpValue))=[];
    tmpValue(isnan(tmpValue))=[];
    if size(tmpValue,1)>5%1
        ToBoxplot = [ToBoxplot {tmpValue}];
        ToBoxplotAreas = [ToBoxplotAreas {currAreaToPlot}];
        
        ToKruskal = [ToKruskal;tmpValue];
%        ToKruskalAreas = [ToKruskalAreas,repmat(currAreaToPlot,size(tmpValue))];
        ToKruskalAreas = [ToKruskalAreas;repmat({currAreaToPlot},size(tmpValue))];
        
        %scatter(log(tmpValueNSpikes),tmpValue)
        %pause()
    end
end

[tmp,SortingIndexes] = sort(cellfun(@median, ToBoxplot));
figure
boxplotGroup(ToBoxplot(SortingIndexes),'primarylabels',regexprep(ToBoxplotAreas(SortingIndexes),'_','-'));
xtickangle(45)

%boxplotGroup([{TablesPerArea.VISp5.PercMatchingSpikes},{TablesPerArea.VISp23.PercMatchingSpikes}])

[p,tbl,stats] = kruskalwallis(ToKruskal,ToKruskalAreas);
c = multcompare(stats);

%% Create table concatenated
MainFolderTables = 'C:\Users\simon\OneDrive\Documenti\Lab\PatchSeqAllen\New_Version_1ms';

BigTable = [];
%DirMainFolderTables = dir(fullfile(MainFolderTables,'Table_Noise*.mat'));
DirMainFolderTables = dir(fullfile(MainFolderTables,'Table_Noise 1_2*.mat'));
for counterTable = 1:size(DirMainFolderTables,1)
    load(fullfile(DirMainFolderTables(counterTable).folder,DirMainFolderTables(counterTable).name))
    
    Table.Group = repmat(counterTable,[size(Table,1) 1]);
    
    if isempty(BigTable)
        BigTable = Table;
        %tmp1 = Table.DistanceBetweenSpikes_Matched_MeanAbs;
    else
        BigTable = [BigTable; Table];
        %tmp1 = [tmp1 Table.DistanceBetweenSpikes_Matched_MeanAbs];
        %tmp1 = [tmp1 Table.nSpikes];
        
    end
end

Table = BigTable;

%% correct all factors with shuffle correction
load('C:\Users\simon\OneDrive\Documenti\Lab\PatchSeqAllen\New_Version_1ms\Permutation.mat')


%SelectedFactorName = {'PercMatchingSpikes','DistanceBetweenSpikes_Matched_MeanAbs','DistanceBetweenSpikes_Matched_Std','DistanceBetweenSpikes_All_Std','DistanceBetweenSpikes_All_MeanAbs'};
SelectedFactorName = {'PercMatchingSpikes'};
%SelectedFactorName = 'DistanceBetweenSpikes_Matched_MeanAbs';
%SelectedFactorName = 'DistanceBetweenSpikes_Matched_Std';
%SelectedFactorName = 'DistanceBetweenSpikes_All_Std';
%SelectedFactorName = 'DistanceBetweenSpikes_All_MeanAbs';
%SelectedFactorName = 'nSpikes';


ToCorrectNames = {'line_name','specimen__name','specimen__hemisphere','structure__name','structure__acronym','tag__apical','tag__dendrite_type','donor__species'}; %structure__layer

Table = Table(find(Table.nSpikes<=500),:);
%Table = Table(find(Table.nSpikes>0),:);
Table = Table(find(round(Table.nSpikes)>0),:);

ShuffleCorrectedTable = Table;

for currSelectedFactorIndex = 1:size(SelectedFactorName,2)
    currSelectedFactorName = SelectedFactorName{currSelectedFactorIndex};
    SelectedFactorValues = nan(size(OutputPermutations,1),size(OutputPermutations,2));
    for counterNSpikes = 1:size(OutputPermutations,1)
        for counterNIterations = 1:size(OutputPermutations,2)
            
            SelectedFactorValues(counterNSpikes,counterNIterations) = nanmean(nanmean(OutputPermutations{counterNSpikes,counterNIterations}.(currSelectedFactorName)));
            
        end
    end
    
    %plot(nanmean(SelectedFactorValues,2))
    %pause()
    
    MeanSelectedFactorValues = nanmean(SelectedFactorValues,2);
    %
    
    
    
    
    ShuffleCorrectedTable.(currSelectedFactorName) = ShuffleCorrectedTable.(currSelectedFactorName) - MeanSelectedFactorValues(round(Table.nSpikes));
    %ShuffleCorrectedTable.(currSelectedFactorName) = MeanSelectedFactorValues(round(Table.nSpikes));
    
end



for currSelectedFactorIndex = 1:size(ToCorrectNames,2)
    currSelectedFactorName = ToCorrectNames{currSelectedFactorIndex};
    currValues = ShuffleCorrectedTable.(currSelectedFactorName);
    if isnumeric(currValues)
        currValues = num2str(currValues);
    end
    
    currValues(find(strcmp(currValues,'_Shuff'))) = {'NO'};
    ShuffleCorrectedTable.(currSelectedFactorName) = strcat(currValues,'_Shuff');
end

%Table = [Table;ShuffleCorrectedTable];
Table = ShuffleCorrectedTable;
%%

%tmp1 = Table(Table.Group==3,:).PercMatchingSpikes;
%tmp2 = Table(Table.Group==6,:).PercMatchingSpikes;
tmp1 = Table(Table.Group==3,:).PercMatchingSpikes;
tmp2 = Table(Table.Group==6,:).PercMatchingSpikes;

boxplot([tmp1;tmp2],[ones(size(tmp1));ones(size(tmp2))*2])
ranksum(tmp1,tmp2)

%%
%tmp1 = Table(Table.Group==1,:).PercMatchingSpikes;
%tmp2 = Table(Table.Group==3,:).PercMatchingSpikes;
tmp1 = Table(Table.Group==1,:).nSpikes-Table(Table.Group==3,:).nSpikes;
tmp2 = Table(Table.Group==1,:).PercMatchingSpikes-Table(Table.Group==3,:).PercMatchingSpikes;


figure, hold on
scatter(tmp1,tmp2)
%hL=plot([0 1],[0 1]);   % 45 line visible at 45



%%

























%%


%%

%tmp1 = Table.DistanceBetweenSpikes_Matched_MeanAbs;
%tmp1 = Table.PercMatchingSpikes;
%tmp2 = Table.nr__max_euclidean_distance;
%tmp2 = Table.nr__number_stems;
%tmp2 = Table.nr__number_bifurcations;
%tmp2 = Table.nr__average_contraction;
%tmp2 = Table.nr__average_parent_daughter_ratio;
%tmp2 = Table.ef__fast_trough_v_long_square;
%tmp2 = Table.ef__upstroke_downstroke_ratio_long_square;
%tmp2 = Table.ef__adaptation;
%tmp2 = Table.ef__f_i_curve_slope;
%tmp2 = Table.ef__threshold_i_long_square;
%tmp2 = Table.ef__tau;
%tmp2 = Table.ef__ri;
%tmp2 = Table.ef__avg_isi;
%tmp2 = Table.ef__avg_firing_rate;
%tmp2 = Table.ef__peak_t_ramp;
%tmp2 = Table.ef__vrest;
%tmp2 = Table.si__height;
%tmp2 = Table.si__width;
%tmp1 = Table.csl__normalized_depth;
tmp1 = Table.nSpikes;
tmp2 = Table.ISI_Mean;
%tmp1 = Table.ISI_Std;




%tmpFilter = find(Table.structure__layer==4);
%tmpFilter = find(strcmp(Table.tag__dendrite_type,'spiny'));
%tmpFilter = intersect(find(strcmp(Table.line_name,'Chrna2-Cre_OE25')),find(strcmp(Table.cell_reporter_status,'positive')));
%tmp1 = Table.nSpikes;
%tmp2 = Table.ISI_Std;

%tmp1 = tmp1(tmpFilter);
%tmp2 = tmp2(tmpFilter);


scatter(tmp1,tmp2)
fitlm(tmp1,tmp2)

%%
figure, hold on
UniqueLineNames = unique(Table.line_name);
for counterLines = 1:size(UniqueLineNames,1)
    tmpFilter = intersect(find(strcmp(Table.line_name,UniqueLineNames{counterLines})),find(strcmp(Table.cell_reporter_status,'positive')));
    tmp1 = Table.nSpikes;
    tmp2 = Table.ISI_Std;
    
    tmp1 = tmp1(tmpFilter);
    tmp2 = tmp2(tmpFilter);
    
    
    scatter(tmp1,tmp2)
    title(UniqueLineNames{counterLines})
    pause()
end
%fitlm(tmp1,tmp2)

%%

% 40-42-44-46
% 41-43-?

figure, hold on
%tmpData = h5read('C:\Users\simon\Downloads\474626524_ephys.nwb','/stimulus/presentation/Sweep_40/data');
tmpData = h5read('C:\Users\simon\Downloads\474626524_ephys.nwb','/acquisition/timeseries/Sweep_40/data');
%tmpName = h5read('C:\Users\simon\Downloads\474626524_ephys.nwb','/stimulus/presentation/Sweep_40/aibs_stimulus_description');
tmpName = h5read('C:\Users\simon\Downloads\474626524_ephys.nwb','/stimulus/presentation/Sweep_40/aibs_stimulus_name');
plot(tmpData)


%tmpData = h5read('C:\Users\simon\Downloads\474626524_ephys.nwb','/stimulus/presentation/Sweep_44/data');
tmpData = h5read('C:\Users\simon\Downloads\474626524_ephys.nwb','/acquisition/timeseries/Sweep_44/data');
%tmpName = h5read('C:\Users\simon\Downloads\474626524_ephys.nwb','/stimulus/presentation/Sweep_44/aibs_stimulus_description');
tmpName = h5read('C:\Users\simon\Downloads\474626524_ephys.nwb','/stimulus/presentation/Sweep_44/aibs_stimulus_name');
plot(tmpData)


%tmpData = h5read('C:\Users\simon\Downloads\474626524_ephys.nwb','/stimulus/presentation/Sweep_46/data');
tmpData = h5read('C:\Users\simon\Downloads\474626524_ephys.nwb','/acquisition/timeseries/Sweep_46/data');
%tmpName = h5read('C:\Users\simon\Downloads\474626524_ephys.nwb','/stimulus/presentation/Sweep_46/aibs_stimulus_description');
tmpName = h5read('C:\Users\simon\Downloads\474626524_ephys.nwb','/stimulus/presentation/Sweep_46/aibs_stimulus_name');
plot(tmpData)





%%

histogram(Table.ef__ri)