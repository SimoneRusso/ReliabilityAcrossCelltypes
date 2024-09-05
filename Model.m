



%% run model - single session

clear

%InputType = 'Excitatory';
InputType = 'Inhibitory';

Jitter = 20;%20;
%Jitter = 0;

Ts_N = 500;
nGaussians = 20;



LAll_N = 1;

LAll_param1 = 0.02;
LAll_param2 = 0.08;
LAll_param3 = 0.2;
LAll_param4 = 0.25;
LAll_param5 = 0.05;




re=ones(LAll_N,1); 
a=[LAll_param1*ones(LAll_N,1)];
b=[LAll_param3*ones(LAll_N,1)];

c=[-65*re.^2; -65*ones(LAll_N,1)];
d=[8-6*re.^2; 2*ones(LAll_N,1)];
%S=[L1_Connectivity.*0.7;L2_Connectivity.*0.7]';

v=-65*ones(LAll_N,1); % Initial values of v
u=b.*v; % Initial values of u
firings=[]; % spike timings

%currS = zeros(LAll_N,LAll_N,SynapseDuration);

Ts_V = nan(LAll_N,100);


x = 0.1:0.1:10;
y = [zeros(1,49) 5 zeros(1,50)];  %y = gaussmf(x,[0.03 5])*5;
if strcmp(InputType,'Inhibitory')
    y = -y;
end


TsI = zeros(1,Ts_N);


RandomJitter = rand(1,nGaussians)*Jitter;
OriginalTimeWindowGaussian = 201:300;

TsPeakGaussians = [];
for counterGaussians = 1:nGaussians
    currTimeWindowGaussian = OriginalTimeWindowGaussian + round(RandomJitter(counterGaussians));
    TsPeakGaussians(counterGaussians) = median(currTimeWindowGaussian)
    TsI(currTimeWindowGaussian) = TsI(currTimeWindowGaussian) + y;
end

for t=1:Ts_N % simulation of 1000 ms??
    
    
    I = TsI(t);
    
    
    %if t== 800
    %    v(1) = 30;
    %   
    %end
    
    
    fired=find(v>=30); % indices of spikes
    v(fired)=30;
    
    Ts_V(:,t) = v;
    
    firings=[firings; t+0*fired,fired];
    v(fired)=c(fired);
    u(fired)=u(fired)+d(fired);
    
    %tmpS = zeros(size(S));
    %tmpS(:,fired) = S(:,fired);
    
    %if t== 800
    %    tmpS = tmpS+0.5;
    %   
    %end
    
    
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
    u=u+a.*(b.*v-u); % stability
    
    
end

fig=figure;
hold on

YLim = [-105 50];
plot(Ts_V)
scatter(TsPeakGaussians,[30+[1:nGaussians]*((max(YLim)-1-30)/nGaussians)],10,'k','filled')
ylim(YLim)
xlim([150 400])

set(gcf,'renderer','painters')

%saveas(fig,fullfile('C:\Users\simon\OneDrive\Documenti\Lab\PatchSeqAllen\Figures\Model\', ['InputType-' InputType '_Jitter-' num2str(Jitter) '_NGaussians-' num2str(nGaussians) '.svg']))



%% run model - exploration evoked frequency as a function of number of inputs 

clear



ParamExplor_NGaussians = [1:50];
NIterations = 500;
NSpikes = nan(size(ParamExplor_NGaussians,2),NIterations);
TimeWindow = nan(size(ParamExplor_NGaussians,2),NIterations);
for counterParamExplor = 1:size(ParamExplor_NGaussians,2)
    for counterIterations = 1:NIterations;
        
        InputType = 'Excitatory';
        %InputType = 'Inhibitory';
        
        %Jitter = 10;%20;
        %Jitter = 0;
        Jitter = 50;%20;
        %Jitter = 200;%20;
        
        Ts_N = 500;
        %Ts_N = 1000;
        nGaussians = ParamExplor_NGaussians(counterParamExplor);
        
        
        
        LAll_N = 1;
        
        LAll_param1 = 0.02;
        LAll_param2 = 0.08;
        LAll_param3 = 0.2;
        LAll_param4 = 0.25;
        LAll_param5 = 0.05;
        
        
        
        
        re=ones(LAll_N,1);
        a=[LAll_param1*ones(LAll_N,1)];
        b=[LAll_param3*ones(LAll_N,1)];
        
        c=[-65*re.^2; -65*ones(LAll_N,1)];
        d=[8-6*re.^2; 2*ones(LAll_N,1)];
        %S=[L1_Connectivity.*0.7;L2_Connectivity.*0.7]';
        
        v=-65*ones(LAll_N,1); % Initial values of v
        u=b.*v; % Initial values of u
        firings=[]; % spike timings
        
        %currS = zeros(LAll_N,LAll_N,SynapseDuration);
        
        Ts_V = nan(LAll_N,100);
        
        
        y = [zeros(1,49) 5 zeros(1,50)]; 
        
        %x = 0.1:0.1:10;
        %y = gaussmf(x,[0.03 5])*5;
        if strcmp(InputType,'Inhibitory')
            y = -y;
        end
        
        
        TsI = zeros(1,Ts_N);
        
        
        RandomJitter = rand(1,nGaussians)*Jitter;
        OriginalTimeWindowGaussian = 201:300;
        
        TsPeakGaussians = [];
        for counterGaussians = 1:nGaussians
            currTimeWindowGaussian = OriginalTimeWindowGaussian + round(RandomJitter(counterGaussians));
            TsPeakGaussians(counterGaussians) = median(currTimeWindowGaussian);
            TsI(currTimeWindowGaussian) = TsI(currTimeWindowGaussian) + y;
        end
        
        for t=1:Ts_N % simulation of 1000 ms??
            
            
            I = TsI(t);
            
            
            %if t== 800
            %    v(1) = 30;
            %
            %end
            
            
            fired=find(v>=30); % indices of spikes
            v(fired)=30;
            
            Ts_V(:,t) = v;
            
            firings=[firings; t+0*fired,fired];
            v(fired)=c(fired);
            u(fired)=u(fired)+d(fired);
            
            %tmpS = zeros(size(S));
            %tmpS(:,fired) = S(:,fired);
            
            %if t== 800
            %    tmpS = tmpS+0.5;
            %
            %end
            
            
            v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
            v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
            u=u+a.*(b.*v-u); % stability
            
            
        end
        
        
        NSpikes(counterParamExplor,counterIterations) = size(firings,1);
        if size(firings,1)>0
            TimeWindow(counterParamExplor,counterIterations) = max(firings(:,1))-min(firings(:,1));
        end
    end
    counterParamExplor
end


TimeWindow(TimeWindow==0) = 25;
TimeWindow = TimeWindow/1000;
%imagesc(NSpikes)
fig = figure
plot(nanmean(NSpikes,2))
%plot(nanmean(NSpikes./TimeWindow,2))
%plot(nanmean(NSpikes,2)./0.025)

%fig=figure;
%hold on
%
%YLim = [-105 50];
%plot(Ts_V)
%scatter(TsPeakGaussians,[30+[1:nGaussians]*((max(YLim)-1-30)/nGaussians)],10,'k','filled')
%ylim(YLim)
%xlim([150 400])



set(gcf,'renderer','painters')

%saveas(fig,fullfile('C:\Users\simon\OneDrive\Documenti\Lab\PatchSeqAllen\Figures\Model\', ['ParamExploration_NGaussians.svg']))
%saveas(fig,fullfile('C:\Users\simon\OneDrive\Documenti\Lab\PatchSeqAllen\Figures\Model\', ['ParamExploration_NGaussians_New.svg']))


MI_x = nanmean(NSpikes,2)./0.025;
MI_y = ParamExplor_NGaussians;
%mi_cont_cont(MI_x,MI_y,3)
mutualinfo(MI_x,MI_y)

%% run model - exploration evoked frequency/ n spikes as a function of jitter 

clear



ParamExplor_Jitter = [0:50];
NIterations = 500;%500;
NSpikes = nan(size(ParamExplor_Jitter,2),NIterations);
TimeWindow = nan(size(ParamExplor_Jitter,2),NIterations);
Integral = nan(size(ParamExplor_Jitter,2),NIterations);

for counterParamExplor = 1:size(ParamExplor_Jitter,2)
    for counterIterations = 1:NIterations;
        
        InputType = 'Excitatory';
        %InputType = 'Inhibitory';
        
        Jitter = ParamExplor_Jitter(counterParamExplor);%20;
        %Jitter = 0;
        
        Ts_N = 500;
        nGaussians = 20; %ParamExplor_NGaussians(counterParamExplor);
        
        
        
        LAll_N = 1;
        
        LAll_param1 = 0.02;
        LAll_param2 = 0.08;
        LAll_param3 = 0.2;
        LAll_param4 = 0.25;
        LAll_param5 = 0.05;
        
        
        
        
        re=ones(LAll_N,1);
        a=[LAll_param1*ones(LAll_N,1)];
        b=[LAll_param3*ones(LAll_N,1)];
        
        c=[-65*re.^2; -65*ones(LAll_N,1)];
        d=[8-6*re.^2; 2*ones(LAll_N,1)];
        %S=[L1_Connectivity.*0.7;L2_Connectivity.*0.7]';
        
        v=-65*ones(LAll_N,1); % Initial values of v
        u=b.*v; % Initial values of u
        firings=[]; % spike timings
        
        %currS = zeros(LAll_N,LAll_N,SynapseDuration);
        
        Ts_V = nan(LAll_N,100);
        
        
        y = [zeros(1,49) 5 zeros(1,50)]; 
        
        %x = 0.1:0.1:10;
        %y = gaussmf(x,[0.03 5])*5;
        if strcmp(InputType,'Inhibitory')
            y = -y;
        end
        
        
        TsI = zeros(1,Ts_N);
        
        
        RandomJitter = rand(1,nGaussians)*Jitter;
        OriginalTimeWindowGaussian = 201:300;
        
        TsPeakGaussians = [];
        for counterGaussians = 1:nGaussians
            currTimeWindowGaussian = OriginalTimeWindowGaussian + round(RandomJitter(counterGaussians));
            TsPeakGaussians(counterGaussians) = median(currTimeWindowGaussian);
            TsI(currTimeWindowGaussian) = TsI(currTimeWindowGaussian) + y;
        end
        
        for t=1:Ts_N % simulation of 1000 ms??
            
            
            I = TsI(t);
            
            
            %if t== 800
            %    v(1) = 30;
            %
            %end
            
            
            fired=find(v>=30); % indices of spikes
            v(fired)=30;
            
            Ts_V(:,t) = v;
            
            firings=[firings; t+0*fired,fired];
            v(fired)=c(fired);
            u(fired)=u(fired)+d(fired);
            
            %tmpS = zeros(size(S));
            %tmpS(:,fired) = S(:,fired);
            
            %if t== 800
            %    tmpS = tmpS+0.5;
            %
            %end
            
            
            v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
            v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
            u=u+a.*(b.*v-u); % stability
            
            
        end
        
        tmp_Ts = Ts_V(200:400);
        if strcmp(InputType, 'Inhibitory')
            %Integral(counterParamExplor,counterIterations) = nansum(tmp_Ts);
            %Integral(counterParamExplor,counterIterations) = nansum(tmp_Ts(tmp_Ts<0));
            %Integral(counterParamExplor,counterIterations) = nanmedian(tmp_Ts(tmp_Ts<0));
            %Integral(counterParamExplor,counterIterations) = sqrt(mean((tmp_Ts.^2)));
            Integral(counterParamExplor,counterIterations) = nanmin(tmp_Ts);
        else
            %Integral(counterParamExplor,counterIterations) = nansum(tmp_Ts);
            %Integral(counterParamExplor,counterIterations) = nansum(tmp_Ts(tmp_Ts>0));
            %Integral(counterParamExplor,counterIterations) = nanmedian(tmp_Ts(tmp_Ts>0));
            %Integral(counterParamExplor,counterIterations) = sqrt(mean((tmp_Ts.^2)));
            Integral(counterParamExplor,counterIterations) = nanmax(tmp_Ts);
        end
        
        NSpikes(counterParamExplor,counterIterations) = size(firings,1);
        if size(firings,1)>0
            TimeWindow(counterParamExplor,counterIterations) = max(firings(:,1))-min(firings(:,1));
        end
    end
    counterParamExplor
end


TimeWindow(TimeWindow==0) = 25;
TimeWindow = TimeWindow/1000;
%imagesc(NSpikes)
fig = figure
%plot(nanmean(NSpikes,2)./0.025)
plot(nanmean(NSpikes,2))

fig2 = figure
plot(ParamExplor_Jitter,nanmean(Integral,2))

%fig=figure;
%hold on
%
%YLim = [-105 50];
%plot(Ts_V)
%scatter(TsPeakGaussians,[30+[1:nGaussians]*((max(YLim)-1-30)/nGaussians)],10,'k','filled')
%ylim(YLim)
%xlim([150 400])



set(gcf,'renderer','painters')

%saveas(fig,fullfile('C:\Users\simon\OneDrive\Documenti\Lab\PatchSeqAllen\Figures\Model\', ['ParamExploration_Variability.svg']))


%% run model - Number of spikes as a function of number of inputs and jittering

clear



ParamExplor_Jitter = [0:50];
ParamExplor_NGaussians = [1:50];
NIterations = 100;
NSpikes = nan(size(ParamExplor_NGaussians,2),NIterations,size(ParamExplor_Jitter,2));

for counterParamExplorJitter = 1:size(ParamExplor_Jitter,2)
    for counterParamExplor = 1:size(ParamExplor_NGaussians,2)
        for counterIterations = 1:NIterations;
            
            InputType = 'Excitatory';
            %InputType = 'Inhibitory';
            
            %Jitter = 20;%20;
            %Jitter = 0;
            Jitter = ParamExplor_Jitter(counterParamExplorJitter);%20;
            
            Ts_N = 500;
            nGaussians = ParamExplor_NGaussians(counterParamExplor);
            
            
            
            LAll_N = 1;
            
            LAll_param1 = 0.02;
            LAll_param2 = 0.08;
            LAll_param3 = 0.2;
            LAll_param4 = 0.25;
            LAll_param5 = 0.05;
            
            
            
            
            re=ones(LAll_N,1);
            a=[LAll_param1*ones(LAll_N,1)];
            b=[LAll_param3*ones(LAll_N,1)];
            
            c=[-65*re.^2; -65*ones(LAll_N,1)];
            d=[8-6*re.^2; 2*ones(LAll_N,1)];
            %S=[L1_Connectivity.*0.7;L2_Connectivity.*0.7]';
            
            v=-65*ones(LAll_N,1); % Initial values of v
            u=b.*v; % Initial values of u
            firings=[]; % spike timings
            
            %currS = zeros(LAll_N,LAll_N,SynapseDuration);
            
            Ts_V = nan(LAll_N,100);
            
            
            y = [zeros(1,49) 5 zeros(1,50)]; 
            
            %x = 0.1:0.1:10;
            %y = gaussmf(x,[0.03 5])*5;
            if strcmp(InputType,'Inhibitory')
                y = -y;
            end
            
            
            TsI = zeros(1,Ts_N);
            
            
            RandomJitter = rand(1,nGaussians)*Jitter;
            OriginalTimeWindowGaussian = 201:300;
            
            TsPeakGaussians = [];
            for counterGaussians = 1:nGaussians
                currTimeWindowGaussian = OriginalTimeWindowGaussian + round(RandomJitter(counterGaussians));
                TsPeakGaussians(counterGaussians) = median(currTimeWindowGaussian);
                TsI(currTimeWindowGaussian) = TsI(currTimeWindowGaussian) + y;
            end
            
            for t=1:Ts_N % simulation of 1000 ms??
                
                
                I = TsI(t);
                
                
                %if t== 800
                %    v(1) = 30;
                %
                %end
                
                
                fired=find(v>=30); % indices of spikes
                v(fired)=30;
                
                Ts_V(:,t) = v;
                
                firings=[firings; t+0*fired,fired];
                v(fired)=c(fired);
                u(fired)=u(fired)+d(fired);
                
                %tmpS = zeros(size(S));
                %tmpS(:,fired) = S(:,fired);
                
                %if t== 800
                %    tmpS = tmpS+0.5;
                %
                %end
                
                
                v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
                v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
                u=u+a.*(b.*v-u); % stability
                
                
            end
            
            
            NSpikes(counterParamExplor,counterIterations,counterParamExplorJitter) = size(firings,1);
        end
    end
    counterParamExplorJitter
end

fig = figure
imagesc(squeeze(nanmean(NSpikes,2)./0.025))

%set(gcf,'renderer','painters')

%saveas(fig,fullfile('C:\Users\simon\OneDrive\Documenti\Lab\PatchSeqAllen\Figures\Model\', ['ParamExploration_NGaussians.svg']))
%saveas(fig,fullfile('C:\Users\simon\OneDrive\Documenti\Lab\PatchSeqAllen\Figures\Model\', ['ParamExploration_NGaussians_New.svg']))

