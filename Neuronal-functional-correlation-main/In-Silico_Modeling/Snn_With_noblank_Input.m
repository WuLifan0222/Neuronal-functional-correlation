clear;
clc;
addpath('.\util\');

% set random seed
seed = 1234; 
rng(seed);    
%% Parameters
opt.seed = seed;
opt.nNeuron = 1e4; % number of neuron
opt.g = 1.2; % The connection strength
opt.dT = 0.25; % The simulation timestep, in sec 
opt.TAll = 1330; % Total simulation time, in sec 

opt.tau = 1; % Decay costant of RNN units
opt.tauWN = 2; % Decay costant of noise

opt.StiDur = 3; % Stimulus duration, in sec 
opt.StiGap = 10; % Stimulation gap, in sec
opt.JitterSig = 0; % The gap jitter range
opt.StiAmp = 3; % StiAmplitude, if a list, the amplitude will be organized in a pseudo-random manner
opt.StiTar= datasample(1:2000, 250, 'Replace', false); % Stimulus target neuron
opt.ampInWN = 0.3;
opt.SimTimes = 1; % The simulation repeat times
opt.Weightdir=[]; % The dir of pre-defined weight matrix, (if not empty, otherwise will be randomly generated)

nN=opt.nNeuron;
g=opt.g;
dT=opt.dT;
outname='visual_stim';
outpath=['output/',outname,'/'];
mkdir(outpath);

%% Distinguish neurons corresponding to different stimuli 
opt.StiNum = 5;
    
num_neurons_per_stim = 200;
sti_index_sequences = cell(5, 1);
all_neurons = opt.StiTar;  

for i = 1:opt.StiNum
    sti_index_sequences{i} = datasample(all_neurons, num_neurons_per_stim, 'Replace', false);
end


%% Simulation parameters
TAll=opt.TAll;
nT=TAll/dT; % The simulation total steps
tau=opt.tau;
StiAmp=opt.StiAmp;
StiTar=opt.StiTar;
nSTy=length(StiAmp); % The number of types of stimulation 

%% External sti
IpS=abs(normrnd(0,0.2,[length(opt.StiTar),1]));
[IpS,~]=sort(IpS,'descend'); % The projection weight
stiSqc=zeros(1,nT);
stiSqc_with_label=zeros(1,nT);
nS=floor(TAll/(opt.StiDur+opt.StiGap))-1;
nSE=floor(nS/nSTy);
StiArray=repmat(StiAmp,nSE);
StiArray=StiArray(randperm(length(StiArray)));
nS=nSE*nSTy;
sedge=zeros(1,nS);
eedge=zeros(1,nS);
E=zeros(nN,nT);

index_loop = 1; 
for iss=1:nS
    if iss==1
        estt=opt.StiGap/dT;
    else
        % estt=round(eend+(opt.StiDur+rand()*opt.JitterSig-opt.JitterSig/2)/dT+opt.StiGap/dT);
        estt=round(eend+opt.StiGap/dT);
    end
    eend=estt+(opt.StiDur)/dT;
    if eend>nT-100 % If too many, stop
        nS = iss-1;
        StiArray=StiArray(1:nS);
        sedge=sedge(1:nS);
        eedge=eedge(1:nS);
        break;        
    end
    sedge(iss)=estt+1;
    eedge(iss)=eend;
    stiSqc(estt+1:eend)=StiArray(iss);
    stiSqc_with_label(estt+1:eend) = index_loop;
    E(sti_index_sequences{index_loop}, estt+1:eend)=IpS(1:num_neurons_per_stim)*stiSqc(estt+1:eend);
    index_loop = index_loop + 1;
    if index_loop > opt.StiNum
        index_loop = 1;
    end
end

T=tabulate(StiArray);
nSE=min(T(:,2));

% E=zeros(nN,nT); 
% E(opt.StiTar,:)=IpS*stiSqc;

figure();
plot(stiSqc); title('Stimulation sequence');


%% Run for multiple times
SimTimes=opt.SimTimes;
FrAll=zeros(SimTimes,nN,nT);
SAll=zeros(SimTimes,nN,nT);
SRAll=zeros(SimTimes,nN,nT);
for sstt=1:SimTimes
    % Randomize a connectivity matrix
    M = zeros(nN, nN);
    local_range = 300; % Local Connection Range
    for i = 1:nN
        neighbors = max(1, i - local_range):min(nN, i + local_range);
        M(i, neighbors) = g * randn(1, length(neighbors)) / sqrt(local_range * 2);
    end
    % save([outpath,'J',num2str(sstt),'.mat'],'J');

    % set up white noise inputs
    tauWN=opt.tauWN;
    ampInWN=opt.ampInWN;
    ampWN = sqrt((tauWN/dT));
    iWN = ampWN*randn(nN, nT);
    inputWN = ones(nN,nT);
    for tt = 2: nT
        inputWN(:, tt) = iWN(:, tt) + (inputWN(:, tt - 1) - iWN(:, tt))*exp(-(dT/tauWN));
    end
    inputWN = ampInWN*inputWN; % input noise

    ini=rand(nN,1)*2-1; % The initial membrane potential
    
    % Run the network
%   [Fr,S]=runSNN_v1(ini,J,nT,tau,E,inputWN,dT);
    [Fr,S,frac]=runSNN_v2_showfrac(ini,M,nT,tau,E,inputWN,dT);

    % Collect result
    FrAll(sstt,:,:)=Fr;
    SAll(sstt,:,:)=S;
    SR=getSpkRateHG(S,20);
    SRAll(sstt,:,:)=SR;
    
    close all;
    fprintf('%d / %d \n',sstt,SimTimes);
end

%% adjust real data, saving result
start_edge = [];
end_edge = [];
n = length(stiSqc);

if n == 0
    return;
end

current_state = stiSqc(1);
if current_state ~= 0
    start = 1;
    start_edge = [start_edge, start];
end

for i = 2:n
    if stiSqc(i) ~= current_state
        if current_state == 0 && stiSqc(i) ~= 0
            start = i;
            start_edge = [start_edge, start];
            current_state = stiSqc(i);
        elseif current_state ~= 0 && stiSqc(i) == 0
            end_edge = [end_edge, i-1];
            current_state = 0;
        end
    end
end

if current_state ~= 0
    end_edge = [end_edge, n];
end
for i = 1:length(start_edge)
    stimuili_label_ind(i) = mod(i-1, 5)+1;
end
stimuli_array = stiSqc;
stimuli_array_with_label = stiSqc_with_label;

save_file_name = fullfile(outpath, sprintf('visual_stimuli_with_label_%d.mat', seed));
save(save_file_name, 'end_edge', 'start_edge','stimuili_label_ind', 'stimuli_array', 'stimuli_array_with_label','Fr','opt');

% figure();
% imagesc(Fr);