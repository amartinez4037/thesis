%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Thesis Work
%    Run after starting bcilab and eeglab
%       Can also run the Thesis_Start script
%    Test to see how basic filtering looks
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set variables
%clc
% Channels used for ICA
%chanid = [2 4 6 9 10 11 12 13];
chanid - [9 11 13];
numchan = length(chanid);

% Trials to be used for analysis
%trial = {'R03', 'R04', 'R07', 'R08', 'R11', 'R12'};
trial = {'R03', 'R07', 'R11'}; % For testing
numtrials = length(trial);  

% Subjects
%subject = {'S001','S002','S003','S004','S005','S006'};
subject = {'S001','S002'};%,'S003','S004','S005','S006'}; % For testing
numsubjects = length(subject);  

% Epochs
numepochs = 7;

% Features : [avg, pwr, ene]
feat2use = [0,0,1];


% Define paths for data location and storage
% Paths for folders were data located and where to store data
homepath = '/home/amartinez/thesis/PhysionetData/EDF/'; % Location of EDF files
filepath = '/home/amartinez/thesis/userscripts2/'; % Location of storage folder

% Prealocate size for features
features = zeros(26, numepochs*numruns*numtrials*numsubjects);
size(features)

% Set each variable depending if process needs to be done
do_import = 1; % If 0 will not edit channels, filter or artifact removal
do_edit_channels = 1;
do_edit_events = 1;
do_filter = 1;
do_artifact_removal = 1;

do_epoch = 1;
do_features_wavelet = 1;
do_NN = 0;

%% Start looping through subjects and trials
for s = 1:numsubjects
    for t = 1:numtrials
        fprintf('\n*****\nProcessing Subject: %s, Trial: %s\n*****\n', subject{s}, trial{t});
        
        % Add to path the subject
        foldername = [homepath subject{s} '/'];
        % Full name of each file for subject and trial together
        savename = [subject{s} trial{t}];
        % Full name of the file to retrieve
        filename = [foldername savename '.edf'];
              
        % Check if the file exists
        if exist(filename, 'file') <=0
            fprintf('Warning file does not exist');
            
        else
            %% Import Data
            if (do_import)
                % Import using BIOSIG - for EDF files
                EEG = pop_biosig(filename);
            else
                fprintf('Skip Org Import  ');
            end
            

            %% Edit Channel Locations
            if (do_edit_channels && do_import)
                % Need to be edited because of extra periods added at the end of each
                % channel name. Every channel has been included. 
                EEG = pop_chanedit(EEG, 'lookup','standard-10-5-cap385.elp',...
                    'changefield',{1 'labels' 'Fc5'},'changefield',{2 'labels' 'Fc3'},...
                    'changefield',{3 'labels' 'Fc1'},'changefield',{4 'labels' 'Fcz'},...
                    'changefield',{5 'labels' 'Fc2'}, 'changefield',{6 'labels' 'Fc4'},...
                    'changefield',{7 'labels' 'Fc6'},'changefield',{8 'labels' 'C5'},...
                    'changefield',{9 'labels' 'C3'},'changefield',{10 'labels' 'C1'},...
                    'changefield',{11 'labels' 'Cz'},'changefield',{12 'labels' 'C2'},...
                    'changefield',{13 'labels' 'C4'},'changefield',{14 'labels' 'C6'},...
                    'changefield',{15 'labels' 'Cp5'},'changefield',{16 'labels' 'Cp3'},...
                    'changefield',{17 'labels' 'Cp1'},'changefield',{18 'labels' 'Cpz'},...
                    'changefield',{19 'labels' 'Cp2'},'changefield',{20 'labels' 'Cp4'},...
                    'changefield',{21 'labels' 'Cp6'},'changefield',{22 'labels' 'Fp1'},...
                    'changefield',{23 'labels' 'Fpz'},'changefield',{24 'labels' 'Fp2'},...
                    'changefield',{25 'labels' 'Af7'},'changefield',{26 'labels' 'Af3'},...
                    'changefield',{27 'labels' 'Afz'},'changefield',{28 'labels' 'Af4'},...
                    'changefield',{29 'labels' 'Af8'},'changefield',{30 'labels' 'F7'},...
                    'changefield',{31 'labels' 'F5'},'changefield',{32 'labels' 'F3'},...
                    'changefield',{33 'labels' 'F1'},'changefield',{34 'labels' 'Fz'},...
                    'changefield',{35 'labels' 'F2'},'changefield',{36 'labels' 'F4'},...
                    'changefield',{37 'labels' 'F6'},'changefield',{38 'labels' 'F8'},...
                    'changefield',{39 'labels' 'Ft7'},'changefield',{40 'labels' 'Ft8'},...
                    'changefield',{41 'labels' 'T7'},'changefield',{42 'labels' 'T8'},...
                    'changefield',{43 'labels' 'T9'},'changefield',{44 'labels' 'T10'},...
                    'changefield',{45 'labels' 'Tp7'},'changefield',{46 'labels' 'Tp8'},...
                    'changefield',{47 'labels' 'P7'},'changefield',{48 'labels' 'P5'},...
                    'changefield',{49 'labels' 'P3'},'changefield',{50 'labels' 'P1'},...
                    'changefield',{51 'labels' 'Pz'},'changefield',{52 'labels' 'P2'},...
                    'changefield',{53 'labels' 'P4'},'changefield',{54 'labels' 'P6'},...
                    'changefield',{55 'labels' 'P8'},'changefield',{56 'labels' 'Po7'},...
                    'changefield',{57 'labels' 'Po3'},'changefield',{58 'labels' 'Poz'},...
                    'changefield',{59 'labels' 'Po4'},'changefield',{60 'labels' 'Po8'},...
                    'changefield',{61 'labels' 'O1'},'changefield',{62 'labels' 'Oz'},...
                    'changefield',{63 'labels' 'O2'},'changefield',{64 'labels' 'Iz'},...
                    'lookup','standard-10-5-cap385.elp');
            else
                fprintf('Skip Channel Locations  ');
            end
            

            %% Edit events
            if (do_import && do_edit_events && do_edit_channels)
             
                duration1 = EEG.event(1).duration;
                duration2 = EEG.event(2).duration;
                sub1 = subject{s};
                trial1 = trial{t};

                num_annotations = length(EEG.event);
                annotations = zeros(num_annotations, 4);
                last = 0;
                
                for kk = 1: num_annotations 
                    annotations(kk,1) = kk; % urevent number
                    annotations(kk,2) = EEG.event(kk).type;
                    annotations(kk,3) = (EEG.event(kk).latency);
                    annotations(kk,4) = (EEG.event(kk).duration);
                    
                    %annotations(kk,5) = last;
                    %annotations(kk,6) = (EEG.urevent(kk).latency - 1);
                    %last = last + annotations(kk,4); % add duration of last event
                end

                % Annotate events based on function
                events_found = event_finder(duration1, duration2, sub1, trial1);
                
                % Do a check to make sure the variables known are true
                for kk = 1: num_annotations
                    if annotations(kk,2) ~= events_found(kk,2)
                        error('Types are not the same');
                    else
                        fprintf('t');
                    end
                    
                    if annotations(kk,4) ~= events_found(kk,4)
                        error('Durations are not the same');
                    else
                        fprintf('d');
                    end   
                end

                fprintf('\nChange info for EEG event and urevent\n');
                
                for kk = 1: 30
                    EEG.event(kk).urevent = events_found(kk,1);  
                    EEG.event(kk).type = events_found(kk,2);
                    EEG.event(kk).latency = events_found(kk,3);
                    EEG.event(kk).duration = events_found(kk,4);
                    EEG.urevent(kk).type = events_found(kk,2);
                    EEG.urevent(kk).latency = events_found(kk,3);
                    EEG.urevent(kk).duration = events_found(kk,4);
                end
                
            else
                fprintf('Skip Event edits  ');
            end
            


            %% Filter data
            if (do_filter && do_import && do_edit_channels)
                % To bandpass filter
                EEG = pop_eegfiltnew(EEG, 0.5, 45, 1056, 0, [], 0);
                
            else
                fprintf('Skip Freq Filtering  ');
            end
            


            %% EOG and EMG filtering and saving data
                % Values set to default
            if (do_artifact_removal && do_import && do_edit_channels && do_filter)
                % EOG removal using BSS - sobi - ...
                EEG = pop_autobsseog( EEG, [123], [123], 'sobi',...
                    {'eigratio', [1000000]}, 'eog_fd', {'range',[2 21]});

                % EMG removal using BSS
                EEG = pop_autobssemg( EEG, [81.9188], [81.9188], 'bsscca',...
                    {'eigratio', [1000000]}, 'emg_psd', {'ratio', [10], 'fs', [160], 'femg', [15],...
                    'estimator', spectrum.welch({'Hamming'}, 80),'range', [0 32]});

                % Save the dataset as it is ready to be epoched
                    % SaveName: Subject{s}Trial{t}_EqochReady.set
                %fprintf('\nSaving dataset as: %s_EpochReady.set\n', savename)
                EEG = pop_editset(EEG, 'setname', [savename '_EpochReady']);
                EEG = pop_saveset(EEG,'filename',[savename '_EpochReady.set'],'filepath',filepath);
            else
                fprintf(' Skip Artifact Removal\n');
            end



            %% Epoch data
            if (do_epoch)
                % Porper dataset will be loaded before epoching
                
                % 2 epoch total T1(L), T2(R)
                %  T1 - 12628 - [0 4.1] 
                %  T2 - 12884 - [0 4.1]

                % Sides information 
                %T = {'12628', '12884'};
                T = {'2', '3'};

                P = [0 4.1];
        
                for si = 1:numsides
                    for pp = 1:numprepost
                        %% Epoch                        
                        % Load Proper Dataset
                        EEG = pop_loadset('filename', [savename '_EpochReady.set'],'filepath',filepath);
              
                        % Epoch
                        EEG = pop_epoch(EEG, { T{si} }, P(pp,:) , 'epochinfo', 'yes');
                        
                        % ICA 
                        %EEG = pop_runica(EEG, 'icatype','runica','dataset',1,'options',{'extended' 1},'chanind',chanid );
                        %EEG.icaact % Not finding icaact ???
                        
                        % Save the dataset for each epoched set
                        fprintf('\nSaving %s %s %s\n\n', savename, sides{si}, prepost{pp});
                        EEG = pop_editset(EEG, 'setname', [savename '_' sides{si} '_' prepost{pp} '']);
                        EEG = pop_saveset(EEG, 'filename',...
                            [savename '_' sides{si} '_' prepost{pp} '_ICA.set'],'filepath',filepath);

                    end
                end

            else
                fprintf('Skipped Epochs and ICA\n');
            end % end epoch extraction


            %% Feature Extraction Wavelets
            if (do_features_wavelet)
                
                % Features for wavelets
                % Features for each channel:
                for si = 1: numsides % T1 and T2
                    for ty = 1: numtypes % ERD, ERS, and MRCP
                        
                        fprintf('\nPrepping wavelet features for %s_%s_%s\n',savename, sides{si}, types{ty});                        
                        
                        % Load Proper Dataset
                        EEG = pop_loadset('filename',...
                            [savename '_' sides{si} '_' types{ty} '.set'],'filepath',filepath);                         
                        
                        EEG.icaact = eeg_getdatact(EEG, 'channel', chanid);
                        
                        %wavchan = 3;
                        for e = 1: numepochs % For the 7 epochs
                            % Preallocate zeros to avoid resizing matrix each iteration
                            wav1 = zeros(1,numchan);
                            wav2 = zeros(1,numchan);

                            % Loop through each channel to find desired value
                            % icaact(channel, data/channel, epoch number)
                            for n = 1:numchan % channel number

                                fprintf('*** channel %d epoch %d', n, e);

                                tempwave = EEG.icaact(n,:,e);
                                
                                [C,L] = wavedec(tempwave,3,'db1');
                                [cD1,cD2,cD3] = detcoef(C,L,[1,2,3]);

                                cD3sort = sort(cD3,'descend');
                                cD3sort(1,1:5)

                                
                                wav1(n) = mean(tempwave);
                                wav2(n) = sum(tempwave.^2);
                            end
                            
                            % Set all values between 0.1 and 0.9 and transpose matrix
                            mi = 0.1;  mx = 0.9;
                            wav1 = mapminmax(wav1, mi, mx)';
                            wav2 = mapminmax(wav2, mi, mx)';


                            % Set side to proper value
                            if strcmp(sides{si},'T1')
                                side = 0;
                            elseif strcmp(sides{si},'T2')
                                side = 1;
                            else
                                side = 5;
                            end

                            %fprintf('\nFeature set for %s %s %s\n', savename, sides{m}, types{u});
                            feature = [wav1; wav2; wav2; type; side];
                            
                            r = ty + numtypes*(si-1); % for indexing
                            ind = e + numepochs*(r - 1) + numruns*numepochs*(t - 1) + numruns*numtrials*numepochs*(s - 1);
                            features(:,ind) = feature;
                               
                        end
                        % uncomment this to make sure features are being
                        % extracted correctly
                        %start = ind - 6 %numepochs*((si - 1)*numtypes + (ty - 1)) + 1
                        %finish = start + numepochs - 1
                        %features(:,start:finish)
                    end
                end
                
            else
                fprintf('Skipped Feature selection wavelets\n');
            end

        end
    end
end


if (do_NN && do_features)
    Avg = (1:8); Pow = (9:16); Ene = (17:24);
    Typ = 25; Tar = 26;
    range = [Avg; Pow; Ene];

    tot = sum(feat2use);
    k = 1;

    for rng = 1 : tot   

        while feat2use(k) ~= 1
            fprintf('Skipped K')
            k = k+1;
        end

        nninputs(range(rng,:),:) = features(range(k,:),:);
        fprintf('range = %d', range(k,1))
        k = k+1;
    end

    [row, col] = size(nninputs);

    nninputs(row+1,:) = features(Typ,:);
    %nninputs(25,:) = mapminmax(nninputs(25,:),mi, mx);
    nntargets = features(Tar,:);
    %nntargets = mapminmax(nntargets,mi,mx)
    %size(nninputs)
    %size(nntargets)
    save('nninoutenergy', 'nninputs', 'nntargets')
else
    fprintf('NO NN Inputs/outputs\n');
end

