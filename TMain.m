%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Thesis Work
%   Run after running TStart script (starts EEGLAB toolbox)
%   Requires pop_biosig and AAR toolboxes to be installed
%   
%   Program consists of the following operations:
%       Set Variables
%           - Sets the variables for what work will be done on (can choose)
%       Import Data (pop_biosid toolbox if working with EDF+)
%           - Imports data and sets it equal to EEG
%       Edit Channels
%           - Edits the names of channels so that EEGLAB can automatically 
%             find their location
%           - Currently only works if 64 channnels imported
%       Edit Events
%           - Edits the event annotations 
%               - Currently events not being imported correctly and are edited
%                 here to show correct location
%       Filter Data
%           - Filter all data in EEG (choose frequencies)
%       EOG and EMG filtering (AAR toolbox)
%           - EOG filtering
%               - removes eye movement/blink artifacts
%           - EMG filtering
%               - removes muscle movement artifacts
%       Epoch and ICA Data
%           - EEG data is epoched based on event annotations
%           - ICA is done depending on channels chosen
%       Rhythm Isolate
%           - Filter data into:
%               - ERD 
%               - ERS
%               - MRCP
%       Find Features
%           - Finds the features and preps them for input into NN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set variables
% Channels used for ICA
chanid = [2 4 6 9 10 11 12 13];
%chanid - [9 11 13];
numchan = length(chanid);

% Trials to be used for analysis
trial = {'R03', 'R07', 'R11'};
numtrials = length(trial);  
%numtrials = 2; % Change for testing

% Subjects
subject = {'S001','S002','S003','S004','S005','S006'};
numsubjects = length(subject);  
%numsubjects = 1; % Change for testing

% Pre or post for epochs
prepost = {'PRE', 'POST'};
numprepost = length(prepost);

% PRE and POST information
PRE = [-2 0]; 
POST = [4.1 6.1];
P = [PRE; POST];

% Types - for use in epoching data
types = {'ERD', 'ERS', 'MRCP'};
numtypes = length(types);
%numtypes = 2;

% Sides - T1=>Left, T2=>Right
sides = {'T1', 'T2'};
numsides = length(sides);

% Sides information - event tag 
%T = {'12628', '12884'};
T = {'2', '3'};

% Runs = types x sides
numruns = numtypes * numsides;

% Epochs - set to 7 even though one will have 8 due to it refrencing out of bounds
%   If set to 1 then only analyzing 1 epoch, unless only checking the icaact
%   This does not effect how many epochs are created only the number to be analyzed
numepochs = 7;

% Define paths for data location and storage
homepath = '~/thesis/PhysionetData/EDF/'; % Location of EDF files
filepath = '~/thesis/userscripts/'; % Location of storage folder

fprintf('\nPerforming on data in: %s',filepath);
fprintf('\nPre = %.1f to %.1f\nPost = %.1f to %.1f\n', PRE, POST);
fprintf('\nNumber of types = %d\n, number of epochs = %d\n', numtypes, numepochs);


%% Set each variable depending if process needs to be done
do_import = 0; % If 0 will not edit channels, filter or artifact removal
do_edit_channels = 1;
do_edit_events = 1;
do_filter = 1;
do_artifact_removal = 1;

do_epoch_ICA = 0;
do_rhythmiso = 0;

% Features - using [avg, pwr, ene]
feat2use = [1,0,0];
do_features = 0;
do_NN = 0;

% Features - using wavelets - either highest coef or avg of x coef
numWavAvg = 10;
numwavcoef = 10;
do_features_wavelet = 1;
do_wave_avg = 1;
do_NN_wave = 1;


%% Prealocate size for features (num of features, number of feature sets)
if (do_features_wavelet && do_wave_avg)
    features = zeros((120/numWavAvg)*numchan + 2, numepochs*numruns*numtrials*numsubjects);
    fprintf('Performing Wave Averging of every %d coefficients\n', numWavAvg);
    featSize = size(features)

elseif (do_features_wavelet)
    features = zeros(numwavcoef*2*numchan + 2, numepochs*numruns*numtrials*numsubjects);
    fprintf('Performing Top %d Wave Coefficients\n', numwavcoef);
    featSize = size(features)

elseif (do_features)
    features = zeros(24 + 2, numepochs*numruns*numtrials*numsubjects);
    fprintf('Performing avg = %d, pwr = %d, and ene = %d\n', feat2use);
    featSize = size(features)
else
    fprintf('No features being determined');
end



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
            %fprintf('Warning file does not exist');
            error('File for %s %s does not exist', subject{s}, trial{t});
            
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
                
                for kk = 1: num_annotations 
                    annotations(kk,1) = kk; % urevent number
                    annotations(kk,2) = EEG.event(kk).type;
                    annotations(kk,3) = (EEG.event(kk).latency);
                    annotations(kk,4) = (EEG.event(kk).duration);
                end

                % Annotate events based on function
                events_found = event_finder(duration1, duration2, sub1, trial1);
                
                % Check to ensure at minimum variables known are true in new events
                for kk = 1: num_annotations
                    if annotations(kk,2) ~= events_found(kk,2)
                        error('Types are not the same');
                    end
                    
                    if annotations(kk,4) ~= events_found(kk,4)
                        error('Durations are not the same');
                    end   
                end

                fprintf('\nChanging info for EEG event and urevent\n');
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
                %EEG = pop_eegfiltnew(EEG, 0.5, 45, 1056, 0, [], 0);
                EEG = pop_eegfiltnew(EEG, 0.5, 45);
                
                % Low edge of frequency filter: 0.5Hz
                %EEG = pop_iirfilt( EEG, 0.5, 0, [], [0]);
                %EEG = pop_eegfiltnew(EEG, [], 0.5, 1056, true, [], 0);

                % High edge of frequency filter: 45Hz 
                    % If over 50 Hz then need to add Notch filter
                %EEG = pop_iirfilt( EEG, 0, 45, [], [0]);
                %EEG = pop_eegfiltnew(EEG, [], 45, 48, 0, [], 0);
            else
                fprintf('Skip Freq Filtering  ');
            end
            

            %% EOG (eye movements/blinks) and EMG (Muscle artifacts) filtering and saving data
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


            %% Epoch data, perform ICA, and save the dataset
            if (do_epoch_ICA)

                for si = 1:numsides
                    for pp = 1:numprepost                        
                        % Load Proper Dataset
                        EEG = pop_loadset('filename', [savename '_EpochReady.set'],'filepath',filepath);
              
                        % Epoch into pre-right, post-right, pre-left, and post-left
                        EEG = pop_epoch(EEG, { T{si} }, P(pp,:) , 'epochinfo', 'yes');
                        
                        % ICA 
                        EEG = pop_runica(EEG, 'icatype','runica','dataset',1,'options',{'extended' 1},'chanind',chanid );

                        $ Check to see if data in EEG.icaact
                        %size(EEG.icaact)
                        %EEG = eeg_checkset(EEG, 'ica');
                        %size(EEG.icaact)
                        %EEG.icaact(1,1:10,1) % Not finding icaact ???
                        
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
            

            %% Rhythm Isolation
            if (do_rhythmiso)
                              
                for si = 1:numsides
                    for ty = 1:numtypes
                        fprintf('\nRhythm Isolation for %s_%s_%s\n',savename, sides{si}, types{ty});
                        
                        % k value for 3 uses PRE as well
                        if ty == 3
                            k = 1;
                        else
                            k = ty;
                        end
                        
                        % Load proper dataset
                        EEG = pop_loadset('filename', [savename '_' sides{si} '_' prepost{k} '_ICA.set'],'filepath',filepath);
                       
                        if ty == 3
                            % Low Pass 3 Hz filter
                            %EEG = pop_eegfiltnew(EEG, [], 3, 264, 0, [], 0);
                            EEG = pop_eegfiltnew(EEG, [], 3);

                            %EEG = pop_iirfilt( EEG, 0, 3, [], [0]);
                        else
                            % bandpass 8-30 Hz filter 
                            %EEG = pop_eegfiltnew(EEG, 8, 30, 264, 0, [], 0);
                            EEG = pop_eegfiltnew(EEG, 8, 30);

                            % Low edge of frequency filter: 8Hz
                            %EEG = pop_iirfilt( EEG, 8, 0, [], [0]);

                            % High edge of frequency filter: 30Hz
                            %EEG = pop_iirfilt( EEG, 0, 30, [], [0]);
                        end

                        % Delete old activation and scalp maps
                        %EEG.icaact = [];
                        %EEG.icawinv = []; % IF DELETED NEED TO FILL IN HOW???

                        % recalculate acts and winv using specified channels
                        %EEG.icaact = eeg_getdatact(EEG, 'channel', chanid);

                        % Save the dataset (as)
                        EEG = pop_editset(EEG, 'setname',  [savename '_' sides{si} '_' types{ty}]);
                        EEG = pop_saveset(EEG, 'filename', [savename '_' sides{si} '_' types{ty} '.set'],'filepath',filepath);
                    end
                end
                
            else
                fprintf('Skipped Rhythm Isolation\n');
            end % end rhythm isolation
            

            %% Feature Extraction
            if (do_features)
                %% Features
                % Features for each channel:
                    % 8 for average
                    % 8 for power
                    % 8 for energy
                    % 1 for type
                    % 1 for side
                
                for si = 1: numsides % T1 and T2
                    for ty = 1: numtypes % ERD, ERS, and MRCP
                        
                        fprintf('\nFinding features for %s_%s_%s\n',savename, sides{si}, types{ty});                        
                        
                        % Load Proper Dataset
                        EEG = pop_loadset('filename',...
                            [savename '_' sides{si} '_' types{ty} '.set'],'filepath',filepath);                         
                        
                        %EEG.icaact = eeg_getdatact(EEG, 'channel', chanid);
                        %act1 = eeg_getdatact(EEG, 'channel', chanid);
                        %size(act1)
                        %act1(1,1:10,1)

                        % As defined in line 977 of EEG_checkset EEG.icaact is:
                        EEG.icaact = (EEG.icaweights * EEG.icasphere) * EEG.data(EEG.icachansind,:);
                        %size(EEG.icaact)
                        %EEG.icaact(1,1:10)

                        act = EEG.icaweights * EEG.icasphere;
                        
                        for e = 1: numepochs  % Loop through the 7 epochs
                            % Preallocate zeros to avoid resizing matrix each iteration
                            avg = zeros(1,numchan);  pwr = zeros(1,numchan);  ene = zeros(1,numchan);
                            
                            for n = 1:numchan  % Loop throudh each channel
                                % icaact(channel, data/channel, epoch number)
                                %avg(n) = mean(mean(EEG.icaact(n,:,:)));
                                %avg(n) = mean(EEG.icaact(n,:,e));
                                avg(n) = mean(EEG.icaact(n,(e-1)*320+1:e*320));
                                    %avg(n) = mean(act(n,:));
                                %pwr(n) = bandpower(EEG.icaact(n,:,1), 320, [8 30]);
                                %pwr(n) = bandpower(EEG.icaact(n,:,e)', length(EEG.icaact), [0 30], 1, 3); %'
                                %pwr(n) = sum(EEG.icaact(n,:,e).^2);
                                pwr(n) = sum(EEG.icaact(n,(e-1)*320+1:e*320).^2);
                                    %pwr(n) = sum(act(n,:).^2);
                                %ene(n) = sum(sum(EEG.icaact(n,:,:).^2));
                                %ene(n) = sum(EEG.icaact(n,:,e).^2);
                                ene(n) = sum(EEG.icaact(n,(e-1)*320+1:e*320).^2);
                                    %ene(n) = sum(act(n,:).^2);
                            end
                            
                            % Set all values between 0.1 and 0.9 and transpose matrix
                            mi = 0.1;  mx = 0.9;
                            avg = mapminmax(avg, mi, mx)'; %'
                            pwr = mapminmax(pwr, mi, mx)'; %'
                            ene = mapminmax(ene, mi, mx)'; %'

                            % Set side to proper value
                            if strcmp(sides{si},'T1')
                                side = 0;
                            elseif strcmp(sides{si},'T2')
                                side = 1;
                            else
                                error('Side for %s %s is incorrect', subject{s}, trial{t}); 
                            end

                            % Set type to proper value
                            if strcmp(types{ty},'ERD')
                                type = 1;
                            elseif strcmp(types{ty},'ERS')
                                type = 2;
                            elseif strcmp(types{ty},'MRCP')
                                type = 3;
                            else
                                error('Type for %s %s is incorrect', subject{s}, trial{t});
                            end

                            feature = [avg; pwr; ene; type; side];
                            
                            r = ty + numtypes*(si-1); % for indexing
                            ind = e + numepochs*(r - 1) + numruns*numepochs*(t - 1) + numruns*numtrials*numepochs*(s - 1);
                            features(:,ind) = feature;
                                                     
                        end
                    end
                end
      
            else
                fprintf('Skipped Feature selection\n');
            end % end feature extraction
            

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
                        
                        %EEG.icaact = eeg_getdatact(EEG, 'channel', chanid);
                        %act1 = eeg_getdatact(EEG, 'channel', chanid);
                        %size(act1)
                        %act1(1,1:10,1)

                        % As defined in line 977 of EEG_checkset EEG.icaact is:
                        EEG.icaact = (EEG.icaweights * EEG.icasphere) * EEG.data(EEG.icachansind,:);
                        %size(EEG.icaact)
                        %EEG.icaact(1,1:10)

                        %act = EEG.icaweights * EEG.icasphere;

                        for e = 1: numepochs % For the 7 epochs
                            % Preallocate zeros to avoid resizing matrix each iteration
                            if (do_wave_avg)
                                numAvgs = 120/numWavAvg;
                                wavT = zeros(1,numAvgs*numchan);
                            else
                                wavT = zeros(1,numwavcoef*numchan*2); % mult by 2 for cD2 and cD3
                            end

                            % Loop through each channel to find desired value
                            for n = 1:numchan % channel number
                                % icaact(channel, data per channel, epoch number)
                                %tempwave = EEG.icaact(n,:,e);
                                tempwave = EEG.icaact(n,(e-1)*320+1:e*320);
                                
                                % Do wavelet analysis to the third level
                                [C,L] = wavedec(tempwave,3,'db1');
                                [cD1,cD2,cD3] = detcoef(C,L,[1,2,3]);

                                if (do_wave_avg)  
                                    lencD2 = length(cD2);
                                    lencD3 = length(cD3);

                                    k = lencD2/numWavAvg;
                                    l = lencD3/numWavAvg;
                                    wav = zeros(1,k+l);

                                    for kk = 1:k  % average of cD2 channels
                                        indWav = (kk-1)*numWavAvg;
                                        wav(kk) = mean(cD2(indWav + 1 : indWav + numWavAvg));
                                    end

                                    for ll = 1:l  % average of cD3 channels
                                        indWav = (ll-1)*numWavAvg;
                                        wav(kk+ll) = mean(cD3(indWav + 1 : indWav + numWavAvg));
                                    end

                                    indWav = (n-1)*numAvgs;
                                    wavT(1, indWav + 1 : indWav + numAvgs) = wav(1, 1 : numAvgs);

                                else  % Take X highest amount of coefficients
                                    cD3sort = sort(cD3,'descend');
                                    cD2sort = sort(cD2,'descend');

                                    wav = zeros(1, numwavcoef*2);

                                    wav(1, 1             : numwavcoef)   = cD3sort(1,1:numwavcoef);
                                    wav(1, numwavcoef + 1: numwavcoef*2) = cD2sort(1,1:numwavcoef);

                                    indWav = (n-1)*numwavcoef*2;
                                    wavT(1, indWav + 1 : indWav + numwavcoef*2) = wav(1, 1 : numwavcoef*2);
                                end

                            end

                            %wavTsee = wavT

                            % Set all values between 0.1 and 0.9 and transpose matrix
                            mi = 0.1;  mx = 0.9;
                            wavT = mapminmax(wavT, mi, mx)'; %'

                            %wavTseeAfter = wavT' % '
                            
                            % Set type to proper value
                            if strcmp(types{ty},'ERD')
                                type = 1;
                            elseif strcmp(types{ty},'ERS')
                                type = 2;
                            elseif strcmp(types{ty},'MRCP')
                                type = 3;
                            else
                                error('Side for %s %s is incorrect', subject{s}, trial{t});
                            end

                            % Set side to proper value
                            if strcmp(sides{si},'T1')
                                side = 0;
                            elseif strcmp(sides{si},'T2')
                                side = 1;
                            else
                                error('Type for %s %s is incorrect', subject{s}, trial{t});
                            end

                            %fprintf('\nFeature set for %s %s %s\n', savename, sides{m}, types{u});
                            feature = [wavT; type; side];
                            
                            r = ty + numtypes*(si-1); % for indexing
                            ind = e + numepochs*(r - 1) + numruns*numepochs*(t - 1) + numruns*numtrials*numepochs*(s - 1);
                            features(:,ind) = feature;
                               
                        end
                    end
                end
                
            else
                fprintf('Skipped Feature selection wavelets\n');
            end

        end
    end
end


%% NN for wavelet features
if (do_NN_wave && do_features_wavelet)

    Tar = featSize(1);
    nninputs = features(1:Tar - 1,:);
    nntargets = features(Tar,:);

else
    fprintf('NO NN Inputs/outputs for wavelets\n');
end


%% NN for avg, pow, and ene 
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
        % fprintf('range = %d', range(k,1))
        k = k+1;
    end

    [row, col] = size(nninputs);
    nninputs(row+1,:) = features(Typ,:);
    nntargets = features(Tar,:);

else
    fprintf('NO NN Inputs/outputs\n');
end


%% Save NN inputs and outputs and start NN training
if (do_NN | do_NN_wave)
    nninputs(1:9,1:5)
    sizenninputs = size(nninputs)
    sizenntargets = size(nntargets)
    save('nninoutenergy', 'nninputs', 'nntargets')

    TNNLoop
end
