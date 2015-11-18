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
%       Epoch Data
%           - EEG data is epoched based on event annotations
%       Find Features
%           - Finds the features and preps them for input into NN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set variables
%chanid = [9 11 13];
chanid = [2 4 6 9 10 11 12 13];
numchan = length(chanid);
 
% Trials to be used for analysis
%trial = {'R04', 'R08', 'R12'};
trial = {'R03', 'R07', 'R11'};
numtrials = length(trial);  
%numtrials = 2; % Change for testing

% Subjects
subject = {'S001','S002','S003','S004','S005','S006','S007','S008','S009','S010',...
    'S011','S012','S013','S014','S015','S016','S017','S018','S019','S020'};
%numsubjects = length(subject);  
numsubjects = 6; % Change for testing

% Set the epoch time and find how many data points there are
epochtime = [0 4.0];
datapoints = 160 * (epochtime(2) - epochtime(1));

% To perform wavelet decomposition to cA3 must be evenly divisible to 3/8 of datapoints
wavecoefnum = datapoints * 3 / 8;

% Sides - T1=>Left, T2=>Right
sides = {'T1', 'T2'};
numsides = length(sides);

% Sides information 
%T = {'12628', '12884'};
T = {'2', '3'};

% Runs = sides
numruns = numsides;

% Define paths for data location and storage
homedir = pwd;
homepath = fullfile(homedir, '/Data/PhysionetData/EDF/'); % Location of EDF files
filepath = fullfile(homedir, '/Data/userscriptswave/'); % Location of storage folder
%homepath = '~/thesis/Data/PhysionetData/EDF/'; % Location of EDF files
%filepath = '~/thesis/Data/userscriptswave/'; % Location of storage folder

fprintf('********************************************************\n');
fprintf('\nPerforming on data in: %s with %d subjects',filepath, numsubjects);
fprintf('\nEpoch time = %.1f to %.1f with %d datapoints', epochtime, datapoints);
fprintf('\nWith that many datapoints there are %d points to work on up to third level\n', wavecoefnum);
fprintf('********************************************************\n\n');

%% Set each variable depending if process needs to be done
do_import = 0; % If 0 will not edit channels, filter or artifact removal
do_edit_channels = 1;
do_edit_events = 1;
do_filter = 1;
do_artifact_removal = 1;

do_epoch = 0;

% Features - using wavelets - either highest coef or avg of x coef
numWavAvg = 10;
numwavcoef = 10;
do_features_wavelet = 1;
do_wave_avg = 0;
do_NN_wave = 1;
%% Prealocate size for features (num of features, number of feature sets)
if (do_features_wavelet && do_wave_avg)
    features = zeros((wavecoefnum/numWavAvg)*numchan + 1, 15*numtrials*numsubjects);
    fprintf('Performing Wave Averging of every %d coefficients\n', numWavAvg);
    fprintf('Feature size is: [%d, %d]\n', size(features));
    featSize = size(features);

elseif (do_features_wavelet)
    features = zeros(numwavcoef*2*numchan + 1, 15*numtrials*numsubjects);
    fprintf('Performing Top %d Wave Coefficients\n', numwavcoef);
    fprintf('Feature size is: [%d, %d]\n', size(features));
    featSize = size(features);

else
    fprintf('No features being determined');
end


%% Start looping through subjects and trials
for s = 1:numsubjects
    for t = 1:numtrials
        fprintf('\n******** Processing Subject: %s, Trial: %s\n\n', subject{s}, trial{t});

        % Finish setting the path to where files located and to where files are written
        foldername = [homepath subject{s} '/']; % Add to path the subject
        savename = [subject{s} trial{t}]; % Full name of each file for subject and trial
        filename = [foldername savename '.edf']; % Full name of the file to retrieve

              
        % Check if the file exists
        if exist(filename, 'file') <=0
            error('File for %s %s does not exist', subject{s}, trial{t});
            
        else % If the file exists then perform processing

            %% Import Data
            if (do_import) % Import using BIOSIG - for EDF files
                EEG = pop_biosig(filename);
            else
                fprintf('Skip Import  ');
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
                        error('\nTypes are not the same');
                    end
                    
                    if annotations(kk,4) ~= events_found(kk,4)
                        error('\nDurations are not the same');
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
                %EEG = pop_eegfiltnew(EEG, 0.5, 50, 1056, 0, [], 0);
                EEG = pop_eegfiltnew(EEG, 0.5, 50);

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
                    % SaveName: Subject{s}Trial{t}_EqochReady4s.set
                EEG = pop_editset(EEG, 'setname', [savename '_EpochReady4s']);
                EEG = pop_saveset(EEG,'filename',[savename '_EpochReady4s.set'],'filepath',filepath);
            else
                fprintf(' Skip Artifact Removal\n');
            end



            %% Epoch data and save the dataset
            if (do_epoch)

                for si = 1:numsides                     
                    % Load Proper Dataset
                    EEG = pop_loadset('filename', [savename '_EpochReady4s.set'],'filepath',filepath);
              
                    % Epoch
                    EEG = pop_epoch(EEG, { T{si} }, epochtime , 'epochinfo', 'yes');

                    % Save the dataset for each epoched set
                    fprintf('\nSaving %s %s\n\n', savename, sides{si});
                    EEG = pop_editset(EEG, 'setname', [savename '_' sides{si} '_FeatReady4s']);
                    EEG = pop_saveset(EEG, 'filename',...
                        [savename '_' sides{si} '_FeatReady4s.set'],'filepath',filepath);

                end

            else
                fprintf('Skipped Epochs\n');
            end % end epoch extraction
            


            %% Feature Extraction Wavelets
            if (do_features_wavelet)
                
                % Features for wavelets
                % Features for each channel:
                r = 0;
                for si = 1: numsides % T1 and T2 
                    fprintf('\nPrepping wavelet features for %s_%s\n',savename, sides{si});                        
                    
                    % Load Proper Dataset
                    EEG = pop_loadset('filename',...
                        [savename '_' sides{si} '_FeatReady4s.set'],'filepath',filepath);   

                    numepochs = EEG.trials;

                    for e = 1: numepochs % For the 7 or 8 epoch
                        % Preallocate zeros to avoid resizing matrix each iteration
                        if (do_wave_avg)
                            numAvgs = wavecoefnum/numWavAvg;
                            wavT = zeros(1,numAvgs*numchan);
                        else
                            wavT = zeros(1,numwavcoef*numchan*2); % mult by 2 for cD2 and cD3
                        end

                        % Loop through each channel to find desired value
                        for n = 1:numchan % channel number

                            %fprintf('\nChannel ID is: %d and epoch is: %d\n', chanid(n), e);
                            tempwave = EEG.data(chanid(n),:,e);
                            
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

                                indWav = (n-1) * numAvgs;
                                wavT(1, indWav + 1 : indWav + numAvgs) = wav(1, 1 : numAvgs);

                            else  % Take X highest amount of coefficients
                                cD3sort = sort(cD3,'descend');
                                cD2sort = sort(cD2,'descend');

                                wav = zeros(1, numwavcoef*2);

                                wav(1, 1             : numwavcoef)   = cD3sort(1,1:numwavcoef);
                                wav(1, numwavcoef + 1: numwavcoef*2) = cD2sort(1,1:numwavcoef);

                                indWav = (n-1) * numwavcoef*2;
                                wavT(1, indWav + 1 : indWav + numwavcoef*2) = wav(1, 1 : numwavcoef*2);
                            end

                        end

                        % Set all values between 0.1 and 0.9 and transpose matrix
                        mi = 0.1;  mx = 0.9;
                        wavT = mapminmax(wavT, mi, mx)'; %'

                        % Set side to proper value
                        if strcmp(sides{si},'T1')
                            side = 0;
                        elseif strcmp(sides{si},'T2')
                            side = 1;
                        else
                            error('Side for %s %s is incorrect', subject{s}, trial{t});
                        end

                        feature = [wavT; side];

                        ind = e + r + 15*(t - 1) + 15*numtrials*(s - 1);
                        features(:,ind) = feature;
                           
                    end
                    r = numepochs;
                end
                
            else
                fprintf('Skipped Feature selection wavelets\n');
            end

        end
    end
end

% Print Procedure
if (do_features_wavelet && do_wave_avg)
    fprintf('\nPerforming Wave Averging of every %d coefficients\n', numWavAvg);

elseif (do_features_wavelet)
    fprintf('\nPerforming Top %d Wave Coefficients\n', numwavcoef);

else
    fprintf('\nNo features being determined');
end

%% NN for wavelet features
if (do_NN_wave && do_features_wavelet)

    Tar = featSize(1);
    nninputs = features(1:Tar - 1,:);
    nntargets = features(Tar,:);

    nninputs(1:9,1:5);
    sizenninputs = size(nninputs);
    sizenntargets = size(nntargets);

    save('nninoutenergy', 'nninputs', 'nntargets')

    TNNLoop

else
    fprintf('NO NN Inputs/outputs for wavelets\n');
end
