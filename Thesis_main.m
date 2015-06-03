%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Thesis Work
%    Run after starting bcilab and eeglab
%    Test to see how basic filtering looks
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set variables
%clc
% Channels used for ICA
chanid = [9 11 13];
numchan = length(chanid);

% Trials to be used for analysis
trial = {'R03', 'R07', 'R11'};
numtrials = length(trial);  
%numtrials = 2; % Change for testing

% Subjects
subject = {'S001','S002','S003','S004','S005','S006'};
numsubjects = length(subject);  
%numsubjects = 1; % Only using 1 for testing

% Types - for use in epoching data
types = {'ERD', 'ERS'}; %, 'MRCP'};
numtypes = length(types);

% Sides - T1=>Left, T2=>Right
sides = {'T1', 'T2'};
numsides = length(sides);

% Runs = types x sides
numruns = numtypes * numsides;

% Epochs
numepochs = 1;

% Pre or post for epochs
prepost = {'PRE'};%, 'POST'};
numprepost = length(prepost);

% Define paths for data location and storage
% Paths for folders were data located and where to store data
homepath = '/home/amart/Physionet Data/EDF/'; % Location of EDF files
filepath = '/home/amart/BCILAB-master/userscripts/'; % Location of storage folder

% Prealocate size for features
features = zeros(26, numepochs*numruns*numtrials*numsubjects);

% Set each variable depending if process needs to be done
do_import = 1; % If 0 will not edit channels, filter or artifact removal
do_edit_channels = 1;
do_filter = 1;
do_artifact_removal = 1;

do_epoch = 1;
do_rhythmiso = 0;
do_features = 0;

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
                fprintf('\nSkipped Import\n');
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
                fprintf('\nSkipped Channel Locations\n');
            end
            
            %% Filter data
            if (do_filter && do_import && do_edit_channels)
                % Low edge of frequency filter: 0.5Hz
                EEG = pop_iirfilt( EEG, 10, 0, [], [0]);

                % High edge of frequency filter: 64Hz
                EEG = pop_iirfilt( EEG, 0, 20, [], [0]);
            else
                fprintf('\nSkipped Filtering\n');
            end
            
            %% EOG and EMG filtering and saving data
            if (do_artifact_removal && do_import && do_edit_channels && do_filter)
                % EOG removal using BSS - sobi - ...
%                 EEG = pop_autobsseog( EEG, [123], [123], 'sobi',...
%                     {'eigratio', [1000000]}, 'eog_fd', {'range',[2 21]});
% 
%                 % EMG removal using BSS
%                 EEG = pop_autobssemg( EEG, [81.9188], [81.9188], 'bsscca',...
%                     {'eigratio', [1000000]}, 'emg_psd', {'ratio', [10], 'fs', [160], 'femg', [15],...
%                     'estimator', spectrum.welch({'Hamming'}, 80),'range', [0 32]});

                % Save the dataset as it is ready to be epoched
                %fprintf('\nSaving dataset as: %s_EpochReady\n', savename)
                EEG = pop_editset(EEG, 'setname', [savename '_EpochReady']);
                EEG = pop_saveset(EEG,'filename',['aaa' savename '_EpochReady.set'],'filepath',filepath);
            else
                fprintf('\nSkipped Artifact Removal\n');
            end

            %% Epoch data
            if (do_epoch)
                % Porper dataset will be loaded before epoching
                % Each epoch will be saved after going through ICA
                
                % 4 epochs total T1(L), T2(R)
                %  T1 - 12628 - ERD/MRCP - [-2 0] 
                %  T2 - 12884 - ERD/MRCP - [-2 0]
                %  T1 - 12628 - ERS - [4.1 6.1] 
                %  T2 - 12884 - ERS - [4.1 6.1]

                % Sides information 
                T = {'12628', '12884'};
                
                % PRE and POST information
                PRE = [-0.5 2.5]; 
                POST = [4.1 6.1];
                P = [PRE; POST];
        
                for si = 1:numsides
                    for pp = 1:numprepost
                        %% Epoch                        
                        % Load Proper Dataset
                        EEG = pop_loadset('filename', ['aaa' savename '_EpochReady.set'],'filepath',filepath);
              
                        % Epoch
                        EEG = pop_epoch(EEG, { T{si} }, P(pp,:) , 'epochinfo', 'yes');
                        
                        % ICA 
                        %EEG = pop_runica(EEG, 'icatype','runica','dataset',1,'options',{'extended' 1},'chanind',chanid );

                        % Save the dataset for each epoched set
                        %fprintf('\nSaving %s %s %s\n\n', savename, sides{m}, prepost{u});
                        EEG = pop_editset(EEG, 'setname', [savename '_' sides{si} '_' prepost{pp} '']);
                        EEG = pop_saveset(EEG, 'filename',...
                            ['aaa' savename '_' sides{si} '_' prepost{pp} '.set'],'filepath',filepath);

                    end
                end

            else
                fprintf('\nSkipped Epochs\n');
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
                            EEG = pop_iirfilt( EEG, 0, 3, [], [0]);
                        else
                            % Low edge of frequency filter: 8Hz
                            EEG = pop_iirfilt( EEG, 8, 0, [], [0]);
                            % High edge of frequency filter: 30Hz
                            EEG = pop_iirfilt( EEG, 0, 30, [], [0]);
                        end

                        % Delete old activation and scalp maps
                        EEG.icaact = [];
                        % EEG.icawinv = []; % IF DELETED NEED TO FILL IN HOW??????????????

                        % recalculate acts and winv using specified channels
                        EEG.icaact = eeg_getdatact(EEG, 'channel', chanid);

                        % Save the dataset (as)
                        EEG = pop_editset(EEG, 'setname', [savename '_' sides{si} '_' types{ty}]);
                        EEG = pop_saveset(EEG, 'filename', [savename '_' sides{si} '_' types{ty} '.set'],'filepath',filepath);
                    end
                end
                
            else
                fprintf('\nSkipped Rhythm Isolation\n');
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
                        
                        for e = 1: numepochs % For the 7 epochs
                            % Preallocate zeros to avoid resizing matrix each iteration
                            avg = zeros(1,numchan);
                            pwr = zeros(1,numchan);
                            ene = zeros(1,numchan);
                            
                            % icaact(channel, data/channel, epoch number)
                            for n = 1:numchan % channel number
                                %avg(n) = mean(mean(EEG.icaact(n,:,:)));
                                avg(n) = mean(EEG.icaact(n,:,e));
                                %pwr(n) = bandpower(EEG.icaact(n,:,1), 320, [8 30]);
                                pwr(n) = bandpower(EEG.icaact(n,:,e), length(EEG.icaact), [0 30]);
                                %ene(n) = sum(sum(EEG.icaact(n,:,:).^2));
                                ene(n) = sum(EEG.icaact(n,:,e).^2);
                            end
                            
                            % Set all values between 0.1 and 0.9 and transpose matrix
                            mi = 0.1;  mx = 0.9;
                            avg = mapminmax(avg, mi, mx)';
                            pwr = mapminmax(pwr, mi, mx)';
                            ene = mapminmax(ene, mi, mx)';

                            % Set side to proper value
                            if strcmp(sides{si},'T1')
                                side = 0;
                            elseif strcmp(sides{si},'T2')
                                side = 1;
                            else
                                side = 5;
                            end

                            % Set type to proper value
                            if strcmp(types{ty},'ERD')
                                type = 1;
                            elseif strcmp(types{ty},'ERS')
                                type = 2;
                            elseif strcmp(types{ty},'MRCP')
                                type = 3;
                            else
                                type = 5;
                            end

                            %fprintf('\nFeature set for %s %s %s\n', savename, sides{m}, types{u});
                            feature = [avg; pwr; ene; type; side];
                            % feature';
                            
                            r = ty + numtypes*(si-1); % for indexing
                            ind = e + numepochs*(r - 1) + numruns*numepochs*(t - 1) + numruns*numtrials*numepochs*(s - 1);
                            features(:,ind) = feature;
                            
                            
                        end
                        % uncomment this to make sure features are being
                        % extrcted correctly
                        %start = ind - 6 %numepochs*((si - 1)*numtypes + (ty - 1)) + 1
                        %finish = start + numepochs - 1
                        %features(:,start:finish)
                    end
                end
      
            else
                fprintf('\nSkipped Feature selection\n');
            end % end feature extraction

        end
    end
end

% Print the features
%features(:,1) = [1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,1];
%features;
nninputs = features(9:25,:);
%nninputs(25,:) = mapminmax(nninputs(25,:),mi, mx);
nntargets = features(26,:)
%nntargets = mapminmax(nntargets,mi,mx)





