%%%%%%%%%%%%
%
%   TStart
%       File for starting EEGLAB and BCILAB if needed
%       Includes changing/adding proper directories
%
%%%%%%%%%%

clear;
clc;

%% Start BCILAB
% addpath(genpath('/home/amart/BCILAB-master')); % For engrX Server
% addpath(genpath('/opt/matlab/opt/BCILAB-1.1')); % For Quantum Server
%cd /opt/matlab/opt/BCILAB-1.1

%bcilab;


%% Start EEGLAB
% cd /opt/matlab/opt/eeglab13_4_3_b % 
%rmpath /opt/matlab/opt/eeglab13_4_4b/ % if needed to remove path of non-working eeglab
%cd ~/eeglab13_4_4b/ % CD or AddPath for eeglab

eeglab; % Start EEGLAB


%% % Bring in AAR toolbox - once brought in and redrawn no need to repeat
% eeglabRoot = fileparts(which('eeglab'));
% url = 'https://github.com/germangh/eeglab_plugin_aar/archive/master.zip';
% unzip(url, [eeglabRoot filesep 'plugins']);
% addpath(genpath(eeglabRoot));
% eeglab redraw;


%% Clear workspace
clc;


%% Change directory into location of scripts to run
cd ~/thesis % Change into home thesis folder


%% Run Proper Script
TMainWavelets
%TNNLoop
