%%%%%%%%%%%%
%
%   Thesis Work
%
%%%%%%%%%%

%% Add path to BCILAB
addpath(genpath('/home/amart/BCILAB-master')); % For engrX Server
% addpath(genpath('/opt/matlab/opt/BCILAB-1.1')); % For Quantum Server

% Bring in AAR toolbox
%eeglabRoot = fileparts(which('eeglab'));
%url = 'https://github.com/germangh/eeglab_plugin_aar/archive/master.zip';
%unzip(url, [eeglabRoot filesep 'plugins']);
%addpath(genpath(eeglabRoot));
%eeglab redraw;

%% Start bcilab and eeglab
bcilab;
eeglab;

%cd /home/amart/
%addpath(genpath('/home/amart/BCILAB-master'));
cd ~
% Clear workspace
clc;