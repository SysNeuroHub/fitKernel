function [saveServer, rootFolder] = getReady()
%[saveServer, rootFolder] = getReady()

switch getenv('COMPUTERNAME')
    
    case 'MU00175834'
        addpath(genpath('C:/Users/dshi0006/git'))
        saveServer = 'Z:\Shared\Daisuke\cuesaccade_data';%'E:/tmp/cuesaccade_data';
        %saveFigFolder = [saveServer, '/20220722'];
        %mkdir(saveFigFolder);
        rootFolder = '//storage.erc.monash.edu.au/shares/R-MNHS-Physio/SysNeuroData/Monash Data/Joanita/';
        
    case 'MU00011697' %quetzal linux machine
        %addpath(genpath('/home/daisuke/Documents/git'));
        rmpath(genpath('/home/daisuke/Documents/git/VariabilityEarlyVisualCortex/matlab/'))
        rmpath(genpath('/home/daisuke/Documents/git/dsbox/chronux_2_11/'));
        %saveServer = '~/Documents/cuesaccade_data';
        %saveServer = '/mnt/syncitium/Daisuke/cuesaccade_data';
        saveServer = '/media/daisuke/cuesaccade_data';
        rootFolder = '/mnt/MBI/Monash Data/Joanita/';
        
    case 'MU00108396'
        addpath(genpath('/home/localadmin/Documents/MATLAB'));
        saveServer = '/mnt/syncitium/Daisuke/cuesaccade_data';
        rootFolder = '/mnt/physio/Monash Data/Joanita/2021/cuesaccade_data/';
end
