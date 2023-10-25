function [saveServer, rootFolder] = getReady()
%[saveServer, rootFolder] = getReady()

switch getenv('COMPUTERNAME')
    
    case 'MU00175834'
        addpath(genpath('C:/Users/dshi0006/git'))
        saveServer = 'Z:\Shared\Daisuke\cuesaccade_data';%'E:/tmp/cuesaccade_data';
        %saveFigFolder = [saveServer, '/20220722'];
        %mkdir(saveFigFolder);
        rootFolder = '//storage.erc.monash.edu.au/shares/R-MNHS-Physio/SysNeuroData/Monash Data/Joanita/';
        
    case 'MU00011697'
        saveServer = '~/Documents/cuesaccade_data';
        rootFolder = '/mnt/MBI/Monash Data/Joanita/';
        addpath(genpath('~/Documents/git'));
        
    case 'MU00108396'
        addpath(genpath('/home/localadmin/Documents/MATLAB'));
        saveFolder = '/mnt/syncitium/Daisuke/cuesaccade_data';
        rootFolder = '/mnt/physio/Monash Data/Joanita/2021/cuesaccade_data/';
end
