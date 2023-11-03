function [prefDir_q, prefBin] = quantizeDir(prefDir, cardinalDir)
%[prefDir_q, prefBin] = quantizeDir(prefDir, cardinalDir)

for ii=1:size(prefDir,1)
    for jj=1:size(prefDir,2)
        %[~,prefBin(ii,jj)] = min(abs(prefDir(ii,jj) - cardinalDir));
            [~,prefBin(ii,jj)] = min(abs(circ_dist(pi/180*prefDir(ii,jj), pi/180*cardinalDir)));%getEyeDirMtx
            prefDir_q(ii,jj) = cardinalDir(prefBin(ii,jj));
    end
end
