function  [pdiam_prctile, pdiam] = getPupilDiameter(eyeData_rmblk_cat)
%         [pdiam_prctile, pdiam] = getPupilDiameter(eyeData_rmblk_cat);

pdiam = sqrt(eyeData_rmblk_cat.parea); %convert area to diameter

[~, I] = sort(pdiam, 'ascend'); %B = pdiam(I)
%I_prctile = 1/length(pdiam)*(1:length(pdiam));

[~,pdiam_order] = sort(I);
pdiam_prctile = 100/length(pdiam)*pdiam_order;

% ax(1)=subplot(211);
% plot(pdiam);
% ax(2)=subplot(212);
% plot(pdiam_prctile);
% linkaxes(ax(:),'x');