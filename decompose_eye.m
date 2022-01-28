function eyeData_tr = decompose_eye(eyeData_cat, t_tr)

%eyeData_rmblk_tr = eyeData;
%t={eyeData.t};
%psth_tr = cell(nTrials, 1);
nTrials = length(t_tr);
headidx = 1;
for itr = 1:nTrials
    nFrames = length(t_tr{itr});
    thisx = eyeData_cat.x(headidx:headidx+nFrames-1);
    thisy = eyeData_cat.y(headidx:headidx+nFrames-1);
    thispwdth = eyeData_cat.pwdth(headidx:headidx+nFrames-1);
    thisphght = eyeData_cat.phght(headidx:headidx+nFrames-1);
    
    %eyeData_cat_tr(itr).t = t{itr};
    %eyeData_tr(itr).x = thisx;
    %eyeData_tr(itr).y = thisy;
    %eyeData_tr(itr).pwdth = thispwdth;
    %eyeData_tr(itr).phght = thisphght;
    eyeData_c = marmodata.eye(t_tr{itr}, thisx, thisy, thispwdth, thisphght);
    eyeData_tr(itr) = eyeData_c;
    
    headidx = headidx+nFrames;
end
