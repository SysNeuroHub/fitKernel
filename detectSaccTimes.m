function msaccTimes = detectSaccTimes(x,y,dt, lambda, minDur)
%msaccTimes = detectSaccTimes(x,y,dt, lambda, minDur)
%detects microsaccade start/end according to Engbert and Kliegel vis res
%
% msaccTimes [#saccade x 2(start/end)]

% x = eyeData_rmotl_cat.x;
% y = eyeData_rmotl_cat.y;
% dt = eyeData_rmotl_cat.dt;
vx = (x(3:end)+x(2:end-1)-x(1:end-2)-x(2:end-1))/6/dt;
vy = (y(3:end)+y(2:end-1)-y(1:end-2)-y(2:end-1))/6/dt;
sigmax = median(vx.^2) - (median(vx)).^2;
sigmay = median(vy.^2) - (median(vy)).^2;
exx = find(abs(vx) > lambda*sigmax)+1;
exy = find(abs(vy) > lambda*sigmay)+1;
exIdx=union(exx,exy);
msaccTraceInit = zeros(length(x),1);
msaccTraceInit(exIdx)=1;
msaccTimesInit = trace2Event(msaccTraceInit);
theseSaccs = find(diff(msaccTimesInit,1,2)>minDur/dt);
msaccTimes = msaccTimesInit(theseSaccs,:);
