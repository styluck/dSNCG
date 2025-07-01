% wsRun1.Me(2,:) = wsRun3.Me(1,:);
% wsRun1.Dev(2,:) = wsRun3.Dev(1,:);
% wsRun1.ObjMe(2,:) = wsRun3.ObjMe(1,:);
% wsRun1.ObjDev(2,:) = wsRun3.ObjDev(1,:);
% 
% wsRun2.Me(2,:) = wsRun4.Me(1,:);
% wsRun2.Dev(2,:) = wsRun4.Dev(1,:);
% wsRun2.ObjMe(2,:) = wsRun4.ObjMe(1,:);
% wsRun2.ObjDev(2,:) = wsRun4.ObjDev(1,:);
slct = [4,   7, 8,9];
algs = algs([1 3 4 5]); % algs = algs(slct);
Dev1 = Dev1(slct,:);
Dev2 = Dev2(slct,:);
Me1 = Me1(slct,:);
Me2 = Me2(slct,:);

ObjDev1 = ObjDev1(slct,:);
ObjDev2 = ObjDev2(slct,:);
ObjMe1 = ObjMe1(slct,:);
ObjMe2 = ObjMe2(slct,:);

