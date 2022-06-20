function [f] = assignFaultPcSimple(f, fault, refPc)
%
% See Zulqarnain et al., IJGGC (2018; equation 3)
%
% NOTE that during simulateScheduleAd, within (My)BlackOilPc, sG for all
% fault cells (used here) is passed in ascending order of global fault cell
% id. Hence, fperm and fporo must also be passed in ascending order of
% global fault cell id, otherwise the evaluation would be wrongly done.
%

fperm = fault.perm;  % maximum perm will usually be in the along (updip) direction
fporo = fault.poro;
f.pcOW{fault.satRegNum} = @(sw)   pcFault(sw, fperm, fporo, refPc);

function v = pcFault(sw, perm, poro, refPc)
    v = refPc.val(sw).*sqrt((refPc.perm .* poro)./(refPc.poro .* perm));
end

end