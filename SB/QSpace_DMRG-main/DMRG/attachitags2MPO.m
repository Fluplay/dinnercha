function MPOwithitags = attachitags2MPO(MPO, sitename, opbondname)

% < Description >
%
% MPOwithitags = attachitags2MPO(MPO, sitename, opbondname)
%
% attaching itags to MPO
%
% < Input >
% MPO : [QSpace array] array of MPO. The leg direction is bottom, top, left, right
% sitename : [string] The itags of physical legs
% opbondname : [string] The itags of virtual bond legs
%
% < Option >
%
% < Output >
% MPOwithitags : [QSpace array] array of MPO with itags. 'sitename site', 'sitename site*', 'opbondname bond', 'opbondname bond*'.
%
% Written by Minsoo Kim (Dec. 21, 2023)

MPOwithitags = MPO;
schar = convertStringsToChars(sitename);
Ochar = convertStringsToChars(opbondname);

MPOwithitags(1).info.itags = {strcat(schar,char(string(1))),strcat(schar,char(string(1)),'*'),...
    strcat(Ochar,'start'),strcat(Ochar,char(string(1)),'*')};
MPOwithitags(end).info.itags = {strcat(schar,char(string(numel(MPO)))),strcat(schar,char(string(numel(MPO))),'*'),...
    strcat(Ochar,char(string(numel(MPO)-1))),strcat(Ochar,'end*')};

if (numel(MPO)>2)
    for itN = 2:numel(MPO)-1
    MPOwithitags(itN).info.itags = {strcat(schar,char(string(itN))),strcat(schar,char(string(itN)),'*'),...
        strcat(Ochar,char(string(itN-1))),strcat(Ochar,char(string(itN)),'*')};
    end
end

end