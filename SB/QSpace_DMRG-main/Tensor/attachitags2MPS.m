function [MPSwithitags] = attachitags2MPS(MPS, sitename, bondname)
MPSwithitags = MPS;
schar = convertStringsToChars(sitename);
bondchar = convertStringsToChars(bondname);

MPSwithitags(1).info.itags = {'startvac', strcat(bondchar,char(string(1)), '*'), strcat(schar,char(string(1)))};
MPSwithitags(end).info.itags = {strcat(bondchar,char(string(numel(MPS)-1))), 'end', strcat(schar,char(string(numel(MPS))))};

if (numel(MPS)>2)
    for itN = 2:numel(MPS)-1
    MPSwithitags(itN).info.itags={strcat(bondchar,char(string(itN-1))), strcat(bondchar,char(string(itN)), '*'), strcat(schar,char(string(numel(MPS))))};
    end
end

end