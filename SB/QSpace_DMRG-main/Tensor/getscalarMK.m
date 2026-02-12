function returnscalar = getscalarMK(scalartensor)

% < Description >
%
% returnscalar = getscalarMK(scalartensor)
%
% For the scalar tensor, which has only one reduced matrix element, get the
% scalar value. It is used to prevent returning empty value. If the tensor
% is empty, it returns 0.
% **NB!**
% It is independent of the CGCs!!
%
% < Input >
% scalartensor : [QSpace tensor] It has only one reduced matrix element.
%
% < Option >
%
% < Output >
% returnscalar : [numeric] The scalar value of the tensor. If the tensor is
% empty, it is 0.
%
% Written by Minsoo Kim (Dec. 21, 2023)


if isempty(scalartensor.data)
    returnscalar=0;
    %sprintf('empty tensor!')
    return
end
returnscalar=scalartensor.data{1}(1);

end