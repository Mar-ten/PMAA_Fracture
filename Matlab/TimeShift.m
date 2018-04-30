function [ shifted ] = TimeShift( data,shift )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
shifted=zeros(size(data));
shift=-shift;
if shift>=0
    
    for i=1:length(shifted)-shift
        shifted(i)=data(i+shift);
    end
else
    for i=-shift+1:length(shifted)
        shifted(i)=data(i+shift);
    end
end

% figure
% plot(1:length(data),data,1:length(data),shifted)
% legend('data','shifted')
end

