% This function is used to calculate the response std
% Normalized to the number of neurons and time
function std=calResNormStd(R)
% R: the response matrix, neuron number x response time x ntrial
[~,~,nTrial]=size(R);
meanSR=squeeze(mean(R,3));
Dis=zeros(size(meanSR));
for sstt=1:nTrial
    Dis=Dis+(squeeze(R(:,:,sstt))-meanSR).^2;
end
std=sum(Dis(:))./sum(meanSR(:).^2)/nTrial;
