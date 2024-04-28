function [] = make_aud_stims(fname,prefix)


AUD = audioread(fname);
INF = audioinfo(fname);

m = mean(abs(AUD));
AUD(:,1) = AUD(:,1) * (m(2)./m(1));

intrvl = INF.SampleRate*.9;

if logical(rem(intrvl,1))
    error('Sample rate precludes dividing audio into 900 ms segments. Please resample the audio file.');
end

numSegs = round(INF.TotalSamples ./ intrvl) - 11;
toskip = floor((INF.TotalSamples - (numSegs .* intrvl)) / 2);

for i = 1:numSegs
    tAUD = AUD(toskip + intrvl*(i-1) + 1: toskip + intrvl*i,:);
    mt = mean(abs(tAUD));
    tAUD(:,1) = tAUD(:,1) .* (m(2)./mt(1));
    tAUD(:,2) = tAUD(:,2) .* (m(2)./mt(2));
    audiowrite([prefix '_' num2str(i,'%03d') '.wav'],tAUD,INF.SampleRate);
end
