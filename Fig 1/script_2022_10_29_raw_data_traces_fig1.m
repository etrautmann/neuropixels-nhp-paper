
%%  2) Load data

% Note: currently figure contains data from 2022-07-03 bank1 stitching
% datasets. Should try replotting with 2021-05-25 dataset for Cousteau
% which was excellent too. 

filename = '/Volumes/emt_ssd_4/data/NeuropixelsNHP/2022-07-03_stitching_dataset/pacman-task_c_220703_neu_stitching_bank1_g0_imec0/pacman-task_c_220703_neu_stitching_bank1_g0_t0.imec0.ap.bin'
% filename = '/Volumes/emt_ssd_6/data/pacman-task/cousteau/2021-05-25/neuropixels/pacman-task_c_210525_neu_g0/pacman-task_c_210525_neu_g0_imec0/pacman-task_c_210525_neu_g0_t0.imec0.ap.bin'

% this dataset is okay, not great: 
% filename = '/Volumes/emt_ssd_4/data/pacman-task/igor/raw/2022-05-26/neuropixels/pacman-task_i_220526_neu_g0/pacman-task_i_220526_neu_g0_imec1/pacman-task_i_220526_neu_g0_t0.imec1.ap.bin'

fid = fopen(filename)

[~,ks_run,~] = fileparts(filename);
[~,ks_run,~] = fileparts(ks_run); % strip extensions
[~,ks_run,~] = fileparts(ks_run) % strip extensions

% fid = fopen('/Volumes/churchland-locker/Jumanji/pacman-task/cousteau/raw/2021-05-25/neuropixels/pacman-task_c_210525_neu_g0/pacman-task_c_210525_neu_g0_imec0/pacman-task_c_210525_neu_g0_t0.imec0.ap.bin')

%%

MICROVOLT_PER_SAMPLE = 2.34  % standard for gain = 500 setting

msOffset = 0;
% secOffset = 1.7;
% secs = .05;

secOffset = 1000;
secsToLoad = 10;

chanStart = 1
chanStop = 384
% chanStart = 147
% chanStop = 149


% chanStart = 1  
% chanStop = 276
% scale = 10/MICROVOLT_PER_SAMPLE;
scale = 1
nChan = 385;
fSamp = 30000;
nTap = 1;
lineWidth = 1

bytesOffset = nChan * fSamp * secOffset;

bytesOffset = 0
fseek(fid, bytesOffset,'bof');
nSamp = secsToLoad*fSamp;

d = double(fread(fid, [nChan nSamp],'*int16'));
d = d';
d(:,end) = []; % strip off sync channel 

% crossChanMean = mean(d,2);
% d = d-int16(repmat(crossChanMean, 1, size(d,2)));
% 
% tMean = mean(d,1);
% d = d-int16(repmat(tMean, size(d,1), 1));


tic


% 1) common mode subtraction
% 1.1) center each channel
d = bsxfun(@minus,d, mean(d,1));

% 1.2) subtract off cross channel median
d = bsxfun(@minus,d, median(d,2));
toc

% 2) bandpass filter
sampRate = 30000;
lowcut = 250;
highcut = 5000;
[b,a] = butter(2, [lowcut highcut]/(sampRate/2));

for i=1:size(d,2)
    d(:,i) = filter(b,a,d(:,i));
end

%%

chanStart = 1
chanStop = 384

% % filter across channels with moving average filter
% for i=1:size(d2,1)
%     d2(i,:) = movmean(d2(i,:),7);
% end
% size(d)
%

traceOffset = 30;

msToPlot =  18;
startTime = 0;
stopSamp = (startTime  + msToPlot)*(fSamp/1000);

startSamp = ceil(startTime*(fSamp/1000)+eps);

t = (startSamp:stopSamp)/fSamp;

% d2 = d(startSamp:stopSamp,:)*MICROVOLT_PER_SAMPLE;
d2 = d(startSamp:stopSamp,:)*MICROVOLT_PER_SAMPLE;
plotOffset = (traceOffset*repmat(1:size(d2,2), size(d2,1), 1));
% plotOffset = int16(scale*repmat(1:size(d,2), size(d,1), 1));
% plotOffset = 0
% dOffset = d2;%+plotOffset;
dOffset = d2+plotOffset;


fig1 = figure(1); clf;
% subplot(121)
p = plot(t', dOffset(:,chanStart+1:2:chanStop),'color',[1 0 .5 .85] ,'linewidth', .2);
hold on
p = plot(t', dOffset(:,chanStart:2:chanStop),'color', [0 .2 .8 .85],'linewidth', .2);
axis tight
axis on

hold on
% subplot(122)

% set(p(badChans),'color',[1 0 0])
% set(p(markChans),'color',[1 0 0])

axis tight
axis off

set(fig1, 'WindowStyle','normal')
% set(fig1,'Position',[0 0 400 600])
% set(fig1,'Position',[0 0 3500 600])
set(fig1,'color',[0 0 0])
set(fig1,'color',[1 1 1])

hold on

yScaleBarLength = 250;
xScaleLength = 5;



plot(min(xlim)*ones(1,2),[0 yScaleBarLength],'k','linewidth',5)
hold on
plot(min(xlim)*ones(1,2) + [0 xScaleLength/1000],[0 0],'k','linewidth',5)


%%  save output

figPath = '/Users/erictrautmann/Dropbox/columbia/manuscripts - columbia/Neuropixels NHP paper/figures/Fig. 1 - Probe/assets/raw_waveforms/';

fname = sprintf('%s_t_%d_dur_%d_off_%d.pdf',ks_run, startTime, msToPlot, traceOffset)
fnamefull = fullfile(figPath, fname);

print(fig1,fnamefull,'-dpdf','-bestfit','-painters')


