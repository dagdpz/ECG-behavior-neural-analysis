function site_lfp = ecg_bna_process_LFP_and_MUA(sites,cfg)

% lfp_tfa_process_LFP - function to read in the trial-wise LFP data for all
% sites recorded in a session, compute the LFP time frequency spectrogram,

ts_LFP=cfg.LFP_ts_LFP;
ts_MUA=cfg.MUA_ts_LFP;

% struct to save data for a site
site_lfp = rmfield(sites,{'LFP','MUA'});
fprintf('=============================================================\n');
fprintf('Processing site, %s\n', sites.site_ID);
site_lfp.session = sites.site_ID(1:12);
site_lfp.recorded_hemisphere = upper(sites.target(end));
% %filtering the LFP file

if cfg.lfp.filter{1}
    rawLFP = eegfilt(sites.LFP,1/ts_LFP,cfg.lfp.filter{2},cfg.lfp.filter{3});
else
    rawLFP = sites.LFP;
end

N_cycles=cfg.lfp.n_cycles;
frequencies = cfg.lfp.foi;
frequency_bands=cfg.lfp.frequency_bands;
morlet_borders=1/min(frequencies)*N_cycles/2;
s = N_cycles./(2*pi*frequencies);
%s = repmat(0.05,size(frequencies));

ts=round(cfg.lfp.timestep/ts_LFP);

% fT parameters (use next-pow-of-2)
time = -morlet_borders:1*ts_LFP:morlet_borders;
n_wavelet     = length(time);


site_lfp.tfs.freq             = frequencies;
site_lfp.tfs.resampling_factor=1/ts;
site_lfp.tfs.sr=1/ts_LFP/ts;
sizepreallocator=[size(frequency_bands,1),floor(numel(sites.LFP)/ts)];
site_lfp.tfs.phabp=NaN(sizepreallocator);
site_lfp.tfs.powbp=NaN(sizepreallocator);
sizepreallocator=[numel(frequencies),floor(numel(sites.LFP)/ts)];
site_lfp.tfs.pha=NaN(sizepreallocator);
site_lfp.tfs.pow=NaN(sizepreallocator);

site_lfp.tfs.lfp=NaN([1,floor(numel(sites.LFP)/ts)]);

trials_block=[sites.block];
blocks_with_LFP=unique(trials_block);
samples_past=0;
samples_past_resampled=0;


for b=1:numel(blocks_with_LFP) 
    B=blocks_with_LFP(b);
    b_samples=trials_block==B;
    
    bs=samples_past_resampled+1;
    be=floor(sum(sites.LFP_samples(b_samples))/ts)+samples_past_resampled;  
    bs_original=samples_past+1;
    be_original=sum(sites.LFP_samples(b_samples))+samples_past;  
    samples_past=be_original;
    samples_past_resampled=be;    
    
    concat_raw = double(sites.LFP(bs_original:be_original))*1000000; % scale here is really bad 
    
    
    n_data        = size(concat_raw,2);
    n_convolution = n_wavelet+n_data;
    n_conv_pow2   = pow2(nextpow2(n_convolution));
    half_wavelet_len = ceil(length(time)/2);
    
    dataft = fft(concat_raw,n_conv_pow2);
    
    site_lfp.tfs.n_samples_per_block(:,b)=[B,be-bs+1];
    dat = nanmean(reshape(concat_raw(1:end-mod(size(concat_raw,2), ts)),ts,[]),1);
        
    for f=1:length(frequencies)
        % create wavelet        
        wavelet = (ts_LFP/(s(f)*sqrt(2*pi))) * exp(2*1i*pi*frequencies(f).*time) .* exp(-time.^2./(2*(s(f)^2)));   
        % convolution
        datconv = ifft(fft(wavelet,n_conv_pow2).*dataft);
        
        dat = datconv(1:n_convolution);
        dat = dat(half_wavelet_len+1:end-half_wavelet_len); % +1 is for one overlap point at the start
        % resample here:
        dat = nanmean(reshape(dat(1:end-mod(size(dat,2), ts)),ts,[]),1);
        
        % extract pha values of reshaped data:
        site_lfp.tfs.pha(f,bs:be)= dat./abs(dat);
        % extracted Power of each trial
        site_lfp.tfs.pow(f,bs:be) = abs(dat).^2;       
    end
    
    for f=1:size(frequency_bands,1)
        %% think about how to use good filters without causing errors for short periods
        %  fltered_data = eegfilt(concat_LFP, round(1/ts),frequency_bands(f,1), []);
        %  fltered_data = eegfilt(fltered_data, round(1/ts), [], frequency_bands(f,2));
%         dat = eegfilt(concat_raw, 1/ts_LFP, frequency_bands(f,1), frequency_bands(f,2));
        [u, v]=butter(3, 2*frequency_bands(f,:)*ts_LFP); % band-pass filter
        dat = filtfilt(u,v,concat_raw);
        H=hilbert(dat);
        H = mean(reshape(H(1:end-mod(size(H,2), ts)),ts,[]),1);
        absH=abs(H);
        site_lfp.tfs.phabp(f,bs:be) = H./absH;
        site_lfp.tfs.powbp(f,bs:be) = absH.^2;               
    end
    
    %% redo for evoked, so it can be filtered independent of the rest
    concat_raw=double(rawLFP(bs_original:be_original))*1000000; % scale here is really bad 
    dat = nanmean(reshape(concat_raw(1:end-mod(size(concat_raw,2), ts)),ts,[]),1);    
    site_lfp.tfs.lfp(1,bs:be)= dat;
end
site_lfp.tfs.pha=site_lfp.tfs.pha(:,1:be);
site_lfp.tfs.pow=site_lfp.tfs.pow(:,1:be);
site_lfp.tfs.phabp=site_lfp.tfs.phabp(:,1:be);
site_lfp.tfs.powbp=site_lfp.tfs.powbp(:,1:be);

site_lfp.tfs.lfp=site_lfp.tfs.lfp(:,1:be);

% Noise rejection - is this even still feasable?
% tic
% site_lfp = ecg_bna_reject_noisy_lfp_trials( site_lfp, lfp_tfa_cfg.noise );
% toc
end

