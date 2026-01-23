function site_out = ecg_bna_resampled_spectra(site,cfg)

% lfp_tfa_process_LFP - function to read in the trial-wise LFP data for all
% site recorded in a session, compute the LFP time frequency spectrogram,
switch cfg.to_trigger
    case 'MUA'
        dospectralanalysis=0;
        zscoreraw=1;
    case 'LFP'
        dospectralanalysis=1;
        zscoreraw=0;
end
rawname=lower(cfg.to_trigger); % mua or lfp
origname=cfg.to_trigger;
samples_fn=[origname '_samples'];
ts_orig=cfg.ts_original;


fprintf('=============================================================\n');
fprintf('Processing site, %s\n', site.site_ID);

%filtering for evoked
if cfg.lfp.filter{1} && strcmp(cfg.to_trigger,'LFP')
    rawstream = eegfilt(site.(origname),1/ts_orig,cfg.lfp.filter{2},cfg.lfp.filter{3});
else
    rawstream = site.(origname);
end

N_cycles=cfg.lfp.n_cycles;
frequencies = cfg.lfp.foi;
frequency_bands=cfg.lfp.frequency_bands;
morlet_borders=1/min(frequencies)*N_cycles/2;
s = N_cycles./(2*pi*frequencies);
%s = repmat(0.05,size(frequencies));

ts=round(cfg.lfp.timestep/ts_orig);

% fT parameters (use next-pow-of-2)
time = -morlet_borders:1*ts_orig:morlet_borders;
n_wavelet     = length(time);


TFS.resampling_factor=1/ts;
TFS.sr=1/ts_orig/ts;

if dospectralanalysis
    TFS.freq             = frequencies;
    sizepreallocator=[size(frequency_bands,1),floor(numel(rawstream)/ts)];
    TFS.phabp=NaN(sizepreallocator);
    TFS.powbp=NaN(sizepreallocator);
    sizepreallocator=[numel(frequencies),floor(numel(rawstream)/ts)];
    TFS.pha=NaN(sizepreallocator);
    TFS.pow=NaN(sizepreallocator);
end

TFS.evoked=NaN([1,floor(numel(rawstream)/ts)]);

trials_block=[site.block];
valid_blocks=unique(trials_block);
samples_past=0;
samples_past_resampled=0;


for b=1:numel(valid_blocks)
    B=valid_blocks(b);
    b_samples=trials_block==B;
    
    bs=samples_past_resampled+1;
    be=floor(sum(site.(samples_fn)(b_samples))/ts)+samples_past_resampled;
    bs_original=samples_past+1;
    be_original=sum(site.(samples_fn)(b_samples))+samples_past;
    samples_past=be_original;
    samples_past_resampled=be;
    
    TFS.n_samples_per_block(:,b)=[B,be-bs+1];
    if dospectralanalysis
        concat_raw    = double(site.(origname)(bs_original:be_original))*1000000; % scale here is really bad
        n_data        = size(concat_raw,2);
        n_convolution = n_wavelet+n_data;
        n_conv_pow2   = pow2(nextpow2(n_convolution));
        half_wavelet_len = ceil(length(time)/2);
        dataft = fft(concat_raw,n_conv_pow2);
        
        
        for f=1:length(frequencies)
            % create wavelet
            wavelet = (ts_orig/(s(f)*sqrt(2*pi))) * exp(2*1i*pi*frequencies(f).*time) .* exp(-time.^2./(2*(s(f)^2)));
            % convolution
            datconv = ifft(fft(wavelet,n_conv_pow2).*dataft);            
            dat = datconv(1:n_convolution);
            dat = dat(half_wavelet_len+1:end-half_wavelet_len); % +1 is for one overlap point at the start
            % resample here:
            dat = nanmean(reshape(dat(1:end-mod(size(dat,2), ts)),ts,[]),1);            
            % extract pha values of reshaped data:
            TFS.pha(f,bs:be) = dat./abs(dat);
            % extracted Power of each trial
            TFS.pow(f,bs:be) = abs(dat).^2;
        end
        
        for f=1:size(frequency_bands,1)
            %% think about how to use good filters without causing errors for short periods
            %  fltered_data = eegfilt(concat_LFP, round(1/ts),frequency_bands(f,1), []);
            %  fltered_data = eegfilt(fltered_data, round(1/ts), [], frequency_bands(f,2));
            %         dat = eegfilt(concat_raw, 1/ts_orig, frequency_bands(f,1), frequency_bands(f,2));
            [u, v]=butter(3, 2*frequency_bands(f,:)*ts_orig); % band-pass filter
            dat = filtfilt(u,v,concat_raw);
            H=hilbert(dat);
            H = mean(reshape(H(1:end-mod(size(H,2), ts)),ts,[]),1);
            absH=abs(H);
            TFS.phabp(f,bs:be) = H./absH;
            TFS.powbp(f,bs:be) = absH.^2;
        end
    end
    
    %% redo for evoked, so it can be filtered independent of the rest
    concat_raw=double(rawstream(bs_original:be_original))*1000000; % scale here is really bad
    dat = nanmean(reshape(concat_raw(1:end-mod(size(concat_raw,2), ts)),ts,[]),1);
    
    if zscoreraw
        dat = zscore(dat);
    end
    
    TFS.evoked(1,bs:be)= dat;
end

if dospectralanalysis
    TFS.pha=TFS.pha(:,1:be);
    TFS.pow=TFS.pow(:,1:be);
    TFS.phabp=TFS.phabp(:,1:be);
    TFS.powbp=TFS.powbp(:,1:be);
end

TFS.evoked=TFS.evoked(:,1:be);

% store output
site_out = rmfield(site,{'LFP','MUA'});
site_out.session = site.site_ID(1:12);
site_out.recorded_hemisphere = upper(site.target(end));
site_out.(rawname) = TFS;

% Noise rejection - is this even still feasable?
% tic
% site_out = ecg_bna_reject_noisy_lfp_trials( site_out, lfp_tfa_cfg.noise );
% toc
end

