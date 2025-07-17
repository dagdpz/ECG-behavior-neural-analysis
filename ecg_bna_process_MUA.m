function site_mua = ecg_bna_process_MUA(sites,cfg,ts_original)
dospectralanalysis=0;

% struct to save data for a site
site_mua = rmfield(sites,{'LFP','MUA'}); %% remove MUA too

fprintf('=============================================================\n');
fprintf('Processing site, %s\n', sites.site_ID);
site_mua.session = sites.site_ID(1:12);
site_mua.recorded_hemisphere = upper(sites.target(end));
% %filtering the LFP file
% sites.LFP = eegfilt(sites.LFP,1/ts_original,4,[]);

N_cycles=cfg.lfp.n_cycles;
frequencies = cfg.lfp.foi;
frequency_bands=cfg.lfp.frequency_bands;
morlet_borders=1/min(frequencies)*N_cycles/2;
s = N_cycles./(2*pi*frequencies);

ts=round(cfg.lfp.timestep/ts_original);

% fT parameters (use next-pow-of-2)
time = -morlet_borders:1*ts_original:morlet_borders;
n_wavelet     = length(time);


site_mua.tfs.resampling_factor=1/ts;
site_mua.tfs.sr=1/ts_original/ts;
if dospectralanalysis
    sizepreallocator=[size(frequency_bands,1),floor(numel(sites.LFP)/ts)];
    site_mua.tfs.phabp=NaN(sizepreallocator);
    site_mua.tfs.powbp=NaN(sizepreallocator);
    sizepreallocator=[numel(frequencies),floor(numel(sites.LFP)/ts)];
    site_mua.tfs.pha=NaN(sizepreallocator);
    site_mua.tfs.pow=NaN(sizepreallocator);
    %site_mua.tfs.lfp=NaN([1,floor(numel(sites.LFP)/ts)]);
end

site_mua.tfs.mua=NaN([1,floor(numel(sites.MUA)/ts)]);

trials_block=[sites.block];
blocks_with_LFP=unique(trials_block);
samples_past=0;
samples_past_resampled=0;


for b=1:numel(blocks_with_LFP) 
    B=blocks_with_LFP(b);
    b_samples=trials_block==B;
    bs=samples_past_resampled+1;
    be=floor(sum(sites.MUA_samples(b_samples))/ts)+samples_past_resampled;  
    bs_original=samples_past+1;
    be_original=sum(sites.MUA_samples(b_samples))+samples_past;  
    samples_past=be_original;
    samples_past_resampled=be;    
    
    concat_raw = real(double(sites.MUA(bs_original:be_original))*1000000); % scale here is really bad 
    %% for some freaking reason, mua can be a complex......
    
    site_mua.tfs.n_samples_per_block(:,b)=[B,be-bs+1];
    dat = nanmean(reshape(concat_raw(1:end-mod(size(concat_raw,2), ts)),ts,[]),1);
    site_mua.tfs.mua(1,bs:be)= dat;
    
    if dospectralanalysis
        
        n_data        = size(concat_raw,2);
        n_convolution = n_wavelet+n_data;
        n_conv_pow2   = pow2(nextpow2(n_convolution));
        half_wavelet_len = ceil(length(time)/2);
        
        % get fT of data
        dataft = fft(concat_raw,n_conv_pow2);
        site_mua.tfs.freq             = frequencies;
        
        for f=1:length(frequencies)
            % create wavelet
            wavelet = (ts_original/(s(f)*sqrt(2*pi))) * exp(2*1i*pi*frequencies(f).*time) .* exp(-time.^2./(2*(s(f)^2)));
            % convolution
            datconv = ifft(fft(wavelet,n_conv_pow2).*dataft);
            
            dat = datconv(1:n_convolution);
            dat = dat(half_wavelet_len+1:end-half_wavelet_len); % +1 is for one overlap point at the start
            % resample here:
            dat = nanmean(reshape(dat(1:end-mod(size(dat,2), ts)),ts,[]),1);
            
            % extract pha values of reshaped data:
            site_mua.tfs.pha(f,bs:be)= dat./abs(dat);
            % extracted Power of each trial
            site_mua.tfs.pow(f,bs:be) = abs(dat).^2;
        end
        
        for f=1:size(frequency_bands,1)
            %% think about how to use good filters without causing errors for short periods
            %  fltered_data = eegfilt(concat_LFP, round(1/ts),frequency_bands(f,1), []);
            %  fltered_data = eegfilt(fltered_data, round(1/ts), [], frequency_bands(f,2));
            %         dat = eegfilt(concat_raw, 1/ts_original, frequency_bands(f,1), frequency_bands(f,2));
            [u, v]=butter(3, 2*frequency_bands(f,:)*ts_original); % band-pass filter
            dat = filtfilt(u,v,concat_raw);
            H=hilbert(dat);
            H = mean(reshape(H(1:end-mod(size(H,2), ts)),ts,[]),1);
            absH=abs(H);
            site_mua.tfs.phabp(f,bs:be) = H./absH;
            site_mua.tfs.powbp(f,bs:be) = absH.^2;
        end
    end
end
if dospectralanalysis
site_mua.tfs.pha=site_mua.tfs.pha(:,1:be);
site_mua.tfs.pow=site_mua.tfs.pow(:,1:be);
site_mua.tfs.phabp=site_mua.tfs.phabp(:,1:be);
site_mua.tfs.powbp=site_mua.tfs.powbp(:,1:be);
end

site_mua.tfs.mua=site_mua.tfs.mua(:,1:be);
end

