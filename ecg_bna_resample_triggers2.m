function  Trigger_samples = ecg_bna_resample_triggers2(Triggers,offset_per_block,Stream_samples_per_block,sr)

Stream_blocks=Stream_samples_per_block(1,:);
Stream_samples=Stream_samples_per_block(2,:);
trigger_blocks=offset_per_block(1,:);
trigger_offset=offset_per_block(2,:);

events=fieldnames(Triggers)';

for f=events
    
    past_blocks=[];
    ts_real=Triggers.(f{:}).ts;
    ts_surrogate=Triggers.(f{:}).shuffled_ts;
    
    Real_ts=[];    
    Surrogate_ts=[];    
    for b=Stream_blocks
        BS=Stream_blocks==b;   %this block's samples in the LFP
        BT=trigger_blocks==b;
        if ~any(BT)
            display(['Block' num2str(b) 'not present in Triggers, check for potential bugs!']);
        end
                
        idx=Triggers.(f{:}).blocks==b;
        
        R=ts_real(idx)-trigger_offset(BT);
        S=round(R*sr);
        S(S<=0)=0;
        S(S>Stream_samples(BS))=0;
        S=S+sum(Stream_samples(ismember(Stream_blocks,past_blocks)));
        S(S==sum(Stream_samples(ismember(Stream_blocks,past_blocks))))=0;
        
        Real_ts=[Real_ts S];
        
        %% not shuffled, but surrogate!
        
        idx=Triggers.(f{:}).shuffled_blocks==b;
        R=ts_surrogate(:,idx)-trigger_offset(BT);
        S=round(R*sr);
        S(S<=0)=0;
        S(S>Stream_samples(BS))=0;
        S=S+sum(Stream_samples(ismember(Stream_blocks,past_blocks)));
        S(S==sum(Stream_samples(ismember(Stream_blocks,past_blocks))))=0;
        
        Surrogate_ts=[Surrogate_ts S];
                
        past_blocks=[past_blocks b];
    end
    Trigger_samples.([f{:} '_real'])=Real_ts;
    Trigger_samples.([f{:} '_shuffled'])=Surrogate_ts;
end