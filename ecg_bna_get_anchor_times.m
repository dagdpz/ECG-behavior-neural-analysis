function blockstart=ecg_bna_get_anchor_times(monkey,session,blocks)

for B=blocks
    BLOCK_PATH=['Y:\Data\TDTtanks\' monkey '_phys\' session '\Block-' num2str(B)];
    OPTIONS={'TYPE', {'EPOCS'}};
    data = TDTbin2mat_working(BLOCK_PATH, OPTIONS{:});
    blockstart(B)=seconds(datetime(data.info.utcStartTime)-datetime('00:00:00'));
end

blockstart= blockstart- blockstart(1);

end