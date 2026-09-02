monkeys={'Magnus','Bacchus'};
signals={'lfp','mua'};
task={'Task','Rest'};
for mon=1:numel(monkeys)
    monkey=monkeys{mon};
    [numericData, textData, rawData] = xlsread('Y:\Projects\Pulv_bodysignal\Figures\significance_table.xls',[monkey '_nsigbins_per_site']);
    nuclei={'VPL','dPul','MD'};
    for n=1:numel(nuclei);
        nuc=nuclei{n};
        targets=rawData(3:end,2);
        rows=~cellfun(@isempty,strfind(targets,nuc));
        
        for s=1:numel(signals)
            for t=1:numel(task)
                sig=signals{s};
                tas=task{t};
                
                label_pre=[sig '_' tas '_pre'];
                label_post=[sig '_' tas '_post'];
                
                col_pre=~cellfun(@isempty,strfind(rawData(2,:),label_pre));
                col_post=~cellfun(@isempty,strfind(rawData(2,:),label_post));
                n_before=numericData(rows,find(col_pre)-2);
                n_after=numericData(rows,find(col_post)-2);
                
                p=signrank(n_after-n_before); %paired
                
                disp([monkey '-' tas '-' nuc '-' sig '-p=' num2str(p)])
                
            end
        end
        
        
        
    end
end


