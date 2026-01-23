function ecg_bna_draw_outlines(matrix, color)
    [rows, cols] = size(matrix);
    if ~exist('color','var') || isempty(color)
        color = [0 0 0];
    end
    hold on;
    x = [];
    y = [];
    for i = 1:rows
        for j = 1:cols
            if matrix(i,j) == 1
                if i == 1 || matrix(i-1,j) == 0
                    x = [x, j-0.5, j+0.5, NaN];
                    y = [y, i-0.5, i-0.5, NaN];
                end
                if i == rows || matrix(i+1,j) == 0
                    x = [x, j-0.5, j+0.5, NaN];
                    y = [y, i+0.5, i+0.5, NaN];
                end
            end
        end
    end
    for j = 1:cols
        for i = 1:rows
            if matrix(i,j) == 1
                if j == 1 || matrix(i,j-1) == 0
                    x = [x, j-0.5, j-0.5, NaN];
                    y = [y, i-0.5, i+0.5, NaN];
                end
                if j == cols || matrix(i,j+1) == 0
                    x = [x, j+0.5, j+0.5, NaN];
                    y = [y, i-0.5, i+0.5, NaN];
                end
            end
        end
    end
    plot(x, y, 'Color', color, 'LineWidth', 1);

% 
% 
% [FX,FY]=gradient(ROI_Matrix);
% patch(ROI_Matrix)
% FX = FX~=0;FY = FY~=0;
% hold on;
% %change in x, horizontal line
% for i = 1:size(FX,1)
%     for j = 1:size(FX,2)
%         if FX(i,j)
%             plot([i,i],[j,j+1],'k-');
%         end
%     end
% end
% %change in y, verticle line
% for i = 1:size(FY,1)
%     for j = 1:size(FY,2)
%         if FY(i,j)
%             plot([i,i+1],[j,j],'k-');
%         end
%     end
% end
% hold off;
%
end