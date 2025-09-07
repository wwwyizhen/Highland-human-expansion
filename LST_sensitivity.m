natural_mask20=1-hua20;
natural_mask20(countries==255)=0;
folder = 'mod11c3\tif_clip'; 
outFolder = fullfile(folder, 'tif');  
files = dir(fullfile(folder, '2020*.tif'));

for k=1:length(files)
    lst=importdata(fullfile(folder, files(k).name));
    lst(lst<0)=NaN;
    lst_stack(:,:,k)=lst;
end
%%
window_tif=zeros(m,n);
for i=6:m
    for j=6:n
        if (hill_change(i,j)>0||plain_change(i,j)>0)&& (countries(i,j)<255) 
           
            ws = window_size(natural_mask20, i, j);
            fprintf('→ %s\n', num2str(ws));
            window_tif(i,j)=ws;
            row_start = max(i - ws, 1);
            row_end   = min(i + ws, m);
            col_start = max(j - ws, 1);
            col_end   = min(j + ws, n);
            natural_mask_window=natural_mask20(row_start:row_end, col_start:col_end);
            for kk=1:length(files)
                lst_window=lst_stack(row_start:row_end, col_start:col_end,kk);
                lst_window_valid=lst_window((natural_mask_window>0.6)&~isnan(lst_window));
                if ~isempty(lst_window_valid)
                    dlst_month(i,j,kk)=lst_window(ws+1,ws+1)-mean(lst_window_valid(:));
                end
            end
        end
    end
end
%%
avg_dlst = nanmean(dlst_month, 3);

tifname=['LST\2020_dlst.tif']
geotiffwrite(tifname, avg_dlst, Ref, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);

%%
dlst2020=importdata('2020_dlst.tif');
dlst2000=importdata('2000_dlst.tif');
dlst=dlst2020-dlst2000;
%%
hill_p=zeros(m,n);
hill_co=zeros(m,n);
plain_r2=zeros(m,n);
valid_count=1;
for i=1:m
    for j=1:n
        if hill_change(i,j)>0
            row_start = max(i - 20, 1);
            row_end   = min(i + 20, m);
            col_start = max(j - 20, 1);
            col_end   = min(j + 20, n);
            hill_change_window=hill_change(row_start:row_end, col_start:col_end);
            dlst_window=dlst(row_start:row_end, col_start:col_end);
            pos=find(hill_change_window>0);
            if length(pos)>50
                disp(valid_count)
                pa = polyfit(hill_change_window(pos), dlst_window(pos), 1);
                [rho,pval]=corr(hill_change_window(pos), dlst_window(pos));
                hill_co(i,j)=pa(1);
               hill_p(i,j)=pval;

                valid_count=valid_count+1;

            end
        end
    end
end
hill_co(hill_p>0.1)=NaN;
% imagesc(hill_co)

  tifname1=['LST\hill_LST_sensitivity.tif']
    

