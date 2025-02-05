time_win_visual = [0.05,0.2];
time_win_baseline = [-0.1,0.05];
x_grid = unique(data_single_bank.xs);
y_grid = unique(data_single_bank.ys);
RF_mat = zeros(length(y_grid),length(x_grid),size(data_single_bank.spiketrain,3));
N_rep_mat = zeros(length(y_grid),length(x_grid),size(data_single_bank.spiketrain,3));
for p=1:length(data_single_bank.xs)
    spk_visual = mean(data_single_bank.spiketrain(p,data_single_bank.ts>=time_win_visual(1)&data_single_bank.ts<time_win_visual(2),:),2)*1000;
    spk_baseline = mean(data_single_bank.spiketrain(p,data_single_bank.ts>=time_win_baseline(1)&data_single_bank.ts<time_win_baseline(2),:),2)*1000;
    x = x_grid==data_single_bank.xs(p);
    y = y_grid==data_single_bank.ys(p);
    RF_mat(y,x,:) = RF_mat(y,x,:)+spk_visual-spk_baseline;
    N_rep_mat(y,x,:) = N_rep_mat(y,x,:)+1;
end
RF_mat = RF_mat./N_rep_mat;


selected_ids = [8, 21, 45, 57, 58, 67, 71, 91, 105, 124, 126, 130, 131, 132, 133, 139, 140, 151, 153, 154, 159, 162, 166, 168, 170, 171, 172, 179, 180, 183, 184, 185, 186, 187, 188, 191, 192, 194, 195, 198, 201, 202, 204, 205, 206, 208, 211, 214, 215, 216, 217, 218, 220, 221, 224, 226, 227, 228, 230, 231, 233, 235, 236, 238, 240, 243, 246, 249, 250, 252, 253, 256, 260, 261, 265, 267, 269, 270, 278, 280, 283, 284, 285, 287, 288, 289, 291, 292, 295, 297, 301, 303, 304, 305, 308, 309, 310, 311, 312, 313, 314, 315, 319, 325, 326, 327, 329, 330, 334, 335, 338, 342, 343, 345, 346, 348, 349, 351, 352, 353, 355, 356, 357, 361, 362, 364, 366, 374, 381, 382, 384, 390, 392, 395, 397, 399, 402, 403, 405, 406, 407, 409, 411, 412, 418, 420, 423, 425, 426, 427, 428, 429, 430, 433, 434, 435, 436, 439, 443, 456, 457, 460, 466, 467, 471, 480, 481, 484, 485, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 502, 503, 504, 505, 506, 507, 508, 509, 511, 512, 513, 514, 515, 517, 518, 520, 521, 522, 523, 525, 526, 529, 532, 533, 535, 536, 540, 542, 544, 545, 546, 548, 549, 550, 600, 602, 606, 609, 612, 613, 624, 626, 627, 628, 629, 667]; 


[depths,sorted_id] = sort(data_single_bank.unit_depths(selected_ids));
selected_ids_sorted = selected_ids(sorted_id);


FigH = figure('Position', get(0, 'Screensize'));
smooth_sigma = 0.5;
resize_factor = 4;
x_grid_large = linspace(x_grid(1), x_grid(end), length(x_grid)*resize_factor);
y_grid_large = linspace(y_grid(1), y_grid(end), length(y_grid)*resize_factor);
RF_smoothed = zeros(size(RF_mat,1)*resize_factor, size(RF_mat,2)*resize_factor, size(RF_mat,3));
for n=1:length(selected_ids_sorted)
    subplot(11,22,n);
    hold off
    RF_smoothed(:,:,selected_ids_sorted(n)) = imresize(imgaussfilt(RF_mat(:,:,selected_ids_sorted(n)),smooth_sigma), resize_factor);
    imagesc(x_grid_large, y_grid_large, RF_smoothed(:,:,selected_ids_sorted(n)));
    hold on
    axis xy;
    colormap('hot');
    title(sprintf('%d,%0.1f',selected_ids_sorted(n),mean(mean(RF_mat(:,:,selected_ids_sorted(n))))))
    colormap('hot');
    xline(0,'w', 'linewidth', 2);
    yline(0,'w', 'linewidth', 2);
    set(gca,'xlim',[-5,15])
    set(gca,'ylim',[-15,15])
    axis off;
end


figure('Renderer','opengl', 'Position', get(0, 'Screensize'))
y_lim = [0, 4000];
y_lim = [0, 20000];
subplot(2,5,[1 6])
hold on
threshold = 0.9;
hFills = cell(1,length(selected_ids_sorted));
for n=1:length(selected_ids_sorted)
    if depths(n)<y_lim(2)
        RF_normalized = RF_smoothed(:,:,selected_ids_sorted(n))/max(max(RF_smoothed(:,:,selected_ids_sorted(n))));

        patch( [0 0 0 0], [max(y_grid) min(y_grid) min(y_grid) max(y_grid)] , [max(depths) max(depths) 0 0], 'k', 'FaceAlpha', 0.1);
        patch( [max(x_grid) min(x_grid) min(x_grid) max(x_grid)] , [0 0 0 0], [max(depths) max(depths) 0 0], 'k', 'FaceAlpha', 0.1);
        [~, h_contourf] = contourf(x_grid_large, y_grid_large, RF_normalized+depths(n), [0.85, 0.85]+depths(n), 'LineColor', 'w', 'ContourZLevel', depths(n));

    end
end
view([10, 5])
set(gca,'zlim',[depths(1),depths(end)])

