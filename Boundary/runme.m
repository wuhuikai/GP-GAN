clear;clc;close all;
addpath('./Mesh2d v24/');

test_list = table2cell(readtable('../random_test_list.txt', 'Delimiter', ';', 'ReadVariableNames', false));
result_folder = '../image_blending_comparation_result';
mkdir(result_folder);

zero_masks = 0;
zero_masks_id = [];
total_size = length(test_list);
for idx = 1:total_size
    fprintf('Processing %d/%d ...\n', idx, total_size);
    
    src = im2double(imread(test_list{idx, 1}));
    ftrg = im2double(imread(test_list{idx, 2}));
    mask = logical(imread(test_list{idx, 3}));
    
    img_name_template = sprintf('%s/%d_%%s.png', result_folder, idx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                    Naive                    %  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\t Copy and Paste\n')
    res_naive = NaiveBlending(ftrg, src, mask);
    imwrite(res_naive, sprintf(img_name_template, 'copy-paste'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %              Poisson                        %  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\t Poisson Image Editing\n')
    res_poisson = PoissonBlending(ftrg, src, mask);
    imwrite(res_poisson, sprintf(img_name_template, 'poisson'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %             Convolution Pyramid             %  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\t Convolution Pyramid\n')
    res_convpyr = ConvPyrBlending(ftrg, src, mask);
    imwrite(res_convpyr, sprintf(img_name_template, 'conv-pyramid'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %            Modified Poisson                %  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\t Modified Poisson Image Editing\n')
    res_modpoisson = ModPoissonBlending(ftrg, src, mask);
    imwrite(res_modpoisson, sprintf(img_name_template, 'modified-poisson'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %            Multi-Splines                    %  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\t Multi Splines\n')    
    if sum(double(mask(:))) == 0
        res_mults = ftrg;
        zero_masks = zero_masks + 1;
        zero_masks_id = [zero_masks_id, idx];
    else
        offset = MultiSplinesBlending(ftrg, ftrg, src, mask, 1, 1);
        res_mults = res_naive + offset;
    end
    imwrite(res_mults, sprintf(img_name_template, 'multi-splines'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %            Mean-Value Coordinates           %  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\t Mean-Value Coordinates\n')
    
    if sum(double(mask(:))) == 0
        res_mvc = ftrg;
    else
        res = MVCBlending(ftrg, src, mask);
        mvc_mask = repmat(double(mask), [1, 1, 3]);
        res_mvc = res_naive .* (1-mvc_mask) + res .* mvc_mask;
    end
    imwrite(res_mvc, sprintf(img_name_template, 'mvc'));
    
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %            MultiBand_Blending               %  
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %multi_band is not proper for object insertion, we just give a example here
%         fprintf('\t MultiBand Blending\n')
%         level = 4;%set the pyramid level
%         res = MultiBandBlending(ftrg, src, mask, level);
%         mb_mask = repmat(double(mask), [1, 1, 3]);
%         res_mbb = res_naive .* (1-mb_mask) + res .* mb_mask;
%         imwrite(res_mbb, sprintf(img_name_template, 'multi-band'));
end

fprintf('Total %d zero masks\n', zero_masks);
save zero_masks.mat zero_masks_id;