function para = prepare_para(para)
fprintf('Loading parameters...');tic

matObj = matfile([para.dir.load_kSpace_dir,para.dir.load_kSpace_name]);
varlist = who(matObj,'kSpace_info');
if ~isempty(varlist)
    load([para.dir.load_kSpace_dir,para.dir.load_kSpace_name],'kSpace_info')
    para.Recon.nSMS = max(kSpace_info.phase_mod) + 1;
    para.kSpace_info = kSpace_info;
end

para.time = datestr(clock,'yymmdd_hhMMSS');

disp('RawData:')
disp([para.dir.load_kSpace_dir, para.dir.load_kSpace_name])

if para.self_gating.flag
    para.dir.save_recon_img_name = [para.dir.load_kSpace_name(1:end-4),'_binned.mat'];
else
    para.dir.save_recon_img_name = [para.dir.load_kSpace_name(1:end-4),'_all_data.mat'];
end

para.dir.save_recon_img_dir = strcat(pwd,'/ReconData/');
if isempty(dir(para.dir.save_recon_img_dir))
    mkdir(para.dir.save_recon_img_dir);
end

if isempty(dir([pwd,'/RawData/']))
    mkdir([pwd,'/RawData/']);
end

para.Recon.epsilon = eps('single');
para.Recon.step_size = 2;
para.Recon.break = 0;

para.CPUtime.load_para_time = toc;toc;fprintf('\n');
end
