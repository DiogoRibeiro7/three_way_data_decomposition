function generate_scr_reference_fixtures(reference_repo, fixture_root)
%GENERATE_SCR_REFERENCE_FIXTURES Generate deterministic SCR MATLAB outputs.
%
% reference_repo must point to a local clone of:
%   https://github.com/moniar412/SCR3waydata
%
% fixture_root must point to inst/extdata/matlab_reference in this repository.
%
% The script deliberately uses fixed numeric inputs rather than random
% generation so cross-language differences cannot be attributed to RNGs.

addpath(reference_repo);

input_dir = fullfile(fixture_root, 'inputs');
output_dir = fullfile(fixture_root, 'outputs');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

X = dlmread(fullfile(input_dir, 'X.csv'), ',');
U0 = dlmread(fullfile(input_dir, 'U.csv'), ',');
TB2 = dlmread(fullfile(input_dir, 'TB_s2.csv'), ',');
TB3 = dlmread(fullfile(input_dir, 'TB_s3.csv'), ',');
TC3 = dlmread(fullfile(input_dir, 'TC_s3.csv'), ',');
SV3 = dlmread(fullfile(input_dir, 'SV_s3.csv'), ',');
SO3 = dlmread(fullfile(input_dir, 'SO_s3.csv'), ',');

eps_value = 1e-8;
display_flag = 0;

% H baseline.
[U_h, Mmu_h, Sig_h, dif_h, like_h, bic_h, it_h] = ...
    mixhom(X, U0, eps_value, display_flag);

write_numeric(fullfile(output_dir, 'H_U.csv'), U_h);
write_numeric(fullfile(output_dir, 'H_Mmu.csv'), Mmu_h);
write_numeric(fullfile(output_dir, 'H_Sig.csv'), Sig_h);
write_numeric(fullfile(output_dir, 'H_scalars.csv'), ...
    [dif_h, like_h, bic_h, it_h]);

% S2 baseline.
[U_s2, TB_s2, SV_s2, Y_s2, like_s2, bic_s2] = ...
    t2mixt(X, U0, TB2, eps_value, display_flag);
M_s2 = Y_s2 * TB_s2';

write_numeric(fullfile(output_dir, 'S2_U.csv'), U_s2);
write_numeric(fullfile(output_dir, 'S2_SV.csv'), SV_s2);
write_numeric(fullfile(output_dir, 'S2_M.csv'), M_s2);
write_numeric(fullfile(output_dir, 'S2_scalars.csv'), ...
    [like_s2, bic_s2]);

% S3 baseline. X has 4 columns = 2 variables x 2 occasions.
[U_s3, TB_s3, TC_s3, SO_s3, SV_s3, Y_s3, like_s3, bic_s3] = ...
    t3mixs(X, U0, TB3, TC3, SV3, SO3, eps_value, display_flag);
M_s3 = Y_s3 * kron(TC_s3', TB_s3');
Sigma_s3 = kron(SO_s3, SV_s3);

write_numeric(fullfile(output_dir, 'S3_U.csv'), U_s3);
write_numeric(fullfile(output_dir, 'S3_M.csv'), M_s3);
write_numeric(fullfile(output_dir, 'S3_Sigma.csv'), Sigma_s3);
write_numeric(fullfile(output_dir, 'S3_scalars.csv'), ...
    [like_s3, bic_s3]);

fid = fopen(fullfile(output_dir, 'REFERENCE_SOURCE.txt'), 'w');
fprintf(fid, 'Repository: moniar412/SCR3waydata\n');
fprintf(fid, 'Generated with mixhom.m, t2mixt.m, t3mixs.m\n');
fprintf(fid, 'Tolerance: %.17g\n', eps_value);
fclose(fid);

end


function write_numeric(path, value)
% Use a representation supported by both MATLAB and GNU Octave.
dlmwrite(path, value, ',', 'precision', '%.17g');
end
