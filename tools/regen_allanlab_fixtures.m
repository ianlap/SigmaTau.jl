function regen_allanlab_fixtures()
%REGEN_ALLANLAB_FIXTURES Regenerate the AllanLab cross-validation fixture.
%
%   Reads `reference/validation/stable32gen.DAT` (10-line header, 8192
%   phase samples), runs every AllanLab deviation kernel that returns
%   EDF/CI on a per-octave m grid, recomputes the CI at p = 0.95 to match
%   Stable32 / SigmaTau convention (AllanLab's default is 1-sigma), and
%   writes a single CSV at
%   `reference/validation/allanlab_out/allanlab_data_full.csv`.
%
%   Schema:
%       Type,AF,Tau,N,Alpha,Sigma,SigmaMin,SigmaMax,EDF
%
%   Type labels match the Stable32 fixture so cross-filtering by Type
%   pulls the same family across all three references.
%
%   Run from the repo root via either of:
%       matlab -batch "run('tools/regen_allanlab_fixtures.m')"
%       (MATLAB GUI) cd('tools'); regen_allanlab_fixtures
%
%   Self-locating: paths are derived from this file's own location, so
%   the working directory at invocation does not matter.

    %-- Locate repo root + AllanLab package via this script's path.
    script_dir   = fileparts(mfilename('fullpath'));
    repo_root    = fileparts(script_dir);
    allanlab_dir = fullfile(repo_root, 'legacy', 'masterclock-kflab', 'AllanLab');
    addpath(allanlab_dir);
    import allanlab.*

    %-- Inputs
    dat_path = fullfile(repo_root, 'reference', 'validation', 'stable32gen.DAT');
    out_dir  = fullfile(repo_root, 'reference', 'validation', 'allanlab_out');
    out_csv  = fullfile(out_dir, 'allanlab_data_full.csv');

    if ~isfolder(out_dir)
        mkdir(out_dir);
    end

    %-- Load phase fixture: skip the 10-line header.
    fid = fopen(dat_path, 'r');
    if fid < 0
        error('Cannot open %s', dat_path);
    end
    cleanup = onCleanup(@() fclose(fid));
    for i = 1:10
        fgetl(fid);
    end
    x = fscanf(fid, '%g');
    clear cleanup;  % triggers fclose

    N = length(x);
    if N ~= 8192
        warning('Expected 8192 phase samples, got %d', N);
    end
    tau0 = 1.0;

    %-- Per-octave m grid: 1, 2, 4, ..., 4096 — matches the Stable32 /
    %   allantools fixtures so the rows align across all three CSVs.
    m_list = 2.^(0:floor(log2(N/2)));

    %-- Function handle / Stable32-compatible Type label per kernel.
    specs = { ...
        @adev,     'Overlapping Allan'; ...
        @mdev,     'Modified Allan'; ...
        @hdev,     'Overlapping Hadamard'; ...
        @mhdev,    'Modified Hadamard'; ...
        @tdev,     'Time'; ...
        @ldev,     'Hadamard Time'; ...
        @totdev,   'Total'; ...
        @mtotdev,  'Modified Total'; ...
        @htotdev,  'Hadamard Total'; ...
        @mhtotdev, 'Modified Hadamard Total' ...
    };

    p_target = 0.683;  % match the Stable32 fixture's "Confidence Factor = 0.683"
                        % header (see reference/validation/stable32out/stable32out.txt).
                        % AllanLab kernels also default to p = 0.683 internally; we
                        % pass it explicitly so the convention is visible at this
                        % seam rather than implicit.

    fid_out = fopen(out_csv, 'w');
    if fid_out < 0
        error('Cannot open %s for writing', out_csv);
    end
    out_cleanup = onCleanup(@() fclose(fid_out));

    fprintf(fid_out, 'Type,AF,Tau,N,Alpha,Sigma,SigmaMin,SigmaMax,EDF\n');

    fprintf('AllanLab fixture regen — N=%d, %d kernels × %d AFs\n', ...
            N, size(specs, 1), length(m_list));

    for k = 1:size(specs, 1)
        fn    = specs{k, 1};
        label = specs{k, 2};

        fprintf('  %s ... ', label);
        t_start = tic;
        try
            [tau, sigma, edf, ~, alpha] = fn(x, tau0, m_list);
        catch ME
            fprintf('FAILED (%s)\n', ME.message);
            continue;
        end

        % Force column orientation.
        tau   = tau(:);
        sigma = sigma(:);
        edf   = edf(:);
        alpha = alpha(:);

        % Recompute CI at p_target. compute_ci falls back to a Gaussian
        % formula when EDF is non-finite, in which case it needs Neff;
        % use a generous (N - L_min) so the fallback bound stays loose
        % rather than spuriously tight. The chi-square branch (used
        % whenever EDF is finite and >= 1) ignores Neff entirely.
        Neff = max(N - 2 * m_list(:), 1);
        Neff = Neff(1:length(sigma));
        ci   = compute_ci(sigma, edf, p_target, alpha, Neff);

        n_written = 0;
        for j = 1:length(tau)
            if isnan(sigma(j))
                continue
            end
            fprintf(fid_out, '%s,%d,%.17e,%d,%g,%.17e,%.17e,%.17e,%.6e\n', ...
                    label, m_list(j), tau(j), N, alpha(j), sigma(j), ...
                    ci(j, 1), ci(j, 2), edf(j));
            n_written = n_written + 1;
        end

        fprintf('%d rows (%.1fs)\n', n_written, toc(t_start));
    end

    fprintf('\nWrote %s\n', out_csv);
end
