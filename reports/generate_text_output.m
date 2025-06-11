function generate_text_output(summary, params, thresholds, output_path, varargin)

    fprintf('Generating results file\n');
    
    assert(isa(summary, 'ExperimentReport'));

    filename = [tempname '.txt'];
    fid = fopen(filename, 'wt');    
    
    sc = startsWith(params.Results.cfg_type, 'sc');
    if (sc)
        tag = 'CB';
        assert(nargin == 8);
    else
        tag = 'UMI';
        assert(nargin == 4);
    end
    
    fprintf(fid, 'INPUT\n\n');
    
    fprintf(fid, '%-55s %s\n', 'CARLIN Amplicon:', params.Results.CARLIN_amplicon);
    fprintf(fid, '%-55s %s\n', 'CfgType:', params.Results.cfg_type);    
    if (iscell(params.Results.fastq_file))
        fastq_files = {params.Results.fastq_file{:}};
        fprintf(fid, '%-55s %s\n', 'Input FastQ File(s):', fastq_files{1});
        for i = 2:length(fastq_files)
            fprintf(fid, '%-55s %s\n', '', fastq_files{i});
        end
    else
        fprintf(fid, '%-55s %s\n', 'Input FastQ File(s):', params.Results.fastq_file);
    end
    if (sc)
        fprintf(fid, '%-55s %s\n', 'Reference barcode file:', params.Results.ref_CB_file);
    end
    fprintf(fid, '%-55s %s\n', 'Output directory:', params.Results.outdir);
    
    fprintf(fid, '\n%-55s %+10s %+6s\n\n', 'READ BREAKDOWN', 'Reads', '%');
    
    for x = fieldnames(summary.reads)'
        fprintf(fid, '%-55s %10d %6.0f\n', [strrep(x{1}, 'tag', tag) ':'], ...
            summary.reads.(x{1}), max(round(double(summary.reads.(x{1}))/double(summary.reads.in_fastq)*100),0));
    end
    
    fprintf(fid, '\nTHRESHOLDS\n\n');
    
    if (sc)
        for x = fieldnames(thresholds.CB)'
            fprintf(fid, '%-55s %10d\n', ['CB_' x{1} ':'], thresholds.CB.(x{1}));
        end
        fprintf(fid, '\n');
        for x = fieldnames(thresholds.UMI)'
            fprintf(fid, '%-55s %10d\n', ['UMI_' x{1} ':'], thresholds.UMI.(x{1}));
        end
    else
        for x = fieldnames(thresholds)'
            fprintf(fid, '%-55s %10d\n', ['UMI_' x{1} ':'], thresholds.(x{1}));
        end
    end
    
    if (sc)
        fprintf(fid, '\nCELL BARCODE BREAKDOWN\n\n');
    else
        fprintf(fid, '\nUMI BREAKDOWN\n\n');        
    end
    
    for x = fieldnames(summary.N)'
        fprintf(fid, '%-55s %10d\n', [erase(x{1}, '_tags') ':'], summary.N.(x{1}));
    end
    
    if (sc)
        % Gross - I didn't design ExperimentSummary to save this
        % information, so I have to recalculate here.
        fprintf(fid, '\nUMI BREAKDOWN\n\n');        
        tag_collection = varargin{1};
        tag_collection_denoised = varargin{2};
        tag_denoise_map = varargin{3};
        tag_called_allele = varargin{4};
        Numi_uncleaned = sum(cellfun(@length, {tag_collection.CBs.UMIs}));
        Numi_matched = sum(cellfun(@length, {tag_collection.CBs(ismember({tag_collection.CBs.CB}, keys(tag_denoise_map.CB))).UMIs}));
        Numi_cleaned = sum(cellfun(@length, {tag_collection_denoised.CBs.UMIs}));
        Numi_common = sum(cellfun(@(x) sum(cellfun(@sum, {x.SEQ_weight})>=thresholds.UMI.chosen), {tag_collection_denoised.CBs.UMIs}));
        Numi_called = sum(cellfun(@(x) size(vertcat(x.allele),1), {tag_called_allele.umi_call_result}));
        Numi_eventful = sum(cellfun(@(x) sum(cellfun(@(y) ~isequal(degap(y.get_seq), summary.CARLIN_def.seq.CARLIN), ...
                                         num2cell(vertcat(x.allele)))), {tag_called_allele.umi_call_result}));
        fprintf(fid, '%-55s %10d\n', 'uncleaned:', Numi_uncleaned);
        fprintf(fid, '%-55s %10d\n', 'matched:', Numi_matched);
        fprintf(fid, '%-55s %10d\n', 'cleaned:', Numi_cleaned);
        fprintf(fid, '%-55s %10d\n', 'common:', Numi_common);
        fprintf(fid, '%-55s %10d\n', 'called:', Numi_called);
        fprintf(fid, '%-55s %10d\n', 'eventful:', Numi_eventful);
    end    
    
    fprintf(fid, '\nPREFERENTIAL AMPLIFICATION\n\n');    
    
    fprintf(fid, '%-55s %10d\n', sprintf('Mean reads per edited %s:', tag), max(round(summary.reads.eventful_tags_total/summary.N.eventful_tags),0));
    fprintf(fid, '%-55s %10d\n', sprintf('Mean reads per unedited %s:', tag), ...
        max(round((summary.reads.called_tags_total-summary.reads.eventful_tags_total)/(summary.N.called_tags-summary.N.eventful_tags)),0));
    
    fprintf(fid, '\nSEQUENCE HETEROGENEITY\n\n');    
    
    fprintf(fid, '%-55s %10.0f\n', sprintf('Mean %% reads in edited %s folded into consensus:', tag), ...
        max(round(summary.reads.eventful_tags_allele/summary.reads.eventful_tags_total*100),0));
    fprintf(fid, '%-55s %10.0f\n', sprintf('Mean %% reads in unedited %s folded into consensus:', tag), ...
        max(round((summary.reads.called_tags_allele-summary.reads.eventful_tags_allele)/...
                  (summary.reads.called_tags_total-summary.reads.eventful_tags_total)*100),0));
          
    fprintf(fid, '\nALLELES\n\n');
    
    fprintf(fid, '%-55s %10d\n', 'Total (including template):', length(summary.alleles));
    fprintf(fid, '%-55s %10d\n', 'Singletons (including template):', sum(summary.allele_freqs==1));
    fprintf(fid, '%-55s %10.0f\n', sprintf('%% %ss edited:', tag), max(round(summary.N.eventful_tags/summary.N.called_tags*100),0));
    fprintf(fid, '%-55s %10d\n', 'Effective Alleles:', round(effective_alleles(summary)));
    fprintf(fid, '%-55s %10.2f\n', 'Diversity Index (normalized by all):', diversity_index(summary, false));
    fprintf(fid, '%-55s %10.2f\n', 'Diversity Index (normalized by edited):', diversity_index(summary, true));
    fprintf(fid, '%-55s %10.2f\n', 'Mean CARLIN potential (by allele):', ...
        max(mean(cellfun(@(x) summary.CARLIN_def.N.segments-length(Mutation.find_modified_sites(summary.CARLIN_def, x)), summary.alleles)),0));
    
    fclose(fid);
    copyfile(filename, [output_path '/Results.txt']);
    
    Mutation.ToFile(summary.CARLIN_def, summary.alleles, output_path, 'AlleleAnnotations.txt');
    
    s = cellfun(@(x) strjoin(x,','), summary.allele_colony, 'un', false);
    fid = fopen([output_path '/AlleleColonies.txt'], 'wt');
    if (~isempty(s))
        fprintf(fid, '%s\n', s{:});
    end
    fclose(fid);

end