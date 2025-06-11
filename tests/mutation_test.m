function tests = mutation_test
    tests = functiontests(localfunctions);
end

function test_mutation_allele_invariance(testCase)
    
    CARLIN_def = CARLIN_amplicon(parse_amplicon_file('OriginalCARLIN'));

    [folder, ~, ~] = fileparts(mfilename('fullpath'));
    
    golden_mut_list = Mutation.FromFile(CARLIN_def, sprintf('%s/data/OriginalCARLIN/Sanger75Annotations.txt', folder));    
    blacklist = cellfun(@(x) isequal(x, 'X') || isequal(x, '?'), golden_mut_list);    
    alleles = cellfun(@(x) Mutation.apply_mutations(CARLIN_def, x), golden_mut_list(~blacklist), 'un', false);
    [~, alleles] = cellfun(@(x) CARLIN_def.cas9_align(degap(x.get_seq)), alleles, 'un', false);
    called_mut_list = cellfun(@(x) Mutation.identify_cas9_events(CARLIN_def, x), alleles, 'un', false);
    verifyEqual(testCase, called_mut_list, golden_mut_list(~blacklist));
        
    Mutation.ToFile(CARLIN_def, called_mut_list, [folder '/Output/OriginalCARLIN'], 'Sanger75Reannotations.txt');    
    golden_strings = splitlines(fileread(sprintf('%s/data/OriginalCARLIN/Sanger75Annotations.txt', folder)));
    if (isempty(golden_strings{end}))
        golden_strings = golden_strings(1:end-1);
    end
    golden_strings = golden_strings(~blacklist);
    called_strings = splitlines(fileread(sprintf('%s/Output/OriginalCARLIN/Sanger75Reannotations.txt', folder)));
    if (isempty(called_strings{end}))
        called_strings = called_strings(1:end-1);
    end
    verifyEqual(testCase, called_strings, golden_strings);
   
end