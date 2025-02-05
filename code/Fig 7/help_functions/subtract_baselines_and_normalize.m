function orange_matrix_minus_baseline_normalized = subtract_baselines_and_normalize(orange_matrix,unit_baselines)
    orange_matrix_minus_baseline = orange_matrix - repmat(unit_baselines,[1,size(orange_matrix,2)]);
    orange_matrix_minus_baseline_normalized = orange_matrix_minus_baseline./repmat((max(abs(orange_matrix_minus_baseline),[],2)),[1,size(orange_matrix,2)]);
