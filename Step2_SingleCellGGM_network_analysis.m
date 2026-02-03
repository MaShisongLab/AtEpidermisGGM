% MATLAB code

% Load gene expression matrix and gene names
expr_matrix = h5read('at.epidermis.merged.h5','/x');
genes = h5read('at.epidermis.merged.h5','/genenames');

% Filter out genes expressed in fewer than 10 cells
idx = sum(expr_matrix > 0, 2) >= 10;
expr_matrix = expr_matrix(idx, :);
expr_matrix = expr_matrix';
genes = genes(idx);

% Perform SingleCellGGM network analysis and false discovery rate control
ggm = SingleCellGGM(expr_matrix, 20000, genes);
fdr = fdr_control(expr_matrix, ggm);

% Save results
save('ggm.epidermis.merged.mat', '-v7.3')
writetable(ggm.SigEdges(:,1:3), ...
           'ggm.epidermis.merged.edges.txt', 'Delimiter', 'tab')
writecell(ggm.gene_name, 'ggm.epidermis.merged.allgenes.txt')

