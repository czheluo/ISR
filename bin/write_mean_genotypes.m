% Write the mean genotypes to a space-delimited text file, in which
% each line of the text file corresponds to a marker (SNP). This
% function is used in MATLAB script geno_pca.m.
function write_mean_genotypes (y, marker, file)
  p = length(y);
  f = fopen(file,'w');
  fprintf(f,'marker mean\n');
  for i = 1:p
    fprintf(f,'%s %0.6f\n',marker{i},y(i));
  end
  fclose(f);
