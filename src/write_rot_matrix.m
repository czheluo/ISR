% Write the rotation matrix R to a space-delimited text file, in which
% each line of the text file corresponds to a marker (SNP), and each
% space-delimited column corresponds to a principal component. This
% function is used in MATLAB script geno_pca.m.
function write_rot_matrix (R, marker, file)
  [p k] = size(R);
  f = fopen(file,'w');
  fprintf(f,'marker');
  fprintf(f,' PC%d',1:k);
  fprintf(f,'\n');
  for i = 1:p
    fprintf(f,'%s',marker{i});
    fprintf(f,' %0.6f',R(i,:));
    fprintf(f,'\n');
  end
  fclose(f);
