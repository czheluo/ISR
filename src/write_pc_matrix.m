% Write the matrix of principal components pc to a space-delimited text
% file, in which each line of the text file corresponds to a sample, and
% each space-delimited column corresponds to a principal component. This
% function is used in MATLAB script geno_pca.m.
function write_rot_matrix (pc, id, file)
  [n k] = size(pc);
  f = fopen(file,'w');
  fprintf(f,'id');
  fprintf(f,' PC%d',1:k);
  fprintf(f,'\n');
  for i = 1:n
    fprintf(f,'%s',id{i});
    fprintf(f,' %0.4f',pc(i,:));
    fprintf(f,'\n');
  end
  fclose(f);

