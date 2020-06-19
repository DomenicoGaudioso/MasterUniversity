%% print_truss_results.m
%
fprintf(fid, '-------------------------------------------------------- \n');
fprintf(fid, ' \n\n\n ******* PRINTING ANALYSIS RESULTS **************\n\n\n');
%
%
%
% Print global force vector
%
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid,'Global force vector f \n');
fprintf(fid,' %g\n',f);
fprintf(fid,'\n');
%
%
% Print Epslon solution vector
%
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid,'epslon deformation:  delta \n');
fprintf(fid,' %8.5f\n',epslon);
fprintf(fid,'\n');
%
% Print Reaction force
%
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid,'Reaction:  realction \n');
fprintf(fid,' %8.5f\n',r);
fprintf(fid,'\n');
%
% Print Members actions
%
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid, 'Members actions \n');
fprintf(fid, 'element    Normal force            action\n');
for i=1:nel
if N(i) > 0
fprintf(fid,  ' %g,         %8.5f,         %s\n',i, N(i), 'Tension');
else
fprintf(fid,  ' %g,         %8.5f,         %s\n',i, N(i), 'Compression');
end
end