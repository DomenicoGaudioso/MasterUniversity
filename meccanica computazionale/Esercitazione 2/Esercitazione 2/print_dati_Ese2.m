%% print_truss_model
fprintf(fid, ' ******* MODEL DATA ESERCITAZIONE 2**************\n\n\n');
% Print Nodal coordinates
%
fprintf(fid, '------------------------------------------------------ \n');
%
fprintf(fid, 'Number of nodes: %g\n', nnd );
fprintf(fid, 'Number of elements: %g\n', nel );
fprintf(fid, 'Number of nodes per element: %g\n', nne );
fprintf(fid, 'Number of degrees of freedom per node: %g\n', nodof);
fprintf(fid, 'Number of degrees of freedom per element: %g\n\n\n', eldof);

fprintf(fid, '------------------------------------------------------ \n');
%
fprintf(fid, 'thermal variation [°C]: %g\n',t);
fprintf(fid,'thermal expansion coefficient: %g\n',alpha);
fprintf(fid, 'Young modulus steel [KN/m^2]: %g\n',E);

fprintf(fid, '------------------------------------------------------ \n');
%
fprintf(fid, 'Node      X           Y \n');
for i=1:nnd
fprintf(fid,' %g,     %07.2f,     %07.2f\n',i, x(i), y(i));
end
fprintf(fid,'\n');
% Print element connectivity
%
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid, 'Element    Node_1     Node_2 \n');
for i=1:nel
fprintf(fid,' %g,         %g,          %g\n',i,conn(i,1),conn(i,2));
end
fprintf(fid,'\n');
%
% Print element property
%
fprintf(fid, '------------------------------------------------------ \n');
%
fprintf(fid, 'Element     Area\n');
for i=1:nel
fprintf(fid,   '%g,        %g\n',i, A(i));
end
fprintf(fid,'\n');
% Print Nodal loads_concentrated
%
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid, 'Node   concentrated_load_X     concentrated_load_Y\n');
for i=1:nnd
fprintf(fid,  '%g,         %07.2f,                  %07.2f\n',i, p(i,1), p(i,2));
end
%
% Print Nodal loads_distributed
%
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid, 'Node   distributed_load_X     distributed_load_Y\n');
for i=1:nnd
fprintf(fid,  '%g,         %07.2f,                  %07.2f\n',i, pq(i,1), pq(i,2));
end
%
% Print thermal forces
%
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid, 'Node     thermal_load_X     thermal_load_Y\n');
for i=1:nnd
fprintf(fid,  '%g,         %07.2f,             %07.2f\n',i, pt(i,1), pt(i,2));
end
%
% Print Total forces
%
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid, 'Node   Total forces_load_X      Total forces_load_Y\n');
for i=1:nnd
fprintf(fid,  '%g,          %07.2f,                  %07.2f\n',i, P(i,1), P(i,2));
end
%
% Print Nodal freedom
%
fprintf(fid, '------------------------------------------------------ \n');
%
fprintf(fid, 'Node      disp_U            disp_V\n');
for i=1:nnd
fprintf(fid,  '%g,        %g,         %g\n',i, vx(i), vy(i));
end
fprintf(fid,'\n');
%
% Print Mass element
%
fprintf(fid, '------------------------------------------------------ \n');
%
fprintf(fid, 'element    Mass [KN]\n');
for i=1:nel
fprintf(fid,  '%g,        %g\n',i,M(i));
end
fprintf(fid,'\n');
% Print Mass element
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid,'\n');
fprintf(fid,'Total Mass, Mtot[KN] = %g\n',Mtot);
fprintf(fid,'\n');
%
%
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid,'\n');
fprintf(fid,'Total number of active degrees of freedom, n = %g\n',c);
fprintf(fid,'\n');
%