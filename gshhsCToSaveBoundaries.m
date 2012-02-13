% example script showing how to produce a 'SaveBoundaries' file

% read in and mangle the gshhs-c dataset
cd ~/code/gshhs
load gshhs_c_t1l0_points_unique.dat
glon = gshhs_c_t1l0_points_unique(:,1);
glat = gshhs_c_t1l0_points_unique(:,2);
glon=mod(glon+360,360);

% open a new file for output, and write in the format required for mpi-scvt
cd ~/code/mpi-scvt/branches/boundary-projections
filename = "SaveBoundaries.gshhs-c";
fid = fopen (filename, "w");
for ii = 1 : length(glon)-1
	fprintf(fid,"%f %f\n",glat(ii),glon(ii));
	fprintf(fid,"%f %f\n",glat(ii+1),glon(ii+1));
end
fclose (fid);
