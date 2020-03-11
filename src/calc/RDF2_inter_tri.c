//
//----------------------------------------------------------------------
#define cmax_length 1000
#define species_length 4
#define pi 3.14159265
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>


/*
    Calculate cos(angle) between vecs
*/
double proj(double Ax,double Ay,double Az,double Bx,double By,double Bz);
double proj(double Ax,double Ay,double Az,double Bx,double By,double Bz)
{
	double dot,norm,res;
	dot = Ax*Bx + Ay*By + Az*Bz;
	norm = Ax*Ax + Ay*Ay + Az*Az;
	res = dot/norm;
	return res;
}

/*
    A type to hold position data
*/
struct Pos {
	double *x, *y, *z;
};

/*
    A type to hold all the data from a lammps file.
*/
struct Lammps_File {
    struct Pos R;
    struct Pos sR;
    int *mol, *type;
    int nx, ny, nz;
    int natom;
    int dt;
    double xlo, xhi, ylo, yhi, zlo, zhi;
    double xlo_bound, xhi_bound, ylo_bound, yhi_bound, zlo_bound, zhi_bound;
    double xy, xz, yz;
};

/*
    Will get the histogram count.

    This will wrap the coordinates first, using the scaled vectors, and
    then use these scaled vectors to find the distances to build the
    histogram.
    
    Inputs:
        * data <Lammps_File> => The data from the snapshot file
        * bin <int> => The number of bins
    Outputs:
        <int*> The bin counts
*/
void get_histogram(struct Lammps_File data, int bins, int *counts, double dr, int atom_type) {
    int diff_mol, correct_atom_type, bh, i, j;
    double dist, rx, ry, rz;
    double svectorx, svectory, svectorz;

    // Init the histogram at the beginning of each step
    for (i = 0; i<bins; ++i){ counts[i] = 0; }

    // Loop over all pairs of atoms
	for(j = 0; j<data.natom - 1; ++j) {
        printf("%d/%d                   \r", j, data.natom - 1);

		for(i = j + 1;i<data.natom; ++i) {

            // If the pair of atoms are the correct type and on different mols then count them
            diff_mol = data.mol[i] != data.mol[j];
            correct_atom_type = data.type[i] == atom_type && data.type[j] == atom_type;
			if(diff_mol && correct_atom_type) {

                // svectors are for the coordinate wrapping
				svectorx = data.sR.x[j] - data.sR.x[i];
				svectory = data.sR.y[j] - data.sR.y[i];
				svectorz = data.sR.z[j] - data.sR.z[i];

				svectorx = svectorx - rint(svectorx);
				svectory = svectory - rint(svectory);
				svectorz = svectorz - rint(svectorz);

                // The actual displacement vectors
				rx = (data.xhi - data.xlo) * svectorx + (data.xy) * svectory + (data.xz) * svectorz;
				ry = (data.yhi - data.ylo) * svectory + (data.yz) * svectorz;
				rz = (data.zhi - data.zlo) * svectorz;

				dist = sqrt( rx*rx + ry*ry + rz*rz );

                // Add to the count
				bh = floor(dist / dr);
				if (bh < bins) { counts[bh] = counts[bh] + 2; }
			}
		}
	}
}


FILE *fp,  *fout;

// fp related
char file_path[cmax_length];
char folder[cmax_length], name[cmax_length];

// dummies
char buffer[cmax_length];
int int_buffer;
double d_buffer;

//
int i, j, k;
int start, step, stop;
int atoms;

// unwrap
int nx, ny, nz;
double xmin, xmax, ymin, ymax, zmin, zmax, Lx, Ly, Lz;

//
int dt;
int bins;
int *hist;
double re, dr;
double *g, *rn;
double xPBC, yPBC, zPBC, r2, r;

int type_N_1, type_N_2, sub_atoms_1 = 0, sub_atoms_2 = 0;
int *type_list_1, *type_list_2, *type, *belongs_to_type_list_1, *belongs_to_type_list_2;
int scale;

//
double V;

//
int flag;
double d_temp;
double S[4];
double MIN_x, MAX_x, MIN_y, MAX_y;

double ax, ay, az, bx, by, bz, cx, cy, cz, norm;
double ax_hat, ay_hat, az_hat, bx_hat, by_hat, bz_hat, cx_hat, cy_hat, cz_hat;
double a0x_hat = 1.0, a0y_hat = 0.0, a0z_hat = 0.0;
double b0x_hat = 0.0, b0y_hat = 1.0, b0z_hat = 0.0;
double c0x_hat = 0.0, c0y_hat = 0.0, c0z_hat = 1.0;

double projaa0;
double projab0;
double projac0;
double projba0;
double projbb0;
double projbc0;
double projca0;
double projcb0;
double projcc0;
double nom1, nom2, nom3, denom;
double a, b, c, a_sign, b_sign, c_sign, a_norm, b_norm, c_norm;
double a_norm_REF, b_norm_REF, c_norm_REF;
//int nx, ny, nz;

double ih_11, ih_12, ih_13, ih_22, ih_23, ih_33, *sx, *sy, *sz;
double svectorx, svectory, svectorz, rx, ry, rz;


int main(int argc, char *argv[])
{
	int bh;

	//double dr,re,*g,*rn;
	//------------------------------------------------------------------
    if(argc == 1)
    {
        printf("\n./RDF2 <folder> <name> <start> <step> <stop> <dr> <re> <type_list>\n\n");exit(-1);
    }

	// read arguments
	sprintf(folder, "%s", argv[1]);
	sprintf(name, "%s", argv[2]);
	start = atoi(argv[3]);
	step = atoi(argv[4]);
	stop = atoi(argv[5]);
    dr = atof(argv[6]);
    re = atof(argv[7]);

    dt = ((stop - start) / step) + 1;

    // rdf initialization
	// number of bins
	bins = floor( (re / dr) + 0.5);
	// bin count array
	// RDF
	g = (double*)  malloc( bins * sizeof(double) );
	// shell radius
	rn = (double*)  malloc( bins*sizeof(double) );
	hist = (int*) malloc( bins * sizeof(int) );
    
	for (i=0; i<bins; ++i) {
        rn[i] = (i + 0.5) * dr;
        g[i] = 0.0;
    }

	// read type config
	sprintf(file_path, "%s", argv[8]);

    // Check if the folder exists
	fp = fopen(file_path, "r");
	if(fp == NULL){
        printf("Cannot locate %s\n", file_path);
        exit(-1);
    }

	// locate group 1
	while( fgets(buffer, cmax_length, fp) != NULL) {
        if(strcmp(buffer, "type_1\n") == 0)
            break;
    }
	fgets(buffer, cmax_length, fp);
    sscanf(buffer,"%d",&type_N_1);
	type_list_1 = (int*) malloc(type_N_1 * sizeof(int));
	for(i=0; i<type_N_1; ++i)
	{
		fgets(buffer, cmax_length, fp);
        sscanf(buffer, "%d", &type_list_1[i]);
	}
	rewind(fp);

	// locate group 2
	while(fgets(buffer, cmax_length, fp) != NULL) {
        if(strcmp(buffer,"type_2\n") == 0)
            break;
    }
	fgets(buffer, cmax_length, fp);
    sscanf(buffer, "%d", &type_N_2);
	type_list_2 = (int*) malloc(type_N_2 * sizeof(int));
	for(i = 0;i<type_N_2;++i)
	{
		fgets(buffer, cmax_length, fp);
        sscanf(buffer, "%d", &type_list_2[i]);
	}
	fclose(fp);

	// scale
	scale = atoi(argv[9]);

	// read each proposed dump file to check if they exist
    for (k=start; k<=stop; k=k+step) {
    	sprintf(file_path, "%s/%s%d.dat", folder, name, k);
    	fp = fopen(file_path, "r");
    	if(fp == NULL){
            printf("Cannot locate %s\n", file_path);
            exit(-1); 
        }
    }

	for(i=0; i<3; ++i) {
        fgets(buffer, cmax_length, fp);
    }
	fgets(buffer, cmax_length, fp);
    sscanf(buffer, "%d", &atoms);
	fclose(fp);

    //////////////////////////////////////////////
    //                                          //
    // Init the Lammps File -this holds all the //
    //                       Lammps File data.  //
    //                                          //
    //////////////////////////////////////////////
    struct Lammps_File data;
    // All positions
    data.R.x = (double*) malloc( atoms * sizeof(double) );
    data.R.y = (double*) malloc( atoms * sizeof(double) );
    data.R.z = (double*) malloc( atoms * sizeof(double) );
    // All scaled positions
    data.sR.x = (double*) malloc( atoms * sizeof(double) );
    data.sR.y = (double*) malloc( atoms * sizeof(double) );
    data.sR.z = (double*) malloc( atoms * sizeof(double) );
    // The molecule number and atom type
    data.mol = (int*) malloc( atoms * sizeof(int) );
    data.type = (int*) malloc( atoms * sizeof(int) );
    data.natom = atoms;

    
    // Whether the atom is of type 1 or 2 -init arrays
    belongs_to_type_list_1 = (int*) malloc( atoms * sizeof(int) );
    for(i = 0; i<atoms; ++i)  belongs_to_type_list_1[i] = 0;
    belongs_to_type_list_2 = (int*) malloc( atoms * sizeof(int) );
    for(i = 0; i<atoms; ++i)  belongs_to_type_list_2[i] = 0;


    

	// start reading the dump files
	for(k=start; k<=stop; k = k + step)
	{
        /*
        Will compare 2 doubles and return an int depending on which is bigger

        Inputs:
            * a <double> => First double to compare
            * b <double> => Second double to compare
        Outputs:
            <-1, 0, 1> -1 if a < b; 0 if a == b; 1 if a > b
        */
        int compare_double(const double a, const double b) {
            if (a == b)  return 0;
            else if (a < b) return -1;
            else return 1;
        }


        /*
        Will read a lammps file from a filepath and return a Lammps_File object

        This will simple iterate over all lines and call fgets on each one to parse
        them.
        The unit cell vectors will also be calculated and stored in the Lammps_File
        struct.

        Inputs:
            * filepath <char[]> => The path pointing towards the lammps file.
        Outputs:
            <Lammps_File> The struct holding all the lammps file data.
        */
        Lammps_File read_lammps_file(char filepath[]) {
            struct Lammps_File data;
            char buffer[cmax_length];
            double S[4];

            
            // Create the filepath
            fp = fopen(filepath, "r");

            //Check the file exists
            if (fp == NULL) {
                printf("Cannot Locate %s\n", filepath);
                exit(-1);
            }
            // Read line by line
            //
            // Get step num and natom
            fgets(buffer, cmax_length, fp);
            fgets(buffer, "%d", &data.dt);
            fgets(buffer, cmax_length, fp);
            fgets(buffer, "%d", &data.natom);
            //
            // Get cell vec parameters
		    fgets(buffer, cmax_length, fp);
            sscanf(buffer, "%lf\t%lf\t%lf", &data.xlo_bound, &data.xhi_bound, &data.xy);
		    fgets(buffer, cmax_length, fp);
            sscanf(buffer, "%lf\t%lf\t%lf", &data.ylo_bound, &data.yhi_bound, &data.xz);
		    fgets(buffer, cmax_length, fp);
            sscanf(buffer, "%lf\t%lf\t%lf", &data.zlo_bound, &data.zhi_bound, &data.yz);

            // Calc vecs from parameters
            S[0] = 0.0;
		    S[1] = data.xy;
		    S[2] = data.xz;
		    S[3] = data.xy + data.xz;

            qsort(S, 4, sizeof(double), compare);

		    MIN_x = S[0];
            MAX_x = S[3];

            // Read atoms




        }
        // Create the filepath string
		sprintf(file_path, "%s/%s%d.dat", folder, name, k);
		// open dump file
		fp = fopen(file_path, "r");

        // Check if the file is there
		if(fp == NULL){
            printf("Cannot locate %s\n", file_path);
            exit(-1);
        }

        // Will read the header
		for(i=0; i<3; ++i){
            fgets(buffer, cmax_length, fp);
        }

		fgets(buffer, cmax_length, fp);
		fgets(buffer, cmax_length, fp);

		fgets(buffer, cmax_length, fp); 
        sscanf(buffer, "%lf\t%lf\t%lf", &data.xlo_bound, &data.xhi_bound, &data.xy);
		fgets(buffer, cmax_length, fp);
        sscanf(buffer, "%lf\t%lf\t%lf", &data.ylo_bound, &data.yhi_bound, &data.xz);
		fgets(buffer, cmax_length, fp);
        sscanf(buffer, "%lf\t%lf\t%lf", &data.zlo_bound, &data.zhi_bound, &data.yz);

		fgets(buffer, cmax_length, fp);

        void function(Lammps_File *data) {
            double S[4];

            // supecell
            //  
            S[0] = 0.0;
            S[1] = data.xy;
            S[2] = data.xz;
            S[3] = data.xy + data.xz;

        }
		// supecell
		//
		S[0] = 0.0;
		S[1] = data.xy;
		S[2] = data.xz;
		S[3] = data.xy + data.xz;

		// bubble sort
		flag = 1;
		while(flag == 1)
		{
            flag = 0;
			for(j=1; j<4; ++j)
			{
                if(S[j] < S[j-1])
				{
                    d_temp = S[j-1];
                    S[j-1] = S[j];
                    S[j] = d_temp;
					flag = 1;
                }
			}
		}
		MIN_x = S[0];
        MAX_x = S[3];

		if(data.yz<0.0){ MIN_y = data.yz; } else { MIN_y = 0.0; }
		if(data.yz>0.0){ MAX_y = data.yz; } else { MAX_y = 0.0; }

		data.xlo = data.xlo_bound - MIN_x;
		data.xhi = data.xhi_bound - MAX_x;
		data.ylo = data.ylo_bound - MIN_y;
		data.yhi = data.yhi_bound - MAX_y;
		data.zlo = data.zlo_bound;
		data.zhi = data.zhi_bound;
		
		// read atoms
		for(i = 0;i<atoms;++i)
		{
			// read
			fgets(buffer, cmax_length, fp);
            sscanf(buffer,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%d\t%d\t%d", &int_buffer, &data.mol[i], &data.type[i], 
                                                                   &data.R.x[i], &data.R.y[i], &data.R.z[i],
                                                                   &data.nx, &data.nz, &nz);
		}

		// check once
		if(k == start)
		{
			for(i=0; i<atoms; ++i)
			{
				for(j=0; j<type_N_1; ++j)
				{
					if(data.type[i] == type_list_1[j]) {
                        belongs_to_type_list_1[i] = 1;
                        ++sub_atoms_1;
                        break;
                    }
				}
				for(j=0; j<type_N_2; ++j)
				{
					if(data.type[i] == type_list_2[j]) {
                        belongs_to_type_list_2[i]=1;
                        ++sub_atoms_2;
                        break;
                    }
				}
			}
		}

		// shift
		for(j = 0;j<atoms;++j)
		{
			data.R.x[j] = data.R.x[j]-data.xlo;data.R.y[j] = data.R.y[j]-data.ylo;data.R.z[j] = data.R.z[j]-data.zlo;
		}

		// oblique supercell vectors
		ax = (data.xhi - data.xlo);// - (0.0);
		ay = (0.0)    ;// - (0.0);
		az = (0.0)    ;// - (0.0);

		bx = (data.xy)     ;// - (0.0);
		by = (data.yhi - data.ylo);// - (0.0);
		bz = (0.0)    ;// - (0.0);

		cx = (data.xz)     ;// - (0.0);
		cy = (data.yz)     ;// - (0.0);
		cz = (data.zhi - data.zlo);// - (0.0);

		V = ax * by * cz;

		// normalize
		norm = sqrt(ax * ax + ay * ay + az * az);
		ax_hat = ax/norm;
		ay_hat = ay/norm;
		az_hat = az/norm;

		norm = sqrt(bx * bx + by * by + bz * bz);
		bx_hat = bx/norm;
		by_hat = by/norm;
		bz_hat = bz/norm;

		norm = sqrt(cx * cx + cy * cy + cz * cz);
		cx_hat = cx/norm;
		cy_hat = cy/norm;
		cz_hat = cz/norm;

		// calculate projections
		projaa0 = proj(ax_hat, ay_hat, az_hat, a0x_hat, a0y_hat, a0z_hat);
		projab0 = proj(ax_hat, ay_hat, az_hat, b0x_hat, b0y_hat, b0z_hat);
		projac0 = proj(ax_hat, ay_hat, az_hat, c0x_hat, c0y_hat, c0z_hat);
		projba0 = proj(bx_hat, by_hat, bz_hat, a0x_hat, a0y_hat, a0z_hat);
		projbb0 = proj(bx_hat, by_hat, bz_hat, b0x_hat, b0y_hat, b0z_hat);
		projbc0 = proj(bx_hat, by_hat, bz_hat, c0x_hat, c0y_hat, c0z_hat);
		projca0 = proj(cx_hat, cy_hat, cz_hat, a0x_hat, a0y_hat, a0z_hat);
		projcb0 = proj(cx_hat, cy_hat, cz_hat, b0x_hat, b0y_hat, b0z_hat);
		projcc0 = proj(cx_hat, cy_hat, cz_hat, c0x_hat, c0y_hat, c0z_hat);

		// projection denominator
		denom = (-projac0 * projbb0 * projca0+projab0 * projbc0 * projca0+projac0 * projba0 * projcb0-projaa0 * projbc0 * projcb0-projab0 * projba0 * projcc0+projaa0 * projbb0 * projcc0);

		// oblique supercell vector lengths
		a_norm_REF = sqrt( (data.xhi - data.xlo)  *  (data.xhi - data.xlo) );
		b_norm_REF = sqrt( data.xy  *  data.xy + (data.yhi - data.ylo)  *  (data.yhi - data.ylo) );
		c_norm_REF = sqrt( data.xz  *  data.xz + data.yz  *  data.yz + (data.zhi - data.zlo)  *  (data.zhi - data.zlo) );

		// inverse box matrix h
		ih_11 = 1.0 / (data.xhi - data.xlo);
        ih_12 = -data.xy / ( (data.xhi - data.xlo)  *  (data.yhi - data.ylo));
        ih_13 = (data.xy  *  data.yz - (data.yhi - data.ylo)  *  data.xz) / ((data.xhi - data.xlo)  *  (data.yhi - data.ylo)  *  (data.zhi - data.zlo));
		ih_22 = 1.0 / (data.yhi - data.ylo);
        ih_23 = -data.yz / ((data.yhi - data.ylo)  *  (data.zhi - data.zlo));
		ih_33 = 1.0/(data.zhi - data.zlo);
		for(j = 0;j<atoms; ++j)
		{

			// projections
			nom1 =  -1.0  *  (projbc0 * projcb0 * data.R.x[j] - projbb0 * projcc0 * data.R.x[j] - projbc0 * projca0 * data.R.y[j] + projba0 * projcc0 * data.R.y[j] + projbb0 * projca0 * data.R.z[j] - projba0 * projcb0 * data.R.z[j]);
			nom2 = (projac0 * projcb0 * data.R.x[j] - projab0 * projcc0 * data.R.x[j] - projac0 * projca0 * data.R.y[j] + projaa0 * projcc0 * data.R.y[j] + projab0 * projca0 * data.R.z[j] - projaa0 * projcb0 * data.R.z[j]);
			nom3 =  - 1.0 * (projac0 * projbb0 * data.R.x[j] - projab0 * projbc0 * data.R.x[j] - projac0 * projba0 * data.R.y[j] + projaa0 * projbc0 * data.R.y[j] + projab0 * projba0 * data.R.z[j] - projaa0 * projbb0 * data.R.z[j]);

			// coordinates with respect to oblique system
			a = nom1/denom;
			b = nom2/denom;
			c = nom3/denom;

			//  + / -  direction
			if((a * ax_hat - 0.0) * (data.xhi - data.xlo - 0.0) + (a * ay_hat - 0.0) * (0.0 - 0.0) + (a * az_hat - 0.0) * (0.0 - 0.0) < 0.0 ){a_sign =  - 1.0;}else{a_sign =  + 1.0;}
			if((b * bx_hat - 0.0) * (data.xy - 0.0) + (b * by_hat - 0.0) * (data.yhi - data.ylo - 0.0) + (b * bz_hat - 0.0) * (0.0 - 0.0) < 0.0 ){b_sign =  - 1.0;}else{b_sign =  + 1.0;}
			if((c * cx_hat - 0.0) * (data.xz - 0.0) + (c * cy_hat - 0.0) * (data.yz - 0.0) + (c * cz_hat - 0.0) * (data.zhi - data.zlo - 0.0) < 0.0 ){c_sign =  - 1.0;}else{c_sign =  + 1.0;}

			// image calculation
			a_norm = (a * ax_hat) * (a * ax_hat) + (a * ay_hat) * (a * ay_hat) + (a * az_hat) * (a * az_hat);a_norm = sqrt(a_norm);
			b_norm = (b * bx_hat) * (b * bx_hat) + (b * by_hat) * (b * by_hat) + (b * bz_hat) * (b * bz_hat);b_norm = sqrt(b_norm);
			c_norm = (c * cx_hat) * (c * cx_hat) + (c * cy_hat) * (c * cy_hat) + (c * cz_hat) * (c * cz_hat);c_norm = sqrt(c_norm);
			data.nx = (int)floor(a_norm/a_norm_REF);if(a_sign<0.0)data.nx =  - data.nx - 1;
			data.nz = (int)floor(b_norm/b_norm_REF);if(b_sign<0.0)data.nz =  - data.nz - 1;
			nz = (int)floor(c_norm/c_norm_REF);if(c_sign<0.0)nz =  - nz - 1;

			// periodic wrapping
			data.R.x[j] = data.R.x[j] - data.nx * (data.xhi - data.xlo) - data.nz * (data.xy) - nz * (data.xz);
			data.R.y[j] = data.R.y[j] - data.nx * (0.0) - data.nz * (data.yhi - data.ylo) - nz * (data.yz);
			data.R.z[j] = data.R.z[j] - data.nx * (0.0) - data.nz * (0.0) - nz * (data.zhi - data.zlo);

			// scaled coordinates
			data.sR.x[j] = ih_11 * data.R.x[j] + ih_12 * data.R.y[j] + ih_13 * data.R.z[j];
			data.sR.y[j] = ih_22 * data.R.y[j] + ih_23 * data.R.z[j];
			data.sR.z[j] = ih_33 * data.R.z[j];
		}

        get_histogram(data, bins, hist, dr, 1);

        // rdf calculation
		for(i = 0;i<bins; ++i)
		{
			g[i] = g[i] + (hist[i] * V)/(4.0 * pi * sub_atoms_1 * sub_atoms_2 * scale * dr * rn[i] * rn[i]);
		}

		fclose(fp);
	}

    fout = fopen("./rdf.dat", "w + ");
    fputs("radius,rdf\n", fout);
    // rdf calculation
	for(i = 0;i<bins; ++i)
	{
		g[i] = g[i]/dt;
		//printf("\r%lf\t%lf\n",rn[i],g[i]);
        fprintf(fout, "%lf,%lf\n",rn[i],g[i]);
	}
    fclose(fout);

    printf("\r                                \nAll Done!");



	//------------------------------------------------------------------
	// free memory

	free(g);free(rn);free(hist);
    free(data.R.x);free(data.R.y);free(data.R.z);

    free(data.mol);

    free(type_list_1);
    free(type_list_2);
    free(data.type);
    free(belongs_to_type_list_1);
    free(belongs_to_type_list_2);

	free(data.sR.x);free(data.sR.y);free(data.sR.z);

	return 0;
}
