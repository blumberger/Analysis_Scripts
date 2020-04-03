

//

//----------------------------------------------------------------------

#define cmax_length 1000

#define species_length 4

#include<stdio.h>

#include<stdlib.h>

#include<unistd.h>

int main(int argc, char *argv[])

{

	FILE *fp,*xyz_fp,*xyz_unwrapped_fp;

	

	// fp related

	char current_folder[cmax_length],file_path[cmax_length];

	char folder[cmax_length],name[cmax_length];

	

	// dummies

	char buffer[cmax_length];

	char species[species_length];

	int int_buffer;

	double d_buffer;

	

	//

	int i,j,k;

	int start,step,stop;

	int atom_types,atoms,atom_ID;	

	char title[cmax_length];

	char **label;

	int *ID,*atomicID;

	double x,y,z;

	

	// unwrap

	int nx,ny,nz;

	double ax,ay,az;

	double bx,by,bz;

	double cx,cy,cz;

	

	//

	int flag;

	double d_temp;

	double xlo_bound,xhi_bound,ylo_bound,yhi_bound,zlo_bound,zhi_bound,xy,xz,yz;

	double S[4],xlo,xhi,ylo,yhi,zlo,zhi;

	double MIN_x,MAX_x,MIN_y,MAX_y;

		

	//------------------------------------------------------------------

	// print exe info

	if(argc<6)

	{

		printf(

		"\n*Reads lammps dump files ./<folder>/<filename>#.dat and creates files suitable for Jmol animation\n"

		"*Can read the species information through file pipelining\n"

		"*Supports custom dump files with the format id mol type x y z ix iy iz\n"

		"*Generates:\n"

		"\t- xyz_animation.dat (coordinates with periodic wrapping)\n"

		"\t -xyz_animation_unwrapped.dat (unwrapped coordinates)\n\n"

		"./animate <folder> <filename> <start> <step> <stop>\n\n"

		);

		exit(-1);

	}

	

	// read arguments

	sprintf(folder,"%s",argv[1]);

	sprintf(name,"%s",argv[2]);

	start=atoi(argv[3]);

	step=atoi(argv[4]);

	stop=atoi(argv[5]);

    sprintf(folder, ".");
    sprintf(name, "snapshot");
    start = 0;
    stop = 0;
    step = 1;

	

	// current folder

	getcwd(current_folder,cmax_length);

	

	// Jmol title

	sprintf(title,"%s","xyz_animation");

	

	// read the first dump file

	sprintf(file_path,"%s/%s/%s%d.dat",current_folder,folder,name,start);

	fp=fopen(file_path,"r");

	if(fp==NULL){printf("Cannot locate %s\n",file_path);exit(-1);}		// file check

	for(i=0;i<3;++i){fgets(buffer,cmax_length,fp);}

	fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&atoms);

	fgets(buffer,cmax_length,fp);

	fgets(buffer,cmax_length,fp);

	fgets(buffer,cmax_length,fp);

	fgets(buffer,cmax_length,fp);

	fgets(buffer,cmax_length,fp);


	atomicID=(int*)malloc(atoms*sizeof(int));

	for(i=0;i<atoms;++i){fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d",&int_buffer,&int_buffer,&atomicID[i]);}

	fclose(fp);

	

	// species check

	atom_types=1;

	for(i=0;i<atoms;++i){if(atomicID[i]>atom_types){atom_types=atomicID[i];}}

	ID=(int*)malloc(atom_types*sizeof(int));

	label=(char**)malloc(atom_types*sizeof(char*));// species

	for (i=0;i<atom_types;++i){label[i]=(char*)malloc(species_length*sizeof(char));}

	// ask for species labels

	printf("Found %d atomic species. Identify:\n",atom_types);

	for(i=0;i<atom_types;++i){printf("[%d]: ",i+1);ID[i]=i+1;scanf("%s",label[i]);}

	

	// create Jmol files

	sprintf(file_path,"%s/xyz_animation.dat",current_folder);

	xyz_fp=fopen(file_path,"w+");

	

	sprintf(file_path,"%s/xyz_animation_unwrapped.dat",current_folder);

	xyz_unwrapped_fp=fopen(file_path,"w+");

	

	// start reading the dump files

	for(k=start;k<=stop;k=k+step)

	{

		// open dump file

		sprintf(file_path,"%s/%s/%s%d.dat",current_folder,folder,name,k);

		fp=fopen(file_path,"r");

		if(fp==NULL){printf("Cannot locate %s\n",file_path);exit(-1);}		// file check

		

		for(i=0;i<3;++i){fgets(buffer,cmax_length,fp);}

		fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&atoms);

		fgets(buffer,cmax_length,fp);

		fgets(buffer,cmax_length,fp);sscanf(buffer,"%lf\t%lf\t%lf",&xlo_bound,&xhi_bound,&xy);

		fgets(buffer,cmax_length,fp);sscanf(buffer,"%lf\t%lf\t%lf",&ylo_bound,&yhi_bound,&xz);

		fgets(buffer,cmax_length,fp);sscanf(buffer,"%lf\t%lf\t%lf",&zlo_bound,&zhi_bound,&yz);

		

		//

		

		S[0]=0.0;

		S[1]=xy;

		S[2]=xz;

		S[3]=xy+xz;

		

		// bubble sort

		flag=1;

		while(flag==1)

		{

			flag=0;

			for(i=1;i<4;++i)

			{

				if(S[i]<S[i-1])

				{


					d_temp=S[i-1];S[i-1]=S[i];S[i]=d_temp;

									

					flag=1;

				}

			}	

		}
        printf("\nBOB %d\n", k);

		MIN_x=S[0];MAX_x=S[3];

		

		if(yz<0.0){MIN_y=yz;}else{MIN_y=0.0;}

		if(yz>0.0){MAX_y=yz;}else{MAX_y=0.0;}

		

		xlo=xlo_bound-MIN_x;

		xhi=xhi_bound-MAX_x;

		ylo=ylo_bound-MIN_y;

		yhi=yhi_bound-MAX_y;

		zlo=zlo_bound;

		zhi=zhi_bound;

		

        printf("\nxlo, xhi: %f   %f\n", xlo, xhi);
        printf("\nylo, yhi: %f   %f\n", ylo, yhi);
        printf("\nzlo, zhi: %f   %f\n", zlo, zhi);
		

		ax=xhi-xlo;

		ay=ylo-ylo;

		az=zlo-zlo;

		

		bx=xy+xlo-xlo;

		by=yhi-ylo;

		bz=zlo-zlo;

		

		cx=xz+xlo-xlo;

		cy=yz+ylo-ylo;

		cz=zhi-zlo;

        printf("\na = (%f   %f   %f)\n", ax, ay, az);
        printf("\nb = (%f   %f   %f)\n", bx, by, bz);
        printf("\nc = (%f   %f   %f)\n", cx, cy, cz);
		

		//

		

		fgets(buffer,cmax_length,fp);

		

		// append in animation.dat and animation_unwrapped.dat

		

		//

		fprintf(xyz_fp,"%d\n%s\n",atoms+4,title);

		

		fprintf(xyz_fp,"%s\t%lf\t%lf\t%lf\n","point",xlo,ylo,zlo);

		fprintf(xyz_fp,"%s\t%lf\t%lf\t%lf\n","point",xhi,ylo,zlo);

		fprintf(xyz_fp,"%s\t%lf\t%lf\t%lf\n","point",xy+xlo,yhi,zlo);

		fprintf(xyz_fp,"%s\t%lf\t%lf\t%lf\n","point",xz+xlo,yz+ylo,zhi);

		

		fprintf(xyz_unwrapped_fp,"%d\n%s\n",atoms,title);

		for(i=0;i<atoms;++i)

		{

			// read

			fgets(buffer,cmax_length,fp);

			sscanf(buffer,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%d\t%d\t%d",&int_buffer,&int_buffer,&atom_ID,&x,&y,&z,&nx,&ny,&nz);

			// figure out the species

			for(j=0;j<atom_types;++j){if(atom_ID==ID[j]){sprintf(species,"%s",label[j]);}}

			// write

			fprintf(xyz_fp,"%s\t%lf\t%lf\t%lf\n",species,x,y,z);

			fprintf(xyz_unwrapped_fp,"%s\t%lf\t%lf\t%lf\n",species,x+nx*ax+ny*bx+nz*cx,y+nx*ay+ny*by+nz*cy,z+nx*az+ny*bz+nz*cz);

		}

	

		fclose(fp);

	}

	

	fclose(xyz_fp);	

	fclose(xyz_unwrapped_fp);	

	printf("\nGenerated %s/xyz_animation.dat\n",current_folder);	

	printf("\nGenerated %s/xyz_animation_unwrapped.dat\n",current_folder);

	

	//------------------------------------------------------------------

	// free memory

	

	for (i=0;i<atom_types;++i){free(label[i]);}free(label);

	free(ID);

	free(atomicID);

	

	return 0;

}
