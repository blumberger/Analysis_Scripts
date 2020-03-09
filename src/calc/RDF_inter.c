//
//----------------------------------------------------------------------
#define cmax_length 1000
#define species_length 4
#define pi 3.14159265
#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include<unistd.h>
#include<math.h>
double proj(double Ax,double Ay,double Az,double Bx,double By,double Bz);
double proj(double Ax,double Ay,double Az,double Bx,double By,double Bz)
{
	double dot,norm,res;
	dot=Ax*Bx+Ay*By+Az*Bz;
	norm=Ax*Ax+Ay*Ay+Az*Az;
	res=dot/norm;
	return res;
}
int main(int argc, char *argv[])
{
	FILE *fp;

	// fp related
	char current_folder[cmax_length],file_path[cmax_length];
	char folder[cmax_length],name[cmax_length];

	// dummies
	char buffer[cmax_length];
	int int_buffer;

	//
	int i,j,k;
	int start,step,stop;
	int atoms;
	double *x,*y,*z;

	// unwrap
	int nx,ny,nz;

	//

    //
    int dt;
	int bins;
	int *hist;
	double re,dr;
	double *g,*rn;
    double r;

    int *mol;

    int type_N_1,type_N_2,sub_atoms_1=0,sub_atoms_2=0;
    int *type_list_1,*type_list_2,*type,*belongs_to_type_list_1,*belongs_to_type_list_2;
    int scale;


	//
	double V;

	//
	int flag;
	double d_temp;
	double xlo_bound,xhi_bound,ylo_bound,yhi_bound,zlo_bound,zhi_bound,xy,xz,yz;
	double S[4],xlo,xhi,ylo,yhi,zlo,zhi;
	double MIN_x,MAX_x,MIN_y,MAX_y;

	double ax,ay,az,bx,by,bz,cx,cy,cz,norm;
	double ax_hat,ay_hat,az_hat,bx_hat,by_hat,bz_hat,cx_hat,cy_hat,cz_hat;
	double a0x_hat=1.0,a0y_hat=0.0,a0z_hat=0.0;
	double b0x_hat=0.0,b0y_hat=1.0,b0z_hat=0.0;
	double c0x_hat=0.0,c0y_hat=0.0,c0z_hat=1.0;

	double projaa0;
	double projab0;
	double projac0;
	double projba0;
	double projbb0;
	double projbc0;
	double projca0;
	double projcb0;
	double projcc0;
	double nom1,nom2,nom3,denom;
	double a,b,c,a_sign,b_sign,c_sign,a_norm,b_norm,c_norm;
	double a_norm_REF,b_norm_REF,c_norm_REF;
	//int nx,ny,nz;

	double ih_11,ih_12,ih_13,ih_22,ih_23,ih_33,*sx,*sy,*sz;
	double svectorx,svectory,svectorz,rx,ry,rz;

	int bh;//,bins,*hist;
	//double dr,re,*g,*rn;


	//------------------------------------------------------------------

    if(argc==1)
    {
        printf("\n./RDF2 <folder> <name> <start> <step> <stop> <dr> <re> <type_list>\n\n");exit(-1);
    }

	// read arguments
	sprintf(folder,"%s",argv[1]);
	sprintf(name,"%s",argv[2]);
	start=atoi(argv[3]);
	step=atoi(argv[4]);
	stop=atoi(argv[5]);

    //
    dr=atof(argv[6]);
    re=atof(argv[7]);
    //

    dt=(stop-start)/step+1;

    // rdf initialization
	// number of bins
	bins=floor(re/dr+0.5);
	// bin count array
	hist=(int*)malloc(bins*sizeof(int));
	// RDF
	g=(double*)malloc(bins*sizeof(double));
	// shell radius
	rn=(double*)malloc(bins*sizeof(double));
	for (i=0;i<bins;++i){rn[i]=((i+1.0)-0.5)*dr;g[i]=0.0;}

	// current folder
	getcwd(current_folder,cmax_length);

	// read type config
	sprintf(file_path,"%s/%s",current_folder,argv[8]);
	fp=fopen(file_path,"r");
	if(fp==NULL){printf("Cannot locate %s\n",file_path);exit(-1);}		// file check
	// locate group 1
	while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"type_1\n")==0)break;
	fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&type_N_1);
	type_list_1=(int*)malloc(type_N_1*sizeof(int));
	for(i=0;i<type_N_1;++i)
	{
		fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&type_list_1[i]);
	}
	rewind(fp);
	// locate group 2
	while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"type_2\n")==0)break;
	fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&type_N_2);
	type_list_2=(int*)malloc(type_N_2*sizeof(int));
	for(i=0;i<type_N_2;++i)
	{
		fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&type_list_2[i]);
	}
	fclose(fp);

	// scale
	scale=atoi(argv[9]);

	// read the first dump file and save atoms
	sprintf(file_path,"%s/%s/%s%d.dat",current_folder,folder,name,start);
	fp=fopen(file_path,"r");
	if(fp==NULL){printf("Cannot locate %s\n",file_path);exit(-1);}		// file check
	for(i=0;i<3;++i){fgets(buffer,cmax_length,fp);}
	fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&atoms);
	fclose(fp);

	// preallocations
    x=(double*)malloc(atoms*sizeof(double));
    y=(double*)malloc(atoms*sizeof(double));
    z=(double*)malloc(atoms*sizeof(double));
    sx=(double*)malloc(atoms*sizeof(double));
    sy=(double*)malloc(atoms*sizeof(double));
    sz=(double*)malloc(atoms*sizeof(double));

    mol=(int*)malloc(atoms*sizeof(int));
    type=(int*)malloc(atoms*sizeof(int));
    belongs_to_type_list_1=(int*)malloc(atoms*sizeof(int));for(i=0;i<atoms;++i)belongs_to_type_list_1[i]=0;
    belongs_to_type_list_2=(int*)malloc(atoms*sizeof(int));for(i=0;i<atoms;++i)belongs_to_type_list_2[i]=0;

	// start reading the dump files
	for(k=start;k<=stop;k=k+step)
	{
        for (i=0;i<bins;++i){hist[i]=0;}
		// open dump file
		sprintf(file_path,"%s/%s/%s%d.dat",current_folder,folder,name,k);
		fp=fopen(file_path,"r");
		if(fp==NULL){printf("Cannot locate %s\n",file_path);exit(-1);}		// file check

		for(i=0;i<3;++i){fgets(buffer,cmax_length,fp);}
		fgets(buffer,cmax_length,fp);//sscanf(buffer,"%d",&atoms);
		fgets(buffer,cmax_length,fp);
		fgets(buffer,cmax_length,fp);sscanf(buffer,"%lf\t%lf\t%lf",&xlo_bound,&xhi_bound,&xy);
		fgets(buffer,cmax_length,fp);sscanf(buffer,"%lf\t%lf\t%lf",&ylo_bound,&yhi_bound,&xz);
		fgets(buffer,cmax_length,fp);sscanf(buffer,"%lf\t%lf\t%lf",&zlo_bound,&zhi_bound,&yz);
		fgets(buffer,cmax_length,fp);
		// supecell

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
			for(j=1;j<4;++j)
			{
				if(S[j]<S[j-1])
				{

					d_temp=S[j-1];S[j-1]=S[j];S[j]=d_temp;

					flag=1;
				}
			}
		}
		MIN_x=S[0];MAX_x=S[3];

		if(yz<0.0){MIN_y=yz;}else{MIN_y=0.0;}
		if(yz>0.0){MAX_y=yz;}else{MAX_y=0.0;}

		xlo=xlo_bound-MIN_x;
		xhi=xhi_bound-MAX_x;
		ylo=ylo_bound-MIN_y;
		yhi=yhi_bound-MAX_y;
		zlo=zlo_bound;
		zhi=zhi_bound;

		//



		// read atoms
		for(i=0;i<atoms;++i)
		{
			// read
			fgets(buffer,cmax_length,fp);

      sscanf(buffer,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%d\t%d\t%d",&int_buffer,&mol[i],&type[i],&x[i],&y[i],&z[i],&nx,&ny,&nz);


		}

		// check once
		if(k==start)
		{
			for(i=0;i<atoms;++i)
			{
				for(j=0;j<type_N_1;++j)
				{
					if(type[i]==type_list_1[j]){belongs_to_type_list_1[i]=1;++sub_atoms_1;break;}
				}
				for(j=0;j<type_N_2;++j)
				{
					if(type[i]==type_list_2[j]){belongs_to_type_list_2[i]=1;++sub_atoms_2;break;}
				}
			}
		}

		// shift
		for(j=0;j<atoms;++j)
		{
			x[j]=x[j]-xlo;y[j]=y[j]-ylo;z[j]=z[j]-zlo;
		}

		// oblique supercell vectors
		ax=(xhi-xlo);// - (0.0);
		ay=(0.0)    ;// - (0.0);
		az=(0.0)    ;// - (0.0);
		bx=(xy)     ;// - (0.0);
		by=(yhi-ylo);// - (0.0);
		bz=(0.0)    ;// - (0.0);
		cx=(xz)     ;// - (0.0);
		cy=(yz)     ;// - (0.0);
		cz=(zhi-zlo);// - (0.0);

		V=ax*by*cz;

		// normalize
		norm=ax*ax+ay*ay+az*az;norm=sqrt(norm);
		ax_hat=ax/norm;
		ay_hat=ay/norm;
		az_hat=az/norm;
		norm=bx*bx+by*by+bz*bz;norm=sqrt(norm);
		bx_hat=bx/norm;
		by_hat=by/norm;
		bz_hat=bz/norm;
		norm=cx*cx+cy*cy+cz*cz;norm=sqrt(norm);
		cx_hat=cx/norm;
		cy_hat=cy/norm;
		cz_hat=cz/norm;

		// calculate projections
		projaa0=proj(ax_hat,ay_hat,az_hat,a0x_hat,a0y_hat,a0z_hat);
		projab0=proj(ax_hat,ay_hat,az_hat,b0x_hat,b0y_hat,b0z_hat);
		projac0=proj(ax_hat,ay_hat,az_hat,c0x_hat,c0y_hat,c0z_hat);
		projba0=proj(bx_hat,by_hat,bz_hat,a0x_hat,a0y_hat,a0z_hat);
		projbb0=proj(bx_hat,by_hat,bz_hat,b0x_hat,b0y_hat,b0z_hat);
		projbc0=proj(bx_hat,by_hat,bz_hat,c0x_hat,c0y_hat,c0z_hat);
		projca0=proj(cx_hat,cy_hat,cz_hat,a0x_hat,a0y_hat,a0z_hat);
		projcb0=proj(cx_hat,cy_hat,cz_hat,b0x_hat,b0y_hat,b0z_hat);
		projcc0=proj(cx_hat,cy_hat,cz_hat,c0x_hat,c0y_hat,c0z_hat);

		// projection denominator
		denom=(-projac0*projbb0*projca0+projab0*projbc0*projca0+projac0*projba0*projcb0-projaa0*projbc0*projcb0-projab0*projba0*projcc0+projaa0*projbb0*projcc0);

		// oblique supercell vector lengths
		a_norm_REF=(xhi-xlo)*(xhi-xlo);a_norm_REF=sqrt(a_norm_REF);
		b_norm_REF=(xy)*(xy)+(yhi-ylo)*(yhi-ylo);b_norm_REF=sqrt(b_norm_REF);
		c_norm_REF=(xz)*(xz)+(yz)*(yz)+(zhi-zlo)*(zhi-zlo);c_norm_REF=sqrt(c_norm_REF);

		// inverse box matrix h
		ih_11=1.0/(xhi-xlo);ih_12=-xy/((xhi-xlo)*(yhi-ylo));ih_13=(xy*yz-(yhi-ylo)*xz)/((xhi-xlo)*(yhi-ylo)*(zhi-zlo));
		ih_22=1.0/(yhi-ylo);ih_23=-yz/((yhi-ylo)*(zhi-zlo));
		ih_33=1.0/(zhi-zlo);
		for(j=0;j<atoms;++j)
		{

			// projections
			nom1=-1.0*(projbc0*projcb0*x[j]-projbb0*projcc0*x[j]-projbc0*projca0*y[j]+projba0*projcc0*y[j]+projbb0*projca0*z[j]-projba0*projcb0*z[j]);
			nom2=(projac0*projcb0*x[j]-projab0*projcc0*x[j]-projac0*projca0*y[j]+projaa0*projcc0*y[j]+projab0*projca0*z[j]-projaa0*projcb0*z[j]);
			nom3=-1.0*(projac0*projbb0*x[j]-projab0*projbc0*x[j]-projac0*projba0*y[j]+projaa0*projbc0*y[j]+projab0*projba0*z[j]-projaa0*projbb0*z[j]);

			// coordinates with respect to oblique system
			a=nom1/denom;
			b=nom2/denom;
			c=nom3/denom;

			// +/- direction
			if((a*ax_hat-0.0)*(xhi-xlo-0.0)+(a*ay_hat-0.0)*(0.0-0.0)+(a*az_hat-0.0)*(0.0-0.0) < 0.0 ){a_sign=-1.0;}else{a_sign=+1.0;}
			if((b*bx_hat-0.0)*(xy-0.0)+(b*by_hat-0.0)*(yhi-ylo-0.0)+(b*bz_hat-0.0)*(0.0-0.0) < 0.0 ){b_sign=-1.0;}else{b_sign=+1.0;}
			if((c*cx_hat-0.0)*(xz-0.0)+(c*cy_hat-0.0)*(yz-0.0)+(c*cz_hat-0.0)*(zhi-zlo-0.0) < 0.0 ){c_sign=-1.0;}else{c_sign=+1.0;}

			// image calculation
			a_norm=(a*ax_hat)*(a*ax_hat)+(a*ay_hat)*(a*ay_hat)+(a*az_hat)*(a*az_hat);a_norm=sqrt(a_norm);
			b_norm=(b*bx_hat)*(b*bx_hat)+(b*by_hat)*(b*by_hat)+(b*bz_hat)*(b*bz_hat);b_norm=sqrt(b_norm);
			c_norm=(c*cx_hat)*(c*cx_hat)+(c*cy_hat)*(c*cy_hat)+(c*cz_hat)*(c*cz_hat);c_norm=sqrt(c_norm);
			nx=(int)floor(a_norm/a_norm_REF);if(a_sign<0.0)nx=-nx-1;
			ny=(int)floor(b_norm/b_norm_REF);if(b_sign<0.0)ny=-ny-1;
			nz=(int)floor(c_norm/c_norm_REF);if(c_sign<0.0)nz=-nz-1;

			// periodic wrapping
			x[j]=x[j]-nx*(xhi-xlo)-ny*(xy)-nz*(xz);
			y[j]=y[j]-nx*(0.0)-ny*(yhi-ylo)-nz*(yz);
			z[j]=z[j]-nx*(0.0)-ny*(0.0)-nz*(zhi-zlo);

			// scaled coordinates
			sx[j]=ih_11*x[j]+ih_12*y[j]+ih_13*z[j];
			sy[j]=           ih_22*y[j]+ih_23*z[j];
			sz[j]=                      ih_33*z[j];
			//if(i==frames-1)printf("%lf\t%lf\t%lf\n",sx[j],sy[j],sz[j]);
			//if(i==frames-1)printf("%lf\t%lf\t%lf\n",rint(sx[j]),rint(sy[j]),rint(sz[j]));
		}

		for(j=0;j<atoms-1;++j)
		{
			for(i=j+1;i<atoms;++i)
			{

				if(mol[i]!=mol[j] && ((belongs_to_type_list_1[i]==1 && belongs_to_type_list_2[j]==1) || (belongs_to_type_list_1[j]==1 && belongs_to_type_list_2[i]==1)))
				{

					//printf("%d-%d\n",j+1,i+1);
					svectorx=sx[j]-sx[i];
					svectory=sy[j]-sy[i];
					svectorz=sz[j]-sz[i];

					svectorx=svectorx-rint(svectorx);
					svectory=svectory-rint(svectory);
					svectorz=svectorz-rint(svectorz);

					rx=(xhi-xlo)*svectorx+(xy)*svectory+(xz)*svectorz;
					ry=(yhi-ylo)*svectory+(yz)*svectorz;
					rz=(zhi-zlo)*svectorz;

					r=rx*rx+ry*ry+rz*rz;r=sqrt(r);

					//printf("%lf\n",r);

					bh=floor(r/dr);
					if (bh<bins)
					{
						hist[bh]=hist[bh]+2;
					}
//					else
//					{
//						printf("%lf > %lf !!!\nAborting...\n\n",r,re);exit(-1);
//					}

				}

			}
		}




        // rdf calculation
		for(i=0;i<bins;++i)
		{
			g[i]=g[i]+(hist[i]*V)/(4.0*pi*sub_atoms_1*sub_atoms_2*scale*dr*rn[i]*rn[i]);
		}

		fclose(fp);
	}

    // rdf calculation
	for(i=0;i<bins;++i)
	{
		g[i]=g[i]/dt;
		printf("%lf\t%lf\n",rn[i],g[i]);
	}



	//------------------------------------------------------------------
	// free memory

	free(g);free(rn);free(hist);
    free(x);free(y);free(z);

    free(mol);

    free(type_list_1);
    free(type_list_2);
    free(type);
    free(belongs_to_type_list_1);
    free(belongs_to_type_list_2);

	free(sx);free(sy);free(sz);

	return 0;
}
