// This code was originated from R package `monoreg`(https://cran.r-project.org/package=monoreg)
#include <stdio.h>
#include <math.h>
#include <stdlib.h> 

void test_code(void){
	printf("-----------------------------\n");
	printf("this is test code in monreg.c\n");
	printf("-----------------------------\n");
}

double fakesum(double x, double y){
	return 123;
}

double max(double x, double y){
	if(x > y) return x;
	else return y;
}

double min(double x, double y){
	if(x < y) return x;
	else return y;
}

double faculty(int n){
	int nfac;
	int i;
	nfac=1;
	for(i=1; i<n; i++){
		nfac=nfac*(i+1);
	}
	return((double) nfac);
}

double gewicht(double *d,int i,int r){
	for(i=0;i<=r;i++){
		d[i] = pow((double) -1,(double) i) * ((faculty(r) / (faculty(i) * faculty(r-i))) / (pow(faculty(2 * r) / pow(faculty(r),2.0),0.5)));
	}
	return(*d);
}

void nadwats(double *x, double *y, double *z, int m, int n, double h, int K, double *p){
	long i,j;
	double zaehl,nenn,arg,kfunkt;


	for(i=0; i<n; ++i){
	    zaehl=0.0;
		nenn=0.0;
		for(j=0; j<m; ++j){
			arg = (x[j] - z[i]) / h ;
			if(fabs(arg) <= 1.0){
				switch (K){
					case 1:
						kfunkt = 3.0*(1.0 - pow(arg,2.0))/4.0;
						break;
					case 2:
						kfunkt=1.0/2.0;
						break;
					case 3:
						kfunkt=15.0/16.0*pow((1-pow(arg,2.0)),2.0);
						break;
					case 4:
						kfunkt=35.0/32.0*pow((1-pow(arg,2.0)),3.0);
						break;
					case 5:
						kfunkt=1.0-fabs(arg);
						break;
					case 6:
						kfunkt=3.1415/4.0*cos(3.1415/2.0*arg);
						break;
				}
			}
			else
				kfunkt=0.0;

			zaehl = zaehl + kfunkt * y[j];
			nenn = nenn + kfunkt;
		}
		if(nenn == 0.0)
			p[i] = 0.0;
		else
			p[i] = zaehl / nenn;
	}

}


void loclinear(double *x, double *y, double *z, int m, int n, double h, int K, double *p){
	long i,j;
	double zaehl,nenn,arg,s0,s1,s2,t0,t1;
	double *kfunkt = malloc(m * sizeof(double));



	for(i=0; i < n; ++i){
		s0 = 0.0;
		s1 = 0.0;
		s2 = 0.0;
		t0 = 0.0;
		t1 = 0.0;
		for(j=0; j < m; ++j){
			arg = (x[j] - z[i]) / h;
        if(fabs(arg) <= 1.0){
            switch (K){
				case 1:
					kfunkt[j] = 3.0*(1.0 - pow(arg,2.0))/4.0;
					break;
				case 2:
					kfunkt[j] = 1.0/2.0;
					break;
				case 3:
					kfunkt[j] = 15.0/16.0*pow((1-pow(arg,2.0)),2.0);
					break;
				case 4:
					kfunkt[j] = 35.0/32.0*pow((1-pow(arg,2.0)),3.0);
					break;
				case 5:
					kfunkt[j] = 1.0-fabs(arg);
					break;
				case 6:
					kfunkt[j] = 3.1415/4.0*cos(3.1415/2.0*arg);
					break;
			}
        }
        else
			kfunkt[j] = 0.0;

		s0 = s0 + kfunkt[j];
        s1 = s1 + (x[j] - z[i]) * kfunkt[j];
        s2 = s2 + (x[j] - z[i]) * (x[j] - z[i]) * kfunkt[j];
        t0 = t0 + kfunkt[j] * y[j];
        t1 = t1 + (x[j] - z[i]) * kfunkt[j] * y[j];
		}
		zaehl = (s2 * t0) - (s1 * t1);
		nenn = (s2 * s0) - (s1 * s1);
		if(nenn == 0.0){	
           p[i] = zaehl / (((m - 1.0) * h * h) / (5.0 * m));
		}
	    else{
			p[i] = zaehl / nenn;
		}
	}
    free(kfunkt);

}


double mdach_i(double t,double *t0, double *minv, int k)
{
	int j;
	double d1,mdachi;

	d1 = fabs(t - minv[0]);
	for(j=1; j<k; j++)
	{
		if(d1 < fabs(t - minv[j]))
			break;
		else
		{
			d1 = fabs(t - minv[j]);
		}
	}
	mdachi = t0[j-1];
	return(mdachi);
}



void mdach_i_inv (double *x, double *y, double *z, double *t,int *m, int *n, int *l, int *tflag, double *hd, double *hr, int *Kd, int *Kr, double *degree, int *ldeg, int *inverse, double *erg){
	long i,j,k;
	double sum,arg,integral,dist,mini,maxi;
	double *minv,*t0;
	double *p = malloc(*n * sizeof(double));



	if(*ldeg == 1){
		if(*degree == 0.0){
			nadwats(x,y,z,*m,*n,*hr,*Kr,p);
		}
		else
			loclinear(x,y,z,*m,*n,*hr,*Kr,p);
	}
	else{
		for(j=0; j<*n; j++){
			p[j] = degree[j];
		}
	}

	if(*inverse == 0){
		t0 = malloc(5001 * sizeof(double));
		mini = p[0];
		maxi = p[0];
		for(i=1; i<*n; i++){
			if(p[i] < mini)
				mini = p[i];
			else if(p[i] > maxi)
				maxi = p[i];
		}
		t0[0] = mini - *hd;
		dist = (maxi - mini + 2.0 * *hd) / ((double) (5000));
		for(i=1; i<5001; i++){
			t0[i]=t0[i-1]+dist;            
		}
		k = 5001;
	}
	else if(*inverse == 1)
		k = *l;


	if(*tflag == 1 && *inverse == 0){
		mini = x[0];
		maxi = x[0];
		for(i=1; i<*m; i++){
			if(x[i] < mini)
				mini = x[i];
			else if(x[i] > maxi)
				maxi = x[i];
		}
		t[0] = mini;
		dist = (maxi - mini) / ((double) (*l-1));
		for(i=1; i<*l; i++){
			t[i] = t[i-1] + dist;
		}
	}

	if(*tflag == 1 && *inverse == 1){
		mini = p[0];
		maxi = p[0];
		for(i=1; i<*n; i++){
			if(p[i] < mini)
				mini = p[i];
			else if(p[i] > maxi)
				maxi = p[i];
		}
		t[0] = mini;
		dist = (maxi - mini) / ((double) (*l - 1));
		for(i=1; i<*l; i++){
			t[i] = t[i-1] + dist;            
		}
		t0 = malloc(*l * sizeof(double));
		t0 = t;
	}


	minv = malloc(k * sizeof(double));

	for(j=0; j<k; j++){
		sum=0.0;
		for(i=0; i<*n; i++){
			arg = (p[i] - t0[j]) / *hd;
			if(arg > 1.0){
				integral=0.0;
			}
			else{
				switch (*Kd){
					case 1:
						integral=*hd * (3.0 / 4.0 * ((1.0 - max(-1.0,arg))-1.0/3.0*(1.0-pow(max(-1.0,arg),3.0))));
						break;
					case 2:
						integral=1.0 / 2.0 * *hd * (1.0 - max(arg,-1.0));
						break;
					case 3:
						integral= *hd * (1.0/16.0*(15.0*(1.0-max(arg,-1.0))-10.0*(1.0-pow(max(arg,-1.0),3.0))+3.0*(1.0-pow(max(arg,-1.0),5.0))));
						break;
					case 4:
						integral= *hd * (1.0/32.0*(35.0*(1.0-max(arg,-1.0))-35.0*(1.0-pow(max(arg,-1.0),3.0))+21.0*(1.0-pow(max(arg,-1.0),5.0))-5.0*(1.0-pow(max(arg,-1.0),7.0))));
						break;
					case 5:
						integral= *hd * (1.0/2.0-max(arg,-1.0)+1.0/2.0*max(arg,-1.0)*fabs(max(arg,-1.0)));
						break;
					case 6:
						integral= *hd * (1.0/2.0*(1-sin(max(arg,-1)*3.1415/2.0)));
						break;
					case 0:
						integral=0.0;
					break;
				}
			}

			sum = sum + integral;
		}
		minv[j] = 1.0 / ((double) *n * *hd) * sum;
	}

	if(*inverse == 0){
		for(j=0; j<*l; j++){
			erg[j] = mdach_i(t[j],t0,minv,k);
			
		}
	}
	else{
		for(j=0; j<k; j++){
			erg[j] = minv[j];
		}
	}

}
			 



void mdach_a_inv (double *x, double *y, double *z, double *t,int *m, int *n, int *l, int *tflag, double *hd, double *hr, int *Kd, int *Kr, double *degree, int *ldeg, int *inverse, double *erg){
	long i,j,k;
	double sum,arg,integral,dist,mini,maxi;
	double *minv,*t0;
	double *p = malloc(*n * sizeof(double));


	if(*ldeg == 1){
		if(*degree == 0.0){
			nadwats(x,y,z,*m,*n,*hr,*Kr,p);
		}
		else
			loclinear(x,y,z,*m,*n,*hr,*Kr,p);
	}
	else{
		for(j=0; j<*n; j++){
			p[j] = degree[j];
		}
	}


	if(*inverse == 0){
		t0 = malloc(5001 * sizeof(double));
		mini = p[0];
		maxi = p[0];
		for(i=1; i<*n; i++){
			if(p[i] < mini)
				mini = p[i];
			else if(p[i] > maxi)
				maxi = p[i];
		}
		t0[0] = mini - *hd;
		dist = (maxi - mini + 2.0 * *hd) / ((double) (5000));
		for(i=1; i<5001; i++){
			t0[i]=t0[i-1]+dist;            
		}
		k = 5001;
	}
	else if(*inverse == 1)
		k = *l;



	if(*tflag == 1 && *inverse == 0){
		mini = x[0];
		maxi = x[0];
		for(i=1; i<*m; i++){
			if(x[i] < mini)
				mini = x[i];
			else if(x[i] > maxi)
				maxi = x[i];
		}
		t[0] = mini;
		dist = (maxi - mini) / ((double) (*l-1));
		for(i=1; i<*l; i++){
			t[i] = t[i-1] + dist;
		}
	}

	if(*tflag == 1 && *inverse == 1){
		mini = p[0];
		maxi = p[0];
		for(i=1; i<*n; i++){
			if(p[i] < mini)
				mini = p[i];
			else if(p[i] > maxi)
				maxi = p[i];
		}
		t[0] = mini;
		dist = (maxi - mini) / ((double) (*l - 1));
		for(i=1; i<*l; i++){
			t[i] = t[i-1] + dist;            
		}
		t0 = malloc(*l * sizeof(double));
		t0 = t;
	}


	minv = malloc(k * sizeof(double));

	for(j=0; j<k; j++){
		sum=0.0;
		for(i=0; i<*n; i++){
			arg = (p[i] - t0[j]) / *hd;
			if(arg <= -1.0){
				integral=0.0;
			}
			else{
				switch (*Kd){
					case 1:
						integral=*hd * (1.0 / 4.0 * ((2.0 + 3.0*min(1.0,arg))-pow(min(1.0,arg),3.0)));
						break;
					case 2:
						integral=1.0 / 2.0 * *hd * (min(arg,1.0) + 1.0);
						break;
					case 3:
						integral= *hd * (1.0/16.0*(15.0*((min(arg,1.0))+1.0)-10.0*(pow(min(arg,1.0),3.0)+1.0)+3.0*(pow(min(arg,1.0),5.0)+1.0)));
						break;
					case 4:
						integral= *hd * (1.0/32.0*(35.0*(min(arg,1.0)+1.0)-35.0*(pow(min(arg,1.0),3.0)+1.0)+21.0*(pow(min(arg,1.0),5.0)+1.0)-5.0*(pow(min(arg,1.0),7.0)+1.0)));
						break;
					case 5:
						integral= *hd * ((-1.0)/2.0+min(arg,1.0)-1.0/2.0*min(arg,1.0)*fabs(min(arg,1.0)));
						break;
					case 6:
						integral= *hd * (1.0/2.0*(1.0+sin(min(arg,1.0)*3.1415/2.0)));
						break;
					case 0:
						integral=0.0;
					break;
				}
			}

			sum = sum + integral;
		}
		minv[j] = 1.0 / ((double) *n * *hd) * sum;
	}

	if(*inverse == 0){
		for(j=0; j<*l; j++){
			erg[j] = mdach_i(t[j],t0,minv,k);
					}
	}
	else{
		for(j=0; j<k; j++){
			erg[j] = minv[j];
		}
	}
}




void sdach_i_inv (double *x, double *y, double *z, double *t,int *m, int *n, int *l, int *tflag, double *h, double *hd, double *hr,int *K, int *Kd, int *Kr, double *mdegree, int *lmdeg, double *sdegree, int *lsdeg, int *inverse, double *erg){

    long i,j,k;
	double sum,arg,integral,dist,mini,maxi;
	double *sinv,*t0;
	double *p = malloc(*m * sizeof(double));
    double *e = malloc(*m * sizeof(double));
    double *s = malloc(*n * sizeof(double));



	
	if(*lsdeg == 1){
		if(*lmdeg == 1){
			if(*mdegree == 0){
				nadwats(x,y,x,*m,*m,*h,*K,p);
			}
			else{
				loclinear(x,y,x,*m,*m,*h,*K,p);
			}
		}
		else{
			for(j=0; j<*m; j++){
				p[j] = mdegree[j];
			}
		}


		for(i=0; i<*m; i++){
			e[i] = pow(y[i] - p[i],2.0);
		}

		if(*sdegree == 0.0){
			nadwats(x,e,z,*m,*n,*hr,*Kr,s);
		}
		else
			loclinear(x,e,z,*m,*n,*hr,*Kr,s);
	}
	else{
		for(j=0; j<*n; j++){
			s[j] = sdegree[j];        
		}
	}


if(*inverse == 0){
		t0 = malloc(5001 * sizeof(double));
		mini = s[0];
		maxi = s[0];
		for(i=1; i<*n; i++){
			if(s[i] < mini)
				mini = s[i];
			else if(s[i] > maxi)
				maxi = s[i];
		}
		t0[0] = mini - *hd;
		dist = (maxi - mini + 2.0 * *hd) / ((double) (5000));
		for(i=1; i<5001; i++){
			t0[i]=t0[i-1]+dist;            
		}
		k = 5001;
	}
	else if(*inverse == 1)
		k = *l;


if(*tflag == 1 && *inverse == 0){
		mini = x[0];
		maxi = x[0];
		for(i=1; i<*m; i++){
			if(x[i] < mini)
				mini = x[i];
			else if(x[i] > maxi)
				maxi = x[i];
		}
		t[0] = mini;
		dist = (maxi - mini) / ((double) (*l-1));
		for(i=1; i<*l; i++){
			t[i] = t[i-1] + dist;
		}
	}



	if(*tflag == 1 && *inverse == 1){
		mini = s[0];
		maxi = s[0];
		for(i=1; i<*n; i++){
			if(s[i] < mini)
				mini = s[i];
			else if(s[i] > maxi)
				maxi = s[i];
		}
		t[0] = mini;
		dist = (maxi - mini) / ((double) (*l - 1));
		for(i=1; i<*l; i++){
			t[i] = t[i-1] + dist;            
		}
		t0 = malloc(*l * sizeof(double));
		t0 = t;
	}

   sinv = malloc(k * sizeof(double));

	for(j=0; j<k; j++){
		sum=0.0;
		for(i=0; i<*n; i++){
			arg = (s[i] - t0[j]) / *hd;
			if(arg > 1.0){
				integral=0.0;
			}
			else{
				switch (*Kd){
					case 1:
						integral=*hd * (3.0 / 4.0 * ((1.0 - max(-1.0,arg))-1.0/3.0*(1.0-pow(max(-1.0,arg),3.0))));
						break;
					case 2:
						integral=1.0 / 2.0 * *hd * (1.0 - max(arg,-1.0));
						break;
					case 3:
						integral= *hd * (1.0/16.0*(15.0*(1.0-max(arg,-1.0))-10.0*(1.0-pow(max(arg,-1.0),3.0))+3.0*(1.0-pow(max(arg,-1.0),5.0))));
						break;
					case 4:
						integral= *hd * (1.0/32.0*(35.0*(1.0-max(arg,-1.0))-35.0*(1.0-pow(max(arg,-1.0),3.0))+21.0*(1.0-pow(max(arg,-1.0),5.0))-5.0*(1.0-pow(max(arg,-1.0),7.0))));
						break;
					case 5:
						integral= *hd * (1.0/2.0-max(arg,-1.0)+1.0/2.0*max(arg,-1.0)*fabs(max(arg,-1.0)));
						break;
					case 6:
						integral= *hd * (1.0/2.0*(1-sin(max(arg,-1)*3.1415/2.0)));
						break;
					case 0:
						integral=0.0;
					break;
				}
			}

	sum = sum + integral;
		}
		sinv[j] = 1.0 / ((double) *n * *hd) * sum;
		
	}


if(*inverse == 0){
		for(j=0; j<*l; j++){
			erg[j] = mdach_i(t[j],t0,sinv,k);
            		}
	}
	else{
		for(j=0; j<k; j++){
			erg[j] = sinv[j];
		}
	}
}




void sdach_a_inv (double *x, double *y, double *z, double *t,int *m, int *n, int *l, int *tflag, double *h, double *hd, double *hr,int *K, int *Kd, int *Kr, double *mdegree, int *lmdeg, double *sdegree, int *lsdeg, int *inverse, double *erg){

    long i,j,k;
	double sum,arg,integral,dist,mini,maxi;
	double *sinv,*t0;
	double *p = malloc(*m * sizeof(double));
    double *e = malloc(*m * sizeof(double));
    double *s = malloc(*n * sizeof(double));



	
	if(*lsdeg == 1){
		if(*lmdeg == 1){
			if(*mdegree == 0){
				nadwats(x,y,x,*m,*m,*h,*K,p);
			}
			else{
				loclinear(x,y,x,*m,*m,*h,*K,p);
			}
		}
		else{
			for(j=0; j<*m; j++){
				p[j] = mdegree[j];
			}
		}



for(i=0; i<*m; i++){
			e[i] = pow(y[i] - p[i],2.0);
		}

		if(*sdegree == 0.0){
			nadwats(x,e,z,*m,*n,*hr,*Kr,s);
		}
		else
			loclinear(x,e,z,*m,*n,*hr,*Kr,s);
	}
	else{
		for(j=0; j<*n; j++){
			s[j] = sdegree[j];        
		}
	}


if(*inverse == 0){
		t0 = malloc(5001 * sizeof(double));
		mini = s[0];
		maxi = s[0];
		for(i=1; i<*n; i++){
			if(s[i] < mini)
				mini = s[i];
			else if(s[i] > maxi)
				maxi = s[i];
		}
		t0[0] = mini - *hd;
		dist = (maxi - mini + 2.0 * *hd) / ((double) (5000));
		for(i=1; i<5001; i++){
			t0[i]=t0[i-1]+dist;            
		}
		k = 5001;
	}
	else if(*inverse == 1)
		k = *l;


if(*tflag == 1 && *inverse == 0){
		mini = x[0];
		maxi = x[0];
		for(i=1; i<*m; i++){
			if(x[i] < mini)
				mini = x[i];
			else if(x[i] > maxi)
				maxi = x[i];
		}
		t[0] = mini;
		dist = (maxi - mini) / ((double) (*l-1));
		for(i=1; i<*l; i++){
			t[i] = t[i-1] + dist;
		}
	}



	if(*tflag == 1 && *inverse == 1){
		mini = s[0];
		maxi = s[0];
		for(i=1; i<*n; i++){
			if(s[i] < mini)
				mini = s[i];
			else if(s[i] > maxi)
				maxi = s[i];
		}
		t[0] = mini;
		dist = (maxi - mini) / ((double) (*l - 1));
		for(i=1; i<*l; i++){
			t[i] = t[i-1] + dist;            
		}
		t0 = malloc(*l * sizeof(double));
		t0 = t;
	}

   sinv = malloc(k * sizeof(double));

	for(j=0; j<k; j++){
		sum=0.0;
		for(i=0; i<*n; i++){
			arg = (s[i] - t0[j]) / *hd;
			if(arg <= -1.0){
				integral=0.0;
			}
			else{
				switch (*Kd){
					case 1:
						integral=*hd * (1.0 / 4.0 * ((2.0 + 3.0*min(1.0,arg))-pow(min(1.0,arg),3.0)));
						break;
					case 2:
						integral=1.0 / 2.0 * *hd * (min(arg,1.0) + 1.0);
						break;
					case 3:
						integral= *hd * (1.0/16.0*(15.0*((min(arg,1.0))+1.0)-10.0*(pow(min(arg,1.0),3.0)+1.0)+3.0*(pow(min(arg,1.0),5.0)+1.0)));
						break;
					case 4:
						integral= *hd * (1.0/32.0*(35.0*(min(arg,1.0)+1.0)-35.0*(pow(min(arg,1.0),3.0)+1.0)+21.0*(pow(min(arg,1.0),5.0)+1.0)-5.0*(pow(min(arg,1.0),7.0)+1.0)));
						break;
					case 5:
						integral= *hd * ((-1.0)/2.0+min(arg,1.0)-1.0/2.0*min(arg,1.0)*fabs(min(arg,1.0)));
						break;
					case 6:
						integral= *hd * (1.0/2.0*(1.0+sin(min(arg,1.0)*3.1415/2.0)));
						break;
					case 0:
						integral=0.0;
					break;
				}
			}

			sum = sum + integral;
		}
		sinv[j] = 1.0 / ((double) *n * *hd) * sum;
	}

	if(*inverse == 0){
		for(j=0; j<*l; j++){
			erg[j] = mdach_i(t[j],t0,sinv,k);
            		}
	}
	else{
		for(j=0; j<k; j++){
			erg[j] = sinv[j];
		}
	}
}




void sdach2_i_inv (double *x, double *y, double *z, double *t, int *m, int *n, int*l, int *r, int *tflag, double *h, double *hd, int *K, int *Kd, double *degree, int *ldeg, int *inverse, double *erg){

long i,j,k;
double sum,arg,integral,dist,mini,maxi;
double *sinv,*t0;
double *d = malloc(*r * sizeof(double));
double *delta = malloc((*m - *r) * sizeof(double));
double *e = malloc((*m - *r) * sizeof(double));
double *s = malloc(*n * sizeof(double));
double *xneu = malloc((*m - *r) * sizeof(double));



	gewicht(d,*r,*r);




	for(i=0;i<(*m - *r);i++){
		delta[i] = 0.0;
		for(j=0;j<=*r;j++){
			delta[i] = delta[i] + d[j] * y[i+j];
		}
		e[i] = pow(delta[i],2.0);
	}

	for(i=0; i<(*m - *r); i++){
		xneu[i] = x[i];
	}

	if(*ldeg==1){
		if(*degree == 0){
			nadwats(xneu,e,z,(*m - *r),*n,*h,*K,s);
		}
		else{
			loclinear(xneu,e,z,(*m - *r),*n,*h,*K,s);
		}
	}
	else{
		for(j=0; j<(*n - *r); j++){
			s[j] = degree[j];  
		}
	}




if(*inverse == 0){
		t0 = malloc(5001 * sizeof(double));
		mini = s[0];
		maxi = s[0];
		for(i=1; i<*n; i++){
			if(s[i] < mini)
				mini = s[i];
			else if(s[i] > maxi)
				maxi = s[i];
		}
		t0[0] = mini - *hd;
		dist = (maxi - mini + 2.0 * *hd) / ((double) (5000));
		for(i=1; i<5001; i++){
			t0[i]=t0[i-1]+dist;            
		}
		k = 5001;
	}
	else if(*inverse == 1)
		k = *l;


if(*tflag == 1 && *inverse == 0){
		mini = x[0];
		maxi = x[0];
		for(i=1; i<*m; i++){
			if(x[i] < mini)
				mini = x[i];
			else if(x[i] > maxi)
				maxi = x[i];
		}
		t[0] = mini;
		dist = (maxi - mini) / ((double) (*l-1));
		for(i=1; i<*l; i++){
			t[i] = t[i-1] + dist;
		}
	}



	if(*tflag == 1 && *inverse == 1){
		mini = s[0];
		maxi = s[0];
		for(i=1; i<*n; i++){
			if(s[i] < mini)
				mini = s[i];
			else if(s[i] > maxi)
				maxi = s[i];
		}
		t[0] = mini;
		dist = (maxi - mini) / ((double) (*l - 1));
		for(i=1; i<*l; i++){
			t[i] = t[i-1] + dist;            
		}
		t0 = malloc(*l * sizeof(double));
		t0 = t;
	}

   sinv = malloc(k * sizeof(double));

	for(j=0; j<k; j++){
		sum=0.0;
		for(i=0; i<*n; i++){
			arg = (s[i] - t0[j]) / *hd;
			if(arg > 1.0){
				integral=0.0;
			}
			else{
				switch (*Kd){
					case 1:
						integral=*hd * (3.0 / 4.0 * ((1.0 - max(-1.0,arg))-1.0/3.0*(1.0-pow(max(-1.0,arg),3.0))));
						break;
					case 2:
						integral=1.0 / 2.0 * *hd * (1.0 - max(arg,-1.0));
						break;
					case 3:
						integral= *hd * (1.0/16.0*(15.0*(1.0-max(arg,-1.0))-10.0*(1.0-pow(max(arg,-1.0),3.0))+3.0*(1.0-pow(max(arg,-1.0),5.0))));
						break;
					case 4:
						integral= *hd * (1.0/32.0*(35.0*(1.0-max(arg,-1.0))-35.0*(1.0-pow(max(arg,-1.0),3.0))+21.0*(1.0-pow(max(arg,-1.0),5.0))-5.0*(1.0-pow(max(arg,-1.0),7.0))));
						break;
					case 5:
						integral= *hd * (1.0/2.0-max(arg,-1.0)+1.0/2.0*max(arg,-1.0)*fabs(max(arg,-1.0)));
						break;
					case 6:
						integral= *hd * (1.0/2.0*(1-sin(max(arg,-1)*3.1415/2.0)));
						break;
					case 0:
						integral=0.0;
					break;
				}
			}

	sum = sum + integral;
		}
		sinv[j] = 1.0 / ((double) *n * *hd) * sum;
		
	}


if(*inverse == 0){
		for(j=0; j<*l; j++){
			erg[j] = mdach_i(t[j],t0,sinv,k);
					}
	}
	else{
		for(j=0; j<k; j++){
			erg[j] = sinv[j];
		}
	}

}






void sdach2_a_inv (double *x, double *y, double *z, double *t, int *m, int *n, int*l, int *r, int *tflag, double *h, double *hd, int *K, int *Kd, double *degree, int *ldeg, int *inverse, double *erg){

long i,j,k;
double sum,arg,integral,dist,mini,maxi;
double *sinv,*t0;
double *d = malloc(*r * sizeof(double));
double *delta = malloc((*m - *r) * sizeof(double));
double *e = malloc((*m - *r) * sizeof(double));
double *s = malloc(*n * sizeof(double));
double *xneu = malloc((*m - *r) * sizeof(double));



	gewicht(d,*r,*r);




	for(i=0;i<(*m - *r);i++){
		delta[i] = 0.0;
		for(j=0;j<=*r;j++){
			delta[i] = delta[i] + d[j] * y[i+j];
		}
		e[i] = pow(delta[i],2.0);
	}

	for(i=0; i<(*m - *r); i++){
		xneu[i] = x[i];
	}

	if(*ldeg==1){
		if(*degree == 0){
			nadwats(xneu,e,z,(*m - *r),*n,*h,*K,s);
		}
		else{
			loclinear(xneu,e,z,(*m - *r),*n,*h,*K,s);
		}
	}
	else{
		for(j=0; j<(*n - *r); j++){
			s[j] = degree[j];  
		}
	}




if(*inverse == 0){
		t0 = malloc(5001 * sizeof(double));
		mini = s[0];
		maxi = s[0];
		for(i=1; i<*n; i++){
			if(s[i] < mini)
				mini = s[i];
			else if(s[i] > maxi)
				maxi = s[i];
		}
		t0[0] = mini - *hd;
		dist = (maxi - mini + 2.0 * *hd) / ((double) (5000));
		for(i=1; i<5001; i++){
			t0[i]=t0[i-1]+dist;            
		}
		k = 5001;
	}
	else if(*inverse == 1)
		k = *l;


if(*tflag == 1 && *inverse == 0){
		mini = x[0];
		maxi = x[0];
		for(i=1; i<*m; i++){
			if(x[i] < mini)
				mini = x[i];
			else if(x[i] > maxi)
				maxi = x[i];
		}
		t[0] = mini;
		dist = (maxi - mini) / ((double) (*l-1));
		for(i=1; i<*l; i++){
			t[i] = t[i-1] + dist;
		}
	}



	if(*tflag == 1 && *inverse == 1){
		mini = s[0];
		maxi = s[0];
		for(i=1; i<*n; i++){
			if(s[i] < mini)
				mini = s[i];
			else if(s[i] > maxi)
				maxi = s[i];
		}
		t[0] = mini;
		dist = (maxi - mini) / ((double) (*l - 1));
		for(i=1; i<*l; i++){
			t[i] = t[i-1] + dist;            
		}
		t0 = malloc(*l * sizeof(double));
		t0 = t;
	}

   sinv = malloc(k * sizeof(double));

	for(j=0; j<k; j++){
		sum=0.0;
		for(i=0; i<*n; i++){
			arg = (s[i] - t0[j]) / *hd;
			if(arg <= -1.0){
				integral=0.0;
			}
			else{
				switch (*Kd){
					case 1:
						integral=*hd * (1.0 / 4.0 * ((2.0 + 3.0*min(1.0,arg))-pow(min(1.0,arg),3.0)));
						break;
					case 2:
						integral=1.0 / 2.0 * *hd * (min(arg,1.0) + 1.0);
						break;
					case 3:
						integral= *hd * (1.0/16.0*(15.0*((min(arg,1.0))+1.0)-10.0*(pow(min(arg,1.0),3.0)+1.0)+3.0*(pow(min(arg,1.0),5.0)+1.0)));
						break;
					case 4:
						integral= *hd * (1.0/32.0*(35.0*(min(arg,1.0)+1.0)-35.0*(pow(min(arg,1.0),3.0)+1.0)+21.0*(pow(min(arg,1.0),5.0)+1.0)-5.0*(pow(min(arg,1.0),7.0)+1.0)));
						break;
					case 5:
						integral= *hd * ((-1.0)/2.0+min(arg,1.0)-1.0/2.0*min(arg,1.0)*fabs(min(arg,1.0)));
						break;
					case 6:
						integral= *hd * (1.0/2.0*(1.0+sin(min(arg,1.0)*3.1415/2.0)));
						break;
					case 0:
						integral=0.0;
					break;
				}
			}

			sum = sum + integral;
		}
		sinv[j] = 1.0 / ((double) *n * *hd) * sum;
	}

	if(*inverse == 0){
		for(j=0; j<*l; j++){
			erg[j] = mdach_i(t[j],t0,sinv,k);
					}
	}
	else{
		for(j=0; j<k; j++){
			erg[j] = sinv[j];
		}
	}
}


