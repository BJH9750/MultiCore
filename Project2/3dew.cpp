#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "sys/time.h"
#include <vector>
#include <omp.h>

using namespace std;

#define PIE 3.1415926
float c0 = -2.927222164;
int mm = 5;

float *u, *v, *w, *up, *up1, *up2, 
          *vp, *vp1, *vp2, *wp, *wp1, *wp2, 
          *us, *us1, *us2, *vs, *vs1, *vs2,
          *ws, *ws1, *ws2, *vpp, *density, *vss;

float c[5][7];
float *wave;
float nshot,t0;
float dtx,dtz,dtxz,dr1,dr2,dtx4,dtz4,dtxz4;
int nx,ny,nz,lt,nedge;
float frequency;
float velmax;
float dt;
int ncx_shot1,ncy_shot1,ncz_shot;
int ishot,ncy_shot,ncx_shot;
float unit;
int nxshot,nyshot,dxshot,dyshot;

void free_all();
void iterate(int iter);

int main(int argc, char **argv)
{
    int i,j,k,kk,kkk,l;

    char infile[80],outfile[80],logfile[80],tmp[80];
    FILE  *fin, *fout, *flog;
    struct timeval start,end;
    float all_time;

    if(argc<4)
    {
        printf("please add 3 parameter: inputfile, outfile, logfile\n");
        exit(0);
    }

    strcpy(infile,argv[1]);
    strcpy(outfile,argv[2]);
    strcpy(logfile,argv[3]);

    strcpy(tmp,"date ");
    strncat(tmp, ">> ",3);
    strncat(tmp, logfile, strlen(logfile));
    flog = fopen(logfile,"w");
    fprintf(flog,"------------start time------------\n");
    fclose(flog);
    system(tmp);
    gettimeofday(&start,NULL);

    fin = fopen(infile,"r");
    if(fin == NULL)
    {
        printf("file %s is  not exist\n",infile);
        exit(0);
    }
    fscanf(fin,"nx=%d\n",&nx);
    fscanf(fin,"ny=%d\n",&ny);
    fscanf(fin,"nz=%d\n",&nz);
    fscanf(fin,"lt=%d\n",&lt);
    fscanf(fin,"nedge=%d\n",&nedge);
    fscanf(fin,"ncx_shot1=%d\n",&ncx_shot1);
    fscanf(fin,"ncy_shot1=%d\n",&ncy_shot1);
    fscanf(fin,"ncz_shot=%d\n",&ncz_shot);
    fscanf(fin,"nxshot=%d\n",&nxshot);
    fscanf(fin,"nyshot=%d\n",&nyshot);
    fscanf(fin,"frequency=%f\n",&frequency);
    fscanf(fin,"velmax=%f\n",&velmax);
    fscanf(fin,"dt=%f\n",&dt);
    fscanf(fin,"unit=%f\n",&unit);
    fscanf(fin,"dxshot=%d\n",&dxshot);
    fscanf(fin,"dyshot=%d\n",&dyshot);
    fclose(fin);

    printf("\n--------workload parameter--------\n");
    printf("nx=%d\n",nx);
    printf("ny=%d\n",ny);
    printf("nz=%d\n",nz);
    printf("lt=%d\n",lt);
    printf("nedge=%d\n",nedge);
    printf("ncx_shot1=%d\n",ncx_shot1);
    printf("ncy_shot1=%d\n",ncy_shot1);
    printf("ncz_shot=%d\n",ncz_shot);
    printf("nxshot=%d\n",nxshot);
    printf("nyshot=%d\n",nyshot);
    printf("frequency=%f\n",frequency);
    printf("velmax=%f\n",velmax);
    printf("dt=%f\n",dt);
    printf("unit=%f\n",unit);
    printf("dxshot=%d\n",dxshot);
    printf("dyshot=%d\n\n",dyshot);
    flog = fopen(logfile,"a");
    fprintf(flog,"\n--------workload parameter--------\n");
    fprintf(flog,"nx=%d\n",nx);
    fprintf(flog,"ny=%d\n",ny);
    fprintf(flog,"nz=%d\n",nz);
    fprintf(flog,"lt=%d\n",lt);
    fprintf(flog,"nedge=%d\n",nedge);
    fprintf(flog,"ncx_shot1=%d\n",ncx_shot1);
    fprintf(flog,"ncy_shot1=%d\n",ncy_shot1);
    fprintf(flog,"ncz_shot=%d\n",ncz_shot);
    fprintf(flog,"nxshot=%d\n",nxshot);
    fprintf(flog,"nyshot=%d\n",nyshot);
    fprintf(flog,"frequency=%f\n",frequency);
    fprintf(flog,"velmax=%f\n",velmax);
    fprintf(flog,"dt=%f\n",dt);
    fprintf(flog,"unit=%f\n",unit);
    fprintf(flog,"dxshot=%d\n",dxshot);
    fprintf(flog,"dyshot=%d\n\n",dyshot);
    fclose(flog);

    u       = (float*)malloc(sizeof(float)*nz*ny*nx);
    v       = (float*)malloc(sizeof(float)*nz*ny*nx);
    w       = (float*)malloc(sizeof(float)*nz*ny*nx);
    up      = (float*)malloc(sizeof(float)*nz*ny*nx);
    up1     = (float*)malloc(sizeof(float)*nz*ny*nx);
    up2     = (float*)malloc(sizeof(float)*nz*ny*nx);
    vp      = (float*)malloc(sizeof(float)*nz*ny*nx);
    vp1     = (float*)malloc(sizeof(float)*nz*ny*nx);
    vp2     = (float*)malloc(sizeof(float)*nz*ny*nx);
    wp      = (float*)malloc(sizeof(float)*nz*ny*nx);
    wp1     = (float*)malloc(sizeof(float)*nz*ny*nx);
    wp2     = (float*)malloc(sizeof(float)*nz*ny*nx);
    us      = (float*)malloc(sizeof(float)*nz*ny*nx);
    us1     = (float*)malloc(sizeof(float)*nz*ny*nx);
    us2     = (float*)malloc(sizeof(float)*nz*ny*nx);
    vs      = (float*)malloc(sizeof(float)*nz*ny*nx);
    vs1     = (float*)malloc(sizeof(float)*nz*ny*nx);
    vs2     = (float*)malloc(sizeof(float)*nz*ny*nx);
    ws      = (float*)malloc(sizeof(float)*nz*ny*nx);
    ws1     = (float*)malloc(sizeof(float)*nz*ny*nx);
    ws2     = (float*)malloc(sizeof(float)*nz*ny*nx);
    vpp     = (float*)malloc(sizeof(float)*nz*ny*nx);
    density = (float*)malloc(sizeof(float)*nz*ny*nx);
    vss     = (float*)malloc(sizeof(float)*nz*ny*nx);
    wave    = (float*)malloc(sizeof(float)*lt);

    nshot=nxshot*nyshot;
    t0=1.0/frequency;
    
    omp_set_nested(1);
    int idx, idx1, idx2, idx3;
    #pragma omp parallel for collapse(3)
    for(i = 0; i < nz; ++i){
        for(j = 0; j < ny; ++j){
            for(k = 0; k < nx; ++k){
                idx = i*ny*nx+j*nx+k;
                if(i < 210){
                    vpp[idx] = 2300.;
                    vss[idx] = 1232.;
                    density[idx] = 1.;
                }else if(i >= 260){
                    vpp[idx] = 3500.;
                    vss[idx] = 1909.;
                    density[idx] = 2.5;   
                }else{
                    vpp[idx] = 2800.;
                    vss[idx] = 1509.;
                    density[idx] = 2.;
                }
            }
        }
    }

    for(l = 0; l < lt; ++l){
        float sp = PIE * frequency * (l * dt - t0);
        float fx = 100000. * exp(-sp*sp) * (1.-2.*sp*sp);
        wave[l] = fx;
    }
    
    c[0][0] = 1.66666665;
    c[1][0] = -0.23809525;
    c[2][0] = 0.03968254;
    c[3][0] = -0.004960318;
    c[4][0] = 0.0003174603;
    c[0][1] = 0.83333;
    c[1][1] = -0.2381;
    c[2][1] = 0.0595;
    c[3][1] = -0.0099;
    c[4][1] = 0.0008;

    for(i = 0; i < 5; ++i){
        for(j = 0; j < 5; ++j){
            c[j][i+2] = c[i][1] * c[j][1];
        }
    }
        

    dtx=dt/unit;
    dtz=dt/unit;

    dr1=dtx*dtx/2.;
    dr2=dtz*dtz/2.;


    float vvp2,drd1,drd2,vvs2,tempux2,tempuy2,tempuz2,tempvx2,tempvy2,tempvz2,
            tempwx2,tempwy2,tempwz2,tempuxz,tempuxy,tempvyz,tempvxy,tempwxz,tempwyz;

    fout=fopen(outfile,"wb");

    char **strs = (char **)malloc(sizeof(char *) * (nshot + 1));
    for(i = 0; i <= nshot; ++i){
        strs[i] = (char *)malloc(sizeof(char) * 32);
    }


    for(ishot = 1; ishot <= nshot; ++ishot){
        printf("shot=%d\n",ishot);
        sprintf(strs[ishot], "shot=%d\n",ishot);

	    // flog = fopen(logfile,"a");
        // fprintf(flog,"shot=%d\n",ishot);
	    // fclose(flog);
        #pragma omp parallel for collapse(3) private (idx1)
        for(i = 0; i < nz; ++i){
            for(j = 0; j < ny; ++j){
                for(k = 0; k < nx; ++k){
                    idx1 = i*ny*nx+j*nx+k;
                    u[idx1] = 0.0f;
                    v[idx1] = 0.0f;
                    w[idx1] = 0.0f;
                    up[idx1] = 0.0f;
                    up1[idx1] = 0.0f;
                    up2[idx1] = 0.0f;
                    vp[idx1] = 0.0f;
                    vp1[idx1] = 0.0f;
                    vp2[idx1] = 0.0f;
                    wp[idx1] = 0.0f;
                    wp1[idx1] = 0.0f;
                    wp2[idx1] = 0.0f;
                    us[idx1] = 0.0f;
                    us1[idx1] = 0.0f;
                    us2[idx1] = 0.0f;
                    vs[idx1] = 0.0f;
                    vs1[idx1] = 0.0f;
                    vs2[idx1] = 0.0f;
                    ws[idx1] = 0.0f;
                    ws1[idx1] = 0.0f;
                    ws2[idx1] = 0.0f;
                }// for(k = 0; k < nx; ++k) end
            }// for(j = 0; j < ny; ++j) end
        }// for(i = 0; i < nz; ++i) end

        ncy_shot = ncy_shot1 + (ishot / nxshot) * dyshot;
        ncx_shot = ncx_shot1 + (ishot % nxshot) * dxshot;
        
        for(l = 1; l <= lt; ++l){
            iterate(l);
        }// for(l = 1; l <= lt; ++l)

        fwrite(up+169*ny*nx,sizeof(float),ny*nx,fout);
       
    
    }//for(ishot=1;ishot<=nshot;ishot++) end
    fclose(fout);

    flog = fopen(logfile,"a");
    for(i = 1; i <= nshot; ++i){
        fputs(strs[i], flog);
    }
    fclose(flog);

    free_all();

    gettimeofday(&end,NULL);
    all_time = (end.tv_sec-start.tv_sec)+(float)(end.tv_usec-start.tv_usec)/1000000.0;
    printf("run time:\t%f s\n",all_time);
    flog = fopen(logfile,"a");
    fprintf(flog,"\nrun time:\t%f s\n\n",all_time);
    fprintf(flog,"------------end time------------\n");
    fclose(flog);
    system(tmp);
    return 1;
}

void free_all(){
    free(u);
    free(v);
    free(w);
    free(up);
    free(up1);
    free(up2);
    free(vp);
    free(vp1);
    free(vp2);
    free(wp);
    free(wp1);
    free(wp2);
    free(us);
    free(us1);
    free(us2);
    free(vs);
    free(vs1);
    free(vs2);
    free(ws);
    free(ws1);
    free(ws2);
    free(vpp);
    free(density);
    free(vss);
    free(wave);
}

void iterate(int iter){
    int i,j,k,kk,kkk;
    float vvp2,drd1,drd2,vvs2,tempux2,tempuy2,tempuz2,tempvx2,tempvy2,tempvz2,
        tempwx2,tempwy2,tempwz2,tempuxz,tempuxy,tempvyz,tempvxy,tempwxz,tempwyz;
    int idx, idx1, idx2, idx3;
    float px,sx;
    int xmax, nleft, nright, nfront, nback, ntop, nbottom;
    xmax = (iter * dt * velmax / unit) + 10;

    nleft = ncx_shot - xmax;
    nright = ncx_shot + xmax;
    nfront = ncy_shot - xmax;
    nback = ncy_shot + xmax;
    ntop = ncz_shot - xmax;
    nbottom = ncz_shot + xmax;
    if(nleft < 5) nleft = 5;
    if(nright > nx - 5) nright = nx - 5;
    if(nfront < 5) nfront = 5;
    if(nback > ny - 5) nback = ny - 5;
    if(ntop < 5) ntop = 5;
    if(nbottom > nz - 5) nbottom = nz - 5;

    ntop = ntop - 1;
    nfront = nfront - 1;
    nleft = nleft - 1;

    #pragma omp parallel for ordered collapse(3) private(vvp2, vvs2, px, sx,tempux2,tempuy2,tempuz2,tempvx2,tempvy2,tempvz2,\
                                                            tempwx2,tempwy2,tempwz2,tempuxz,tempuxy,tempvyz,tempvxy,tempwxz,tempwyz, idx2)
    for(k = ntop; k < nbottom; ++k){
        for(j = nfront; j < nback; ++j){
            for(i = nleft; i < nright; ++i){
                idx2 = k*ny*nx+j*nx+i;
                if(i==ncx_shot-1 && j==ncy_shot-1 && k==ncz_shot-1){
                    px=1.;
                    sx=0.;
                }else{
                    px=0.;
                    sx=0.;
                }

                vvp2=vpp[idx2]*vpp[idx2];
                vvs2=vss[idx2]*vss[idx2];

                tempux2=0.0f; tempuy2=0.0f; tempuz2=0.0f; tempvx2=0.0f; tempvy2=0.0f; tempvz2=0.0f; tempwx2=0.0f; tempwy2=0.0f; tempwz2=0.0f; 
                tempuxz=0.0f; tempuxy=0.0f; tempvyz=0.0f; tempvxy=0.0f; tempwxz=0.0f; tempwyz=0.0f;

                for(kk = 1; kk <= mm; ++kk){
                    tempux2 = tempux2 + c[kk-1][0] * (u[k*ny*nx+j*nx+(i+kk)] + u[k*ny*nx+j*nx+(i-kk)]);
                    tempuy2 = tempuy2 + c[kk-1][0] * (u[k*ny*nx+(j+kk)*nx+i] + u[k*ny*nx+(j-kk)*nx+i]);
                    tempuz2 = tempuz2 + c[kk-1][0] * (u[(k+kk)*ny*nx+j*nx+i] + u[(k-kk)*ny*nx+j*nx+i]);

                    tempvx2 = tempvx2 + c[kk-1][0] * (v[k*ny*nx+j*nx+(i+kk)] + v[k*ny*nx+j*nx+(i-kk)]);
                    tempvy2 = tempvy2 + c[kk-1][0] * (v[k*ny*nx+(j+kk)*nx+i] + v[k*ny*nx+(j-kk)*nx+i]);
                    tempvz2 = tempvz2 + c[kk-1][0] * (v[(k+kk)*ny*nx+j*nx+i] + v[(k-kk)*ny*nx+j*nx+i]);

                    tempwx2 = tempwx2 + c[kk-1][0] * (w[k*ny*nx+j*nx+(i+kk)] + w[k*ny*nx+j*nx+(i-kk)]);
                    tempwy2 = tempwy2 + c[kk-1][0] * (w[k*ny*nx+(j+kk)*nx+i] + w[k*ny*nx+(j-kk)*nx+i]);
                    tempwz2  =tempwz2 + c[kk-1][0] * (w[(k+kk)*ny*nx+j*nx+i] + w[(k-kk)*ny*nx+j*nx+i]);

                } //for(kk=1;kk<=mm;kk++) end

                tempux2 = (tempux2+c0*u[idx2])*vvp2*dtx*dtx;
                tempuy2 = (tempuy2+c0*u[idx2])*vvs2*dtx*dtx;
                tempuz2 = (tempuz2+c0*u[idx2])*vvs2*dtz*dtz;

                tempvx2 = (tempvx2+c0*v[idx2])*vvs2*dtx*dtx;
                tempvy2 = (tempvy2+c0*v[idx2])*vvp2*dtx*dtx;
                tempvz2 = (tempvz2+c0*v[idx2])*vvs2*dtz*dtz;

                tempwx2 = (tempwx2+c0*w[idx2])*vvs2*dtx*dtx;
                tempwy2 = (tempwy2+c0*w[idx2])*vvs2*dtx*dtx;
                tempwz2 = (tempwz2+c0*w[idx2])*vvp2*dtz*dtz;

                for(kk=1;kk<=mm;kk++){
                    for(kkk=1;kkk<=mm;kkk++){
                        tempuxz=tempuxz+c[kkk-1][1+kk]*(u[(k+kkk)*ny*nx+j*nx+(i+kk)]
                                        -u[(k-kkk)*ny*nx+j*nx+(i+kk)]
                                        +u[(k-kkk)*ny*nx+j*nx+(i-kk)]
                                        -u[(k+kkk)*ny*nx+j*nx+(i-kk)]);
                        tempuxy=tempuxy+c[kkk-1][1+kk]*(u[k*ny*nx+(j+kkk)*nx+(i+kk)]
                                        -u[k*ny*nx+(j-kkk)*nx+(i+kk)]
                                        +u[k*ny*nx+(j-kkk)*nx+(i-kk)]
                                        -u[k*ny*nx+(j+kkk)*nx+(i-kk)]);

                        tempvyz=tempvyz+c[kkk-1][1+kk]*(v[(k+kkk)*ny*nx+(j+kk)*nx+i]
                                        -v[(k-kkk)*ny*nx+(j+kk)*nx+i]
                                        +v[(k-kkk)*ny*nx+(j-kk)*nx+i]
                                        -v[(k+kkk)*ny*nx+(j-kk)*nx+i]);
                        tempvxy=tempvxy+c[kkk-1][1+kk]*(v[k*ny*nx+(j+kkk)*nx+(i+kk)]
                                        -v[k*ny*nx+(j-kkk)*nx+(i+kk)]
                                        +v[k*ny*nx+(j-kkk)*nx+(i-kk)]
                                        -v[k*ny*nx+(j+kkk)*nx+(i-kk)]);

                        tempwyz=tempwyz+c[kkk-1][1+kk]*(w[(k+kkk)*ny*nx+(j+kk)*nx+i]
                                        -w[(k-kkk)*ny*nx+(j+kk)*nx+i]
                                        +w[(k-kkk)*ny*nx+(j-kk)*nx+i]
                                        -w[(k+kkk)*ny*nx+(j-kk)*nx+i]);
                        tempwxz=tempwxz+c[kkk-1][1+kk]*(w[(k+kkk)*ny*nx+j*nx+(i+kk)]
                                        -w[(k-kkk)*ny*nx+j*nx+(i+kk)]
                                        +w[(k-kkk)*ny*nx+j*nx+(i-kk)]
                                        -w[(k+kkk)*ny*nx+j*nx+(i-kk)]);
                    } // for(kkk=1;kkk<=mm;kkk++) end
                } //for(kk=1;kk<=mm;kk++) end

                up[idx2] = 2.*up1[idx2]-up2[idx2]
                                +tempux2+tempwxz*vvp2*dtz*dtx
                                +tempvxy*vvp2*dtz*dtx;
                vp[idx2] = 2.*vp1[idx2]-vp2[idx2]
                                +tempvy2+tempuxy*vvp2*dtz*dtx
                                +tempwyz*vvp2*dtz*dtx;
                wp[idx2] = 2.*wp1[idx2]-wp2[idx2]
                                +tempwz2+tempuxz*vvp2*dtz*dtx
                                +tempvyz*vvp2*dtz*dtx
                                +px*wave[iter-1];
                us[idx2] = 2.*us1[idx2]-us2[idx2] + tempuy2+tempuz2
                                -tempvxy*vvs2*dtz*dtx-tempwxz*vvs2*dtz*dtx;
                vs[idx2] = 2.*vs1[idx2]-vs2[idx2] + tempvx2+tempvz2
                                -tempuxy*vvs2*dtz*dtx-tempwyz*vvs2*dtz*dtx;
                ws[idx2] = 2.*ws1[idx2]-ws2[idx2] + tempwx2+tempwy2
                                -tempuxz*vvs2*dtz*dtx-tempvyz*vvs2*dtz*dtx;
            }//for(i = nleft; i < nright; i++) end
        }// for(j = nfront; j < nback; j++) end        
    }// for(k = ntop; k < nbottom; k++) end

    #pragma omp parallel for collapse(3) private(idx3)
    for(k = ntop; k < nbottom; ++k){
        for(j = nfront; j < nback; ++j){
            for(i = nleft; i < nright; ++i){
                idx3 = k*ny*nx+j*nx+i;
                u[idx3] = up[idx3] + us[idx3];
                v[idx3] = vp[idx3] + vs[idx3];
                w[idx3] = wp[idx3] + ws[idx3];

                up2[idx3] = up1[idx3];
                up1[idx3] = up[idx3];
                us2[idx3] = us1[idx3];
                us1[idx3] = us[idx3];
                vp2[idx3] = vp1[idx3];
                vp1[idx3] = vp[idx3];
                vs2[idx3] = vs1[idx3];
                vs1[idx3] = vs[idx3];
                wp2[idx3] = wp1[idx3];
                wp1[idx3] = wp[idx3];
                ws2[idx3] = ws1[idx3];
                ws1[idx3] = ws[idx3];
            }// for(i = nleft; i < nright; ++i)
        }// for(j = nfront; j < nback; ++j)
    }// for(k = ntop; k < nbottom; ++k)

}
