/* Rate-State Friction Simulator, based on 'shsq.c'
        This program simulates two types of tests:
1) slide-hold-slide
2) velocity jump

  The simulator uses 5th order Runge-Kutta with adaptive step size control.

  Various friction state evolution laws are possible */

#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	"rsp.h"
#define		SR_2   1.41421356
#define		TRUE 1
#define		FALSE 0
#define		tinyAmt 1e-10			/*to deal w/ rounding summation for time*/
 
/*see get_mu_xt.c for more extensive comments*/
/*Last modified:
24/4/96, written based on foq.c and qi_look.c. Uses adaptive step size control
25/4/96: added code to calculate del_mu and a few other things --so that I can look
	at systematics of strengthening rate.
26/4/96: changed V_REF and MU_REF from definitions to "variables" defined in the rsp structure
		I'll still use a fixed reference velocity, but I'm going to fix mu_ref for 
	a given run, so that mu at 1mic/sec =0.6				   
26/4/96: changed op to only give 9 digits. For std velocities and mu's the first 8 or so 
	are independent of v_ref and mu_ref.

--lots of other changes...

30/3/99: modified to calculate porosity change, via Segall and Rice (95) and Sleep (95) models
6/7/99: modified so ref. friction level is based on 0.6 at 1 mic/s, rather than 0.6 at any velocity

11/2/02: minor modifications of command line, etc. so that it can be used
         by others, geosc508 class in particular.

2021.12.19 set up for multiple velocity steps. Change main op array to time rather than disp.
2022.01.03  
2023.12.03   output max slip velocity --useful for looking at posthold 'overshoot' slip events
cjm
*/

char	stop_at_peak = 'n';
char	*progname;

int main (int ac, char *av[])
{

char	outfile[64], string[65536], string1[65536], string2[65536];
/*char	*strcpy(), *index(), *strcat();*/
char	tim_out ='n', vel_out ='n', state_out ='n', porosity_time_out = 'n', sd_out = 'n';
char	segall_rice_law = 'n', doing_hold = 'n', doing_ramp = 'n';
char	set_ref_fric = 'n', do_mvs = 'n';
FILE	*tim_file, *vel_file, *state_file, *dist_file, *x_phi_file, *t_phi_file, *sd_file, *fopen();
int	ndata,i,j,kk,nsteps, vstep, data_count;
int get_mu_at_xt( double *x, double *mod_mu, double *v_slider, int *ndata, struct rs_parameters *rsp, char *outfile, double *state, double *sd, double *phi, int nsteps, char doing_hold);
double  log(), exp();
double	ceil();
double	max_x, op_inc_disp, op_inc_hold;
double  *v_slider,*state,*sd, *phi, hold_time=0;
double	calc_time, temp_time;
double  temp, temp1, temp2, step_size, step_duration, lpDisp;

struct  rs_parameters rsp;

	progname = av[0];
	if(ac < 15)
	{
        fprintf(stderr,"\t\nThis program solves the rate and state friction laws with 1-d elastic interaction for quasi-static conditions (neither inertia nor radiation damping are included). One or two state variables can be used.\n\n");
        fprintf(stderr,"Usage: %s op_file law v_init v2 hold_time op_inc_hold op_inc_disp max_disp k a b1 Dc1 b2 dc2 -svtdm \n",progname);
	fprintf(stderr,"where law is: d, r, p, or s \n\td=Dieterich, slowness law, r=Runia, slip law,\n\tp=Perrin-Rice, quadratic law, s=Segall and Rice eq'n 17 \n");
        fprintf(stderr,"-t = output time vs mu\n");
        fprintf(stderr,"-v = output ln(v_s/v_init) vs mu\n");
        fprintf(stderr,"-s = output time vs state(theta)\n");
        fprintf(stderr,"-d = output Slider_Displacement vs mu\n");
	fprintf(stderr,"-f = output time vs porosity\n");
        fprintf(stderr,"-p = stop calculation at peak friction on reload, write a data point there\n");
        fprintf(stderr,"-r = set ref fric to 0.6 at 1 mic/s, rather than to arb. value to give 0.6 as initial value\n");
	fprintf(stderr,"-m = do multiple velocity steps, for each step list velocity and time duration \n\t If you use multiple options m has to be last option.\n\t If you want a ramp make v < 0 and give the ramp duration\n\t The velocity v2 is ignored for this option\n\t You can do a hold. Your velocity sequence will be used for the reload after the hold\n");
	fprintf(stderr,"\nSlip vs time is always output for this version\n");
	fprintf(stderr,"Command line definitions:\n");
	fprintf(stderr,"\tv_init is the initial, steady-state velocity \n");
        fprintf(stderr,"\tv2 is the reload velocity following the hold or the 2nd velocity for a simple v-step\n");
        fprintf(stderr,"\thold_time is the hold time in seconds; enter a negative number if you don't want a hold\n");
        fprintf(stderr,"\top_inc_hold is the time increment for output during the hold\n");
        fprintf(stderr,"\top_inc_disp is the time increment for output following the hold\n");
        fprintf(stderr,"\tdisp is the final disp for the calc for a simple v-step, a simple hold, or for a vstep with a ramp\n");
        fprintf(stderr,"\tTo do a 1 state variable model, set dc2 < 0\n");
        fprintf(stderr,"\tNote: k should  have dimensions of 1/length\n");
        fprintf(stderr,"\n");
 
        fprintf(stderr,"\t\nOne or more ASCII text files are written as output. Each row has two columns (time, displacement, friction, etc --as chosen by command line options summarized above) and the output frequency is specified by the command line parameters op_inc_hold and op_inc_disp \n\n");
        fprintf(stderr,"\tExample command lines:\n\n");
        fprintf(stderr,"1: %s junk r 1 1 100  0.1  0.1 200 1e-3 0.008 0.01 10 1 -10 -t\n\n",av[0]);
        fprintf(stderr,"2: %s junk r 1 10 -1  1    0.1 300 1e-3 0.008 0.01 10 1 -10 -t\n\n", av[0]);
        fprintf(stderr,"3: %s junk r 1 10 -1  1    0.05 300 1e-3 0.008 0.01 10 1 -10 -t -v -m 2.5 200 10 31 30 15 1 200 10 40\n\n", av[0]);
        fprintf(stderr,"4: %s junk r 1 10 100 1    0.05 300 1e-3 0.008 0.01 10 1 -10 -t -v -m 2.5 200 10 31 30 15 1 200 10 40\n\n", av[0]);
        fprintf(stderr,"5: %s junk r 1 10 -1  0.1  0.1 300 1e-3 0.008 0.01 10 1 -10 -t -v -m -20 2.3 \n\n", av[0]);
        fprintf(stderr,"6: %s junk r 1 10 30  0.1  0.1 300 1e-3 0.008 0.01 10 1 -10 -t -v -m -20 2.3 \n\n", av[0]);
        fprintf(stderr,"\tThe first example specifies a slide-hold-slide test with one state variable. The output file name is junk, the dieterich state evolution law is used, the initial and re-load velocities are 1 micron/s, hold time is 100 sec, output is written every 0.1 sec during the hold and 0.1 micron after the hold up to a maximum of 200 microns, stiffness is 0.001 (friction) per micron, and constitutive parameters are a=0.008, b=0.01 and Dc of 10 microns. Note that there is nothing intrinsic about the units; the length unit could just as well be taken as meters and the time as hours. The only requirement is that all units be consistent. The files junk.dis and junk.tim are output. Each file contains a one line header and then x,y pairs of: junk.dis: load point displacement and friction junk.tim time and friction. \n\n");
        fprintf(stderr,"\tThe second example specifies a velocity step test from 1 mic/s to 10 mic/s with the Runia law. Data are written every 0.1 micron up to a maximum of 300 microns. The friction parameters and output files are the same as the first example.\n");
        fprintf(stderr,"\tExample 3 is a series of velocity steps from 1 to (v2 is ignored w/ the -m option) 2.5 (200 sec) to 10 (31 sec) to 30 (15 sec) to 1 (200 sec) to 10 (40 sec). Data are written every 0.05 s until the last step is done. \n");
        fprintf(stderr,"\tExample 4 is a 100 sec hold after sliding at 1 µm/s with a post-hold reload that is a series of velocity steps from end-of-hold-v to 2.5 (200 sec) to 10 (31 sec) to 30 (15 sec) to 1 (200 sec) to 10 (40 sec). Data are written every 0.05 s until the last step is done. \n");
        fprintf(stderr,"\tExample 5 is like Example 2 except the velocity jump is from 1 to 20 (v2 is ignored) and happens as a ramp of 100 steps over a duration of 2.3 sec \n");
        fprintf(stderr,"\tExample 6 is like Example 4 except the reload velocity after the hold (v2 is ignored) is ramp up to 20 µm/s in 100 steps over a duration of 2.3 sec \n");

	exit(2);
	}

	switch(*av[2])
	{
		case 'd':
		case 'D':
			rsp.law = 'd';
			break;
		case 'r':
		case 'R':
			rsp.law = 'r';
			break;
		case 'p':
		case 'P':
			rsp.law = 'p';
			break;
		case 's':
		case 'S':
			rsp.law = 'd';		/*use Dieterich law, with eq'n 17*/
			segall_rice_law = 'y';
			break;
		default :
			fprintf(stderr,"Please choose a law by using either d, r, p, s, or S\n");
			exit(1);
	}
	rsp.vel_list[0] = atof(av[3]);	/*init vel*/
	rsp.vel_list[2] = atof(av[4]);	/*v2 or reload vel*/
	hold_time = atof(av[5]);	/*this value is negative if we're not doing a hold */
	op_inc_hold = atof(av[6]);		/*time increment for ouput during the hold*/	
	op_inc_disp = atof(av[7]);		/*time increment for ouput during slip*/
	max_x = atof(av[8]);

	rsp.vel_dur_list[0] = 0;			/*the first step happens immediately*/

 	if(hold_time > 0)
	{
		doing_hold = 'y';	/*default is 'n'*/
		rsp.vel_list[1] = 0;		/*vel_list[1] is the hold*/
		rsp.vel_dur_list[1] = hold_time;
		rsp.vel_dur_list[2] = max_x/rsp.vel_list[2];	/*this is a default if no mvs option*/
		nsteps = 2; /*step 2 is the reload, which we may replace w/ mvs */
		vstep = 2; /*vstep 0 is from v0 to 0, vstep 1 is reload w/ v2 or mvs */
	}
	else	/*simple shs, or mvs*/
	{
		hold_time = 0.0;
		rsp.vel_list[1] = rsp.vel_list[2]; 		/*vel_list[1] is the step from v0 */
		rsp.vel_dur_list[1] = max_x/rsp.vel_list[1];	/*this is a default if no mvs option*/
		nsteps = 1; 
		vstep = 1; /*vstep 0 is from v0 to v-reload or from mvs */
	}

	rsp.stiff = atof(av[9]);
	rsp.a = atof(av[10]);
	rsp.b1 = atof(av[11]);
	rsp.dc1 = atof(av[12]);
	rsp.b2 = atof(av[13]);
	rsp.dc2 = atof(av[14]);

	if(max_x <= 0 || op_inc_disp <= 0)
	{
		fprintf(stderr,"max displacement and output frequency (op_inc_disp) need to be > 0\n");
		exit(1);
	}


        for(i=15;i<ac;i++)
        {
          if(*(av[i]) == '-')
          {

                switch(*(av[i]+1))
                {
                        case 't':
				tim_out = 'y';
                                strcpy(outfile,av[1]); strcat(outfile,".tim"); tim_file = fopen(outfile, "w");
                                break;

                        case 'v':
                                vel_out = 'y';
                                strcpy(outfile,av[1]); strcat(outfile,".vel"); vel_file = fopen(outfile, "w");
                                break;

                        case 's':
                                state_out = 'y';
                                strcpy(outfile,av[1]); strcat(outfile,".sta"); state_file = fopen(outfile, "w");
                                break;
 
                        case 'd':
                                sd_out = 'y';
                                strcpy(outfile,av[1]); strcat(outfile,".sd"); sd_file = fopen(outfile, "w");
                                break;

                        case 'f':
                                porosity_time_out = 'y';
                                strcpy(outfile,av[1]); strcat(outfile,".t_phi"); t_phi_file = fopen(outfile, "w");
                                break;

                        case 'p':
                                stop_at_peak = 'y';
                                break;

                        case 'r':
                                set_ref_fric = 'y';
                                break;

                        case 'm':
                                do_mvs = 'y';
				i++;	/*so that we can start processing immediately*/
				if( (ac - i) % 2 )	/*check to make sure we have an even number of tokens left*/
				{
					fprintf(stderr,"m option requires an even number of tokens after it --pairs of vels and times \n");
					exit(1);
				}
				if(atof(av[i]) < 0)
					doing_ramp = 'y';
                                break;
                }
          }
	  if(do_mvs == 'y')	
	  {
				/*if there's a hold and also multiple velocity steps I assume the post-hold re-load velocity is the first mvs; otherwise we'd need a duration for the std. reload..*/
	   rsp.vel_list[vstep] = atof(av[i]);	/*new loading vel*/
	   rsp.vel_dur_list[vstep] = atof(av[i+1]);	/*duration in seconds*/
           /*fprintf(stderr,"i=%d, ac=%d, vstep %d, nsteps = %d, v=%g, t=%g\n",i,ac,vstep,nsteps, rsp.vel_list[vstep], rsp.vel_dur_list[vstep]);*/
	   i+2 == ac ?  i=i+2 : i++;	/*because we start right away after case m is read*/
	   vstep++;
	   nsteps++;
  	  }	

        }
	if(do_mvs == 'y')		/*to account for last increment when command line is finished */
		nsteps--;

/*work out the six cases		--vel_list[1] is always zero, for a hold, but we'll skip it if there's no hold
1: simple shs: nsteps will be 2, rsp.vel_list will be v0-0-v2 and durations will be 0-hold_time-max_x/v2
2: simple v step no hold: nsteps will be 1, rsp.vel_list will be v0-v2 and durations will be 0-0-max_x/v2
3: no hold but multiple velocity steps: nsteps will be N, rsp.vel_list will be v0-v1-v2-vN and durations will be 0-dur1-dur2-durN
4: hold followed by mvs: nsteps will be N+1, rsp.vel_list will be v0-0-v2-vN and durations will be 0-hold_time-dur2-durN
5: no hold but velocity step as a ramp: nsteps will be 100, rsp.vel_list will be v0-ramp to V100 over ramp duration 
6: hold followed by velocity step as a ramp: nsteps will be 101, rsp.vel_list will be v0-hold_time-ramp to V100 over ramp duration

duration of the first vel is always zero --the step happens immediately.
*/

/* rsp.vel_list[0] is init vel and rsp.vel_list[1] is the 'final' vel if no hold, else [2] is the reload vel. after the hold */


	if(doing_ramp == 'n')  /*Cases 1, 2, 3 & 4 --note that rsp.vel_list[] has been set above*/
	{
		for(i = 0; i<=nsteps; i++)
		{
        		sprintf(string,"v[%d]=%g, %g s;  ",i, rsp.vel_list[i], rsp.vel_dur_list[i]);
        		strcat(string2,string);
		
			if(rsp.dc1/2 < op_inc_disp*rsp.vel_list[i])
				fprintf(stderr,"Note that output time increment is only %g, so for vel=%g µm/s you write data every %g µm, which is < Dc/2. Consider reducing your op time increment\n\n",op_inc_disp, rsp.vel_list[i], op_inc_disp*rsp.vel_list[i]); 
		}
	}
	else 			/* 5 & 6; set up ramp and print a few things*/
	{
		if(doing_hold == 'y')		
		{
			vstep=2;
			rsp.vel_list[2] *= -1.0;	/*change the sign*/
			temp = rsp.vel_dur_list[2];
			temp2 = rsp.vel_list[0] * (1/(30*log(hold_time)));	/*crude guess at final velocity after hold*/
			step_size = (rsp.vel_list[2] - temp2)/100 ;	
			rsp.vel_list[2] =  temp2+step_size; 
		}
		else				/*v step from initial velocity */
		{
			vstep=1;
			rsp.vel_list[1] *= -1.0;	/*change the sign*/
			temp = rsp.vel_dur_list[1];
			step_size = (rsp.vel_list[1] - rsp.vel_list[0])/100 ;	/*step size for ramp*/
			rsp.vel_list[1] =  rsp.vel_list[0]+step_size; 
		}
		step_duration = rsp.vel_dur_list[vstep]/100 ;	/*duration of each increment*/
		rsp.vel_dur_list[vstep] =  step_duration;  
	        /*fprintf(stderr,"nsteps=%d, step_duration is %g, step_size is %g, v[0] is %g, v[1] is %g doing a vel 'step' as a ramp to %g over %g sec  \n",nsteps,step_duration, step_size, rsp.vel_list[0],rsp.vel_list[1],rsp.vel_list[nsteps], temp);*/

		nsteps += 99;
		for(i = vstep+1; i<=nsteps; i++)
		{
			rsp.vel_list[i] =  rsp.vel_list[i-1] + step_size;  
			rsp.vel_dur_list[i] =  step_duration;  
		}
								/*for a ramp, take the total calc time from max_x */
		nsteps++;					/*add one more segment*/
		rsp.vel_list[nsteps] =  rsp.vel_list[nsteps-1] ;  
 		rsp.vel_dur_list[nsteps] = max_x/rsp.vel_list[nsteps-1];   		/*use the final velocity*/

/*fprintf(stderr,"max_x = %g, last vel = %g (%g s) and 2nd last =%g (%g s)\n",max_x,rsp.vel_list[nsteps], rsp.vel_dur_list[nsteps], rsp.vel_list[nsteps-1],rsp.vel_dur_list[nsteps-1]);*/

		if(doing_hold == 'y')		
	        	sprintf(string2,"v[0] is %g, hold for %g s followed by a v_step as a ramp from %g (end-of-hold_v) to %g over %g sec  ",rsp.vel_list[0],hold_time,temp2, rsp.vel_list[nsteps], temp);
		else
	        	sprintf(string2,"v[0] is %g, v[1] is %g doing a vel 'step' as a ramp to %g over %g sec  ",rsp.vel_list[0],rsp.vel_list[1],rsp.vel_list[nsteps-1], temp);
	}
	fprintf(stderr,"%s\n",string2);

/*set up array space and times for op. Get array index values for each v step*/

	calc_time = 0;  /*if there's a hold then hold_time is rsp.vel_dur_list[1], otherwise that's the dur of the first v step*/

	for(i = 1; i<=nsteps; i++)			/*this could start at 0 but vel_dur_list[0] is always zero*/
	{
		calc_time += rsp.vel_dur_list[i];
                if(i==1 && doing_hold == 'y')
			ndata = (int)(ceil(hold_time/op_inc_hold));		
		else	
			ndata += (int)(ceil(rsp.vel_dur_list[i]/op_inc_disp));

	}
/*	ndata++;*/
	if(doing_ramp == 'y' && doing_hold == 'y')	
		ndata--;

					/*array space*/
	rsp.disp_data = (double *)calloc((unsigned)(ndata), 8) ;	
        rsp.model_mu  = (double *)calloc((unsigned)(ndata), 8) ;
        v_slider  = (double *)calloc((unsigned)(ndata), 8) ;
        state  = (double *)calloc((unsigned)(ndata), 8) ;
        sd  = (double *)calloc((unsigned)(ndata), 8) ;
	phi = (double *)calloc((unsigned)(ndata), 8) ;

	rsp.model_mu[0] = rsp.muo;

	/*set up times for op. this is called 'disp_data' for historical reasons... should be time*/

	rsp.vs_row_list[0] = 0;		/*first step is at row 0*/
	rsp.disp_data[0] = 0;	/*first op point is at t=0*/

	sprintf(string1,"");	/*clear string1*/

	data_count=0;
	i=1;
	calc_time = 0;
	for(j = 1; j<=nsteps; j++)
	{
		temp_time = 0;
		if(j==1 && doing_hold == 'y')
		{
				/*for a hold of 100 sec and op_inc=0.3 there are 333.33 steps, so step is at 334*/
			while( (temp_time+op_inc_hold+tinyAmt) < hold_time)	/*hold_time is zero if there's no hold*/
                	{
                        	temp_time += op_inc_hold;
				calc_time  += op_inc_hold;
                        	rsp.disp_data[i++] = calc_time;
                	}
			temp_time += (hold_time - rsp.disp_data[i-1]);		/*deal with any remainder*/
			calc_time += (hold_time - rsp.disp_data[i-1]);		/*deal with any remainder*/
		}
		else
		{
			while( (temp_time+op_inc_disp+tinyAmt) < rsp.vel_dur_list[j])	/*hold_time is zero if there's no hold*/
			{
				temp_time += op_inc_disp;
				calc_time  += op_inc_disp;
                        	rsp.disp_data[i++] = calc_time;
/*fprintf(stderr,"time=%g at i=%d\n",rsp.disp_data[i-1],i-1);*/
			}
			calc_time += rsp.vel_dur_list[j] - temp_time;	/*deal with any remainder*/
		}
		data_count = i;
                rsp.disp_data[i++] = calc_time;
		rsp.vs_row_list[j] = data_count;		/*save the point for each step */
		/*if(j<nsteps)
		{	
        		sprintf(string,"nsteps=%d; At array point %d step to: v[%d]=%g, %g (s);  ",nsteps,rsp.vs_row_list[j],j+1, rsp.vel_list[j+1], rsp.vel_dur_list[j+1]);
        		strcat(string1,string);
		}*/
	}

/*	if(doing_ramp == 'n')
		fprintf(stderr,"%s\n",string1);
	fprintf(stderr,"end: nsteps=%d, ndata = %d; data_count=%d, calc_time = %g, temp_time = %g, i=%d, j=%d, rsp.disp_data[2]=%g, rsp.disp_data[ndata]=%g\n",nsteps,ndata,data_count,calc_time,temp_time,i,j,rsp.disp_data[2],rsp.disp_data[ndata-1]);
*/

/*j=1;
for(i=0;i<ndata;i++)*/

/*
i=0;
for(j = 1; j<=nsteps; j++)
{
  while(i<rsp.vs_row_list[j])
  {
  	i++;
  	fprintf(stderr,"i=%d time=%f mu=%f vel=%g vel_dur=%g step=%d j=%d\n",i,rsp.disp_data[i],rsp.model_mu[i],rsp.vel_list[j],rsp.vel_dur_list[j], rsp.vs_row_list[j],j);
  }
}
*/

/*set up output files, ref. values and get ready to do the calculations*/

        strcpy(outfile,av[1]);
	strcat(outfile,".dis");
	dist_file = fopen(outfile, "w");
	if(segall_rice_law == 'y')	/*write file for lp_disp vs. porosity if using Segall-Rice laws*/
	{
        	strcpy(outfile,av[1]);
		strcat(outfile,".phi");
		x_phi_file = fopen(outfile, "w");
	}

	if(rsp.dc2 < 0)		/*one sv fit*/
        {
		rsp.one_sv_flag=TRUE;
                rsp.b2 = 0.000000000000000000000000000000000000000000000;
                rsp.dc2 = 1e90;
        	rsp.vr_dc2 = 1e-90;
        }
	else
		rsp.one_sv_flag=FALSE;

	rsp.amb = rsp.a - (rsp.b1 + rsp.b2);
	
				/*set V_REF, MU_REF*/
	if(set_ref_fric == 'n')
	{
		rsp.v_ref = 1e-1;	/* Reference mu_ref to mu = 0.6 at .1 mic/sec */
		rsp.muo = 0.6;
		rsp.mu_ref = rsp.muo - rsp.amb*log(rsp.vel_list[0]/rsp.v_ref);
        	fprintf(stderr,"\n%s\n\tmu_ref=%f, V_REF=%g, mu_init=%g, (mu_ref set to give mu=0.6 at v_init)\n",av[1],rsp.mu_ref,rsp.v_ref,rsp.muo);
	}
	else
	{
		rsp.v_ref = 1.0;	/* Reference mu_ref to mu = 0.6 at 1 mic/sec */
		rsp.mu_ref = 0.6;
		rsp.muo = rsp.mu_ref + rsp.amb*log(rsp.vel_list[0]/rsp.v_ref);
        	fprintf(stderr,"\n%s\n\tmu_ref=%f, V_REF=%g, mu_init=%g, (mu_ref set to give mu=0.6 at 1mic/s)\n",av[1],rsp.mu_ref,rsp.v_ref,rsp.muo);
	}


	rsp.lin_term=0;
	rsp.vr_dc1 = rsp.v_ref / rsp.dc1;
	if(!rsp.one_sv_flag)
	{
		rsp.vr_dc2 = rsp.v_ref / rsp.dc2;
                sprintf(string, "a=%g  b1=%g  dc1=%g  b2=%g dc2=%g law=%c k=%g Vo=%g\n",rsp.a, rsp.b1, rsp.dc1, rsp.b2, rsp.dc2, rsp.law, rsp.stiff, rsp.vel_list[0]); 
	}
	else
                sprintf(string, "a=%g  b1=%g  Dc=%g law=%c k=%g Vo=%g\n",rsp.a, rsp.b1, rsp.dc1, rsp.law, rsp.stiff, rsp.vel_list[0]);

	strcat(string,string2);
        fprintf(dist_file,"Load Point Displacement (µm), Friction, %s\n",string);
	if(segall_rice_law == 'y') fprintf(x_phi_file,"%s\n",string);	
        if(tim_out == 'y') fprintf(tim_file,"Time (s), Friction, %s\n",string);
        if(vel_out == 'y') fprintf(vel_file,"ln(V/Vo), Friction, %s\n",string);
        if(state_out == 'y') fprintf(state_file,"Time (s), Friction State Var,%s\n",string);
        if(sd_out == 'y') fprintf(sd_file,"Fault Slip (µm), Friction,%s\n",string);
        if(porosity_time_out == 'y') fprintf(t_phi_file,"Time (s), Porosity%s\n",string);

        if( (get_mu_at_xt(rsp.disp_data, rsp.model_mu, v_slider, &ndata, &rsp, outfile, state, sd, phi, nsteps, doing_hold)) != 0) 
        {
                fprintf(stderr,"%cerror in calc., op file is not complete, ndata=%i\n\n",0x7,ndata);  
        }

					/*write data to op files*/
	if(tim_out == 'y') 		
        {
		fprintf(tim_file,"-10.000\t%.9f\n",rsp.muo);
		fprintf(tim_file,"0.0000\t%.9f\n",rsp.muo);
	}
	if(vel_out == 'y')
		fprintf(vel_file,"0.0000\t%.9f\n",rsp.muo);

	if(state_out == 'y')
	{
		fprintf(state_file,"-10.000\t%.9f\n",state[1]);
		fprintf(state_file,"0.0000\t%.9f\n",state[1]);
	}
	if(sd_out == 'y')
	{
		fprintf(sd_file,"-10.000\t%.9f\n",rsp.model_mu[0]);
		fprintf(sd_file,"%.9f\t%.9f\n",sd[0],rsp.model_mu[0]);
	}
	if(porosity_time_out == 'y')
	{
		fprintf(t_phi_file,"-10.000\t%.9f\n",phi[1]);
		fprintf(t_phi_file,"0.0000\t%.9f\n",phi[1]);
	}

	if( segall_rice_law == 'y')
	{
		fprintf(x_phi_file,"-10.000\t%.9f\n",phi[1]);
		fprintf(x_phi_file,"%.9f\t%.9f\n",rsp.disp_data[0],phi[1]);
	}

	fprintf(dist_file,"-10.000\t%.9f\n",rsp.model_mu[0]);
	fprintf(dist_file,"%.9f\t%.9f\n",rsp.disp_data[0],rsp.model_mu[0]);

j=0;
for(i=0;i<ndata;i++)
{
  if(i>=rsp.vs_row_list[j])
	j++;
  /*fprintf(stderr,"i=%d,x=%f mu=%f vel=%g vel_dur=%g step=%d j=%d\n",i,rsp.disp_data[i],rsp.model_mu[i],rsp.vel_list[j],rsp.vel_dur_list[j], rsp.vs_row_list[j],j);*/
}
	lpDisp=0;
	for(i=1, j=1;i<ndata;i++)
	{
		if(i == rsp.vs_row_list[j])	/*increment velocity*/
			j++;
		lpDisp += (rsp.disp_data[i]-rsp.disp_data[i-1])*rsp.vel_list[j];
		fprintf(dist_file,"%g\t%g\n",lpDisp, rsp.model_mu[i]);
		/*fprintf(dist_file,"%g\t%g\t%d\n",rsp.vel_list[j], rsp.model_mu[i],rsp.vs_row_list[j]);*/
		if(segall_rice_law == 'y')		/*lp_disp vs. phi*/
			fprintf(x_phi_file,"0.000\t%.9f\n",phi[i]);
		if(tim_out == 'y') 
			fprintf(tim_file,"%.9f\t%.9f\n",rsp.disp_data[i],rsp.model_mu[i]);
		if(state_out == 'y') 
			fprintf(state_file,"%.9f\t%.9f\n",rsp.disp_data[i],state[i]);
		if(porosity_time_out == 'y') 
			fprintf(t_phi_file,"%.9f\t%.9f\n",rsp.disp_data[i],phi[i]);

		if(vel_out == 'y')
			fprintf(vel_file,"%.9f\t%.9f\n",log(v_slider[i]/rsp.vel_list[0]),rsp.model_mu[i]);
		if(sd_out == 'y') 
			fprintf(sd_file,"%.9f\t%.9f\n",sd[i],rsp.model_mu[i]);
	}
	
fprintf(stderr,"\n");
exit(0);
}

/* get_mu_func */

/**************************************************************************/

#define DISP_TOL 1e-3   /* disp changed to time in Jan 2022 update; I assume this is in sec */	
#define TWO      2.0
#define SMALL_NUM 1e-90
#define EPSILON  1e-13
#define PGROW   -0.200
#define PSHRINK -0.250
#define FCOR    1.0/15.0
#define SAFETY  0.9
#define ERRCON  6.0e-4
#define VPEAK_WARN_THRES	0.01
#define MUPEAK_WARN_THRES	1e-4

/* This function calculates friction at the displacements in x[] for a few
Rate/State friction law . The result is put in mod_mu[].

  constit. law
                (1)     mu = mu_ref + a ln(v/v_ref) + b ln(v_ref*theta/Dc)
                (2)     d_theta/dt = -(v*theta/Dc) *ln(v*theta/Dc)
  elastic coupling
                (3)     d_mu/dt = k (v_lp - v)

 calculation involves solving (1) for v and subing this into (2), so that:

                (4)     d_theta/dt = (alpha*v_ref*theta/Dc) * ln(alpha*v_ref*theta/Dc)
 where we define:
        alpha = exp[(mu - mu_ref - b*ln(v_ref*theta/Dc)/a]
 or for two sv's alpha = exp[(mu - mu_ref - b1*ln(v_ref*theta1/Dc1)-b2*ln(v_ref*theta2/Dc1))/a]

 Then the solution involves solving (4) with a rewritten version of (3):

                (5)     d_mu/dt = K/sigma_n * (v_lp - alpha/v_ref)

porosity:  (via eq'n 17 or Segall and Rice)
		 	phi = phi_ref - epsilon ln(v_ref*theta/Dc)
for the "s" law, this is just used as a rule to calculate phi from a given theta

        x is the ndata vector of times (used to be displacements) at which mod_mu[i] is needed

        rsp contains things like stiffness etc. It's defined in rsp.h

15/7/94: modified to calculate other rate_state laws
24/4/96: modified to handle multiple velocity steps and to use a reference velocity
29/5/96: option added to compute 2 state variable response 

30/3/99: added porosity calculation, using eq'n 17 of Segall and Rice, 1995
2022.01.03 modified to handle multiple velocity steps
*/

int get_mu_at_xt( double *x, double *mod_mu, double *v_slider, int *ndata, struct rs_parameters *rsp, char *outfile, double *state, double *sd, double *phi, int nsteps, char doing_hold)
/*note that x[] is an array of times; based on Jan 2022 update*/
{
 char past_first_fric_peak = FALSE, hold_over = FALSE;
 double h, hh, next_h, hdid, tnow, ttnd, v;
 double last_mus, mu_s, p_mu, old_mu, psi1, p_psi1, old_psi1, psi2, p_psi2, old_psi2;
 double slip_at_end_hold, slip_at_peak, lpd_at_peak, disp_slider=0, disp=0;	/*disp is load point disp*/
 double mu_err_scale, psi_err_scale, err, max_err;
 double mu_err, psi1_err, psi2_err; 
 double fabs(), log(), exp();
 double my_const[5];			/*for rk calc*/
 double v_s, p_v_s, time=0, calc_time, temp_time;
 double maxSlipVel = 0;
 int i, j, h_err;
 int do_rk( double *mu, double *psi1, double *psi2, double *H, double *v, struct rs_parameters *rsp, double *my_const);
 double	mu_min, mu_max, mu_peak, del_mu;		/* del_mu is peak - base friction */
 double	theta_hold, theta_peak;		/* get from "inside" calc. loop, assume they correspond to min (max) fric*/
 double v_peak;
 double epsilon = 1.7e-4;	/*this is Segall and Rice's value*/ 
 double phi_ref = 0.2;		/*assume 20% porosity, but this is arbitrary*/

 my_const[0]=my_const[1]=0.0;
 my_const[2]=my_const[3]=0.5;
 my_const[4]=1.0;

 v_s = rsp->vel_list[0];			/*initial velocity*/
 mu_peak = mu_max = 0;
 mu_min=mu_s = rsp->muo; 		/* initial mu, from caller*/

 mu_err_scale = rsp->muo;		/* use constant fractional errors*/
 psi_err_scale = rsp->dc1/v_s;        


 switch(rsp->law)				/* set up initial value of state (note I use psi for theta)*/
 {						/* 	assume steady state for initial v*/
        case 'o':
                psi1 = psi2 = -log(v_s/rsp->v_ref);
 		psi_err_scale = fabs(log(v_s/rsp->v_ref));
                break;
        case 'r':
        case 'd':
                psi1 = rsp->dc1/v_s;        /*this will be the same for theta version of r and d, note that  */
                psi2 = rsp->dc2/v_s;	    /*  s and S will also use same state */
                break;
        case 'p':
 		rsp->vr_dc1 = rsp->v_ref/(TWO*rsp->dc1);		/*modify for Perrin-Rice law*/
 		rsp->vr_dc2 = rsp->v_ref/(TWO*rsp->dc2);		/*note: this is set above as rsp.vr_dc..*/
                psi1 = (TWO*rsp->dc1)/v_s;
                psi2 = (TWO*rsp->dc2)/v_s;
                break;
        case 'j':                                       /*check this later, just a guess now 15/7/94 */
		fprintf(stderr,"rice law not implemented\n");
		return(-1);
                break;
 }
 if(rsp->one_sv_flag)
 	rsp->vr_dc2 = SMALL_NUM;  /*b2 is zero for 1sv model, set dc2 to very large number*/

 for(i=0;i<2;i++)	/*fill with prejump mu*/
 {
	mod_mu[i] = mu_s;
	v_slider[i] = v_s;
	state[i] = psi1;
	sd[i]=0.0;
	phi[i] = phi_ref - epsilon*log(psi1*rsp->vr_dc1);
 }

/*start main calc loop*/
 	/*there are nstep lp velocities in rsp->vel_dur_list[], each one has a duration: rsp->vel_dur_list[]*/
	
 		/* first point in x-array is the vel step. we'll write  the initial state for reference */ 
 next_h  = x[1];

 time=0;

 if(doing_hold == 'n')
	hold_over = TRUE;

 j=1;		/*this is the counter for rsp->vel_list[]*/

 i=1;
 rsp->v_lp = rsp->vel_list[1];
 /*fprintf(stderr,"v_lp=%g, vel_list=%g rsp->row_list[1]=%d rsp->row_list[2]=%d\n",rsp->v_lp,rsp->vel_list[1], rsp->vs_row_list[1], rsp->vs_row_list[2]);*/

 while(i < *ndata)
 {
   if(hold_over == FALSE && time >= rsp->vel_dur_list[1])		/*if we're doing a hold*/
   {
	hold_over = TRUE;
 	mu_peak  = 0;					/*reset, so that we get max mu post-hold*/
	slip_at_end_hold = disp_slider;
	fprintf(stderr,"end of hold: mu is %f theta is %f v_slider is %f time is %f\n",mu_s, psi1, v_slider[i-1], time);
	next_h /=100;		/*force it to start small for reload*/
   }

   if(i >= rsp->vs_row_list[j])            
   {
   	rsp->v_lp=rsp->vel_list[++j];            /*new load point velocity*/
 	/*fprintf(stderr,"hold_flag=%c i=%d, time=%g,  v_lp=%g, vel_list=%g\n",hold_over,i,x[i],rsp->v_lp,rsp->vel_list[j]);*/
   }

   ttnd = x[i]-x[i-1];		/* rel time to next data*/

   tnow=0; 					/* time since last data*/ 

   while(tnow < ttnd)
   {
     h_err=TRUE;

     if(tnow+next_h >= ttnd)			/*set time*/
	h = ttnd-tnow;
     else
	h = next_h;

     while(h_err)
     {
			/* save old values */
	old_mu = p_mu = mu_s;
	old_psi1 = p_psi1 = psi1;
	old_psi2 = p_psi2 = psi2;

	hh = h/TWO;				/* half-steps*/

	
        if(     (do_rk(&mu_s, &psi1, &psi2, &hh, &v_s, rsp, my_const) == -1) ||  /*1st half step*/
                (do_rk(&mu_s, &psi1, &psi2, &hh, &v_s, rsp, my_const) == -1) ||  /*2nd half step, using updated vals*/
                (do_rk(&p_mu, &p_psi1, &p_psi2, &h, &p_v_s, rsp, my_const) == -1)    /* full step, updated vals in p_ */
          )
	{
		fprintf(stderr,"\terror from do_rk, row=%d, exiting early\n\n",i);
		fprintf(stderr,"hold_time, del_mu, v_init, v_rl, del_mu_min, mu_peak, mu_max, theta_hold, theta_peak, v_peak, a, b1, dc1, b2, dc2, k, file, law, lp_disp_peak, slip(peak-eo_hold), slip_at_end_hold\n");
		fprintf(stderr,"%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %s, %c, %g, %g, %g\n\n",rsp->vel_dur_list[1], mu_peak-rsp->muo,rsp->vel_list[0],rsp->vel_list[2],rsp->muo-mu_min,mu_peak,mu_max,theta_hold, theta_peak, v_peak, rsp->a,rsp->b1,rsp->dc1,(rsp->one_sv_flag) ? 0 : rsp->b2,(rsp->one_sv_flag) ? 0 : rsp->dc2,rsp->stiff,outfile,rsp->law, lpd_at_peak, slip_at_peak-slip_at_end_hold, slip_at_end_hold);
		printf("\twarning: error from do_rk, row=%d, we exited early, hold was %s and %s peak friction, peak vel = %f.  Check op\n",i-1, hold_over == TRUE ? "over" : "not over", past_first_fric_peak == TRUE ? "past" : "possibly not past", v_peak);
		printf("%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %s, %c, %g, %g, %g\n",rsp->vel_dur_list[1], mu_peak-rsp->muo,rsp->vel_list[0],rsp->vel_list[2],rsp->muo-mu_min,mu_peak,mu_max,theta_hold, theta_peak, v_peak, rsp->a,rsp->b1,rsp->dc1,(rsp->one_sv_flag) ? 0 : rsp->b2,(rsp->one_sv_flag) ? 0 : rsp->dc2,rsp->stiff,outfile,rsp->law, lpd_at_peak, slip_at_peak-slip_at_end_hold, slip_at_end_hold);
		*ndata = i-1;
		return(-1);
	}	
	if(!h)					/*time step too small? */
	{
		fprintf(stderr,"\terror from do_rk, time_step = 0, row=%d, exiting early\t",i);
		return(-1);
	}

	max_err=0;					 /* evaluate error*/
	mu_err = mu_s-p_mu; 
	err = fabs(mu_err/mu_err_scale);
	max_err = (err > max_err) ? err : max_err;

	psi1_err = psi1-p_psi1;
	err = fabs(psi1_err/psi_err_scale);
	max_err = (err > max_err) ? err : max_err;

	psi2_err = psi2-p_psi2;
	err = fabs(psi2_err/psi_err_scale);
	max_err = (err > max_err) ? err : max_err;

	max_err /= EPSILON;				/* scale relative to tolerance*/
	
	if(max_err <=1.0)				/*step succeeded, compute size of next step*/
	{
		hdid = h;
		next_h	= (max_err > ERRCON) ? SAFETY*h*exp(PGROW*log(max_err)) : 4.0*h;
		h_err=FALSE;
	}
	else
	{
		h = SAFETY*h*exp(PSHRINK*log(max_err));	/*truncation error too large, reduce step size*/
		mu_s   = old_mu;
                psi1 = old_psi1;
                psi2 = old_psi2;
	}
     }			/* loop for adaptive step size control*/
     tnow += hdid;
     time += hdid;	/*cumulative time*/
     mu_s += mu_err*FCOR;					/* fifth order bit */
     psi1 += psi1_err*FCOR;
     psi2 += psi2_err*FCOR;
     disp += rsp->v_lp*hdid;

     switch(rsp->law)		/*update slider velocity*/
     {
	case 'o':					/*old, original Runia slip law*/
   	 v_s = rsp->v_ref * exp((mu_s - rsp->mu_ref - (rsp->b1*psi1) - (rsp->b2*psi2))/rsp->a);
	 if (v_s > maxSlipVel) maxSlipVel = v_s;
	 break;

	default:
   	 v_s = rsp->v_ref * exp((mu_s - rsp->mu_ref - (rsp->b1*log(rsp->vr_dc1 * psi1)) - (rsp->b2*log(rsp->vr_dc2 * psi2)))/rsp->a);
	 if (v_s > maxSlipVel) maxSlipVel = v_s;
	 break;
     }
     disp_slider += v_s*hdid;			/*slider disp*/

     if(!hold_over) 				/*save mu min (occurs at the end of the hold, */
     {						/* use flag to indicate when the hold is finished --to avoid later oscillations*/
	mu_min = mu_s; 
	theta_hold = psi1;
     }
     if(mu_s > mu_max)			/*save highest mu, even if not "mu_peak"*/
	mu_max = mu_s;
     if(hold_over && !past_first_fric_peak) 
     {
        /*if(v_s >= rsp->vel_list[0] || mu_s <= last_mus)*/
        if(mu_s < last_mus)
	{
		fprintf(stderr,"past peak friction, v_peak=%g  del_mu=%f mu_peak=%f theta_peak is %f\n",
						v_peak, mu_peak-rsp->muo, mu_peak, theta_peak);
		past_first_fric_peak = TRUE;		

		if(stop_at_peak == 'y')	/*set to special value to stop calc*/
			ttnd = -1;
	}
	else
	{
		mu_peak = mu_s; 
		theta_peak = psi1;
		v_peak = v_s;
		lpd_at_peak = disp;
		slip_at_peak = disp_slider;
	}
     }
     last_mus = mu_s;		/*save for comparison to find peak*/
   }			/* loop for time to next data  */

/* fprintf(stderr,"i=%d, time=%g,  v_lp=%g, vel_list=%g, mu=%g\n",i,x[i],rsp->v_lp,rsp->vel_list[j],mu_s);*/
   if(fabs(time-x[i]) > DISP_TOL )			/*time section error checking*/
   {
   	fprintf(stderr,"\tget_mu_at_xt: error, time not correct, row=%d, array_time=%f, numer_time =%f, vel=%g, j=%dexiting early\t",i,x[i],time,rsp->vel_list[j],j);
   	return(-1);
   }
								/*save (write) new values for op files*/
   mod_mu[i] = mu_s;
   state[i] = psi1;
   sd[i]=disp_slider;
   phi[i] = phi_ref - epsilon*log(rsp->vr_dc1 * psi1);

/*fprintf(stderr,"i=%d,mu=%f\n",i,mod_mu[i]);*/

				/*use updated v_s based on 5th order est. of other params*/
   switch(rsp->law)
   {
	case 'o':					/*old, original Runia slip law*/
   	 v_slider[i] = rsp->v_ref * exp((mu_s - rsp->mu_ref - (rsp->b1*psi1) - (rsp->b2*psi2))/rsp->a);
	 break;

	default:
   	 v_slider[i] = rsp->v_ref * exp((mu_s - rsp->mu_ref - (rsp->b1*log(rsp->vr_dc1 * psi1)) - (rsp->b2*log(rsp->vr_dc2 * psi2)))/rsp->a);
	 break;
   }
   i++;		/*increment data counter*/

   if(ttnd < 0)	/*always >= 0, unless the stop_at_peak flag is on and we've just reached the peak*/
   {
	*ndata = i-1;	
	break;   /*stop calc*/
   }

 }		/* end of loop for data points */

						/*write some things to stdout (maybe a log file)*/
/*hold time = rsp->vel_dur_list[1]*/

fprintf(stderr,"hold_time, del_mu, v_init, v_rl, del_mu_min, mu_peak, mu_max, theta_hold, theta_peak, v_peak, a, b1, dc1, b2, dc2, k, file, law, lp_disp_peak, slip(peak-eo_hold), slip_at_end_hold, maxSlipvelocity\n");
fprintf(stderr,"%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %s, %c, %g, %g, %g, %g\n",rsp->vel_dur_list[1], mu_peak-rsp->muo,rsp->vel_list[0],rsp->vel_list[2],rsp->muo-mu_min,mu_peak,mu_max,theta_hold, theta_peak, v_peak, rsp->a,rsp->b1,rsp->dc1,(rsp->one_sv_flag) ? 0 : rsp->b2,(rsp->one_sv_flag) ? 0 : rsp->dc2,rsp->stiff,outfile,rsp->law, lpd_at_peak, slip_at_peak-slip_at_end_hold, slip_at_end_hold, maxSlipVel);
printf("%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %s, %c, %g, %g, %g, %g\n",rsp->vel_dur_list[1], mu_peak-rsp->muo,rsp->vel_list[0],rsp->vel_list[2],rsp->muo-mu_min,mu_peak,mu_max,theta_hold, theta_peak, v_peak, rsp->a,rsp->b1,rsp->dc1,(rsp->one_sv_flag) ? 0 : rsp->b2,(rsp->one_sv_flag) ? 0 : rsp->dc2,rsp->stiff,outfile,rsp->law, lpd_at_peak, slip_at_peak-slip_at_end_hold, slip_at_end_hold, maxSlipVel);

return(0);	/* success! */
}
#undef SMALL_NUM 
#undef ONE
#undef TWO
#undef PGROW	
#undef PSHRINK
#undef FCOR
#undef SAFETY
#undef ERRCON
#undef EPSILON  
#undef DISP_TOL 

/*************************************************************************/

/*	This func. does most of the work of solving the coupled
equations for rate/state friction and elastic interaction. It's called 
from get_mu_at_xt(). 
*/
#define	BIGNUM 1e30
#define ONE	1.0
#define TWO	2.0
#define FIVE	5.0
#define TEN	10.0
#define P16	ONE/6.0

/*      these are defined in get_mu_at_xt.c
 my_const[0]=my_const[1]=0.0;
 my_const[2]=my_const[3]=0.5;
 my_const[4]=1.0;
*/

int do_rk( double *mu, double *psi1, double *psi2, double *H, double *v, struct rs_parameters *rsp, double *my_const)
{

int	calc_bombed=1;
int 	i;
/*int	isnan();*/
double	alpha,arg;
double	J[5], K[5], M[5];
double	w_psi1, w_psi2, w_mu;
double	old_psi1, old_psi2, old_mu, old_H, old_v;
double	exp(), log();

 while(calc_bombed)
 {
	old_psi1 = *psi1;             /* save last vals in case calc bombs ...*/
	old_psi2 = *psi2;
	old_mu = *mu;
	old_H = *H; 
	old_v = *v; 
	J[0]=K[0]=M[0]=0;	
							/* Runga Kutta calc */
	switch(rsp->law)
	{
	  case 'o':					/*old, original Runia slip law*/	
		for(i=1;i<5;++i)	
		{
			w_psi1 =	*psi1 + J[i-1]*my_const[i];
			w_psi2 =	*psi2 + K[i-1]*my_const[i];
			w_mu = 		*mu   + M[i-1]*my_const[i];
			alpha = exp((w_mu - rsp->mu_ref - rsp->b1*w_psi1 - rsp->b2*w_psi2)/rsp->a);
           		J[i] = (*H) * -rsp->vr_dc1 * alpha * (w_psi1 + log(alpha));
           		K[i] = (*H) * -rsp->vr_dc2 * alpha * (w_psi2 + log(alpha));
           		M[i] = (*H) * rsp->stiff * (rsp->v_lp - rsp->v_ref*alpha)  ;
		} 
		*psi1 +=  (J[1] + TWO*J[2] + TWO*J[3] + J[4])*P16; 
		*psi2 +=  (K[1] + TWO*K[2] + TWO*K[3] + K[4])*P16;
		*mu   +=  (M[1] + TWO*M[2] + TWO*M[3] + M[4])*P16;

		*v = rsp->v_ref * exp((*mu - rsp->mu_ref - (rsp->b1*(*psi1)) - (rsp->b2*(*psi2)))/rsp->a);
		break;
	
	  case 'r':			/*Ruina, slip law, using theta as sv*/
		for(i=1;i<5;++i)	
		{
			w_psi1 =	*psi1 + J[i-1]*my_const[i];
			w_psi2 =	*psi2 + K[i-1]*my_const[i];
			w_mu = 		*mu   + M[i-1]*my_const[i];
			alpha = exp((w_mu - rsp->mu_ref - rsp->b1*log(rsp->vr_dc1*w_psi1) - rsp->b2*log(rsp->vr_dc2*w_psi2) )/rsp->a);
			arg = alpha * rsp->vr_dc1 * w_psi1;
           		J[i] = (*H) * -arg * log(arg);
			arg = alpha * rsp->vr_dc2 * w_psi2;
           		K[i] = (*H) * -arg * log(arg);
           		M[i] = (*H) * rsp->stiff * (rsp->v_lp - rsp->v_ref*alpha);
		} 
		*psi1 +=  (J[1] + TWO*J[2] + TWO*J[3] + J[4])*P16; 
		*psi2 +=  (K[1] + TWO*K[2] + TWO*K[3] + K[4])*P16;
		*mu   +=  (M[1] + TWO*M[2] + TWO*M[3] + M[4])*P16;

		*v = ( rsp->v_ref * exp((*mu - rsp->mu_ref - (rsp->b1*log(rsp->vr_dc1 * *psi1)) - (rsp->b2*log(rsp->vr_dc2 * *psi2)))/rsp->a));
		break;

	  case 'd':			/*Dieterich, slowness law*/
		for(i=1;i<5;++i)	
		{
			w_psi1 =	*psi1 + J[i-1]*my_const[i];
			w_psi2 =	*psi2 + K[i-1]*my_const[i];
			w_mu = 		*mu   + M[i-1]*my_const[i];
			alpha = exp((w_mu - rsp->mu_ref - rsp->b1*log(rsp->vr_dc1*w_psi1) - rsp->b2*log(rsp->vr_dc2*w_psi2) )/rsp->a);
			arg = alpha * rsp->vr_dc1 * w_psi1;
           		J[i] = (*H) * (ONE - arg);
			arg = alpha * rsp->vr_dc2 * w_psi2;
           		K[i] = (*H) * (ONE - arg);
           		M[i] = (*H) * rsp->stiff * (rsp->v_lp - rsp->v_ref*alpha);
		} 
		*psi1 +=  (J[1] + TWO*J[2] + TWO*J[3] + J[4])*P16; 
		*psi2 +=  (K[1] + TWO*K[2] + TWO*K[3] + K[4])*P16;
		*mu   +=  (M[1] + TWO*M[2] + TWO*M[3] + M[4])*P16;

		*v = ( rsp->v_ref * exp((*mu - rsp->mu_ref - (rsp->b1*log(rsp->vr_dc1 * *psi1)) - (rsp->b2*log(rsp->vr_dc2 * *psi2)))/rsp->a));
		break;

	  case 'j':			/*Rice Law*/
		fprintf(stderr,"no rice law yet \n");
        	return(-1);         /* inform caller of problem*/
		break;

	  case 'p':		/*Perrin-Rice Quadratic law*/
		for(i=1;i<5;++i)	
		{
			w_psi1 =	*psi1 + J[i-1]*my_const[i];
			w_psi2 =	*psi2 + K[i-1]*my_const[i];
			w_mu = 		*mu   + M[i-1]*my_const[i];
			alpha = exp((w_mu - rsp->mu_ref - rsp->b1*log(rsp->vr_dc1*w_psi1) - rsp->b2*log(rsp->vr_dc2*w_psi2) )/rsp->a);
			arg = alpha * rsp->vr_dc1 * w_psi1;
           		J[i] = (*H) * (ONE - arg*arg);
			arg = alpha * rsp->vr_dc2 * w_psi2;
           		K[i] = (*H) * (ONE - arg*arg);
           		M[i] = (*H) * rsp->stiff * (rsp->v_lp - rsp->v_ref*alpha);
		} 
		*psi1 +=  (J[1] + TWO*J[2] + TWO*J[3] + J[4])*P16; 
		*psi2 +=  (K[1] + TWO*K[2] + TWO*K[3] + K[4])*P16;
		*mu   +=  (M[1] + TWO*M[2] + TWO*M[3] + M[4])*P16;

		*v = ( rsp->v_ref * exp((*mu - rsp->mu_ref - (rsp->b1*log(rsp->vr_dc1 * *psi1)) - (rsp->b2*log(rsp->vr_dc2 * *psi2)))/rsp->a));
		break;

	   default : 
		fprintf(stderr,"no rs law chosen. --how could this have happened?\n");
        	return(-1);         /* inform caller of problem*/
	}
	if( isnan(*v) || isnan(*mu) || (fabs(*v) > BIGNUM) )
        {
		fprintf(stderr,"do_rk: calculation bombed, retry #%d,\tcurrent parameters are:\na=%g, b1=%g, dc1=%g, b2=%g, dc2=%g,\nmu_o=%g, mu_f=%g, next time step=%g, v=%g old_H=%g, old_v=%g\n",calc_bombed,rsp->a,rsp->b1,rsp->dc1,rsp->b2,rsp->dc2,rsp->muo,rsp->mu_ref,*H/TEN,*v,old_H, old_v);
		if(++calc_bombed == 3 )          /* give up if problem persists */
		{
			fprintf(stderr,"do_rk: calculation bombed out\tcurrent parameters are:\na=%g, b1=%g, dc1=%g, b2=%g, dc2=%g,\nmu_o=%g, mu_f=%g, time step=%g, v=%g old_H=%g, old_v=%g\n",rsp->a,rsp->b1,rsp->dc1,rsp->b2,rsp->dc2,rsp->muo,rsp->mu_ref,*H,*v,old_H, old_v);
               		return(-1);         /* inform caller of problem*/
		}
                else
                {
                        *psi1 = old_psi1;        /* install old vals */
                        *psi2 = old_psi2;
                        *mu   = old_mu;
                        *H    = old_H/TEN;      /* reduce H -stablize calc? */
			*v    = old_v;
                }
	}
        else
        	calc_bombed = FALSE;               	/* all OK, stop and return */ 

 }	/* end of calc_bombed loop */
 return(0);
}
#undef ONE
#undef FIVE
#undef TEN
#undef TWO
#undef P16
#undef BIGNUM

