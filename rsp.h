/* 21 12 21  --changed velocity list to 500 and added vel_dur_list[] */

#include	<string.h>
#include        <math.h>
#include	<stdio.h>

#define	TRUE 1
#define FALSE 0
#define MAX_COL   17
#define MAX_ROW 100000
#define    BELL     '\007'
#define ALPHA   1.0             /*  reflection factor   */
#define BETA    0.5             /*  contraction factor  */
#define GAMMA   2.0             /*  expansion factor    */
#define MAX_PARAM 10            /*  max # of parameters to be optimised */
#define TWOPI 6.283185308
#define PI 3.141592654
#define TRUE 1
#define FALSE 0
#define SMALL 1e-10
#define SQR(a) (a*a)
#define WELCH_WINDOW(j,a,b) (1.0-SQR((((j)-1)-(a))*(b)))   /*see p 445 of numerical recipes in c*/
#define         HALF  0.5000
#define         ROUND(A)  ( ( (A-floor(A) ) > HALF ) ? ceil(A)  : floor(A) )

struct rs_parameters
     {
	char	law;			/*d=Dieterich, r=Runia, j=Rice, p=Perrin*/
	int	one_sv_flag;		/*true if doing a 1sv case*/
	int	op_file_flag;	/* use bits, see lookv4.c for key*/
	int	hold_flag;	/*TRUE if doing slide-hold-slide*/
	int	vs_row;
	int	first_row;
	int	last_row;
	int	mu_col;
	int	disp_col;
	int	mu_fit_col;
	int	weight_row;
	int	end_weight_row;
	int	weight_pts;
	int	weight_control;
	int	peak_row;
	int     added_pts;
	double  *disp_data;
        double  *mu_data;
        double  *model_mu;
	double	*model_tau;
	double	stiff;
	double	sig_n;
	double	v_s;
	double	v_lp;
	double	v_ref;
	double	mu_ref;
	double	tau_ref;
	double	tauo;
	double	muo;
	double	muf;
	double	amb;
	double  AMB;
	double  A;
	double	B1;
	double	B2;
	double	a;
	double	a_er;
	double	a_step;
	double	b1;
	double	b1_er;
	double	b1_step;
	double	b2;
	double	b2_er;
	double	b2_step;
	double	dc1;
	double	dc1_er;
	double	dc1_step;
	double	dc2;
	double	dc2_er;
	double	dc2_step;
	double	total_er;
	double	weight;
	double  vr_dc1;
	double  vr_dc2;
	double	mass;	/* (mass per area) --for inertial calc*/
	double  k_vf;
	double  lin_term;
	double  vel_list[512];
	double  vel_dur_list[512];
	int	n_vsteps;
	int	vs_row_list[512];
     };
struct rs_parameters rs_param ;

