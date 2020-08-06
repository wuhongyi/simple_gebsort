#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include "Rtypes.h"
#include "TROOT.h"
#include "TFile.h"
#include "TRandom.h"
#include "TH1.h"
#include "TH2.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TKey.h"
#include "TSystem.h"
#include "TCutG.h"
#include "TTree.h"
#include "gdecomp.h"
#include "GEBSort.h"
#include "GTMerge.h"
#define NGSGE NGE
#include "gsang.h"

#define TAPE_MOVED 10
#define EBIS_CLOCK 11
#define BETA_FIRED 12

#define ALL2DS 0
#define BASICSONLY 0

/* parameters */
extern wuDGSEVENT wuxaevent[MAXCOINEV];//wuhongyi
extern DGSEVENT XAEvent[MAXCOINEV];
extern int XAng;
extern PARS Pars;
int XAtlkup[NCHANNELS];
int XAtid[NCHANNELS];

// to be determined by the experimentalists ...
int gb_dt_lim = -10, gb_dt_lim2 = 40;   // Gamma-Beta coin gate (EhiBeta)
int gg_lim1 = -20, gg_lim2 = 20;        // Clover-Clover coin gate (hClovClov_DT)
int tg_lim1 = 620, tg_lim2 = 1220;      // Gate on the tape (ClovClov1, ...)

// X-array - histogram definitions - FGK - 062019
// TH1D *hBetaCounter;
// TH2F *hEhiBeta, *hEhiClo;
// TH2F *hGeBeta_DT, *hGeTape_DT, *hGeGe_DT;
TH2F *hClovID, *hClovBetaID, *hClovMult, *hClovBetaMult;
TH2F *hClovTape_DT, *hClovClov_DT, *hClovClov, *hClovClov1;

/* pointers to ROOT spectra */
TH1D *XAhEventCounter;
TH2F *XAhGeCounter, *XAhBGOCounter;
// TH2F *XAhEhiRaw, *XAhEhiRawRaw, *XAhEhiCln, *XAhEhiDrty, *XAhEhiCln_nodop;
TH2F *XAhGeBGO_DT;
// TH2F *XApzraw;
TH2F *XAbase1_diff;
TH2F *XAbase2_diff;
TH1D *XAhrr;
//TH2F *xa_gg;

#if(ALL2DS)
TH2F *XAe2_e1vse1[NGE + 1];
TH2F *XASZe2_e1vse1[NGE + 1];
TH2F *XAhbase1, *XAhbase2;
#endif

/* Gain Calibration */

float XAehigain[NGE + 1];
float XAehioffset[NGE + 1];
float XAehibase[NGE + 1];
float XAbl_array[5000];
float XAPZ[NGE + 1];
float XAave_base[NGE + 1];
long long int xaHeaderID[20];

float XAbase1, XAbase2, XAbase1_av = 0;
long long int ntapemoves;

/*---------------------------------------------------------------------*/
void
XAgetcal (char *file)
{
  int i, ret = 0;
  float b, c, d;
  char mystring[1000];
  FILE *fp;


  /* get pol zero */

  fp = fopen (Pars.xa_PZfn, "r");       // read mode
  if (fp == NULL)
    {
      printf ("could not open the XA PZ file: %s; use default == 1\n", Pars.xa_PZfn);
    }
  else
    {

      // read file and parse

      while (fgets (mystring, 100, fp) != NULL)
        {
          ret = sscanf (mystring, "%i %f ", &i, &b);
          XAPZ[i] = b;
	  //      printf ("ge %3i has pz of %8.4f\n", i, PZ[i]);
        }
      fclose (fp);
    };


  /* get energy cal file */

  fp = fopen (Pars.xa_ecalfn, "r");     // read mode
  if (fp == NULL)
    {
      printf ("could not open the XA cal file: %s; set 0 offset and 1.0 gain\n", Pars.xa_ecalfn);
      for (i = 1; i <= 110; i++)
        {
          XAehigain[i] = 1;
          XAehioffset[i] = 0;
        };
    }
  else
    {

      // read file and parse

      while (fgets (mystring, 100, fp) != NULL)
        {
          ret = sscanf (mystring, "%i %f %f ", &i, &c, &d);
          XAehigain[i] = d;
          XAehioffset[i] = c;
          printf ("ge %3i has offset and gain of: %8.2f %8.4f and a PZ of %8.4f ", i, XAehioffset[i], XAehigain[i],
                  XAPZ[i]);
          if (Pars.enabled[i])
            printf ("enabled\n");
          else
            printf ("DISABLED!\n");
        }
      fclose (fp);
    };

  /* done */

  return;

}


/*-----------------------------------------------------*/
int
sup_XA ()
{

  /* declarations */
  char str[256];
  FILE *fp;
  int i, i1, i2, i7, i8;
  int imod, ichan;

  // functions for making root histograms 
  TH2F *make2D (const char *, int, int, int, int, int, int);
  TH1D *make1D (const char *, int, int, int);

  //#if(BASICSONLY==0)
  // 2-D's for PZ and Baseline
  //  XApzraw = make2D ("XApzraw", NGE + 1, 1, NGE + 1, 2000, 0, 2.0);

  // 2-D's for Energy

  //  XAhEhiRaw = make2D ("XAEhiRaw", LENSP, 0, LENSP + 1, NGE + 1, 1, NGE + 1);
  //  XAhEhiRawRaw = make2D ("XAEhiRawRaw", LENSP, 0, LENSP + 1, NGE + 1, 1, NGE + 1);
  //  XAhEhiCln = make2D ("XAEhiCln", LENSP, 0, LENSP + 1, NGE + 1, 1, NGE + 1);
  //  XAhEhiCln_nodop = make2D ("XAEhiCln_nodop", LENSP, 0, LENSP + 1, NGE + 1, 1, NGE + 1);
  //  XAhEhiDrty = make2D ("XAEhiDrty", LENSP, 0, LENSP + 1, NGE + 1, 1, NGE + 1);

  /* simple gg matix */

  //  xa_gg = make2D ("xa_gg", Pars.GGMAX, 1, Pars.GGMAX, Pars.GGMAX, 1, Pars.GGMAX);
  //  xa_gg->SetXTitle ("g1");
  //  xa_gg->SetYTitle ("g2");

  //  hBetaCounter = make1D ("BetaCounter", 14400, 0, 14400);

  // X-array - 2-D's  
  //  hEhiBeta = make2D ("EhiBeta", LENSP, 0, LENSP + 1, NGE + 1, 1, NGE + 1);
  //  hEhiClo = make2D ("EhiClo", LENSP, 0, LENSP + 1, NGE + 1, 1, NGE + 1);
  //  hGeGe_DT = make2D ("GeGe_DT", 2048, -1024, 1024, NGE + 1, 1, NGE + 1);
  //  hGeBeta_DT = make2D ("GeBeta_DT", 400, -200, 200, NGE + 1, 1, NGE + 1);
  //  hGeTape_DT = make2D ("GeTape_DT", 4000, 0, 4001, 4000, 0, 4001);

  //#endif
  
  hClovID = make2D ("ClovID", LENSP, 0, LENSP + 1, 10, 0, 10);
  hClovBetaID = make2D ("ClovBetaID", LENSP, 0, LENSP + 1, 10, 0, 10);
  hClovMult = make2D ("ClovMult", LENSP, 0, LENSP + 1, 10, 0, 10);
  hClovBetaMult = make2D ("ClovBetaMult", LENSP, 0, LENSP + 1, 10, 0, 10);
  hClovTape_DT = make2D ("ClovTape_DT", 4000, 0, 4001, 4000, 0, 4001);
  hClovClov_DT = make2D ("ClovClov_DT", 2048, -1024, 1024, NGE + 1, 1, NGE + 1);
  hClovClov = make2D ("ClovClov", 4096, 0, 4096, 4096, 0, 4096);
  hClovClov1 = make2D ("ClovClov1", 4096, 0, 4096, 4096, 0, 4096);
  ntapemoves = 0;


#if(ALL2DS)//not used
  XAbase1_diff = make2D ("XAbase1_diff", NGE + 1, 1, NGE + 1, 2048, -1024, 1024);
  XAbase2_diff = make2D ("XAbase2_diff", NGE + 1, 1, NGE + 1, 2048, -1024, 1024);
  XAhbase1 = make2D ("XAhbase1", NGE + 1, 1, NGE + 1, 4096, 0, 16384);
  XAhbase2 = make2D ("XAhbase2", NGE + 1, 1, NGE + 1, 4096, 0, 16384);

  for (i = 1; i <= NGE; i++)
    {
      sprintf (str, "XAe2_e1vse1_%3.3i", i);
      XAe2_e1vse1[i] = make2D (str, 2048, 1, 10000, 1024, 1, 10000);
      sprintf (str, "XASZe2_e1vse1_%3.3i", i);
      XASZe2_e1vse1[i] = make2D (str, 2048, 1, 10000, 1024, 1, 10000);
    };
#endif


  // 2-D's for Rate
  XAhEventCounter = make1D ("XAEvntCounter", 14400, 0, 14400);  // Good for 4 hours if Counts/sec
  XAhrr = make1D ("XArr", 1012, 0, 2);
  XAhGeCounter = make2D ("XAGeCounter", 14400, 0, 14400, NGE + 1, 1, NGE + 1);
  XAhBGOCounter = make2D ("XABGOCounter", 14400, 0, 14400, NGE + 1, 1, NGE + 1);

  // 2-D's for Tacs
  XAhGeBGO_DT = make2D ("XAGeBGO_DT", 400, -200, 200, NGE + 1, 1, NGE + 1);



  /* -------------------- */
  /* read in the map file */
  /* -------------------- */
  for (i = 0; i < NCHANNELS; i++)
    {
      XAtlkup[i] = NOTHING;
      XAtid[i] = NOTHING;
    };

  fp = fopen ("XA_map.dat", "r");
  if (fp == NULL)
    {
      printf ("need an \"XA_map.dat\" file to run\n");
      exit (1);
    };

  printf ("\nXAmapping - started\n");

  i2 = fscanf (fp, "\n%i %i %i %s", &i1, &i7, &i8, str);
  while (i2 == 4)
    {
      XAtlkup[i1] = i7;
      XAtid[i1] = i8;
      i2 = fscanf (fp, "\n%i %i %i %s", &i1, &i7, &i8, str);
      ichan = i1 % 10;
      imod = (i1 - ichan) / 10;
      printf ("map: %4i ( %3i %1i ) %3i %3i %7s -- (official) ", i1, imod, ichan, i7, i8, str);
      switch (i7)
        {
        case GE:
          printf ("GE\n");
          break;
        case BGO:
          printf ("BGO\n");
          break;
        case SIDE:
          printf ("SIDE\n");
          break;
        case AUX:
          printf ("AUX\n");
          break;
        case DSSD:
          printf ("DSSD\n");
          break;
        case FP:
          printf ("FP\n");
          break;
        case XARRAY:
          printf ("XARRAY\n");
          break;
        case CHICO2:
          printf ("CHICO2\n");
          break;
        case SSD:
          printf ("SSD\n");
          break;
        case CLOVER:
          printf ("CLOVER\n");
          break;
        case SPARE:
          printf ("SPARE\n");
          break;
        case SIBOX:
          printf ("SIBOX\n");
          break;
        default:
          printf ("dont know what this is\n");
          break;

        };
    };
  fclose (fp);

  printf ("\nmapping - complete!!\n");

  /* Set Default Calibration and PZ parameters */

  for (i = 0; i <= NGE + 1; i++)
    {
      XAehigain[i] = 1.0;
      XAehioffset[i] = 0.0;
      XAPZ[i] = 1.0;
      XAehibase[i] = 0.0;
    };

  for (i = 0; i < 5000; i++)
    {
      XAbl_array[i] = 0;
    }

  for (i = 0; i < NGE + 1; i++)
    {
      XAave_base[i] = 0;
    }

  /* get the DGS calibration file */

  XAgetcal (Pars.xa_ecalfn);

  /* dgs header ids */

  for (i = 0; i < 20; i++)
    xaHeaderID[i] = 0;

  /* done */
  return (0);
};


/*-----------------------------------------------------*/
int
exit_XA ()
{

  /* declarations */
  printf ("\nbegin exit_XA\n");
  printf ("tape moved %lli times in this dataset\n", ntapemoves);

  /* done */
  printf ("done exit_XA\n");
  return (0);
};


/* ----------------------------------------------------------------- */
int
bin_XA (GEB_EVENT * GEB_event)
{

  /* declarations */
  int i, j, gsid, e, tdiff, tdiff1, tdiff2, st;
  unsigned long long int EvTimeStam0;
  int is_beta = 0;
  static long long int tape_ts, ebis_ts, beta_ts;
  char str[128];
  int RelEvT = 0;
  float Energy, Energy_nodop, top, bot, r1, rr;
  double d1, d2, erawraw[NGE + 1];

  /* prototypes */

  int GebTypeStr (int type, char str[]);
  int DGSEvDecompose_v3 (unsigned int *ev, int len, DGSEVENT * XAEvent, int tlkup[], int tid[]);
  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("\nentered bin_XA:\n");

#if(1)
  /* temp fix of header 24-->XA for testing */

  for (i = 0; i < GEB_event->mult; i++)
    if (GEB_event->ptgd[i]->type == 24)
      GEB_event->ptgd[i]->type = GEB_TYPE_XA;
  //  for (i = 0; i < GEB_event->mult; i++) 
  //    printf("%i is %i\n", i,  GEB_event->ptgd[i]->type);
#endif

  //  for (i=1000;i<1030;i++)
  //    printf("debug1: %i: tlkup %i tid %i\n", i, XAtlkup[i], XAtid[i]);

  /* loop through the coincidence event and fish out XA data */
  /* (gamma rays) count in XAng */

  XAng = 0;
  xaevent_vec.clear();
  for (i = 0; i < GEB_event->mult; i++)
    {
      if (GEB_event->ptgd[i]->type == GEB_TYPE_XA)
        {
          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              GebTypeStr (GEB_event->ptgd[i]->type, str);
              printf ("bin_XA (header): %2i> %2i, %s, TS=%lli\n", i, GEB_event->ptgd[i]->type, str,
                      GEB_event->ptgd[i]->timestamp);
            }

          st =
            DGSEvDecompose_v3 ((unsigned int *) GEB_event->ptinp[i], GEB_event->ptgd[i]->length / sizeof (unsigned int),
                               &XAEvent[XAng], XAtlkup, XAtid);//wuhongyi jta.c
          if (st != 0) return (0);
          /* xaevent_vec.push_back(XAEvent[XAng]);//wuhongyi 这里填充的只是最原始的数据，而且数据存在问题 */

	  wuxaevent[XAng].tpe = XAEvent[XAng].tpe;
	  wuxaevent[XAng].tid = XAEvent[XAng].tid;
	  wuxaevent[XAng].event_timestamp = XAEvent[XAng].event_timestamp;
	  wuxaevent[XAng].last_disc_timestamp = XAEvent[XAng].last_disc_timestamp;
	  wuxaevent[XAng].peak_timestamp = XAEvent[XAng].peak_timestamp;
	  wuxaevent[XAng].timestamp_match_flag = XAEvent[XAng].timestamp_match_flag;
	  wuxaevent[XAng].cfd_valid_flag = XAEvent[XAng].cfd_valid_flag;
	  wuxaevent[XAng].peak_valid_flag = XAEvent[XAng].peak_valid_flag;
	  wuxaevent[XAng].pileup_flag = XAEvent[XAng].pileup_flag;
	  wuxaevent[XAng].sampled_baseline = XAEvent[XAng].sampled_baseline;
	  wuxaevent[XAng].cfd_sample_0 = XAEvent[XAng].cfd_sample_0;
	  wuxaevent[XAng].cfd_sample_1 = XAEvent[XAng].cfd_sample_1;
	  wuxaevent[XAng].cfd_sample_2 = XAEvent[XAng].cfd_sample_2;
	  wuxaevent[XAng].sum1 = XAEvent[XAng].sum1;
	  wuxaevent[XAng].sum2 = XAEvent[XAng].sum2;
	  wuxaevent[XAng].m2begin = XAEvent[XAng].m2begin;
	  wuxaevent[XAng].m2end = XAEvent[XAng].m2end;
	  wuxaevent[XAng].m1begin = XAEvent[XAng].m1begin;
	  wuxaevent[XAng].m1end = XAEvent[XAng].m1end;
	  wuxaevent[XAng].peak_sample = XAEvent[XAng].peak_sample;
	  wuxaevent[XAng].base_sample = XAEvent[XAng].base_sample;
	  xaevent_vec.push_back(wuxaevent[XAng]);
	  
          XAng++;
        }
    }

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("XAng=%i\n", XAng);
      for (i = 0; i < XAng; i++)
        {
          printf ("gsid=XAEvent[i].tid=%i ", XAEvent[i].tid);
          printf ("tpe=%i ", XAEvent[i].tpe);
          printf ("ts %lli\n", XAEvent[i].event_timestamp);
          if (XAEvent[i].tpe == TAPE_MOVED)
            printf ("10 tape moved\n");
          if (XAEvent[i].tpe == EBIS_CLOCK)
            printf ("11 EBIS clock\n");
          if (XAEvent[i].tpe == BETA_FIRED)
            printf ("12 beta fired\n");
        };
    };
  EvTimeStam0 = XAEvent[0].event_timestamp;
  if (EvTimeStam0 == 0)
    EvTimeStam0 = XAEvent[0].event_timestamp;
  RelEvT = (int) ((XAEvent[0].event_timestamp - EvTimeStam0) / 100000000);      // overflow?
  XAhEventCounter->Fill (RelEvT);

  /* Loop */

  for (i = 0; i < XAng; i++)
    {
      gsid = XAEvent[i].tid;
      if (XAEvent[i].tpe == GE)
        {

          /* this should be unnecessary to do; */
          /* but we sometimes crash if we don't */
          /* needs to be looked at */

          if (gsid < 1 || gsid > NGE)
            {
              printf ("bad gsid= %i\n", gsid);
              fflush (stdout);
              gsid = 0;
            };

          XAhGeCounter->Fill ((int) ((XAEvent[0].event_timestamp - EvTimeStam0) / 100000000), gsid);

	  /*           current baseline; SZ's Method1 implementation */

          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              printf ("\n");
              printf ("m1begin= XAEvent[%i].m1begin  = %8i w11,lobit\n", i, XAEvent[i].m1begin);
              printf ("m1end  = XAEvent[%i].m1end    = %8i w11,upbit ", i, XAEvent[i].m1end);
              printf ("diff %i\n", XAEvent[i].m1begin - XAEvent[i].m1end);
              printf ("m2begin= XAEvent[%i].m2begin  = %8i w10,upbit\n", i, XAEvent[i].m2begin);
              printf ("m2end  = XAEvent[%i].m2end    = %8i w10,lobit ", i, XAEvent[i].m2end);
              printf ("diff %i\n", XAEvent[i].m2begin - XAEvent[i].m2end);
            };


          top = (float) (XAEvent[i].sum2 * XAEvent[i].m1begin - XAEvent[i].sum1 * XAEvent[i].m2begin);
          bot = XAEvent[i].sum2 - XAEvent[i].sum1;
          bot -= Pars.xa_MM * (XAEvent[i].m2begin - XAEvent[i].m1begin);
          XAbase1 = top / bot;

          top = (float) (XAEvent[i].m1end * XAEvent[i].m2begin - XAEvent[i].m1begin * XAEvent[i].m2end);
          bot = (XAEvent[i].m2begin - XAEvent[i].m2end) - (XAEvent[i].m1begin - XAEvent[i].m1end);
          XAbase2 = top / bot;

#if(ALL2DS)
          if (Pars.enabled[gsid])
            {
              XAbase1_diff->Fill ((double) gsid, XAEvent[i].m1begin - XAEvent[i].m1end);
              XAbase2_diff->Fill ((double) gsid, XAEvent[i].m2begin - XAEvent[i].m2end);
              XAhbase1->Fill ((double) gsid, (double) XAbase1);
              XAhbase2->Fill ((double) gsid, (double) XAbase2);
            };
#endif

          /* keep a running average of the baseline */
          /* avoid extreme values (per Darek) */
          /* add constraints from SZ as well */

          rr = XAbase1 / XAbase2;
	  	      

          if (rr > 0.0 && rr < 2.0)
            XAhrr->Fill (double (rr));

          if ((XAbase1 < 40000) && (XAbase1 > 0))
            if ((XAEvent[i].sum2 - XAEvent[i].sum1) > 200)
              if (rr > 0.95 && rr < 1.05 && i == 0){	
		if (XAave_base[gsid] > 0){
		  XAave_base[gsid] = XAave_base[gsid] * 0.99 + XAbase1 * 0.01;
		  //                    XAave_base[gsid] = 958.12245237;
		}
                else
                  XAave_base[gsid] = XAbase1;}
	  //	      else
	  //		      XAave_base[gsid] = XAbase1;
	 


          /* this is SZ's baseline restore correction, see */
          /* soffice /home/tl/d6/keep/2018_zhu_baseline.ppt */
          /* ave_base: the running base. note: Pars.xa_MM is a float */

	  //#if(BASICSONLY ==0)

          /* use average base */

          Energy = XAEvent[i].sum2 / Pars.xa_MM - XAEvent[i].sum1 * XAPZ[gsid] / Pars.xa_MM;
          Energy -= (1. - XAPZ[gsid]) * XAave_base[gsid];
	  //	  XAave_base[gsid] = 0;
	  //#endif
	  //#if(0)
          /* use event base (you get noise) */

	  //          Energy = XAEvent[i].sum2 / Pars.xa_MM - XAEvent[i].sum1 * XAPZ[gsid] / Pars.xa_MM;
	  //          Energy -= (1. - XAPZ[gsid]) * base2;
	  //#endif

          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              printf ("\n");
              printf ("XAEvent[i].sum2= %i\n", XAEvent[i].sum2);
              printf ("XAEvent[i].sum1= %i\n", XAEvent[i].sum1);
              printf ("Pars.xa_MM= %f\n", Pars.xa_MM);
              printf ("PZ[%i]= %f\n", gsid, XAPZ[gsid]);
              printf ("base1= %8.1f base2= %8.1f, diff= %8.1f rat= %5.2f\n", XAbase1, XAbase2, XAbase1 - XAbase2,
                      XAbase1 / XAbase2);
              printf ("ave_base[%i]= %f\n", gsid, XAave_base[gsid]);
              printf ("Energy = %f (uncalibrated)\n", Energy);
            };

          /* gain match */

	  //          erawraw[i] = Energy;
	  XAEvent[i].ehiraw = Energy;
          Energy = Energy * XAehigain[gsid] + XAehioffset[gsid];
          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              printf ("Energy = %f (calibrated)\n", Energy);
            };


          /* dopler correct the energy */

          Energy_nodop = Energy;
          if (Pars.vc_xa != 0)
            {
              d1 = angtheta[gsid - 1] / 57.29577951;
              Energy = Energy * (1 - Pars.vc_xa * cos (d1)) / sqrt (1 - Pars.vc_xa * Pars.vc_xa);
            }
          if (Pars.enabled[gsid])
            {
              XAEvent[i].ehi = Energy;
              XAEvent[i].ehi_nodop = Energy_nodop;
	      if(XAEvent[i].ehi>100 && XAEvent[i].ehi<=180){
		XAEvent[i].event_timestamp = XAEvent[i].event_timestamp -(int)((180-XAEvent[i].ehi)/20);
	      }
	  
	      if(XAEvent[i].ehi>0 && XAEvent[i].ehi<= 100){
		XAEvent[i].event_timestamp = XAEvent[i].event_timestamp -(int)((100-XAEvent[i].ehi)/7 + 6);
	      }

            }
          else
            {
              /* mark bad */
              XAEvent[i].ehi = -1;
              XAEvent[i].ehi_nodop = -1;
            }
          XAEvent[i].id = gsid;
	  if(XAEvent[i].event_timestamp > 5334995920870 && XAEvent[i].event_timestamp < 5334995920890){
	    printf("zwq:XAEvent[i].ehi=%f\n",XAEvent[i].ehi);
	  } 

	  //#if(BASICSONLY)

	  //          /* abandon the rest */

	  //          return (0);

	  //#endif

#if(ALL2DS)
          if (Pars.enabled[gsid])
            {
              d1 = (XAEvent[i].sum2 - XAEvent[i].sum1) / Pars.xa_MM;    /* uncorrected energy */
              d2 = XAEvent[i].sum1 / Pars.xa_MM;        /* ~ baseline */
              if (gsid <= NGE)
                if (d1 < (double) 8192)
                  if (d2 < (double) 8192)
                    {
                      XAe2_e1vse1[gsid]->Fill (d1, d2, 1);
                      d1 = XAEvent[i].ehi;      /* SZ corrected */
                      XASZe2_e1vse1[gsid]->Fill (d1, d2, 1);
                    };
            };
#endif

          /* fill SZ's pol zero (raw) spectrum */

          top = XAEvent[i].m2end - XAEvent[i].m1end;
          bot = XAEvent[i].m2begin - XAEvent[i].m1begin;
          r1 = top / bot;
	  //          if (r1 > 0 && r1 < 2.0)
	  //            if (Pars.enabled[gsid])
	  //              XApzraw->Fill (gsid, r1);

        }


    }                           /* for (i = 0; i < ng; i++) */

  /* Energy Histogram loop */

  for (i = 0; i < XAng; i++)
    {
      if (XAEvent[i].tpe == GE)
        {
          if (erawraw[i] > 0 && erawraw[i] < LENSP)
	    //            XAhEhiRawRaw->Fill (erawraw[i], gsid);
	    e = (int) XAEvent[i].ehi;
          if (e > 0 && e < LENSP)
            {
              gsid = XAEvent[i].tid;
	      //              XAhEhiRaw->Fill (e, gsid);
              if (XAEvent[i].flag == 0)
                {
		  //                  XAhEhiCln->Fill (e, gsid);
		  //                  XAhEhiCln_nodop->Fill (XAEvent[i].ehi_nodop, gsid);
                };
              if (XAEvent[i].flag == 1)
		//                XAhEhiDrty->Fill (e, gsid);
		;};
        };
    }

  /* gg matrix */

  for (i = 0; i < XAng; i++)
    if (Pars.enabled[XAEvent[i].tid])
      if (XAEvent[i].tpe == GE)
        if (XAEvent[i].flag == 0)
          if (XAEvent[i].ehi > 0 && XAEvent[i].ehi < Pars.GGMAX)
            for (j = i + 1; j < XAng; j++)
              if (Pars.enabled[XAEvent[j].tid])
                if (XAEvent[j].tpe == GE)
                  if (XAEvent[j].flag == 0)
                    if (XAEvent[j].ehi > 0 && XAEvent[j].ehi < Pars.GGMAX)
                      {
			//                        xa_gg->Fill ((double) XAEvent[i].ehi, (double) XAEvent[j].ehi, 1);
                        if (Pars.CurEvNo <= Pars.NumToPrint)
                          printf ("filled dgs_gg with %f %f\n", XAEvent[i].ehi, XAEvent[j].ehi);
                      }

  /* debug list the gamma rays we found */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("final list: we have %i gamma rays\n", XAng);
      for (i = 0; i < XAng; i++)
        {
          printf ("%2i> ", i);
          printf ("chan_id=%i; ", XAEvent[i].chan_id);
          printf ("board_id=%i; ", XAEvent[i].board_id);
          printf ("id =%i; ", XAEvent[i].id);
          printf ("tpe=%i; ", XAEvent[i].tpe);
          printf ("tid=%i; ", XAEvent[i].tid);
          printf ("EventTS=%llu; ", XAEvent[i].event_timestamp);
          printf ("ehi=%8f ", XAEvent[i].ehi);
          if (XAEvent[i].flag == 1)
            printf ("(dirty)");
          else
            printf ("(clean)");
          if (Pars.enabled[gsid])
            printf ("(enabled)");
          else
            printf ("(disabled)");
          printf ("\n");
          fflush (stdout);
        };
    };

  // initial loop through the channels
  for (i = 0; i < XAng; i++)

    {
      if (XAEvent[i].tpe == TAPE_MOVED)

        {

          /* we can usually afford to print this out every time */
          if (Pars.CurEvNo <= Pars.NumToPrint || 1)

            {
              d1 = XAEvent[i].event_timestamp - tape_ts;
              d1 /= 100000000;
              printf ("tape moved @TS=%lli; time since last tape move: %5.2f sec\n", XAEvent[i].event_timestamp, d1);
            };
          ntapemoves++;
          tape_ts = XAEvent[i].event_timestamp;
        }
      if (XAEvent[i].tpe == EBIS_CLOCK)

        {
          ebis_ts = XAEvent[i].event_timestamp;

          //printf("EBIS Clock %i \n",XAEvent[i].tpe);
        }
      if (XAEvent[i].tpe == BETA_FIRED)

        {
          is_beta = 1;
          beta_ts = XAEvent[i].event_timestamp;
	  //          hBetaCounter->Fill ((int) ((beta_ts - EvTimeStam0) / 100000000));

	  //      printf("Beta Detected %llu %llu \n",XAEvent[i].event_timestamp,beta_ts); 
        }
    };

  /* Energy Histograms loop */
  int cbc[21];
  for (i = 0; i < XAng; i++)

    {

      // time difference between clover crystal 1 and any other crystals 
      gsid = XAEvent[i].tid;
      if (gsid != 1)

        {
          tdiff = (int) (XAEvent[1].event_timestamp - XAEvent[i].event_timestamp);
	  //          hGeGe_DT->Fill (tdiff, gsid);
          tdiff = (int) (XAEvent[i].event_timestamp - XAEvent[1].event_timestamp);
	  //          hGeGe_DT->Fill (tdiff, gsid);
        }
      cbc[i] = 0;
      if (XAEvent[i].tpe == 1)  // select Ge
        {
          e = (int) XAEvent[i].ehi;
          if (e > 10 && e < LENSP)      // select E_gamma > 10 channel
            {
              gsid = XAEvent[i].tid;
	      //              hEhiClo->Fill (e, gsid);
              if (is_beta == 1) //  is_beta = 1 - beta is in the event
                {
                  tdiff = (int) (XAEvent[i].event_timestamp - beta_ts);
		  //                  hGeBeta_DT->Fill (tdiff, gsid);
                  if (tdiff > gb_dt_lim && tdiff < gb_dt_lim2)  // beta-gama coin window  
                    {
		      //                      hEhiBeta->Fill (e, gsid); // Generate beta-gated e_gamma vs id
                      cbc[i] = 1;
                    }
                  if (tape_ts != 0)     // Generate tape histogram for beta-gated e_gamma 
                    {
                      tdiff = (int) ((XAEvent[i].event_timestamp - tape_ts) / 10000000);
		      //                      hGeTape_DT->Fill (e, tdiff);
                    }
                }
            }
        }
    }
  /* Make summed clover. */
  int cloverE[6];
  int ii, clov_pu[6], clov_beta_coin[6], clov_mult[6];
  unsigned long long int clov_ts[6];
  unsigned long long int clov_ts0 = 0;
  for (i = 1; i < 6; i++)

    {
      cloverE[i] = 0;
      clov_pu[i] = 0;
      clov_ts[i] = 0;
      clov_beta_coin[i] = 0;
      clov_mult[i] = 0;
    };
  for (i = 0; i < XAng; i++)

    {
      e = (int) XAEvent[i].ehi;
      if (XAEvent[i].tpe == 1 && e > 10)        // select Ge and E_gamma > 10 ch
        {
          ii = 1;
          for (j = 1; j < 20; j += 4)

            {
              if (XAEvent[i].tid > j - 1 && XAEvent[i].tid < j + 4)

                {
                  clov_mult[ii] = clov_mult[ii] + 1;
                  cloverE[ii] = cloverE[ii] + e;
                  if (clov_ts[ii] == clov_ts0)
                    clov_ts[ii] = XAEvent[i].event_timestamp;
                  if (clov_beta_coin[ii] == 0)
                    clov_beta_coin[ii] = cbc[i];
                  if (XAEvent[i].pileup_flag == 1)
                    clov_pu[ii] = XAEvent[i].pileup_flag;       // pileup
                }
              ii++;
            }
        }
    };
  for (i = 1; i < 6; i++)

    {
      tdiff = (int) ((clov_ts[i] - tape_ts) / 10000000);
      if (clov_beta_coin[i] == 1)
        hClovTape_DT->Fill (cloverE[i], tdiff);
      hClovID->Fill (cloverE[i], i);
      if (clov_beta_coin[i] == 1)
        hClovBetaID->Fill (cloverE[i], i);
      hClovMult->Fill (cloverE[i], clov_mult[i]);
      if (clov_beta_coin[i] == 1)
        hClovBetaMult->Fill (cloverE[i], clov_mult[i]);
    }

  // g-g coincidence histograms
  for (i = 1; i < 6; i++)

    {
      if (cloverE[i] > 10 && clov_beta_coin[i] == 1)

        {
          tdiff1 = (int) ((clov_ts[i] - tape_ts) / 10000000);   // 100 msec clover time vs tape_0
          for (j = i + 1; j < 6; j++)

            {
              if (cloverE[j] > 10 && clov_beta_coin[j] == 1)

                {
                  tdiff2 = (int) ((clov_ts[j] - tape_ts) / 10000000);   // 100 msec clover time vs tape_0
                  tdiff = (int) (clov_ts[j] - clov_ts[i]);
                  hClovClov_DT->Fill (tdiff, i);
                  if (tdiff > gg_lim1 && tdiff < gg_lim2)

                    {
                      hClovClov->Fill (cloverE[i], cloverE[j]);
                      hClovClov->Fill (cloverE[j], cloverE[i]);
                      if (tdiff1 > tg_lim1 && tdiff1 > tg_lim1 && tdiff2 > tg_lim1 && tdiff2 > tg_lim1)

                        {
                          hClovClov1->Fill (cloverE[i], cloverE[j]);
                          hClovClov1->Fill (cloverE[j], cloverE[i]);
                        }
                    }
                }
            }
        }
    }

  /* done */
  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("exit bin_XA\n");
  return (0);
}
