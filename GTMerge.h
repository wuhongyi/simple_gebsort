#ifndef GTMERGE_H 
#define GTMERGE_H
#include <TROOT.h>
#include <vector>
#define PMODE 0644
#define LENSP 16384

#define USEBREAD 0
#define BREAD_BUFSIZE 16384
#define USEZLIB 0

#define MAXFILES 1000
#define MAXCOINEV 1000
#define DEBUG1 0
#define NCHANNELS 1000*11
#define MAXVALIDHIE 32000
#define NGE 110
#define MAXBOARDID 100000

#define TRUE 1
#define FALSE 0
#define NOTDEF -1


#define EOE          0xaaaaaaaa
#define HDRLENBYTES  28
//#define HDRLENINTS    7
#define HDRLENINTS  13  
#define HDRLENWORDS  2*HDRLENINTS
#define MAXLENINTS  519

#define LENEOVWORDS  2

#define MAXTRACELEN 8192
#define EOE          0xaaaaaaaa

#define CC_ID1 10      /* instead of 29 as default */
#define CC_ID2 11      /* instead of 39 as default */

/* instrument types */

#define NOTHING 0
#define GE 1
#define BGO 2
#define SIDE 3
#define AUX 4
#define DSSD 5
#define FP 6
#define XARRAY 7
#define CHICO2 8
#define SSD 9
#define CLOVER 10
#define SPARE 11
#define SIBOX 12

#define MAXTPE 12
#define MAXTID 1000

/*---------------*/
/* single events */
/*---------------*/

//add wuhongyi

typedef struct wuDGSEVENT
{
  unsigned short int      tpe;//实际的探测器类型 tpe[1Ge/2BGO/3SIDE] 
  unsigned short int      tid;//实际的探测器编号 tid[1-110] 

  unsigned long long int  event_timestamp;//时间戳
  unsigned long long int  last_disc_timestamp;//Timestamp of previous discriminator 上一个触发的时间戳的后30位
  unsigned long long int  peak_timestamp;//Timestamp of peak detect

  bool                    timestamp_match_flag;//TSM 有用 
  bool                    cfd_valid_flag;//CV  有用	                   
  bool                    peak_valid_flag;//PV  有用
  bool                    pileup_flag ;//PU   有用

  int                     sampled_baseline;// sampled baseline   有用
  int                     cfd_sample_0;//CFD Sample 0   有用
  int                     cfd_sample_1;//CFD Sample 1   有用
  int                     cfd_sample_2;//CFD Sample 2   有用
  
  int                     sum1;//pre-rise sum   有用
  int                     sum2;//post-rise sum   有用
  unsigned short int      m2begin;//post-rise end sample   有用
  unsigned short int      m2end;//post-rise begin sample   有用
  unsigned short int      m1begin;//pre-rise begin sample   有用
  unsigned short int      m1end;//pre-rise end sample   有用
  unsigned short int      peak_sample;//peak sample   有用
  unsigned short int      base_sample;//base sample   有用
  unsigned short int      m2last_begin_sample;
  unsigned short int      m2last_end_sample;
} wuDGSEVENT;


typedef struct wuDFMAEVENT
{
  int                     ch;    // 计算得到 (postrisesum-prerisesum)/10
  unsigned short int      tpe, tid;//通过变量 id 转换而来。实际的探测器编号
  unsigned long long int  ts;

  int wheel;//extract wheel
  unsigned long long int  prets;//上一个事件的ts
  
  int               baseline;//基线
  int		    postrisesum;//
  int		    prerisesum;//
  int               m2_last_beg;//
  int               m2_last_end;//
  int               prerisebeg;//
  int               preriseend;//
  int               postrisebeg;//
  int               postriseend;//not used
  int               peaksample;//
  int               basesample;//

  unsigned short int      traceLen;
  short int               trace[MAXTRACELEN];
}  wuDFMAEVENT;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
				
typedef struct DGSEVENT
{
public:
  float                   ehi;//修正多普勒的能量
  float                   ehiraw;//计算的原始能量道址
  float                   ehi_nodop;//未修正多普勒的能量
  short int               id;//通过 board_id 和 chan_id 运算得到
  //short int               id_pw[MAXPW];
  unsigned short int      tpe, tid;//通过变量 id 转换而来。实际的探测器编号 tpe[1Ge/2BGO/3SIDE]  tid[1-110] 
  unsigned short int      flag;//反康标记
  unsigned short int      board_id;//硬件信息
  unsigned short int      chan_id;//硬件信息
  unsigned short int      geo_addr;//硬件信息
  unsigned short int      packet_length;//数据信息

  unsigned short int      header_type;//数据信息  4
  unsigned short int      event_type;//数据信息   0
  unsigned short int      header_length;//数据信息 28

  unsigned long long int  event_timestamp;//时间戳，记得输出没有time walk 修正的
  unsigned long long int  last_disc_timestamp;//Timestamp of previous discriminator 上一个触发的时间戳的后30位
  unsigned long long int  peak_timestamp;//Timestamp of peak detect

  unsigned short int      timestamp_match_flag;//TSM 有用 
  unsigned short int      external_disc_flag;//ED    0 应该没用
  unsigned short int      cfd_valid_flag;//CV  有用
  unsigned short int      pileup_only_flag;//PO     0 应该没用
  unsigned short int      offset_flag;//OF     0 应该没用
  unsigned short int      sync_error_flag;//SE     0 应该没用
  unsigned short int      general_error_flag;//GE     0 应该没用

  unsigned short int      peak_valid_flag;//PV  有用
  unsigned short int      pileup_flag ;//PU   有用

  int                     sampled_baseline;// sampled baseline   有用
  int                     cfd_sample_0;//CFD Sample 0   有用
  int                     cfd_sample_1;//CFD Sample 1   有用
  int                     cfd_sample_2;//CFD Sample 2   有用
  int                     sum1;//pre-rise sum   有用
  int                     sum2;//post-rise sum   有用
  
  unsigned short int      m2end;//post-rise begin sample   有用
  float                   cfd_interpolate; // should be added to TS for additional accuracy - in 10ns units.
  unsigned short int      m2begin;//post-rise end sample   有用
  unsigned short int      m2last_begin_sample;
  unsigned short int      m2last_end_sample;
  unsigned short int      m1begin;//pre-rise begin sample   有用
  unsigned short int      m1end;//pre-rise end sample   有用
  unsigned short int      peak_sample;//peak sample   有用
  unsigned short int      base_sample;//base sample   有用
  
  int                     baseline;//计算得到  (m2begin+m2end)/2
  
  unsigned short int      traceLen;//波形长度 4 
  short int               trace[MAXTRACELEN];  


  //unsigned short int      A[MAXPW];
  //unsigned short int      B[MAXPW];
  //unsigned short int      C[MAXPW];
  //unsigned short int      T[MAXPW];


  unsigned long long int  LEDts;
  unsigned long long int  CFDts;
  //unsigned long long int  PEAKts;
  //char                    flag;  
  //short int               baseline;
  //unsigned short int      traceLen;
  //short int               trace[MAXTRACELEN];
} DGSEVENT ;

#define DGSEVENT_BASELEN sizeof(DGSEVENT)-MAXTRACELEN*sizeof(short int)

typedef struct DFMAEVENT
{
  // tpe 5=DSSD 6:FP  12:SIBOX
  // tid  1-160 front   161-320 back
  int               ehi;    // WAS SHORT INT    计算得到 (postrisesum-prerisesum)/10
  short int               id;//通过 board_id 和 chan_id 运算得到
  unsigned short int      tpe, tid;//通过变量 id 转换而来。实际的探测器编号
  unsigned short int      board_id;//硬件信息
  unsigned short int      chan_id;//硬件信息
  unsigned long long int  LEDts;//extract the LED time stamp
  unsigned long long int  CFDts;
  unsigned long long int  PEAKts;
  char                    flag;  
  char 			  pu;
  int			  d2t0;
  int			  d2t1;
  int                     d2t2;
  int wheel;//extract wheel
  unsigned long long int  prevTS;//上一个事件的ts
  int               baseline;//基线

  int		    postrisesum;//
  int		    prerisesum;//
  
  int               m2_last_beg;//
  int               m2_last_end;//

  int               prerisebeg;//
  int               preriseend;//
  
  int               postrisebeg;//
  int               postriseend;//not used
  
  int               peaksample;//
  int               basesample;//

  int               header_type;//数据信息

  unsigned short int      traceLen;//波形长度
  short int               trace[MAXTRACELEN];//波形
}  DFMAEVENT;

#define DFMAEVENT_BASELEN sizeof(DFMAEVENT)-MAXTRACELEN*sizeof(short int)

typedef struct GTEVENT
{
  /* raw */

  unsigned short int      len;
  short int               ehi;
  short int               id;
  short int               module;
  unsigned short int      tpe, tid;
  unsigned short int      board_id;
  unsigned short int      chan_id;
  unsigned long long int  LEDts;
  unsigned long long int  CFDts;
  unsigned long long int  PEAKts;
  char                    flag;  
  short int               baseline;
  int                     rawE;
  unsigned int            hdr[HDRLENINTS];
  unsigned short int      traceLen;
  short int               trace[MAXTRACELEN];
}  GTEVENT;




/*--------------------*/
/* coincidence events */
/*--------------------*/

typedef struct COINEV_struct
  {
  unsigned short int      len;
  unsigned char      lenclean;
  unsigned char      lendirty;
  unsigned char      lenbgo;
  unsigned char      lenaux;
  GTEVENT  GTEvent[MAXCOINEV];
  DGSEVENT DGSEvent[MAXCOINEV];
  DFMAEVENT DFMAEvent[MAXCOINEV];
  } COINEV;


typedef struct CONTROL_struct
  {
  int nOpenFiles;
  int fileActive[MAXFILES];
  int fileEventsRead[MAXFILES];
  int filesiz[MAXFILES];
  int nwritten;
  int nwritten_chunk;
  int chunkno;
  long long int chunksiz;
  int nread;
  int nzeroehi;
  int minGE;
  int minFP;
  int minDSSD;
  int minmult;
  int TSlistelen;
  int noverflowehi;
  int CurEvNo;
  long long int dts_min;
  long long int dts_max;
  int suppressBadAtRead;
  int zzipout;
  int dtsfabort;
  int dtsbabort;
  long long int startTS_lo;
  long long int startTS_hi;
  int startTS;
  int TSlist_lo;
  int TSlist_hi;
  unsigned int waitfordataseconds;
  } CONTROL;


typedef struct stat_struct
  {
  long long int inbytes;
  long long int badid;
  long long int outbytes;
  unsigned int ge_hit[NGE+1];
  unsigned int ge_cs[NGE+1];
  int nTSjumprecover_f[MAXTID];
  int nTSjumprecover_b[MAXTID];
  unsigned int in_hit[MAXTPE][MAXTID];
  unsigned int out_hit[MAXTPE][MAXTID];
  int nbigbufreads;
  int nswaps;
  unsigned int id_hit[NCHANNELS];
  unsigned int GEBIds[30];
  long long int GEBlen[30];
  } MSTAT;


/*-----------------*/
/* header for file */
/*-----------------*/

typedef struct DGSHEADER_struct
  {
  unsigned int      id;
  char              date[19];
  unsigned int      nfiles;
  char              dummy;
  } DGSHEADER;



extern std::vector<wuDGSEVENT> dgsevent_vec;//wuhongyi
extern std::vector<wuDFMAEVENT> dfmaevent_vec;//wuhongyi
extern std::vector<wuDGSEVENT> xaevent_vec;//wuhongyi

#endif
