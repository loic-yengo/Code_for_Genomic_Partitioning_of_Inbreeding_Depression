#include <string.h>
#include <math.h>
#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <random>
#include <ctime>
#include <unordered_map>
#include "Eigen/Core"
#include "Eigen/Dense"

#define PACK_DENSITY 4
#define MASK0   3 // 3 << 2 * 0
#define MASK1  12 // 3 << 2 * 1
#define MASK2  48 // 3 << 2 * 2
#define MASK3 192 // 3 << 2 * 3

using namespace std;
using namespace Eigen;
unordered_map<string,int> index_of_snps;

void decode_plink(char *output, const char *input, const int lengthInput){
  int i, k;
  char tmp, geno;
  int a1, a2;
  
  for(i=0;i<lengthInput;++i){
    tmp = input[i];
    k   = PACK_DENSITY * i;
    geno      = (tmp & MASK0);
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;
    
    geno      = (tmp & MASK1) >> 2; 
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;
    
    geno      = (tmp & MASK2) >> 4; 
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;
    
    geno      = (tmp & MASK3) >> 6; 
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
  }
}


// Main funtion used on cov_files.txt with a sequential access
int main(int argc, char *argv[]){
  
  if(argc==1){
    cerr<<"\tMust specify which chromosome to analyse (e.g. ./fannot 22)."<<endl;
    exit(1);
  }
  
  int chrom_num = atoi(argv[1]);
  int chrom     = chrom_num - 1;
  
  string i_start_chr = argv[2];
  string i_end_chr   = argv[3];
  
  int i_start   = atoi(i_start_chr.c_str());
  int i_end     = atoi(i_end_chr.c_str());
  
  cout<<"# Processing chromosome "<<chrom_num<<"...\n";
  cout<<"# Individual #"<<i_start_chr<<" to #"<<i_end_chr<<".\n";  
  
  string output_g[22] = {
    "annots/ibc_fast/g.chrom1",
    "annots/ibc_fast/g.chrom2",
    "annots/ibc_fast/g.chrom3",
    "annots/ibc_fast/g.chrom4",
    "annots/ibc_fast/g.chrom5",
    "annots/ibc_fast/g.chrom6",
    "annots/ibc_fast/g.chrom7",
    "annots/ibc_fast/g.chrom8",
    "annots/ibc_fast/g.chrom9",
    "annots/ibc_fast/g.chrom10",
    "annots/ibc_fast/g.chrom11",
    "annots/ibc_fast/g.chrom12",
    "annots/ibc_fast/g.chrom13",
    "annots/ibc_fast/g.chrom14",
    "annots/ibc_fast/g.chrom15",
    "annots/ibc_fast/g.chrom16",
    "annots/ibc_fast/g.chrom17",
    "annots/ibc_fast/g.chrom18",
    "annots/ibc_fast/g.chrom19",
    "annots/ibc_fast/g.chrom20",
    "annots/ibc_fast/g.chrom21",
    "annots/ibc_fast/g.chrom22"
  };
  
  // Sum of F per chromsomes
  string output_mf[22] = {
    "annots/ibc_fast/m.chrom1.f",
    "annots/ibc_fast/m.chrom2.f",
    "annots/ibc_fast/m.chrom3.f",
    "annots/ibc_fast/m.chrom4.f",
    "annots/ibc_fast/m.chrom5.f",
    "annots/ibc_fast/m.chrom6.f",
    "annots/ibc_fast/m.chrom7.f",
    "annots/ibc_fast/m.chrom8.f",
    "annots/ibc_fast/m.chrom9.f",
    "annots/ibc_fast/m.chrom10.f",
    "annots/ibc_fast/m.chrom11.f",
    "annots/ibc_fast/m.chrom12.f",
    "annots/ibc_fast/m.chrom13.f",
    "annots/ibc_fast/m.chrom14.f",
    "annots/ibc_fast/m.chrom15.f",
    "annots/ibc_fast/m.chrom16.f",
    "annots/ibc_fast/m.chrom17.f",
    "annots/ibc_fast/m.chrom18.f",
    "annots/ibc_fast/m.chrom19.f",
    "annots/ibc_fast/m.chrom20.f",
    "annots/ibc_fast/m.chrom21.f",
    "annots/ibc_fast/m.chrom22.f"
  };
  
  // Sum of weights per chromsomes
  string output_mn[22] = {
    "annots/ibc_fast/m.chrom1.n",
    "annots/ibc_fast/m.chrom2.n",
    "annots/ibc_fast/m.chrom3.n",
    "annots/ibc_fast/m.chrom4.n",
    "annots/ibc_fast/m.chrom5.n",
    "annots/ibc_fast/m.chrom6.n",
    "annots/ibc_fast/m.chrom7.n",
    "annots/ibc_fast/m.chrom8.n",
    "annots/ibc_fast/m.chrom9.n",
    "annots/ibc_fast/m.chrom10.n",
    "annots/ibc_fast/m.chrom11.n",
    "annots/ibc_fast/m.chrom12.n",
    "annots/ibc_fast/m.chrom13.n",
    "annots/ibc_fast/m.chrom14.n",
    "annots/ibc_fast/m.chrom15.n",
    "annots/ibc_fast/m.chrom16.n",
    "annots/ibc_fast/m.chrom17.n",
    "annots/ibc_fast/m.chrom18.n",
    "annots/ibc_fast/m.chrom19.n",
    "annots/ibc_fast/m.chrom20.n",
    "annots/ibc_fast/m.chrom21.n",
    "annots/ibc_fast/m.chrom22.n"
  };
  
  
  string mbfiles[22] = {
    "annots/ukb-geno/chrom1",
    "annots/ukb-geno/chrom2",
    "annots/ukb-geno/chrom3",
    "annots/ukb-geno/chrom4",
    "annots/ukb-geno/chrom5",
    "annots/ukb-geno/chrom6",
    "annots/ukb-geno/chrom7",
    "annots/ukb-geno/chrom8",
    "annots/ukb-geno/chrom9",
    "annots/ukb-geno/chrom10",
    "annots/ukb-geno/chrom11",
    "annots/ukb-geno/chrom12",
    "annots/ukb-geno/chrom13",
    "annots/ukb-geno/chrom14",
    "annots/ukb-geno/chrom15",
    "annots/ukb-geno/chrom16",
    "annots/ukb-geno/chrom17",
    "annots/ukb-geno/chrom18",
    "annots/ukb-geno/chrom19",
    "annots/ukb-geno/chrom20",
    "annots/ukb-geno/chrom21",
    "annots/ukb-geno/chrom22"
  };
  
  string qcedfreqs[22] = {
    "annots/freq/chrom1.afreq",
    "annots/freq/chrom2.afreq",
    "annots/freq/chrom3.afreq",
    "annots/freq/chrom4.afreq",
    "annots/freq/chrom5.afreq",
    "annots/freq/chrom6.afreq",
    "annots/freq/chrom7.afreq",
    "annots/freq/chrom8.afreq",
    "annots/freq/chrom9.afreq",
    "annots/freq/chrom10.afreq",
    "annots/freq/chrom11.afreq",
    "annots/freq/chrom12.afreq",
    "annots/freq/chrom13.afreq",
    "annots/freq/chrom14.afreq",
    "annots/freq/chrom15.afreq",
    "annots/freq/chrom16.afreq",
    "annots/freq/chrom17.afreq",
    "annots/freq/chrom18.afreq",
    "annots/freq/chrom19.afreq",
    "annots/freq/chrom20.afreq",
    "annots/freq/chrom21.afreq",
    "annots/freq/chrom22.afreq"
  };
  
  string qcedannotfiles[22] = {
    "annots/annot-files/baselineLF2.2.UKB.1.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.2.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.3.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.4.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.5.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.6.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.7.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.8.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.9.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.10.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.11.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.12.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.13.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.14.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.15.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.16.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.17.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.18.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.19.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.20.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.21.qced.annot",
    "annots/annot-files/baselineLF2.2.UKB.22.qced.annot"
  };
  
  string bfile, bedfile, bimfile, frqfile, token;
  bfile      = mbfiles[chrom];
  bimfile    = bfile+".bim";
  bedfile    = bfile+".bed";
  frqfile    = qcedfreqs[chrom];
  
  string line = "";
  ifstream bimStream;
  int p = -1;
  bimStream.open(bimfile.c_str());
  while(bimStream){
    getline(bimStream,line);
    p++;
  }
  bimStream.close();
  cout<<"# "<<p<<" SNPs detected."<<endl;
  
  string fileannot = qcedannotfiles[chrom];
  
  int i, j, k;
  int nsnps_annot   = p; //nsnpbimfiles[chrom];
  
  string annotNames = "annots/baselineLF2.2.UKB.annot-names.txt";
  ifstream annotNamesStream;
  int n_annot       = -1;
  annotNamesStream.open(annotNames.c_str());
  while(annotNamesStream){
    getline(annotNamesStream,line);
    n_annot++;
  }
  annotNamesStream.close();
  cout<<"# "<<n_annot<<" annotations detected from "<<annotNames<<endl;
  
  MatrixXf weights = MatrixXf::Zero(nsnps_annot,n_annot);
  
  string bp, wgts_chr;
  
  cout<<"# Reading annotation file ("<<fileannot<<")...";
  clock_t tic_readannot = clock();
  ifstream tmpStream;
  tmpStream.open(fileannot.c_str());
  for(j=0;j<nsnps_annot;j++){
    tmpStream >> bp;
    index_of_snps.insert({bp,j});
    for(k=0;k<n_annot;k++){
      tmpStream >> wgts_chr;
      weights(j,k) = atof(wgts_chr.c_str());
    }
  }
  tmpStream.close();
  clock_t toc_readannot = clock();
  cout<<".[done].\n";
  printf("# (time elapsed: %f seconds).\n", (float)(toc_readannot - tic_readannot) / CLOCKS_PER_SEC);    
  
  int N = 348501;
  int n = i_end - i_start + 1;  
  cout<<"# Processing genotypes of "<<n<<" individuals.\n";
  VectorXf Fg = VectorXf::Zero(n);
  VectorXf Ng = VectorXf::Zero(n);
  MatrixXf Fm = MatrixXf::Zero(n,n_annot);
  MatrixXf Nm = MatrixXf::Zero(n,n_annot);
  
  int numBytes   = (int)ceil((float)N / PACK_DENSITY);
  char* packed   = new char[numBytes];
  char* unpacked = new char[numBytes * PACK_DENSITY];
  
  // string bfile, bedfile, bimfile, frqfile, token;
  int g_i;
  float p_j, q_j;
  float F0 = 0., F2 = 0., F1 = -1.;
  
  VectorXf F = VectorXf::Zero(n);
  VectorXf S = VectorXf::Zero(n);
  
  /*int p      = nsnpbimfiles[chrom];
   bfile      = mbfiles[chrom];
   bimfile    = bfile+".bim";
   bedfile    = bfile+".bed";
   frqfile    = qcedfreqs[chrom];
   */ 
  
  ifstream influx;  
  influx.open(bedfile.c_str(), std::ios::in | std::ios::binary);
  if(!influx){
    cerr << "[readGenotypes] Error reading file "<<bedfile<<endl;
    exit(1);
  }
  influx.seekg(0, ifstream::end);
  influx.seekg(3, ifstream::beg);
  
  clock_t tic = clock();
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  cout <<"# Starting reading genotypes : ";
  cout << (now->tm_year + 1900) << '-' 
       << (now->tm_mon + 1) << '-'
       <<  now->tm_mday << " at "
       <<  now->tm_hour <<":"<<now->tm_min<<":"<<now->tm_sec
       << endl;
  
  //ifstream bimStream;
  bimStream.open(bimfile.c_str());
  ifstream freqStream;
  string freqFileheader;
  freqStream.open(frqfile.c_str());
  for(k=0;k<6;k++){
    freqStream >> token;
  }
  for(j=0;j<p;j++){
    // Read bim file to get chr and pos
    for(k=0;k<6;k++){
      bimStream >> token;
      if(k==3){
        bp = token;
      }
    }
    int j_in = index_of_snps[bp];
    
    // Read frequency file
    p_j = 0.;
    for(k=0;k<6;k++){
      freqStream >> token;
      if(k==4){
        p_j = atof(token.c_str());
      }
    }
    
    influx.read((char*)packed, sizeof(char) * numBytes);
    decode_plink(unpacked, packed, numBytes);
    
    //cout<<snp<<"\t"<<p_j<<endl;
    if(p_j>=0.01 and p_j<=0.99){
      q_j = 1.0-p_j;
      F0  = p_j/q_j;
      F2  = q_j/p_j;
      int ii = 0;
      for(i=0;i<N;i++){
        g_i  = (int) unpacked[i];
        if(i>=i_start and i<=i_end){
          if(g_i==3){
            F(ii) = S(ii) = 0.;
          }else{
            S(ii) = 1.0;
            if(g_i==0) F(ii) = F0;
            if(g_i==1) F(ii) = F1;
            if(g_i==2) F(ii) = F2;
          }
          ii++;
        }
      } // for i
      Ng += S;
      Fg += F;
      Fm += F * weights.row(j_in);
      Nm += S * weights.row(j_in);
    }
  }// for j
  freqStream.close();
  bimStream.close();
  influx.close();
  
  string tmp_g  = output_g[chrom]+"_"+i_start_chr+"_"+i_end_chr;
  string tmp_mf = output_mf[chrom]+"_"+i_start_chr+"_"+i_end_chr;
  string tmp_mn = output_mn[chrom]+"_"+i_start_chr+"_"+i_end_chr;
  
  cout<<"# Output files:\n";
  cout<<"# 1) "<<tmp_g<<endl;
  cout<<"# 2) "<<tmp_mf<<endl;
  cout<<"# 3) "<<tmp_mn<<endl;
  
  ofstream fileOut_g(tmp_g.c_str());
  ofstream fileOut_mf(tmp_mf.c_str());
  ofstream fileOut_mn(tmp_mn.c_str());
  
  for(i=0;i<n;i++){
    fileOut_g<<Fg(i)<<"\t"<<Ng(i)<<endl;
    fileOut_mf<<Fm(i,0);
    fileOut_mn<<Nm(i,0);
    for(k=1;k<n_annot;k++){
      fileOut_mf<<"\t"<<Fm(i,k);
      fileOut_mn<<"\t"<<Nm(i,k);
    }
    fileOut_mf<<endl;
    fileOut_mn<<endl;
  }
  fileOut_mn.close();
  fileOut_mf.close();
  fileOut_g.close();
  
  
  time_t t2 = time(0);   // get time now
  struct tm * now2 = localtime( & t2 );
  cout <<"\n\tAnalysis ends: ";
  cout << (now2->tm_year + 1900) << '-' 
       << (now2->tm_mon + 1) << '-'
       <<  now2->tm_mday << " at "
       <<  now2->tm_hour <<":"<<now2->tm_min<<":"<<now2->tm_sec
       << endl;
  clock_t toc = clock();
  
  printf("\tTime elapsed: %f seconds\n", (float)(toc - tic) / CLOCKS_PER_SEC);
  
  return EXIT_SUCCESS;
}
