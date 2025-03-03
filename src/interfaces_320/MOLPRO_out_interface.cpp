//  Copyright 2013-2021, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
//  Guillermo Ramirez, Anmol Kumar, Sachin D. Yeole, Shridhar R. Gadre
// 
//  This file is part of DAM320.
// 
//  DAM320 is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
// 
//  DAM320 is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with DAM320.  If not, see <http://www.gnu.org/licenses/>.
//
//------------------------------------------------------------------------
//
//  Version of September 2017
//
//------------------------------------------------------------------------
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <algorithm>
#include <ctime>
#include <functional> 
#include <cctype>
#include <locale>
#include <vector>
using namespace std;

#define abs(n) ((n>0) ? n : -n)
#define minimum(n,m) ((n>m) ? m : n)

string limpia(string);
int min(int,int);

const int MXCEN = 10000;            // Maximum number of centers
const int MXSHELL = 20000;        // Maximum total number of contractions
const int MXSHELLAT = 50;        // Maximum number of contractions per atom
const int MXPRIMCENT = 300;       // Maximum number of primitives per center
const int MXFUN = 30000;            // Maximum number of contracted basis functions
const int MXPRIMCNTR = 20;        // Maximum number of primitives per contraction
const int MXCONTR = 10;            // Maximum number of contractions sharing the same primitives
const int MXREPIRRED = 8;        // Maximum number of irreducible representations
const int MXLBAS = 7;            //    Maximum value of l quantum number in the basis set (for the moment up to "h" functions)
const int MXSIMCENT = 8;        // Maximum number of centers per symmetry function

const double PI = 3.141592653589793;
static int i,j,k,ii,jj,ki;
static int len, ncen, nbasis, izn;
static double zn[MXCEN], x[MXCEN], y[MXCEN], z[MXCEN];
static bool lMOLPRO = false, lzn = false, lcoord = false, lbasis = false, lden = false, lorba = false,
    lorbb = false, lorbenera = false, lorbenerb = false;
static double facts[50];

static string lvals[9]={"s","p","d","f","g","h","i","j","k"};
static int posbasis[6];// Indices of columns (minus 1) of .out file with: Nr, Sym, Nuc, Type, Exponents, Contraction coefs
static int kntrepirred[MXREPIRRED];
static string ptype[3]={"y","z","x"};
static int *indbases;
static int nfun = 0, indrepir, indstate, indstateaMO, indstatebMO, ipos, numrepir, state, statesim;
static string repirvec[MXREPIRRED];
static int *kntshellat;
static int *kntprimit;
static int *pntprimit;
static int *lvec;
static double *primexp;
static double *cfcontr;

// Function prototypes
void readbasisset(ifstream * inputfile, ofstream * outimportfile, ofstream * ggbsfile, int * indbases, string MOLPROfiles);
void readcoordinates(ifstream * inputfile, ofstream * outimportfile, ofstream * ggbsfile, string MOLPROfiles);
void readoptimizedgeometry(ifstream * inputfile, ofstream * outimportfile, ofstream * ggbsfile, string MOLPROfiles);
void writebasisset(ofstream * ggbsfile);
int seekl(string s);
int seekm(string s);
int coefsim(string s);
bool pairCompare(const std::pair<int, double>& firstElem, const std::pair<int, double>& secondElem);

// Structures
struct basesimstr{
    int numcen;
    int centers[MXSIMCENT];
    int icontr;
    int lval;
    int mval;
    int sign[4];
    string repirred;
};
//static struct basesimstr basesim[MXFUN];

basesimstr *basesim = new basesimstr [MXFUN];

#if defined(_WIN32) || defined(WIN32) 
    static const std::string slash = "\\";    
#else
    static const std::string slash = "/";
#endif

#ifdef __MINGW32__
inline
char *strtok_r(char *str, const char *delim, char **save)
{
    char *res, *last;

    if( !save )
        return strtok(str, delim);
    if( !str && !(str = *save) )
        return NULL;
    last = str + strlen(str);
    if( (*save = res = strtok(str, delim)) )
    {
        *save += strlen(res);
        if( *save < last )
            (*save)++;
        else
            *save = NULL;
    }
    return res;
}
#endif

#ifndef _MSC_VER
inline
char* strtok_s(char* s, const char* delim, char** context)
{
        return strtok_r(s, delim, context);
}
#else
inline
int round(double a)
{
        return floor(a+0.5);
}
#endif
    
// trim from start
static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}

int main(int argc,char *argv[])
{
    string ggbsfilename, MOLPROfiles, ss, s, newprojectname, projectname, inpath, outpath, repirred, typeden;


    clock_t startTime = clock();
    int mden = 0;
    bool seekdenmat = false;
    double *OM, *dmat;
    
    if (argc<2) {
        cout << "Filename (with or without the extension  .out): ";
        cin >> projectname ;
        inpath = "." + slash;
        outpath = "." + slash;
        mden = 1;
    }
    else{
        projectname=argv[1];
        if (argc==2){
            mden=1;
            inpath = "." + slash;
            outpath = "." + slash;
        }

        else {
            if (strncmp(argv[2],"?",1)==0){
                seekdenmat = true;
                cout << "To choose one of the following available density matrices, "
                    << "run the interface with the index of the selected matrix as second entry" 
                    << endl ;
            }
            else{
                mden=atoi(argv[2]);
                if (argc >= 4) {
                    inpath = argv[3];
                    if(inpath.find_last_of(slash) != (inpath.length()-1) ) inpath += slash;
                }
                else inpath = "." + slash;
                if (argc >= 5) {
                    outpath = argv[4];
                    if(outpath.find_last_of(slash) != (outpath.length()-1) ) outpath += slash;
                    if (argc == 6)
                        newprojectname = argv[5];
                }
                else outpath = "." + slash;
            }
        }
    }
    if (projectname.size() > 3 && projectname.substr(projectname.size()-4,4) == ".out"){
        projectname = projectname.substr(0,projectname.size()-4);
    }
    cout << "mden = " << mden << endl;
    if (argc < 5)
        newprojectname = projectname;
    if (mden == 0) 
        seekdenmat = true;
    s = outpath + newprojectname + "-MOLPRO_out_interface.out";
    ofstream outimportfile(s.c_str(),ios::out);
    if (!outimportfile){
        cerr << "In MOLPRO_out_interface: unable to open file " << s << endl ;
        exit(1);
    }

    MOLPROfiles = inpath + projectname;
    s = MOLPROfiles + ".out";
    ifstream inputfile;
    inputfile.open(s.c_str(),ios::in);
    if (!inputfile) {
        cerr << "In MOLPRO_out_interface: unable to open file " << s <<  endl ;
        outimportfile << "In MOLPRO_out_interface: unable to open file " << s <<  endl ;
        exit(1);
    }
    

    s = outpath + newprojectname + ".ggbs";
    ggbsfilename = s;
    ofstream ggbsfile(s.c_str(),ios::out);
    if (!ggbsfile) {
        cerr << "In MOLPRO_out_interface: unable to open file " << s << endl ;
        outimportfile << "In MOLPRO_out_interface: unable to open file " << s << endl ;
        exit(1);
    }
    

    s = outpath + newprojectname + ".den";
    ofstream denfile(s.c_str(),ios::out);
    if (!denfile) {
        cerr << "In MOLPRO_out_interface: unable to open file " << s << endl ;
        outimportfile << "In MOLPRO_out_interface: unable to open file " << s << endl ;
        exit(1);
    }

    while(!lMOLPRO && getline(inputfile,s)){
        if( !(s.find("PROGRAM SYSTEM MOLPRO")==string::npos) ) {
            lMOLPRO = true;
        }
    }
    
    if (!lMOLPRO){
        cerr << "File " << MOLPROfiles + ".out is not a MOLPRO output file" << endl ;
        outimportfile  << "File " << MOLPROfiles + ".out is not a MOLPRO output file" << endl ;
        exit(1);
    }
    
    if(seekdenmat){    // Seeks for the available density matrices
        string dmatrices[10];
        int nummatrices;
        nummatrices = 0;
        while(getline(inputfile,s)){
            transform(s.begin(), s.end(),s.begin(),::tolower);
            if( (!(s.find("density")==string::npos)) && (!(s.find("read from record")==string::npos))){
                string ssaux = s;
                getline(inputfile,s);
                getline(inputfile,s);
                if( (s.find("MATRIX")==string::npos) )
                    continue;
                len = ssaux.length();
                char *tokenPrt, *ptr = new char [len+1], *newtoken;
                ssaux.copy(ptr,len,0);
                ptr[len] = 0;        
                strtok_s(ptr,"=",&newtoken);
                tokenPrt = strtok_s(NULL,")",&newtoken);
                dmatrices[nummatrices++] = tokenPrt;
                if (nummatrices >= 10)
                    outimportfile  << "more than 10 density matrices in " << MOLPROfiles + ".out" << endl
                        << "Only the first 10 available in this interface." << endl ;
            }
            
        }
        cout << nummatrices << " Available matrices in " << MOLPROfiles + ".out" << endl ;
        outimportfile << nummatrices << "  Available matrices in " << MOLPROfiles + ".out" << endl ;
        for (i = 0 ; i < nummatrices ; i++){
                cout << "     " << i+1 << ": " << dmatrices[i] << ")" << endl;
                outimportfile  << "     " << i+1 << ": " << dmatrices[i] << ")" << endl;
        }
        exit(0);
    }

//     SEEK AND READ COORDINATES

    readcoordinates(&inputfile, &outimportfile, &ggbsfile, MOLPROfiles);
    indbases = new int[ncen*MXSHELLAT];

//     SEEK AND READ BASIS SET

    readbasisset(&inputfile, &outimportfile, &ggbsfile, indbases, MOLPROfiles);

    ggbsfile.close();
    bool existorbenera = false;
    bool existorbenerb = false; 
    double (*enerorba) = new double[nbasis];
    double (*enerorbb) = NULL;
    int (*norbair) = new int[numrepir];
    int (*norbbir) = new int[numrepir];
    int norba = 0;
    int norbb = 0;
    int kntden = 1;
    for (i = 0 ; i < numrepir ; i++){
        norbair[i] = 0;
        norbbir[i] = 0;
    }
    cout << endl;
    outimportfile << flush ;
    while(getline(inputfile,s)){
        transform(s.begin(), s.end(),s.begin(),::tolower);
        if (!(s.find("geometry changed")==string::npos)){
            ofstream ggbsfile(ggbsfilename.c_str(),ios::out);
            if (!ggbsfile) {
                cerr << "In MOLPRO_out_interface: unable to open file " << ggbsfilename << endl ;
                outimportfile << "In MOLPRO_out_interface: unable to open file " << ggbsfilename << endl ;
                exit(1);
            }
            lcoord = false;
            lzn = false;
            readcoordinates(&inputfile, &outimportfile, &ggbsfile, MOLPROfiles);
            writebasisset(&ggbsfile);
            ggbsfile.close();
        }
        if (!(s.find("end of geometry optimization")==string::npos)){
            ofstream ggbsfile(ggbsfilename.c_str(),ios::out);
            if (!ggbsfile) {
                cerr << "In MOLPRO_out_interface: unable to open file " << ggbsfilename << endl ;
                outimportfile << "In MOLPRO_out_interface: unable to open file " << ggbsfilename << endl ;
                exit(1);
            }
            lcoord = false;
            lzn = false;
            readoptimizedgeometry(&inputfile, &outimportfile, &ggbsfile, MOLPROfiles);
            writebasisset(&ggbsfile);
            ggbsfile.close();
        }
        if (!lorbenera && !(s.find("orbital energies")==string::npos)){
            if( (s.find("orbital energies for positive spin")==string::npos) ) {
                lorbenerb = true;
                existorbenera = true;
            }
            else{
                enerorbb = new double[nbasis];
                existorbenera = true;
                existorbenerb = true;
            }
            int knt = 0;
            int kntorb = 0;
            int irant = 0;
            while(!lorbenera){
                getline(inputfile,s);
                len = s.length();
                char *tokenPrt=NULL, *ptr = new char [len+1], *tokenPrt2=NULL, *newtoken;
                s.copy(ptr,len,0);
                ptr[len] = 0;
                tokenPrt = strtok_s(ptr," ",&newtoken);
                if (tokenPrt == NULL) continue;     // Loops until a non-empty line is found
                if (atof(tokenPrt) == 0.0){     // Checks that the line contains a number as first element
                    norba = kntorb;
                    if (norba == 0) existorbenera = false;
                    lorbenera = true;
                    break;
                }
                char *pospnt = strchr(tokenPrt,'.');
                if (pospnt == NULL){
                    cerr << "Error reading (alpha) orbital energies" << endl ;
                    outimportfile << "Error reading (alpha) orbital energies" << endl ;
                    existorbenera = false;
                    lorbenera = true;
                    break;
                }
                int ir = atoi(pospnt+1);
                if (ir != irant){
                    if (irant > 0) {
                        norbair[irant-1] = knt;
                    }
                    irant++;
                    knt = 0;
                }
                getline(inputfile,s);    // This should be a line with orbital energies
                s.copy(ptr,len,0);
                ptr[len] = 0;
                tokenPrt = strtok_s(ptr," ",&newtoken);
                if (tokenPrt == NULL || atof(tokenPrt) == 0.0){
                    cerr << "Error reading (alpha) orbital energies" << endl ;
                    outimportfile << "Error reading (alpha) orbital energies" << endl ;
                    exit(1);
                }
                while (tokenPrt != NULL){
                    enerorba[kntorb] = atof(tokenPrt);
                    kntorb++;
                    knt++;
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                }
            }
            if (existorbenera)
                norbair[irant-1] = knt;
            else
                norbair[0] = 0;
            knt = 0;
            kntorb = 0;
            irant = 0;
            while(!lorbenerb){
                getline(inputfile,s);
                len = s.length();
                char *tokenPrt=NULL, *ptr = new char [len+1], *tokenPrt2=NULL, *newtoken;
                s.copy(ptr,len,0);
                ptr[len] = 0;
                tokenPrt = strtok_s(ptr," ",&newtoken);
                if (tokenPrt == NULL) continue;
                if (atof(tokenPrt) == 0.0){
                    norbb = kntorb;
                    if (norbb == 0) existorbenerb = false;
                    lorbenerb = true;
                    break;
                }
                char *pospnt = strchr(tokenPrt,'.');
                if (pospnt == NULL){
                    cerr << "Error reading (beta) orbital energies" << endl ;
                    outimportfile << "Error reading (beta) orbital energies" << endl ;
                    existorbenerb = false;
                    lorbenerb = true;
                }
                int ir = atoi(pospnt+1);
                if (ir != irant){
                    if (irant > 0) {
                        norbbir[irant-1] = knt;
                    }
                    irant++;
                    knt = 0;
                }
                getline(inputfile,s);    // This should be a line with orbital energies
                s.copy(ptr,len,0);
                ptr[len] = 0;
                tokenPrt = strtok_s(ptr," ",&newtoken);
                if (tokenPrt == NULL || atof(tokenPrt) == 0.0){
                    cerr << "Error reading (beta) orbital energies" << endl ;
                    outimportfile << "Error reading (beta) orbital energies" << endl ;
                    exit(1);
                }
                while (tokenPrt != NULL){
                    enerorbb[kntorb] = atof(tokenPrt);
                    kntorb++;
                    knt++;
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                }
            }
            if (existorbenerb)
                norbbir[irant-1] = knt;
            else
                norbbir[0] = 0;
        }
        int *iorbaord = new int[norba];

        if (norba > 0){
            std::vector<std::pair<int,double> > orben;
            int ishft = 0;
            int k = 0;
            for (i = 0 ; i < numrepir ; i++){
                for (j = 0 ; j < norbair[i] ; j++){
                    orben.push_back(std::make_pair(j+ishft,enerorba[k]));
                    k++;   
                }
                ishft += kntrepirred[i];
            }
            std::sort(orben.begin(), orben.end(), pairCompare);    // Sorts the alpha orbitals on ascending energy
            
            i = 0;
            for (std::vector<std::pair<int,double> >::iterator it = orben.begin() ; it != orben.end(); ++it, i++){
                iorbaord[i] = (*it).first;
            }
        }
        int *iorbbord = new int[norbb];

        if (norbb > 0){
            std::vector<std::pair<int,double> > orben;
            int ishft = 0;
            int k = 0;
            for (i = 0 ; i < numrepir ; i++){
                for (j = 0 ; j < norbbir[i] ; j++){
                    orben.push_back(std::make_pair(j+ishft,enerorbb[k]));
                    k++;    
                }
                ishft += kntrepirred[i];
            }
            std::sort(orben.begin(), orben.end(), pairCompare);    // Sorts the beta orbitals on ascending energy
            
            i = 0;
            for (std::vector<std::pair<int,double> >::iterator it = orben.begin() ; it != orben.end(); ++it, i++){
                iorbbord[i] = (*it).first;
            }
        }
        
//    ORBITALS
//         if (!lorba && (!(s.find("orbitals")==string::npos)) && (!(s.find("read from record")==string::npos)) &&  ((s.find("beta")==string::npos))){
        if (!lorba && (!(s.find("orbitals orb read from record")==string::npos)) ||  (!(s.find("orbitals orba read from record")==string::npos))
                ||  (!(s.find("orbitals gaorba read from record")==string::npos))){
/*        seeks for alpha or natural orbitals */
            double orbmet;
            int orbasim;
            OM = new double[nfun*nfun];
            int icen, kntir, offset;
            int kntorb = 0;
            ipos = s.find("state");
            s = s.substr(ipos);
            ipos = s.find(".");
            if (ipos > 0){
                indstateaMO = atoi(s.substr(ipos-1,1).c_str());
                orbasim = atoi(s.substr(ipos+1,1).c_str());
            }
            else{
                indstateaMO = -1;
                orbasim = -1;
            }
            for (i = 0 ; i < nfun*nfun ; i++){
                OM[i] = 0.;
            }
            offset = 0;
            while (!lorba){
                double (*orbsim) = NULL;
                int indrepirold = -1;
                indrepir = 0;
                while(getline(inputfile,s)){
//cerr << rtrim(s).length() << " s = " << s << endl;
                    if (s.length() == 0) continue;
                    if (!(s.find("***")==string::npos)){ 
                        // Stores the block of the previous IR
                        if (kntir != kntrepirred[indrepir]*kntrepirred[indrepir]){
                            cerr << "Error reading " << indrepir << " block of " << indstateaMO << "." << orbasim
                                << " alpha MO matrix" << endl ;
                            lorba = true;
                            break;      // desists reading alpha orbitals
                        }
                        if (!(indrepir == numrepir-1)){
                            cerr << "Number of IR blocks read in alpha orbitals matrix " << indrepir+1 << " wrong, it should be " << numrepir << endl;
                            outimportfile << "Number of IR blocks read in alpha orbitals matrix " << indrepir+1 << 
                                    " wrong, it should be " << numrepir << endl;
                            lorba = true;
                            break;      // desists reading alpha orbitals
                        }
                        kntir = 0;
                        for(i = offset ; i < offset+kntrepirred[indrepir] ; i++){
                            for (j = offset ; j < offset+kntrepirred[indrepir] ; j++){
                                orbmet = orbsim[kntir++] / sqrt((double) basesim[i].numcen);
                                for(ii = 0 ; ii < basesim[i].numcen ; ii++){
                                    OM[(indbases[basesim[i].centers[ii]*MXSHELLAT+basesim[i].icontr]
                                        + basesim[i].mval + basesim[i].lval)*nfun+j] += orbmet * basesim[i].sign[ii];
                                }
                            }
                        }
                        offset += kntrepirred[indrepir];
                        lorba = true;
                        break;
                    }
                    else if (!(s.find("SYMMETRY BLOCK")==string::npos)){
                        indrepirold = indrepir;
                        ipos = s.find(".");
                        indrepir = atoi(s.substr(ipos-1,1).c_str())-1;
                        // Stores the block of the previous IR
                        if (indrepir > 0){
                            if (kntir != kntrepirred[indrepirold]*kntrepirred[indrepirold]){
                                cerr << "Error reading block no. " << indrepir << " of MO matrix of state " << indstateaMO << "." << orbasim
                                    << endl ;
                                lorba = true;
                                break;      // desists reading alpha orbitals
                            }
                            kntir = 0;
                            for(i = offset ; i < offset+kntrepirred[indrepirold] ; i++){
                                for (j = offset ; j < offset+kntrepirred[indrepirold] ; j++){
                                    orbmet = orbsim[kntir++] / sqrt((double) basesim[i].numcen);
                                    for(ii = 0 ; ii < basesim[i].numcen ; ii++){
                                        OM[(indbases[basesim[i].centers[ii]*MXSHELLAT+basesim[i].icontr]
                                            + basesim[i].mval + basesim[i].lval)*nfun+j] += orbmet * basesim[i].sign[ii];
                                    }
                                }
                            }
                            offset += kntrepirred[indrepirold];
                        }
                        state = atoi(s.substr(ipos+1,1).c_str());
                        kntir = 0;
                        if (!orbsim){ 
                            delete orbsim;
                            orbsim = NULL;
                        }
                        orbsim = new double[kntrepirred[indrepir]*kntrepirred[indrepir]];
                    }
                    else if (rtrim(s).length() > 12){      // Reads the block of the current IR
                        len = s.length();
                        char *tokenPrt, *ptr = new char [len+1], *newtoken;
                        s.copy(ptr,len,0);
                        ptr[len] = 0;
                        tokenPrt = strtok_s(ptr," ",&newtoken);
                        orbsim[kntir++] = atof(tokenPrt);
                        while((tokenPrt = strtok_s(NULL," ",&newtoken)) != NULL){
                            orbsim[kntir++] = atof(tokenPrt);
                        }
                        delete [] ptr;
                    }
                    else{
                        continue;
                    }
                }
            }
            kntorb = offset;
            ofstream aorbfile;
            s = outpath + newprojectname + ".GAorba";
            aorbfile.open(s.c_str(),ios::out);
            if (!aorbfile) {
                cerr << "In MOLPRO_out_interface: unable to open file " << s << endl ;
                outimportfile << "In MOLPRO_out_interface: unable to open file " << s << endl ;
                lorba = true;
                continue;
            }
            aorbfile << nbasis << "  " << kntorb << "  " << kntorb << endl;
            aorbfile << setprecision(15) ;
            aorbfile.setf(ios::scientific,ios::floatfield);
            for(j=0 ; j < norba ; j++){    // Stores the alpha orbitals whose energies are supplied (if any) sorted on increasing energy
                aorbfile << endl ;
                for (i=0 ; i < nbasis ; i++){
                    aorbfile << OM[i*kntorb+iorbaord[j]] << " ";
                    if ((i+1)%5 == 0)
                        aorbfile << endl;
                }
                aorbfile << endl;
            }
            int knt = 0;
            int korbair = 0;
            for (k = 0 ; k < numrepir ; k++){    // Stores the remaining alpha orbitals ordered by IR
                if (norba > 0) 
                    korbair = norbair[k];
                knt += korbair;
                for(j=0 ; j < kntrepirred[k]-korbair ; j++){
                    aorbfile << endl ;
                    for (i=0 ; i < nbasis ; i++){
                        aorbfile << OM[i*kntorb+knt] << " ";
                        if ((i+1)%5 == 0)
                            aorbfile << endl;
                    }
                    aorbfile << endl;
                    knt++;
                }
            }
            aorbfile.close();
        }    
        if (!lorbb && (!(s.find("Orbitals orbb read from record")==string::npos))
                || (!(s.find("orbitals gaorbb read from record")==string::npos))){
/*    seeks for beta orbitals */
            double orbmet;
            int orbbsim;
            OM = new double[nfun*nfun];
            int kntorb = 0;
            int icen, kntir, offset;
            ipos = s.find("state");
            s = s.substr(ipos);
            ipos = s.find(".");
            if (ipos > 0){
                indstatebMO = atoi(s.substr(ipos-1,1).c_str());
                orbbsim = atoi(s.substr(ipos+1,1).c_str());
            }
            else{
                indstatebMO = -1;
                orbbsim = -1;
            }
            for (i = 0 ; i < nfun*nfun ; i++){
                OM[i] = 0.;
            }
            offset = 0;
            
            while (!lorbb){
                double (*orbsim) = NULL;
                int indrepirold = -1;
                indrepir = 0;
                while(getline(inputfile,s)){
                    if (s.length() == 0) continue;
                    if (!(s.find("***")==string::npos)){ 
                        // Stores the block of the previous IR
                        if (kntir != kntrepirred[indrepir]*kntrepirred[indrepir]){
                            cerr << "Error reading " << indrepir << " block of " << indstatebMO << "." << orbbsim
                                << " beta MO matrix" << endl ;
                            lorbb = true;
                            break;      // desists reading beta orbitals
                        }
                        if (!(indrepir == numrepir-1)){
                            cerr << "Number of IR blocks read in alpha orbitals matrix " << indrepir+1 << " wrong, it should be " << numrepir << endl;
                            outimportfile << "Number of IR blocks read in alpha orbitals matrix " << indrepir+1 << 
                                    " wrong, it should be " << numrepir << endl;
                            lorbb = true;
                            break;      // desists reading beta orbitals
                        }
                        kntir = 0;
                        for(i = offset ; i < offset+kntrepirred[indrepir] ; i++){
                            for (j = offset ; j < offset+kntrepirred[indrepir] ; j++){
                                orbmet = orbsim[kntir++] / sqrt((double) basesim[i].numcen);
                                for(ii = 0 ; ii < basesim[i].numcen ; ii++){
                                    OM[(indbases[basesim[i].centers[ii]*MXSHELLAT+basesim[i].icontr]
                                        + basesim[i].mval + basesim[i].lval)*nfun+j] += orbmet * basesim[i].sign[ii];
                                }
                            }
                        }
                        offset += kntrepirred[indrepir];
                        lorbb = true;
                        break;
                    }
                    else if (!(s.find("SYMMETRY BLOCK")==string::npos)){
                        indrepirold = indrepir;
                        ipos = s.find(".");
                        indrepir = atoi(s.substr(ipos-1,1).c_str())-1;
                        // Stores the block of the previous IR
                        if (indrepir > 0){
                            if (kntir != kntrepirred[indrepirold]*kntrepirred[indrepirold]){
                                cerr << "Error reading " << indrepirold << " block of " << indstatebMO << "." << orbbsim
                                    << " MO matrix" << endl ;
                                lorbb = true;
                                break;      // desists reading beta orbitals
                            }
                            kntir = 0;
                            for(i = offset ; i < offset+kntrepirred[indrepirold] ; i++){
                                for (j = offset ; j < offset+kntrepirred[indrepirold] ; j++){
                                    orbmet = orbsim[kntir++] / sqrt((double) basesim[i].numcen);
                                    for(ii = 0 ; ii < basesim[i].numcen ; ii++){
                                        OM[(indbases[basesim[i].centers[ii]*MXSHELLAT+basesim[i].icontr]
                                            + basesim[i].mval + basesim[i].lval)*nfun+j] += orbmet * basesim[i].sign[ii];
                                    }
                                }
                            }
                            offset += kntrepirred[indrepirold];
                        }
                        state = atoi(s.substr(ipos+1,1).c_str());
                        kntir = 0;
                        if (!orbsim){ 
                            delete orbsim;
                            orbsim = NULL;
                        }
                        orbsim = new double[kntrepirred[indrepir]*kntrepirred[indrepir]];
                    }
                    else if (rtrim(s).length() > 12){      // Reads the block of the current IR
                        len = s.length();
                        char *tokenPrt, *ptr = new char [len+1], *newtoken;
                        s.copy(ptr,len,0);
                        ptr[len] = 0;
                        tokenPrt = strtok_s(ptr," ",&newtoken);
                        orbsim[kntir++] = atof(tokenPrt);
                        while((tokenPrt = strtok_s(NULL," ",&newtoken)) != NULL){
                            orbsim[kntir++] = atof(tokenPrt);
                        }
                        delete [] ptr;
                    }
                    else{
                        continue;
                    }
                }
            }                
            kntorb = offset;
            ofstream aorbfile;
            s = outpath + newprojectname + ".GAorbb";
            aorbfile.open(s.c_str(),ios::out);
            if (!aorbfile) {
                cerr << "In MOLPRO_out_interface: unable to open file " << s << endl ;
                outimportfile << "In MOLPRO_out_interface: unable to open file " << s << endl ;
                lorbb = true;
                continue;
            }
            aorbfile << nbasis << "  " << kntorb << "  " << kntorb << endl;
            aorbfile << setprecision(15) ;
            aorbfile.setf(ios::scientific,ios::floatfield);
            for(j=0 ; j < norbb ; j++){    // Stores the beta orbitals whose energies are supplied (if any) sorted on increasing energy
                aorbfile << endl ;
                for (i=0 ; i < nbasis ; i++){
                    aorbfile << OM[i*kntorb+iorbbord[j]] << " ";
                    if ((i+1)%5 == 0)
                        aorbfile << endl;
                }
                aorbfile << endl;
            }
            int knt = 0;
            int korbbir = 0;
            for (k = 0 ; k < numrepir ; k++){    // Stores the remaining beta orbitals ordered by IR
                if (norbb > 0) 
                    korbbir = norbbir[k];
                knt += korbbir;
                for(j=0 ; j < kntrepirred[k]-korbbir ; j++){
                    aorbfile << endl ;
                    for (i=0 ; i < nbasis ; i++){
                        aorbfile << OM[i*kntorb+knt] << " ";
                        if ((i+1)%5 == 0)
                            aorbfile << endl;
                    }
                    aorbfile << endl;
                    knt++;
                }
            }
            aorbfile.close();
        }
                
//    DENSITY MATRIX
        if (!lden && (!(s.find("density")==string::npos)) && (!(s.find("read from record")==string::npos))){
            string ssaux = s;
            getline(inputfile,s);
            getline(inputfile,s);
            if (s.find("MATRIX")==string::npos)
                continue;
/*    seeks for the density matrix */
            double densimet;
            dmat = new double[nfun*nfun];
            int icen, kntir, offset;
            if (kntden < mden){
                kntden++;
                continue;
            }
            len = ssaux.length();
            char *tokenPrt, *ptr = new char [len+1], *newtoken;
            ssaux.copy(ptr,len,0);
            ptr[len] = 0;        
            strtok_s(ptr,"=",&newtoken);
            tokenPrt = strtok_s(NULL," ",&newtoken);
            typeden = tokenPrt;
            tokenPrt = strtok_s(NULL,")",&newtoken);
            ss = tokenPrt;
            ipos = ss.find(".");
            indstate = atoi(ss.substr(ipos-1,1).c_str());
            statesim = atoi(ss.substr(ipos+1,1).c_str());
            for (i = 0 ; i < nfun*nfun ; i++){
                dmat[i] = 0.;
            }
            offset = 0;
            delete [] ptr;
            while (!lden){
                double (*dsim) = NULL;
                int indrepirold = -1;
                indrepir = 0;
                while(getline(inputfile,s)){ 
//cerr << rtrim(s).length() << " s = " << s << endl;
                    if (s.length() == 0) continue;
                    if (!(s.find("***")==string::npos)){ 
                        // Stores the block of the last IR
                        if (kntir != kntrepirred[indrepir]*kntrepirred[indrepir]){
                            cerr << "Error reading " << indrepir << " block of " << indstate << "." << statesim
                                << " density matrix" << endl ;
                            exit(1);
                        }
                        if (!(indrepir == numrepir-1)){
                            cerr << "Number of IR blocks read in density matrix " << indrepir+1 << " wrong, it should be " << numrepir << endl;
                            outimportfile << "Number of IR blocks read in density matrix " << indrepir+1 << " wrong, it should be " << numrepir << endl;
                            exit(1);
                        }
                        kntir = 0;
                        for(i = offset ; i < offset+kntrepirred[indrepir] ; i++){
                            for (j = offset ; j < offset+kntrepirred[indrepir] ; j++){
                                densimet = dsim[kntir++] / sqrt((double) basesim[i].numcen * basesim[j].numcen );
                                for(ii = 0 ; ii < basesim[i].numcen ; ii++){
                                    for(jj = 0 ; jj < basesim[j].numcen ; jj++){
                                        dmat[(indbases[basesim[i].centers[ii]*MXSHELLAT+basesim[i].icontr]
                                            + basesim[i].mval + basesim[i].lval)
                                            +(indbases[basesim[j].centers[jj]*MXSHELLAT+basesim[j].icontr]
                                            + basesim[j].mval + basesim[j].lval)*nfun] += 
                                                densimet * basesim[i].sign[ii] * basesim[j].sign[jj];
                                    }
                                }
                            }
                        }
                        lden = true;
                        break;
                    }
                    else if (!(s.find("SYMMETRY BLOCK")==string::npos)){
                        indrepirold = indrepir;
                        ipos = s.find(".");
                        indrepir = atoi(s.substr(ipos-1,1).c_str())-1;
                        
                        // Stores the block of the previous IR
                        if (indrepir > 0){
                            if (kntir != kntrepirred[indrepirold]*kntrepirred[indrepirold]){
                                cerr << "Error reading " << indrepirold << " block of " << indstate << "." << statesim
                                    << " density matrix" << endl ;
                                exit(1);
                            } 
                            kntir = 0;
                            for(i = offset ; i < offset+kntrepirred[indrepirold] ; i++){
                                for (j = offset ; j < offset+kntrepirred[indrepirold] ; j++){
                                    densimet = dsim[kntir++] / sqrt((double) basesim[i].numcen * basesim[j].numcen );
                                    for(ii = 0 ; ii < basesim[i].numcen ; ii++){
                                        for(jj = 0 ; jj < basesim[j].numcen ; jj++){
                                            dmat[(indbases[basesim[i].centers[ii]*MXSHELLAT+basesim[i].icontr]
                                                + basesim[i].mval + basesim[i].lval)
                                                +(indbases[basesim[j].centers[jj]*MXSHELLAT+basesim[j].icontr]
                                                + basesim[j].mval + basesim[j].lval)*nfun] += 
                                                    densimet * basesim[i].sign[ii] * basesim[j].sign[jj];
                                        }
                                    }
                                }
                            }
                            offset += kntrepirred[indrepirold];
                        }
                        state = atoi(s.substr(ipos+1,1).c_str());
                        kntir = 0;
                        if (!dsim){ 
                            delete dsim;
                            dsim = NULL;
                        }
                        dsim = new double[kntrepirred[indrepir]*kntrepirred[indrepir]];
                    }
                    else if (rtrim(s).length() > 12){      // Reads the block of the current IR
                        len = s.length();
                        char *tokenPrt, *ptr = new char [len+1], *newtoken;
                        s.copy(ptr,len,0);
                        ptr[len] = 0;
                        tokenPrt = strtok_s(ptr," ",&newtoken);
                        dsim[kntir++] = atof(tokenPrt);
                        while((tokenPrt = strtok_s(NULL," ",&newtoken)) != NULL){
                            dsim[kntir++] = atof(tokenPrt);
                        }
                        delete [] ptr;
                    }
                    else{
                        continue;
                    }
                }
            }
            denfile << " " << nfun << endl;
            denfile << setprecision(15) ;
            denfile.setf(ios::scientific,ios::floatfield);
            int knt = 0;
            for (i = 0 ; i < nfun ; i++){
                for (j = 0 ; j <= i ; j++){
                    denfile << dmat[j+i*nfun] << " " ;
                    knt++;
                    if ((knt)%5 == 0) denfile << endl;
                }
//                     denfile << endl ;
            }
            denfile.unsetf(ios::scientific);
            denfile.close();
        }
    }
    if (!lbasis){
        cerr << "Basis set data not included in " <<  MOLPROfiles + ".out file " << endl;
        outimportfile << "Basis set data not included in " <<  MOLPROfiles + ".out file "  << endl;
        exit(1);
    }
    if (!lden){
        cerr << "Density matrix number " << mden << " not available" << endl;
        outimportfile << "Density matrix number " << mden << " not available" << endl;
        exit(1);
    }
    outimportfile << flush ;
    clock_t endTime = clock();
    clock_t clockTicksTaken = endTime - startTime;
    double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;
    
/*   Prints out some statistics  */
    outimportfile << "MOLPRO_interface output" << endl;
    outimportfile << flush ;
    outimportfile << "==========================" << endl<< endl;
    outimportfile << "inpath = " << inpath << endl;
    outimportfile << "outpath = " << outpath << endl << endl;
    
    outimportfile << "PROJECT NAME: " <<  newprojectname << endl << endl;
    outimportfile << "Number of centers: " << ncen << endl;
    outimportfile << "Number of contractions: " << nfun << " (";
    outimportfile << flush ;
    int knt = 0;
    for (i = 0 ; i < numrepir-1 ; i++){
        outimportfile << kntrepirred[i] << " " << trim(basesim[knt].repirred) << " + "; 
        knt += kntrepirred[i];
    }

    outimportfile << flush ;
    outimportfile << kntrepirred[numrepir-1] << " " << trim(basesim[knt].repirred) << ")" << endl;
    outimportfile << "Number of basis functions: " << nfun << endl;
    outimportfile << "Density type: " << typeden << endl;
    outimportfile << "State symmetry: " << repirvec[statesim-1] << endl;
    outimportfile << "Number of state: " << indstate << endl;
    outimportfile << "Time (in secs) = " << timeInSeconds << endl;
    outimportfile << flush ;
    outimportfile << endl;
    outimportfile.close();
}


//     SEEK AND READ COORDINATES

void readcoordinates(ifstream * inputfile, ofstream * outimportfile, ofstream * ggbsfile, string MOLPROfiles){
    string s;
    while(getline(*inputfile,s) && !lcoord){
/*    seek for atomic numbers and coordinates of the centers */
        if( !(s.find("ATOMIC COORDINATES")==string::npos) ) {
            lcoord = true;
            ncen = 0;
            double xc=0., yc=0., zc=0., zntot=0.;
            while(!lzn){
                getline(*inputfile,s);
                len = s.length();
                if (len == 0) continue;
                char *tokenPrt, *ptr = new char [len+1], *newtoken;
                s.copy(ptr,len,0);
                ptr[len] = 0;
                tokenPrt = strtok_s(ptr," ",&newtoken);
                if (atoi(tokenPrt) == 0) continue;
                else{
                    if (atoi(tokenPrt) != (ncen+1) ) {
                        cerr << "Error reading geometry" << endl ;
                        *outimportfile << "Error reading geometry" << endl ;
                        exit(1);
                    }
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    zn[ncen] = atof(tokenPrt);
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    x[ncen] = atof(tokenPrt);
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    y[ncen] = atof(tokenPrt);
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    z[ncen] = atof(tokenPrt);
                    zntot += zn[ncen];
                    xc += zn[ncen] * x[ncen];
                    yc += zn[ncen] * y[ncen];
                    zc += zn[ncen] * z[ncen];
                    ncen++;
                    while(getline(*inputfile,s)&&!lzn){
                        len = s.length();
                        if (len == 0) {
                            lzn = true;
                            delete [] ptr;
                            break;
                        }
                        s.copy(ptr,len,0);
                        ptr[len] = 0;
                        tokenPrt = strtok_s(ptr," ",&newtoken);
                        if (atoi(tokenPrt) == 0) {
                            lzn = true;
                            delete [] ptr;
                            break;
                        }
                        if (atoi(tokenPrt) != (ncen+1)) {
                            cerr << "Error reading geometry" << endl ;
                            *outimportfile << "Error reading geometry" << endl ;
                            exit(1);
                        }
                        tokenPrt = strtok_s(NULL," ",&newtoken);
                        tokenPrt = strtok_s(NULL," ",&newtoken);
                        zn[ncen] = atof(tokenPrt);
                        tokenPrt = strtok_s(NULL," ",&newtoken);
                        x[ncen] = atof(tokenPrt);
                        tokenPrt = strtok_s(NULL," ",&newtoken);
                        y[ncen] = atof(tokenPrt);
                        tokenPrt = strtok_s(NULL," ",&newtoken);
                        z[ncen] = atof(tokenPrt);
                        zntot += zn[ncen];
                        xc += zn[ncen] * x[ncen];
                        yc += zn[ncen] * y[ncen];
                        zc += zn[ncen] * z[ncen];
                        ncen++;
                    }
                }
            }
//            redefines the coordinates with respect to the center of positive charges and writes them to file .ggbs
             xc = xc / zntot; yc = yc / zntot; zc = zc / zntot;
//            xc = 0.; yc = 0.; zc = 0.;  // Keeps original coordinates without translation of origin to center of positive charges
            *ggbsfile << ncen << endl ;
            (*ggbsfile).setf(ios::scientific,ios::floatfield);
            for(i=0; i<ncen;i++){
                izn = int(round(zn[i]));
                *ggbsfile << setprecision(15) << x[i]-xc << " " << y[i]-yc << " " << z[i]-zc << " " ;
                *ggbsfile << izn << endl ;
            }
        }
    }
    if (!lcoord){
        cerr << "Geometry data not included in " <<  MOLPROfiles + ".out file " << endl;
        *outimportfile << "Geometry data not included in " <<  MOLPROfiles + ".out file "  << endl;
        exit(1);
    }
}


//     SEEK AND READ OPTIMIZED GEOMETRY

void readoptimizedgeometry(ifstream * inputfile, ofstream * outimportfile, ofstream * ggbsfile, string MOLPROfiles){
    string s;
    float unitsconversion = 1.;
    while(getline(*inputfile,s) && !lcoord){
/*    seek for atomic numbers and coordinates of the centers */
        if( !(s.find("Current geometry")==string::npos) ) {
            lcoord = true;
            ncen = 0;
            double xc=0., yc=0., zc=0., zntot=0.;
            if (!(s.find("Angstrom")==string::npos)) unitsconversion = 1.88973;
            getline(*inputfile,s);  // Reads blank line
            getline(*inputfile,s);  // Reads number of centers
            getline(*inputfile,s);  // Reads comment
            while(!lzn){
                getline(*inputfile,s);
                len = s.length();
                if (len == 0) {
                    continue;
                }
                else{
                    char *tokenPrt, *ptr = new char [len+1], *newtoken;
                    s.copy(ptr,len,0);
                    ptr[len] = 0;
                    tokenPrt = strtok_s(ptr," ",&newtoken);
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    x[ncen] = atof(tokenPrt) * unitsconversion;
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    y[ncen] = atof(tokenPrt) * unitsconversion;
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    z[ncen] = atof(tokenPrt) * unitsconversion;
                    zntot += zn[ncen];
                    xc += zn[ncen] * x[ncen];
                    yc += zn[ncen] * y[ncen];
                    zc += zn[ncen] * z[ncen];
                    ncen++;
                    while(getline(*inputfile,s)&&!lzn){
                        len = s.length();
                        if (len == 0) {
                            lzn = true;
                            delete [] ptr;
                            break;
                        }
                        s.copy(ptr,len,0);
                        ptr[len] = 0;
                        tokenPrt = strtok_s(ptr," ",&newtoken);
                        tokenPrt = strtok_s(NULL," ",&newtoken);
                        x[ncen] = atof(tokenPrt) * unitsconversion;
                        tokenPrt = strtok_s(NULL," ",&newtoken);
                        y[ncen] = atof(tokenPrt) * unitsconversion;
                        tokenPrt = strtok_s(NULL," ",&newtoken);
                        z[ncen] = atof(tokenPrt) * unitsconversion;
                        zntot += zn[ncen];
                        xc += zn[ncen] * x[ncen];
                        yc += zn[ncen] * y[ncen];
                        zc += zn[ncen] * z[ncen];
                        ncen++;
                    }
                }
            }
//            redefines the coordinates with respect to the center of positive charges and writes them to file .ggbs
             xc = xc / zntot; yc = yc / zntot; zc = zc / zntot;
//            xc = 0.; yc = 0.; zc = 0.;  // Keeps original coordinates without translation of origin to center of positive charges
            *ggbsfile << ncen << endl ;
            (*ggbsfile).setf(ios::scientific,ios::floatfield);
            for(i=0; i<ncen;i++){
                izn = int(round(zn[i]));
                *ggbsfile << setprecision(15) << x[i]-xc << " " << y[i]-yc << " " << z[i]-zc << " " ;
                *ggbsfile << izn << endl ;
            }
        }
    }
    if (!lcoord){
        cerr << "Geometry data not included in " <<  MOLPROfiles + ".out file " << endl;
        *outimportfile << "Geometry data not included in " <<  MOLPROfiles + ".out file "  << endl;
        exit(1);
    }
}

//     SEEK AND READ BASIS SET

void readbasisset(ifstream * inputfile, ofstream * outimportfile, ofstream * ggbsfile, int * indbases, string MOLPROfiles){
    string oldrepirred, repirred, s, ss;
    bool skipgetline = false;
    while(getline(*inputfile,s) && !lbasis){
/*    seek for the basis set */
        if( !(s.find("BASIS DATA")==string::npos) ) {
            oldrepirred = "";
            for (i = 0 ; i < MXREPIRRED ; i++){
                kntrepirred[i] = 0;
            }
            indrepir = -1;
            while(getline(*inputfile,s)&&(s.find("1.1")==string::npos)) continue;
            size_t posicion= s.find_first_of("AB");
            posbasis[0] = 0;
            posbasis[1] = posicion;
            posbasis[2] = posbasis[1]+3;
            posicion = s.find_first_of("spd",posbasis[2]);
            posbasis[3] = posicion-2;
            posbasis[4] = posbasis[3]+7;
            posbasis[5] = posbasis[4]+16;
            int kntline = 1;
            getline(*inputfile,ss);  // Intended to check whether next line contains expansion coefficients of large basis sets
            if (!(ss.substr(posbasis[4],posbasis[5]-posbasis[4]).find(".")==string::npos)){  // more than one line for contraction coefficients
                skipgetline = true;
            }
            else{
                skipgetline = false;
                kntline++;
                s.append(ss);
            }

            nbasis = 0;
            int kprim = 0, kcntbas = 1;
            primexp = new double [ncen*MXPRIMCENT];
            cfcontr = new double [ncen*MXPRIMCENT];
            pntprimit = new int [ncen*(MXSHELLAT+1)];
            lvec = new int [ncen*MXSHELLAT];
            kntshellat = new int [ncen];
            kntprimit = new int[ncen];
            double primexpaux[MXPRIMCNTR];
            double (*cfcontraux)[MXPRIMCNTR] = new double [10][MXPRIMCNTR];
            int (*vcen) = new int[ncen];
            int (*sgncen) = new int[ncen];
            int icen, knt, kntcen, kntprim, lvalue, mvalue, nfn;
            bool nextfnt, lstore, lexist;
            len = s.length();
            char *tokenPrt[MXCONTR+1], *ptr = new char [len+1], *newtoken;
            s.copy(ptr,posbasis[4]-2,0);
            ptr[posbasis[4]-2] = 0;
            strtok_s(ptr," ",&newtoken);
            knt = 1;
            while (strtok_s(NULL," ",&newtoken) != NULL) {
                knt++;
            }

            if (knt != 4){
                cerr << "Error reading basis set in line " << kntline << endl ;
                *outimportfile << "Error reading basis set in line " << kntline  << endl ;
                exit(1);
            }
            for (i = 0 ; i < ncen ; i++){
                kntshellat[i] = 0;
                kntprimit[i] = 0;
                pntprimit[i*(MXSHELLAT+1)+0] = 0;
            }
            lstore = false;
            bool last = false;
            while (!lbasis){
                if (s.length() < 4 || !(s.substr(0,4).find("N")==string::npos)){
                    lbasis = true;
                    break;
                }
                char *tokenPrt[MXCONTR+1], *ptr = new char [len+1], *newtoken;
                for (i = 0 ; i < ncen ; i++){
                    vcen[i] = -1;
                }
                kntcen = 0, kntprim = 0;
                icen = atoi(s.substr(posbasis[2],posbasis[3]-posbasis[2]).c_str());
                if (icen == 0){
                    cerr << "Error reading basis: center index lacking" << endl;
                    *outimportfile << "Error reading basis: center index lacking" << endl;
                    exit(1);
                }
                vcen[kntcen] = icen-1;
                lvalue = seekl(s.substr(posbasis[3],posbasis[4]-posbasis[3]));
                if (!(s.substr(posbasis[3],posbasis[4]-posbasis[3]).find_first_of("sz0") == string::npos)){
                    lstore = true;
                }
                if (lvalue < 0){
                    cerr << "Error in seekl: " << s.substr(posbasis[3],posbasis[4]-posbasis[3])
                        << " does not correspond to any allowable value of l" << endl ;
                    *outimportfile    << "Error in seekl: " << s.substr(posbasis[3],posbasis[4]-posbasis[3])
                        << " does not correspond to any allowable value of l" << endl ;
                    exit(1);
                }
                s.copy(ptr,len-posbasis[4]+2,posbasis[4]-2);
                ptr[len-posbasis[4]+2] = 0;
                knt = 0;

                tokenPrt[knt] = strtok_s(ptr," ",&newtoken);

                while (tokenPrt[knt] != NULL) {
                    knt++;
                    if (knt > MXCONTR+1){
                        cerr << "Error reading exponents and coefficients: number of contractions in line "
                            << kntline << " higher than maximum allowed: MXCONTR = " << MXCONTR << endl ;
                        *outimportfile << "Error reading exponents and coefficients: number of contractions in line "
                            << kntline << " higher than maximum allowed: MXCONTR = " << MXCONTR << endl ;
                        exit(1);
                    }
                    tokenPrt[knt] = strtok_s(NULL," ",&newtoken);
                }
                if (knt < 1){
                    cerr << "Error reading exponents and coefficients: insufficient number of contractions in line "
                        << kntline << endl ;
                    *outimportfile << "Error reading exponents and coefficients: insufficient number of contractions"
                            << " in line " << kntline << endl ;
                    exit(1);
                }
                primexpaux[kntprim] = atof(tokenPrt[0]);
                nfn = knt-1;
                for (i = 0 ; i < nfn ; i++){
                    cfcontraux[i][kntprim] = atof(tokenPrt[1+i]);
                }
                kntprim++;
                sgncen[kntcen] = coefsim(s.substr(posbasis[3],1));
                mvalue = seekm(s.substr(posbasis[3]+3,2));
                repirred = s.substr(posbasis[1],posbasis[2]-posbasis[1]);
                if (strcmp(repirred.c_str(),oldrepirred.c_str()) != 0){
                    indrepir = atoi(s.substr(posbasis[1]-2,2).c_str()) - 1;
                    kntrepirred[indrepir]  = 0;
                    repirvec[indrepir] = repirred;
                    oldrepirred = repirred;
                }

                kntcen++;
                if (kntcen > MXSIMCENT){
                    cerr << "Number of centers in symmetry function = " << kntcen
                         << " greater than maximum allowed = " << MXSIMCENT << endl;
                    cerr << "Increase value of parameter MXSIMCENT and recompile " << endl;
                    *outimportfile << "Number of centers in symmetry function = " << kntcen
                                   << " greater than maximum allowed = " << MXSIMCENT << endl;
                    *outimportfile << "Increase value of parameter MXSIMCENT and recompile " << endl;
                    exit(1);
                }
                for (j = 0 ; j < nfn-1 ; j++){
                    if (skipgetline){
                        s = ss;
                    }
                    else{
                        getline(*inputfile,s);
                        kntline++;
                    }
                    getline(*inputfile,ss);
                    kntline++;
                    if ((ss.length() < posbasis[4]) || !(ss.substr(posbasis[4],posbasis[5]-posbasis[4]).find(".")==string::npos)){
                        skipgetline = true;
                    }
                    else{
                        skipgetline = false;
                        s.append(ss);
                    }
                    len = s.length();
                    if (len == 0 || (s.substr(posbasis[2]).find(".")==string::npos)){
                        cerr << "Error reading basis set: no tokens in line " << kntline  << endl ;
                        *outimportfile << "Error reading basis set: no tokens in line " << kntline << endl ;
                        exit(1);
                    }
                    else{
                        icen = atoi(s.substr(posbasis[2],posbasis[3]-posbasis[2]).c_str());
                        if (icen > 0){
                            vcen[kntcen] = icen-1;
                            lvalue = seekl(s.substr(posbasis[3],posbasis[4]-posbasis[3]));
                            if (!(s.substr(posbasis[3],posbasis[4]-posbasis[3]).find_first_of("sz0") == string::npos)){
                                lstore = true;
                            }
                            if (lvalue < 0){
                                cerr << "Error in seekl: " << s.substr(posbasis[3],posbasis[4]-posbasis[3])
                                    << " does not correspond to any allowable value of l" << endl ;
                                *outimportfile    << "Error in seekl: " << s.substr(posbasis[3],posbasis[4]-posbasis[3])
                                    << " does not correspond to any allowable value of l" << endl ;
                                exit(1);
                            }
                            sgncen[kntcen] = coefsim(s.substr(posbasis[3],1));
                            kntcen++;
                            if (kntcen > MXSIMCENT){
                                cerr << "Number of centers in symmetry function = " << kntcen
                                     << " greater than maximum allowed = " << MXSIMCENT << endl;
                                cerr << "Increase value of parameter MXSIMCENT and recompile " << endl;
                                *outimportfile << "Number of centers in symmetry function = " << kntcen
                                               << " greater than maximum allowed = " << MXSIMCENT << endl;
                                *outimportfile << "Increase value of parameter MXSIMCENT and recompile " << endl;
                                exit(1);
                            }
                        }
                        s.copy(ptr,len-posbasis[4]+2,posbasis[4]-2);
                        ptr[len-posbasis[4]+2] = 0;
                        knt = 0;
                        tokenPrt[knt] = strtok_s(ptr," ",&newtoken);
                        while (tokenPrt[knt] != NULL) {
                            knt++;
                            if (knt > MXCONTR+1){
                                cerr << "Error reading exponents and coefficients: number of contractions in line "
                                    << kntline << " higher than maximum allowed: MXCONTR = " << MXCONTR << endl ;
                                *outimportfile << "Error reading exponents and coefficients: number of contractions in line "
                                    << kntline << " higher than maximum allowed: MXCONTR = " << MXCONTR << endl ;
                                exit(1);
                            }
                            tokenPrt[knt] = strtok_s(NULL," ",&newtoken);
                        }
                        if (knt != (nfn+1)){
                            cerr << "Error reading basis set: number of tokens in line " << kntline << " knt = " << knt
                                << " different than number of contractions + 1 = " << nfn << endl ;
                            *outimportfile << "Error reading basis set: number of tokens in line " << kntline
                                << " different than number of contractions + 1 = " << nfn << endl ;
                            exit(1);
                        }
                        primexpaux[kntprim] = atof(tokenPrt[0]);
                        for (i = 0 ; i < nfn ; i++){
                            cfcontraux[i][kntprim] = atof(tokenPrt[1+i]);
                        }
                        kntprim++;
                    }
                }
                nextfnt = false;
                while (!nextfnt){
                    if (skipgetline){
                        s = ss;
                    }
                    else{
                        getline(*inputfile,s);
                    }
                    if (s.length() < 4 || !(s.substr(0,4).find("N")==string::npos)){
                        nextfnt = true;
                        lbasis = true;
                        break;
                    }
                    getline(*inputfile,ss);
                    kntline++;
                    if (((ss.length() > posbasis[5]) && !(ss.substr(posbasis[4],posbasis[5]-posbasis[4]).find(".")==string::npos))
                        || ((ss.length() > posbasis[1]) && (ss.length() < posbasis[5]))){
                        skipgetline = true;
                    }
                    else{
                        skipgetline = false;
                        s.append(ss);
                    }
                    len = s.length();
                    if (len != 0 && (s.substr(0,posbasis[2]).find(".")==string::npos) && s.length() > posbasis[2]){
                        icen = atoi(s.substr(posbasis[2],posbasis[3]-posbasis[2]).c_str());
                        if (icen > 0){
                            vcen[kntcen] = icen-1;
                            lvalue = seekl(s.substr(posbasis[3],posbasis[4]-posbasis[3]));
                            if (!(s.substr(posbasis[3],posbasis[4]-posbasis[3]).find_first_of("sz0") == string::npos)){
                                lstore = true;
                            }
                            if (lvalue < 0){
                                cerr << "Error in seekl: " << s.substr(posbasis[3],posbasis[4]-posbasis[3])
                                    << " does not correspond to any allowable value of l" << endl ;
                                *outimportfile    << "Error in seekl: " << s.substr(posbasis[3],posbasis[4]-posbasis[3])
                                    << " does not correspond to any allowable value of l" << endl ;
                                exit(1);
                            }
                            sgncen[kntcen] = coefsim(s.substr(posbasis[3],1));
                            kntcen++;
                            if (kntcen > MXSIMCENT){
                                cerr << "Number of centers in symmetry function = " << kntcen
                                     << " greater than maximum allowed = " << MXSIMCENT << endl;
                                cerr << "Increase value of parameter MXSIMCENT and recompile " << endl;
                                *outimportfile << "Number of centers in symmetry function = " << kntcen
                                               << " greater than maximum allowed = " << MXSIMCENT << endl;
                                *outimportfile << "Increase value of parameter MXSIMCENT and recompile " << endl;
                                exit(1);
                            }
                        }
                        s.copy(ptr,len-posbasis[4]+2,posbasis[4]-2);
                        ptr[len-posbasis[4]+2] = 0;
                        knt = 0;
                        tokenPrt[knt] = strtok_s(ptr," ",&newtoken);
                        while (tokenPrt[knt] != NULL) {
                            if (strlen(tokenPrt[knt]) > 2) knt++;
                            if (knt > MXCONTR+1){
                                cerr << "Error reading exponents and coefficients: number of contractions in line "
                                    << kntline << " higher than maximum allowed: MXCONTR = " << MXCONTR << endl ;
                                *outimportfile << "Error reading exponents and coefficients: number of contractions in line "
                                    << kntline << " higher than maximum allowed: MXCONTR = " << MXCONTR << endl ;
                                exit(1);
                            }
                            tokenPrt[knt] = strtok_s(NULL," ",&newtoken);
                        }
                        if (knt == 0){        // Center with same single contraction as the previous one
                            continue;
                        }
                        if (knt != (nfn+1)){
                            cerr << "Error reading basis set: number of tokens in line " << kntline
                                << " different than number of contractions + 1 = " << nfn << endl ;
                            *outimportfile << "Error reading basis set: number of tokens in line " << kntline
                                << " different than number of contractions + 1 = " << nfn << endl ;
                            exit(1);
                        }
                        primexpaux[kntprim] = atof(tokenPrt[0]);
                        for (i = 0 ; i < nfn ; i++){
                            cfcontraux[i][kntprim] = atof(tokenPrt[1+i]);
                        }
                        kntprim++;
                    }
                    else{
                        nextfnt = true;
                    }
                }
                for (i = 0 ; i < nfn ; i++){
                    basesim[nfun+i].numcen = kntcen;
                    basesim[nfun+i].lval = lvalue;
                    basesim[nfun+i].mval = mvalue;
                    basesim[nfun+i].repirred = repirred;
                    for (j = 0 ; j < kntcen ; j++){
                        basesim[nfun+i].centers[j] = vcen[j];
                        basesim[nfun+i].sign[j] = sgncen[j];
                    }
                }
                for (i = 0 ; i < nfn ; i++){
                    if (vcen[0] > -1){
                        lexist = false;
                        for (j = 0 ; j < kntshellat[vcen[0]] ; j++){
                            if (lvalue != lvec[vcen[0]*MXSHELLAT+j]) continue;
                            lexist = true;
                            for (k = 0, ki = 0; k < kntprim ; k++){
                                if (cfcontraux[i][k] == 0.) continue;
                                if (primexpaux[k] != primexp[vcen[0]*MXPRIMCENT+pntprimit[vcen[0]*(MXSHELLAT+1)+j]+ki] ||
                                    cfcontraux[i][k] != cfcontr[vcen[0]*MXPRIMCENT+pntprimit[vcen[0]*(MXSHELLAT+1)+j]+ki]){
                                    lexist = false;
                                    break;
                                }
                                ki++;
                            }
                            if (lexist) {
                                basesim[nfun+i].icontr = j;
                                break;
                            }
                        }
                        if (!lexist){
                            basesim[nfun+i].icontr = kntshellat[vcen[0]];
                            for(k = 0 ; k < kntcen ; k++){
                                lvec[vcen[k]*MXSHELLAT+kntshellat[vcen[k]]] = lvalue;
                                kntshellat[vcen[k]] += 1;
                                for(j = 0 ; j < kntprim ; j++){
                                    if (cfcontraux[i][j] != 0.){
                                        primexp[vcen[k]*MXPRIMCENT+kntprimit[vcen[k]]] = primexpaux[j];
                                        cfcontr[vcen[k]*MXPRIMCENT+kntprimit[vcen[k]]] = cfcontraux[i][j];
                                        kntprimit[vcen[k]] += 1;
                                        if (kntprimit[vcen[k]] > MXPRIMCENT){
                                             cerr << "Highest number of primitives per center (" << MXPRIMCENT << ") exceeded " << endl ;
                                             cerr << "Increase parameter MXPRIMCENT in Molpro_out_interface.cpp and recompile " << endl ;
                                             *outimportfile << "Highest number of primitives per center (" << MXPRIMCENT << ") exceeded " << endl ;
                                             *outimportfile << "Increase parameter MXPRIMCENT in Molpro_out_interface.cpp and recompile " << endl ;
                                             exit(1);

                                        }
                                    }
                                }
                                pntprimit[vcen[k]*(MXSHELLAT+1)+kntshellat[vcen[k]]] = kntprimit[vcen[k]];
                            }
                        }
                    }
                }
                nfun = nfun + nfn;
                kntrepirred[indrepir] += nfn;
                if (!lbasis){
                    s.copy(ptr,posbasis[4]-2,0);
                    ptr[posbasis[4]-2] = 0;
                    strtok_s(ptr," ",&newtoken);
                    knt = 1;
                    while (strtok_s(NULL," ",&newtoken) != NULL) {
                        knt++;
                    }
                    if (knt != 4){
                        cerr << "Error reading basis set in line " << kntline << endl ;
                        *outimportfile << "Error reading basis set in line " << kntline  << endl ;
                        exit(1);
                    }
                }
                lstore = false;
            }
            int kntbasis = 0;
            for (i = 0 ; i < ncen ; i++){
                for (j = 0 ; j < kntshellat[i] ; j++){
                    indbases[i*MXSHELLAT+j] = kntbasis;
                    kntbasis += 2*lvec[i*MXSHELLAT+j]+1;
                }
            }
            nbasis = kntbasis;
            numrepir = indrepir+1;
//    Writes the basis set to ggbs file
            facts[0] = .5 * sqrt(PI);
            for(i = 1; i < 50 ; i++){ facts[i] = facts[i-1] * .5 * (2*i+1);}
            writebasisset(ggbsfile);
        }
    }

}

void writebasisset(ofstream * ggbsfile){
    for (i = 0 ; i < ncen ; i++){
        *ggbsfile << endl;
        *ggbsfile << kntshellat[i] << endl;
        for (j = 0 ; j < kntshellat[i] ; j++){
            *ggbsfile << pntprimit[i*(MXSHELLAT+1)+j+1]-pntprimit[i*(MXSHELLAT+1)+j] << " " << lvec[i*MXSHELLAT+j] << endl;
            (*ggbsfile).setf(ios::scientific,ios::floatfield);
            *ggbsfile << setprecision(15);
            for (k = pntprimit[i*(MXSHELLAT+1)+j] ; k < pntprimit[i*(MXSHELLAT+1)+j+1] ; k++){
                *ggbsfile << primexp[i*MXPRIMCENT+k] << " " ;
                if ((k-pntprimit[i*(MXSHELLAT+1)+j]+1)%5==0) *ggbsfile << endl;
            }
            *ggbsfile << endl;
            for (k = pntprimit[i*(MXSHELLAT+1)+j] ; k < pntprimit[i*(MXSHELLAT+1)+j+1] ; k++){
                *ggbsfile << cfcontr[i*MXPRIMCENT+k] * sqrt(sqrt(4.L*pow(2.L*primexp[i*MXPRIMCENT+k],2*lvec[i*MXSHELLAT+j]+3))
                    / facts[lvec[i*MXSHELLAT+j]]) << " " ;
                if ((k-pntprimit[i*(MXSHELLAT+1)+j]+1)%5==0) *ggbsfile << endl;
            }
            *ggbsfile << endl;
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------

// Function seekl
int seekl(string s){
    int l = -1;
    transform(s.begin(), s.end(),s.begin(),::tolower);
    for (int i=0 ; i < 9 ; i++){
        if( !(s.find(lvals[i])==string::npos) ){
            l = i;
            break;
        }
    }
    return l;
}
// End of function seekl

// Function seekm
int seekm(string s){
    int aux;
    if (s.find("  ")!=string::npos){
        return 0;
    }
    else if (s.find("x ")!=string::npos){
        return 1;
    }
    else if (s.find("y ")!=string::npos){
        return -1;
    }
    else if (s.find("z ")!=string::npos){
        return 0;
    }
    else {
        aux = atoi(s.substr(0,1).c_str());
        if (s.substr(1,1).find("-")!=string::npos){
            return -aux;
        }
        else{
            return aux;
        }
    }
}
// End of function seekm

// Function coefsim
int coefsim(string s){
    if (s.find("-")!=string::npos)
        return -1;
    else
        return 1;
}
// End of function coefsim

// Function pairCompare
bool pairCompare(const std::pair<int, double>& firstElem, const std::pair<int, double>& secondElem) {
  return firstElem.second < secondElem.second;

}
// End of function pairCompare

