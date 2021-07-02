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
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <algorithm>
#include <ctime>
using namespace std;

#define abs(n) ((n>0) ? n : -n)
#define minimum(n,m) ((n>m) ? m : n)

string limpia(string);

const int MXCEN = 10000;            // Maximum number of centers
const int MXATDIF = 20;            // Maximum number of different atoms
const int MXSHELLAT = 50;        // Maximum number of contractions per atom
const int MXPRIMCENT = 200;       // Maximum number of primitives per center
const int MXFUN = 60000;            // Maximum number of contracted basis functions
const int MXSHELL = 10000;        // Maximum number of shells
const double PI = 3.141592653589793L;
const int MXLBAS = 6;        //    Maximum value of l quantum number in the basis set (for the moment up to "g" functions)
                        //    For higher values of "l" in the basis set, it is mandatory to know the way the functions
                        //    differing in "m" quantum number are stored in TURBOMOLE and their definitions,
                        //    mainly in the signs included in the definitions. This information must be taken into 
                        //    account when storing the density matrix according to DAM order and definitions
int i, j, k, ii, jj, ki, kj, i1, j11, klin, znatdif[MXATDIF];
double zn[MXCEN];

bool lzn = false, lcoord = false, lbasis = false, lclosedsh = false, luhf = false;
string atomnms[104]={"  ","H","HE","LI","BE","B","C","N","O","F"
            ,"NE","NA","MG","AL","SI","P","S","CL","AR","K"
            ,"CA","SC","TI","V","CR","MN","FE","CO","NI","CU"
            ,"ZN","GA","GE","AS","SE","BR","KR","RB","SR","Y"
            ,"ZR","NB","MO","TC","RU","RH","PD","AG","CD","IN"
            ,"SN","SB","TE","I","XE","CS","BA","LA","CE","PR"
            ,"ND","PM","SM","EU","GD","TB","DY","HO","ER","TM"
            ,"YB","LU","HF","TA","W","RE","OS","IR","PT","AU"
            ,"HG","TL","PB","BI","PO","AT","RN","FR","RA","AC"
            ,"TH","PA","U","NP","PU","AM","CM","BK","CF","ES"
            ,"FM","MD","NO","LW"};
string lvals[9]={"s","p","d","f","g","h","i","j","k"};

// Function prototypes
int seekzn(string s);
int seekl(string s);
int seekiatdif(int natdif, int n);
string editinputalpha(string s);

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
#endif


// Main

int main(int argc,char *argv[]) {
    int firstalphaocc, firstbetaocc, lastalphaocc, lastbetaocc, len, nbasis, nalphaorb, nbetaorb, ncen, 
        kntshell = 0, izn;
    double alphaocc, betaocc;
    double x[MXCEN], y[MXCEN], z[MXCEN], facts[MXLBAS+1], fact2l1[MXLBAS+1];
    double *OM;
    double *dmat;
    bool leecontrol;
    string TMfiles, inpath, newprojectname, outpath, projectname, s, s2;
    
    clock_t startTime = clock();
    // Default values for inpath and outpath
    inpath = "." + slash;
    outpath = "." + slash;
    
    if (argc<2) {
        cout << "Please, supply a project name: ";
        cin >> projectname ;
        
    } else{
        projectname = argv[1]; // First argument (if available): project name (name of files .basis, ... without extension) 
        if (argc > 2) {
            inpath = argv[2]; // Second argument (if available): full path to TURBOMOLE .basis, .coord and .mos files
            if (argc > 3) {
                outpath = argv[3]; // Third argument (if available): full path to output files directory
                if (argc > 4){
                    newprojectname = argv[4]; // Third argument (if available): New project name
                }
            }
            else outpath = "." + slash;
        }
        else {
            inpath = "." + slash;
            outpath = "." + slash;
        }
    }
    if (argc < 5)
        newprojectname = projectname;
    
    s = outpath + newprojectname + "-TURBOMOLE_interface.out";
    ofstream outimportfile(s.c_str(),ios::out);
    if (!outimportfile){
        cerr << "Unable to open file " << s << endl ;
        exit(1);
    }
    
    outimportfile << "argv = " << argv[0] ;
    for (i = 1 ; i < argc ; i++){
        outimportfile << " " << argv[i] ;
    }
    outimportfile << endl;

    TMfiles = inpath + projectname;
    s = TMfiles + ".coord";
    ifstream inputcoord;
    inputcoord.open(s.c_str(),ios::in);
    if (!inputcoord) {
        inputcoord.open("coord",ios::in);
    }
    if (!inputcoord) {
        cerr << "Unable to open file " << s <<  endl ;
        outimportfile << "Unable to open file " << s <<  endl ;
        exit(1);
    }

    s = TMfiles + ".basis";
    ifstream inputbasis;
    inputbasis.open(s.c_str(),ios::in);
    if (!inputbasis) {
        inputbasis.open("basis",ios::in);
    }
    if (!inputbasis) {
        cerr << "Unable to open file " << s << endl ;
        outimportfile << "Unable to open file " << s << endl ;
        exit(1);
    }

    s = TMfiles + ".mos";
    ifstream inputalpha, inputbeta;
    inputalpha.open(s.c_str(),ios::in);
    if (!inputalpha) {
        inputalpha.open("mos",ios::in);
    }

    if (inputalpha) lclosedsh = true;
    else lclosedsh = false;

    s = TMfiles + ".control";
    ifstream inputcontrol;
    inputcontrol.open(s.c_str(),ios::in);
    if (!inputcontrol) {
        inputcontrol.open("control",ios::in);
    }
    
    if (!inputcontrol) {
        leecontrol = false;
    }
    else{
        leecontrol = true;
    }

    if (!lclosedsh){    // If the calculation is not closed shell seek for alpha and beta files 
        s = TMfiles + ".alpha";
        inputalpha.open(s.c_str(),ios::in);
        if (!inputalpha) inputalpha.open("alpha",ios::in);
        if (!inputalpha) luhf = false;
        else{
            s = TMfiles + ".beta";
            inputbeta.open(s.c_str(),ios::in);
            if (!inputbeta) inputbeta.open("beta",ios::in);
            if (!inputbeta) luhf = false;
            else luhf = true;
            if (!leecontrol){
                cerr << "Occupancies are required for uhf calculations and control file is not found." << endl
                    << "Please make control file available. Stop." << endl ;
            }
        }
    }
    if (!( lclosedsh || luhf )) {
        cerr << "Unable to open file(s) with molecular orbitals coefficients. Stop " << endl ;
        outimportfile << "Unable to open file(s) with molecular orbitals coefficients. Stop " << endl ;
        exit(1);
    }

    s = outpath + newprojectname + ".ggbs";
    ofstream ggbsfile(s.c_str(),ios::out);
    if (!ggbsfile) {
        cerr << "Unable to open file " << s << endl ;
        outimportfile << "Unable to open file " << s << endl ;
        exit(1);
    }

    s = outpath + newprojectname + ".den";
    ofstream denfile(s.c_str(),ios::out);
    if (!denfile) {
        cerr << "Unable to open file " << s << endl ;
        outimportfile << "Unable to open file " << s << endl ;
        exit(1);
    }
    
//    Scans controlfile if available
    if (leecontrol){
        while(getline(inputcontrol,s)){
            if( !(s.find("$closed shell")==string::npos) ) {
                getline(inputcontrol,s);
                len = s.length();
                char *tokenPrt, *ptr = new char [len+1], *newtoken;
                s.copy(ptr,len,0);
                ptr[len] = 0;
                tokenPrt = strtok_s(ptr," ",&newtoken);
                tokenPrt = strtok_s(NULL,"-",&newtoken);
                firstalphaocc = atoi(tokenPrt);
                tokenPrt = strtok_s(NULL,"(",&newtoken);
                lastalphaocc = atoi(tokenPrt);
                tokenPrt = strtok_s(NULL,")",&newtoken);
                alphaocc = atof(tokenPrt);
                break;
            }
            if( !(s.find("$alpha shell")==string::npos) ) {
                getline(inputcontrol,s);
                len = s.length();
                char *tokenPrt, *ptr = new char [len+1], *newtoken;
                s.copy(ptr,len,0);
                ptr[len] = 0;
                tokenPrt = strtok_s(ptr," ",&newtoken);
                tokenPrt = strtok_s(NULL,"-",&newtoken);
                firstalphaocc = atoi(tokenPrt);
                tokenPrt = strtok_s(NULL,"(",&newtoken);
                lastalphaocc = atoi(tokenPrt);
                tokenPrt = strtok_s(NULL,")",&newtoken);
                alphaocc = atof(tokenPrt);
            }
            if( !(s.find("$beta shell")==string::npos) ) {
                getline(inputcontrol,s);
                len = s.length();
                char *tokenPrt, *ptr = new char [len+1], *newtoken;
                s.copy(ptr,len,0);
                ptr[len] = 0;
                tokenPrt = strtok_s(ptr," ",&newtoken);
                tokenPrt = strtok_s(NULL,"-",&newtoken);
                firstbetaocc = atoi(tokenPrt);
                tokenPrt = strtok_s(NULL,"(",&newtoken);
                lastbetaocc = atoi(tokenPrt);
                tokenPrt = strtok_s(NULL,")",&newtoken);
                betaocc = atof(tokenPrt);
                break;
            }
        }
    }
//     Scans coord file for center coordinates and atom symbols
    ncen = 0;
    double xc=0.L, yc=0.L, zc=0.L, zntot=0.;
    int natdif = 0;
    bool latdif = true;
    while(getline(inputcoord,s)&&!lcoord){
        len = s.length();
        if( !(s.find("$coord")==string::npos) ) {
            lcoord = true;
            
            while(getline(inputcoord,s)&&!lzn){
                if ( !(s.find("$")==string::npos) ) 
                    lzn = true;
                else {
                    len = s.length();
                    char *tokenPrt, *ptr = new char [len+1], *newtoken;
                    s.copy(ptr,len,0);
                    ptr[len] = 0;
                    tokenPrt = strtok_s(ptr," ",&newtoken);
                    x[ncen] = atof(tokenPrt);
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    y[ncen] = atof(tokenPrt);
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    z[ncen] = atof(tokenPrt);
                    tokenPrt = strtok_s(NULL," ",&newtoken); 
                    s2 = string(tokenPrt);
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    zn[ncen] = seekzn(s2);
                    latdif = true;
                    for (j = 0 ; j < ncen ; j++){
                        if (zn[ncen] == zn[j]){
                            latdif = false;
                            break;
                        }
                    }
                    if (latdif) natdif++;
                    zntot += zn[ncen];
                    xc += zn[ncen] * x[ncen];
                    yc += zn[ncen] * y[ncen];
                    zc += zn[ncen] * z[ncen];
                    ncen++;
                    if (ncen > MXCEN){
                        outimportfile << "Number of centers higher than maximum allowable " << MXCEN << endl ;
                        outimportfile << "Redefine parameter MXCEN and recompile " << endl ;
                        exit(1);
                    }
                    delete [] ptr;
                }
            }
        }
    }
    if (!leecontrol){    // Takes a number of electrons equal to the nuclear charge and closhed shell occupancies
        int nelec = int(zntot+0.5L);
        if (nelec%2 != 0){
            cerr << "Odd number of electrons. Occupancies are required and control file is not found." << endl
                    << "Please make control file available. Stop." << endl ;
            exit(1);
        }
        firstalphaocc = 1;
        lastalphaocc = nelec/2;
        alphaocc = 2.L;
    }
    if (natdif > MXATDIF){
        outimportfile << "Highest number of different atoms exceeded" << endl;
        outimportfile << "Current value = " << natdif << endl;
        outimportfile << "Maximum allowed = " << MXATDIF << endl;
        outimportfile << "Change the value of parameter MXATDIF and recompile" << endl;
    }
//    redefines the coordinates with respect to the center of positive charges and writes the to file .ggbs
    xc = xc / zntot; yc = yc / zntot; zc = zc / zntot;
//    xc = 0.; yc = 0.; zc = 0.;  // Keeps original coordinates without translation of origin to center of positive charges
    ggbsfile << ncen << endl ;
    ggbsfile.setf(ios::scientific,ios::floatfield);
    for(i=0; i<ncen;i++){
        izn = int(zn[i]+0.5L);
        ggbsfile << setprecision(15) << x[i]-xc << " " << y[i]-yc << " " << z[i]-zc << " " ;
        ggbsfile << izn << endl ;
    }
//     Reads the basis set from basis file 
    facts[0] = .5L * sqrt(PI);
    fact2l1[0] = 1.L;
    for(i = 1; i < MXLBAS+1 ; i++){ facts[i] = facts[i-1] * .5L * (2*i+1); fact2l1[i] = fact2l1[i-1]*(2*i+1);}
    int *ncontr = new int[natdif];
    double (*primexp)[MXPRIMCENT] = new double [MXATDIF][MXPRIMCENT]; 
    double (*cfcontr)[MXPRIMCENT] = new double [MXATDIF][MXPRIMCENT];
    int (*nprimit)[MXSHELLAT] = new int [MXATDIF][MXSHELLAT];
    int (*lvec)[MXSHELLAT] = new int [MXATDIF][MXSHELLAT];
    int icontr = 0, kntprim = 0;
    int *lmax = new int[natdif];
    bool first_asterisk = false, second_asterisk = false;
    int kntatdif = -1;
    while(getline(inputbasis,s)&&!first_asterisk){
        if (strncmp(s.c_str(),"*",1) == 0){
            first_asterisk = true; 
            break;
        }
    }
    while(getline(inputbasis,s)){
        if (strncmp(s.c_str(),"$end",4) == 0) break;
        if (strncmp(s.c_str(),"*",1) == 0){
            if (second_asterisk){
                second_asterisk = false;
                first_asterisk = true;
                continue;
            } else{
                first_asterisk = false;
                second_asterisk = true;
                icontr = 0; 
                kntprim = 0;
                continue;
            }
        }
        if (first_asterisk) {
            len = s.length();
            char *tokenPrt, *ptr = new char [len+1], *newtoken;
            s.copy(ptr,len,0);
            ptr[len] = 0;
            tokenPrt = strtok_s(ptr," ",&newtoken);
            s2 = string(tokenPrt);
            kntatdif++;
            lmax[kntatdif] = 0;
            znatdif[kntatdif] = seekzn(s2);
            delete [] ptr;
            first_asterisk = false;     // Only reads one record after the first asterisk of each basis
        } 
        else if(second_asterisk) {
            len = s.length();
            char *tokenPrt, *ptr = new char [len+1], *newtoken;
            s.copy(ptr,len,0);
            ptr[len] = 0;
            tokenPrt = strtok_s(ptr," ",&newtoken);
            nprimit[kntatdif][icontr] = atoi(tokenPrt);
            tokenPrt = strtok_s(NULL," ",&newtoken);
            s2 = string(tokenPrt);
            lvec[kntatdif][icontr] = seekl(s2);
            if (lvec[kntatdif][icontr] > lmax[kntatdif])
                lmax[kntatdif] = lvec[kntatdif][icontr];
            delete [] ptr;
            int k0 = kntprim;
            for (int j=0 ; j < nprimit[kntatdif][icontr] ; j++){
                getline(inputbasis,s);
                len = s.length();
                char *tokenPrt2, *ptr2 = new char [len+1];
                s.copy(ptr2,len,0);
                ptr2[len] = 0;
                tokenPrt2 = strtok_s(ptr2," ",&newtoken);
                primexp[kntatdif][kntprim] = atof(tokenPrt2);
                tokenPrt2 = strtok_s(NULL," ",&newtoken);
                if (tokenPrt2 != NULL){
                    if (atof(tokenPrt2) != 0.)
                        cfcontr[kntatdif][kntprim] = atof(tokenPrt2) ;
                    else
                        cfcontr[kntatdif][kntprim] = 1.;
                }
                else{
                    cfcontr[kntatdif][kntprim] = 1.;
                }
                cfcontr[kntatdif][kntprim] = cfcontr[kntatdif][kntprim]* 
                    sqrt(sqrt(4.L*pow(2.L*primexp[kntatdif][kntprim],2*lvec[kntatdif][icontr]+3))
                        / facts[lvec[kntatdif][icontr]]);
                kntprim++;
                delete [] ptr2;
            }
            double aux = 0.L;
            for (int i=0 ; i < nprimit[kntatdif][icontr] ; i++){
                for (int j=0 ; j < nprimit[kntatdif][icontr] ; j++){
                    aux += cfcontr[kntatdif][k0+i] * cfcontr[kntatdif][k0+j]  
                      / ( 2.L * pow(primexp[kntatdif][k0+i]+primexp[kntatdif][k0+j],lvec[kntatdif][icontr]+1.5L));
                }
            }
            aux = 1.L / sqrt(facts[lvec[kntatdif][icontr]] * aux);
            for (int i=0 ; i < nprimit[kntatdif][icontr] ; i++){
                cfcontr[kntatdif][k0+i] *= aux;
            }
            icontr++;
            ncontr[kntatdif] = icontr;
        }
    }
    int iatdif, ncontrtot = 0, lvectot[MXSHELL];
    nbasis = 0;
    for (int icen = 0 ; icen < ncen ; icen++){
        ggbsfile.unsetf(ios::scientific);
        iatdif = seekiatdif(natdif, icen);  
        ggbsfile << endl << ncontr[iatdif] << endl;
        for(int l = 0 ; l <= lmax[iatdif] ; l++){
            int kntprimexp = 0, kntcfcontr = 0;
            for (int icontr = 0 ; icontr < ncontr[iatdif] ; icontr++){
                if (lvec[iatdif][icontr] == l){
                    lvectot[ncontrtot] = lvec[iatdif][icontr];
                    ncontrtot++;
                    if (ncontrtot > MXSHELL){
                        outimportfile << "Highest allowable number of shells exceeded. Increase parameter MXSHELL and recompile" << endl;
                        exit(1);
                    }
                    nbasis += 2 * lvec[iatdif][icontr] + 1; 
                    if (nbasis > MXFUN){
                        outimportfile << "Highest allowable number of basis functions exceeded. Increase parameter MXFUN and recompile" 
                            << endl; 
                        exit(1);
                    }
                    ggbsfile << nprimit[iatdif][icontr] << " " << lvec[iatdif][icontr] << endl;
                    ggbsfile.setf(ios::scientific,ios::floatfield);
                    ggbsfile << setprecision(15);
                    for (int j = 0 ; j < nprimit[iatdif][icontr] ; j++){ 
                        ggbsfile << primexp[iatdif][kntprimexp++] << " ";
                        if ((j+1)%5==0) ggbsfile << endl;
                    }
                    if (nprimit[iatdif][icontr]%5 != 0)    ggbsfile << endl;
                    for (int j = 0 ; j < nprimit[iatdif][icontr] ; j++){ 
                        ggbsfile << cfcontr[iatdif][kntcfcontr++] << " ";
                        if ((j+1)%5==0) ggbsfile << endl;
                    }
                    if (nprimit[iatdif][icontr]%5 != 0)    ggbsfile << endl;
                    ggbsfile.unsetf(ios::scientific);
                }
                else{
                    kntprimexp += nprimit[iatdif][icontr];
                    kntcfcontr += nprimit[iatdif][icontr];
                }
            }
        }
        ggbsfile << endl;
    }
    ofstream lvecfile;
    s = outpath + newprojectname + ".lvec";
    lvecfile.open(s.c_str(),ios::out);
    lvecfile.setf(ios::dec);
    for (i = 0 ; i < ncontrtot ; i++){
        lvecfile << lvectot[i] << " " ;
        if ((i+1)%50 == 0) lvecfile << endl;
    }
    if (i%50 != 0) lvecfile << endl;
//     Reads the molecular orbitals from mos file (alpha component for uhf case)
    OM = new double[nbasis*nbasis];
    while(getline(inputalpha,s) && !(strpbrk(s.c_str(),"$#") == NULL)){    };    // skips header records
    int j = 0;
    nalphaorb = 0;
    while(getline(inputalpha,s) && (s.find(string("$"))==string::npos)){
        if (!(s.find(string("nsaos"))==string::npos)) {
            nalphaorb++;
            j = 0;
            getline(inputalpha,s);
        }
        s = editinputalpha(s);
        len = s.length();
        if (len == 0)    // To prevent errors from reading blank records
            continue;
        char *tokenPrt, *ptr = new char [len+1], *newtoken;
        s.copy(ptr,len,0);
        ptr[len] = 0;
        tokenPrt = strtok_s(ptr," ",&newtoken);
        OM[(j++)*nbasis+nalphaorb] = atof(tokenPrt);
        while (((tokenPrt = strtok_s(NULL," ",&newtoken)) != NULL) && (string(tokenPrt).length() != 1)){
            OM[(j++)*nbasis+nalphaorb] = atof(tokenPrt);
        }
    }
    nalphaorb++;
//    Writes TURBOMOLE alpha orbitals to file .TMaorb
//     ofstream TMaorbfile;
//     s = outpath + projectname + ".TMaorb";
//     TMaorbfile.open(s.c_str(),ios::out);
//     if (!TMaorbfile) {
//         cerr << "Unable to open file " << s << endl ;
//         outimportfile << "Unable to open file " << s << endl ;
//         exit(1);
//     }
//     TMaorbfile << nbasis << "  " << lastalphaocc-firstalphaocc+1 << "  " << nalphaorb << endl;
//     TMaorbfile << setprecision(15) ;
//      TMaorbfile.setf(ios::scientific,ios::floatfield);
//     for(j=0 ; j < nalphaorb ; j++){
//         int knt = 0;
//         TMaorbfile << endl ;
//         for (i=0 ; i < ncontrtot ; i++){
//             for (int m = -lvectot[i] ; m <=  lvectot[i] ; m++){
//                 TMaorbfile << OM[(knt++)*nbasis+j] << " ";
//                 if (knt%5 == 0)
//                     TMaorbfile << endl;
//             }
//         }
//         TMaorbfile << endl;
//     }
//     TMaorbfile.close();
//    Reorders the MO from TURBOMOLE to DAM 
    double px, py, pz, zlmv[2*MXLBAS+1]; 
    double (*mvecbasis) = new double [nbasis];
    for(j=0 ; j < nbasis ; j++){
        int knt=0;
        for (i=0 ; i < ncontrtot ; i++){
            if(lvectot[i] == 1){
                py=OM[(knt+1)*nbasis+j];
                pz=OM[(knt+2)*nbasis+j];
                px=OM[knt*nbasis+j];
                OM[knt*nbasis+j]=py;
                OM[(knt+1)*nbasis+j]=pz;
                OM[(knt+2)*nbasis+j]=px;
            }
            else{
                int sgn = ((lvectot[i]%2) == 0 ?  1 : -1);
                for (k=0; k<= 2*lvectot[i]; k++){zlmv[k] = OM[(knt+k)*nbasis+j];}
                for (k = lvectot[i] ; k > 0 ; k--){
                    OM[(knt + sgn*k + lvectot[i])*nbasis+j] = zlmv[2*k];
                    OM[(knt - sgn*k + lvectot[i])*nbasis+j] = zlmv[2*k-1];
                    sgn = -sgn;
                }
                OM[(knt+lvectot[i])*nbasis+j] = zlmv[0];
            }
            for(k = 0; k < 2*lvectot[i]+1; k++){ mvecbasis[knt+k] = 1.L;}
//            Sign changes caused by the different sign conventions in TURBOMOLE harmonics with respect to the standard convention
//            in cases of l = 3 and l = 4. Tested up to l = 4
            if (lvectot[i] == 3) mvecbasis[knt] = -1.L;
            if (lvectot[i] == 4) {mvecbasis[knt+1] = -1.L;  mvecbasis[knt+6] = -1.L; }
            knt += 2*lvectot[i]+1;
        }
    }
    ofstream aorbfile;
    s = outpath + newprojectname + ".GAorba";
    aorbfile.open(s.c_str(),ios::out);
    if (!aorbfile) {
        cerr << "Unable to open file " << s << endl ;
        outimportfile << "Unable to open file " << s << endl ;
        exit(1);
    }
    aorbfile << nbasis << "  " << lastalphaocc-firstalphaocc+1 << "  " << nalphaorb << endl;
    aorbfile << setprecision(15) ;
     aorbfile.setf(ios::scientific,ios::floatfield);
    for(j=0 ; j < nbasis ; j++){
        int knt = 0;
        aorbfile << endl ;
        for (i=0 ; i < ncontrtot ; i++){
            for (int m = -lvectot[i] ; m <=  lvectot[i] ; m++){
                aorbfile << mvecbasis[knt] * OM[(knt++)*nbasis+j] << " ";
                if (knt%5 == 0)
                    aorbfile << endl;
            }
        }
        aorbfile << endl;
    }
    aorbfile.close();
//    Builds the density matrix for the closed shell case and the alpha component for the open shell case
    dmat = new double[nbasis*nbasis];
    double sum;
    for (i=0 ; i < nbasis ; i++){
        for (j=0 ; j < nbasis ; j++){
            sum = 0;
            for (k = firstalphaocc-1; k < lastalphaocc; k++){
                sum += OM[i*nbasis+k] * OM[j*nbasis+k];
            }
            dmat[i*nbasis+j] = sum * alphaocc * mvecbasis[i] * mvecbasis[j];
        }
    }

//     Reads the beta molecular orbitals component for uhf case
    if (luhf){
        while(getline(inputbeta,s) && !(strpbrk(s.c_str(),"$#") == NULL)){    };    // skips header records
        int j = 0;
        nbetaorb = 0;
        while(getline(inputbeta,s) && (s.find(string("$"))==string::npos)){
            if (!(s.find(string("nsaos"))==string::npos)) {
                nbetaorb++;
                j = 0;
                getline(inputbeta,s);
            }
            s = editinputalpha(s);
            len = s.length();
            char *tokenPrt, *ptr = new char [len+1], *newtoken;
            s.copy(ptr,len,0);
            ptr[len] = 0;
            tokenPrt = strtok_s(ptr," ",&newtoken);
            OM[(j++)*nbasis+nbetaorb] = atof(tokenPrt);
            while (((tokenPrt = strtok_s(NULL," ",&newtoken)) != NULL) && (string(tokenPrt).length() != 1)){
                OM[(j++)*nbasis+nbetaorb] = atof(tokenPrt);
            }
        }
        nbetaorb++;
//    Writes TURBOMOLE beta orbitals to file .TMborb
//         ofstream TMborbfile;
//         s = outpath + projectname + ".TMborb";
//         TMborbfile.open(s.c_str(),ios::out);
//         if (!TMborbfile) {
//             cerr << "Unable to open file " << s << endl ;
//             outimportfile << "Unable to open file " << s << endl ;
//             exit(1);
//         }
//         TMborbfile << nbasis << "  " << lastbetaocc-firstbetaocc+1 << "  " << nbetaorb << endl;
//         TMborbfile << setprecision(15) ;
//         TMborbfile.setf(ios::scientific,ios::floatfield);
//         for(j=0 ; j < nbetaorb ; j++){
//             int knt = 0;
//             TMborbfile << endl ;
//             for (i=0 ; i < ncontrtot ; i++){
//                 for (int m = -lvectot[i] ; m <=  lvectot[i] ; m++){
//                     TMborbfile << OM[(knt++)*nbasis+j] << " ";
//                     if (knt%5 == 0)
//                         TMborbfile << endl;
//                 }
//             }
//             TMborbfile << endl;
//         }
//         TMborbfile.close();
//    Reorders the MO from TURBOMOLE to DAM 
        for(j=0; j<nbasis;j++){
            int knt=0;
            for (i=0 ; i < ncontrtot ; i++){
                if(lvectot[i] == 1){
                    py=OM[(knt+1)*nbasis+j];
                    pz=OM[(knt+2)*nbasis+j];
                    px=OM[knt*nbasis+j];
                    OM[knt*nbasis+j]=py;
                    OM[(knt+1)*nbasis+j]=pz;
                    OM[(knt+2)*nbasis+j]=px;
                }
                else {
                    int sgn = ((lvectot[i]%2) == 0 ?  1 : -1);
                    for (k=0; k<= 2*lvectot[i]; k++){zlmv[k] = OM[(knt+k)*nbasis+j];}
                    for (k = lvectot[i] ; k > 0 ; k--){
                        OM[(knt + sgn*k + lvectot[i])*nbasis+j] = zlmv[2*k];
                        OM[(knt - sgn*k + lvectot[i])*nbasis+j] = zlmv[2*k-1];
                        sgn = -sgn;
                    }
                    OM[(knt+lvectot[i])*nbasis+j] = zlmv[0];
                }
                for(k = 0; k < 2*lvectot[i]+1; k++){ mvecbasis[knt+k] = 1.L;}
//            Sign changes caused by the different sign conventions in TURBOMOLE harmonics with respect to the standard convention
                if (lvectot[i] == 3) mvecbasis[knt] = -1.L;
                if (lvectot[i] == 4) {mvecbasis[knt+1] = -1.L; mvecbasis[knt+6] = -1.L; }
                knt += 2*lvectot[i]+1;
            }
        }
        ofstream borbfile;
        s = outpath + newprojectname + ".GAorbb";
        borbfile.open(s.c_str(),ios::out);
        if (!borbfile) {
            cerr << "Unable to open file " << s << endl ;
            outimportfile << "Unable to open file " << s << endl ;
            exit(1);
        }
        borbfile << nbasis << "  " << lastbetaocc-firstbetaocc+1 << "  " << nbetaorb << endl;
        borbfile << setprecision(15) ;
        borbfile.setf(ios::scientific,ios::floatfield);
        for(j=0 ; j < nbasis ; j++){
            int knt = 0;
            borbfile << endl ;
            for (i=0 ; i < ncontrtot ; i++){
                for (int m = -lvectot[i] ; m <=  lvectot[i] ; m++){
                    borbfile << mvecbasis[knt] * OM[knt*nbasis+j] << " ";
                    knt++;
                    if (knt%5 == 0)
                        borbfile << endl;
                }
            }
            borbfile << endl;
        }
        borbfile.close();
    
//    Adds the beta component to the the density matrix for the open shell case
        double sum;
        for (i=0 ; i < nbasis ; i++){
            for (j=0 ; j < nbasis ; j++){
                sum = 0;
                for (k = firstbetaocc-1; k < lastbetaocc; k++){
                    sum +=  betaocc * OM[i*nbasis+k] * OM[j*nbasis+k] * mvecbasis[i] * mvecbasis[j];
                }
                dmat[i*nbasis+j] += sum;
            }
        }
    }
    
/*   Writes the density matrix to file denfile (.den)   */

    denfile << " " << nbasis << endl;

    denfile << setprecision(15) ;
    denfile.setf(ios::scientific,ios::floatfield);
    for(i=0; i<nbasis;i++){
        for(j=0,klin=1; j<=i; j++, klin++){
            denfile << dmat[i*nbasis+j] << " ";
            if (klin == 5){
                denfile << endl ;
                klin = 0;
            }
        }
        if (klin != 1) denfile << endl ;
    }
    denfile.unsetf(ios::scientific);
    clock_t endTime = clock();
    clock_t clockTicksTaken = endTime - startTime;
    double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;

/*   Prints out some statistics  */
    outimportfile << "TURBOMOLE_interface output" << endl;
    outimportfile << "==========================" << endl<< endl;
    outimportfile << "inpath = " << inpath << endl;
    outimportfile << "outpath = " << outpath << endl << endl;
    
    outimportfile << "PROJECT NAME: " <<  newprojectname << endl << endl;
    outimportfile << "Number of centers: " << ncen << endl;
    if (lclosedsh){ 
        outimportfile << "Number of electrons: " << 2*(lastalphaocc-firstalphaocc+1) << endl;
    }
    if (luhf){
        outimportfile << "Number of alpha electrons: " << lastalphaocc-firstalphaocc+1 << endl;
        outimportfile << "Number of beta electrons: " << lastbetaocc-firstbetaocc+1 << endl;
    }
    outimportfile << "Number of contractions: " << ncontrtot << endl;
    outimportfile << "Number of basis functions: " << nbasis << endl;
    if (lclosedsh){ 
        outimportfile << "Total number of MO: " << nalphaorb << endl ;
        outimportfile << "Number of occupied MO: " << lastalphaocc-firstalphaocc+1 << endl ;
        outimportfile << "Last occupied MO: " << lastalphaocc << endl ;
    }
    else if (luhf){
        outimportfile << "Total number of alpha MO: " << nalphaorb << endl ;
        outimportfile << "Number of occupied alpha MO: " << lastalphaocc-firstalphaocc+1 << endl ;
        outimportfile << "Last occupied alpha MO: " << lastalphaocc << endl ;
        outimportfile << "Total number of beta MO: " << nbetaorb << endl ;
        outimportfile << "Number of occupied beta MO: " << lastbetaocc-firstbetaocc+1 << endl ;
        outimportfile << "Last occupied beta MO: " << lastbetaocc << endl ;
    }
    outimportfile << "Time (in secs) = " << timeInSeconds << endl;
    outimportfile << endl;
}
// End of main

//-----------------------------------------------------------------------------------------------------------------------------

// Function seekzn
int seekzn(string s){
    int zn = 0;
    // trim trailing spaces of s
    size_t endpos = s.find_first_of(" \r\t\0\n");
    if( string::npos != endpos )
    {
        s = s.substr( 0, endpos );
    }
    // set s to uppercase and seek the atom
    transform(s.begin(), s.end(),s.begin(),::toupper);
    for (int i=0 ; i < 104 ; i++){
        if (s.compare(atomnms[i]) == 0){    // only compares the first two characters of s (s may contains trailing blanks)
            zn = i;
            break;
        }
    }
    return zn;
}
// End of function seekzn

//-----------------------------------------------------------------------------------------------------------------------------

// Function seekl
int seekl(string s){
    int l = 0;
    transform(s.begin(), s.end(),s.begin(),::tolower);
    for (int i=0 ; i < 9 ; i++){
        if (s.compare(0,1,string(lvals[i])) == 0){ // only compares the first character of s (s may contains trailing blanks)
            l = i;
            break;
        }
    }
    return l;
}
// End of function seekl

//-----------------------------------------------------------------------------------------------------------------------------

// Function seekiatdif
int seekiatdif(int natdif, int n){
    int j = 0;
    for (int i = 0 ; i < natdif ; i++){
        if (zn[n] == znatdif[i]){
            j = i;
            break;
        }
    }
    return j;
}
// End of function seekiatdif

//-----------------------------------------------------------------------------------------------------------------------------

// Function editinputalpha
string editinputalpha(string s){
    size_t found = 0;
    string s2 = "E";
    while (!((found=s.find(string("D"),found+1))==string::npos)){
        s.replace(found,1,string("E"));
        found += 4;
        s.insert(found," ");
    }
    return s;
}
// End of function editinputalpha

