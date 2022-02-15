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
//  Version of January 2017
//
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstddef>         // std::size_t
#include <stdlib.h>
#include <math.h>
#include <string>
#include <algorithm>
#include <ctime>
#include <sys/stat.h>
using namespace std;

#define abs(n) ((n>0) ? n : -n)
#define minimum(n,m) ((n>m) ? m : n)

inline bool exists_test3 (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

string limpia(string);

const int MXCEN = 10000;            // Maximum number of centers
const int MXATDIF = 20;            // Maximum number of different atoms
const int MXSHELLAT = 50;        // Maximum number of contractions per atom
const int MXPRIMCENT = 200;       // Maximum number of primitives per center
const int MXFUN = 30000;            // Maximum number of contracted basis functions
const int MXSHELL = 20000;        // Maximum number of shells
const double PI = 3.141592653589793L;
const int MXLBAS = 6;        //    Maximum value of l quantum number in the basis set (for the moment up to "g" functions)
                        //    For higher values of "l" in the basis set, it is mandatory to know the way the functions
                        //    differing in "m" quantum number are stored in NWChem and their definitions,
                        //    mainly in the signs included in the definitions. This information must be taken into
                        //    account when storing the density matrix according to DAM order and definitions
int znatdif[MXATDIF];
double zn[MXCEN];

bool lzn = false, lcoord = false, lbasis = false, lclosedsh = true, luhf = false;
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
string editinputMO(string s);

#if defined(_WIN32) || defined(WIN32) 
    static const std::string slash = "\\";
#else
    static const std::string slash = "/";
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
    int fisrtalphaocc, firstbetaocc, lastalphaocc, lastbetaocc, len, nbasis, nalphaorb, norb, nbetaorb, ncen,
        kntshell = 0, izn;
    double x[MXCEN], y[MXCEN], z[MXCEN], facts[MXLBAS+1], fact2l1[MXLBAS+1];
    double *OM;
    double *dmat;
    bool leecontrol;
    string NWCfile, inpath, newprojectname, outpath, projectname, s, s2, movecs;

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
            inpath = argv[2]; // Second argument (if available): full path to NWChem .out and .movecs files
            if (argc > 3) {
                outpath = argv[3]; // Third argument (if available): full path to output files directory
                if (argc > 4){
                    newprojectname = argv[4]; // Fourth argument (if available): New project name
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

    s = outpath + newprojectname + "-NWChem_interface.out";
    ofstream outimportfile(s.c_str(),ios::out);
    if (!outimportfile){
        cerr << "Unable to open file " << s << endl ;
        exit(1);
    }

//cerr << "inpath = " << inpath << endl;

//    outimportfile << "argv = " << argv[0] ;
//    for (int i = 1 ; i < argc ; i++){
//        outimportfile << " " << argv[i] ;
//    }
//    outimportfile << endl;

    NWCfile = inpath + projectname;

//cerr << "NWCfile = " << NWCfile << endl;
    s = NWCfile + ".nwcout";
    ifstream inputfile;
    inputfile.open(s.c_str(),ios::in);
    if (!inputfile) {
        cerr << "Unable to open file " << s <<  endl ;
        outimportfile << "Unable to open file " << s <<  endl ;
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

//    Reads the starting geometry
    double xc=0.L, yc=0.L, zc=0.L, zntot=0.;
    int kntgeom = 0;
    while(getline(inputfile,s)){
        if (strstr(s.c_str(),"Output coordinates")){
//             Scans for center coordinates and atom symbols
            ncen = 0;
            zntot = 0.;
            kntgeom++;
            double unit_transf = 1.;
            if (strstr(s.c_str(),"angstrom")){
                unit_transf = 1.8897259886;
            }
            else if (strstr(s.c_str(),"nanometer")){
                unit_transf = 18.897259886;
            }
            else if (strstr(s.c_str(),"picometer")){
                unit_transf = 0.018897259886;
            }
            getline(inputfile,s);    // Skip one line
            getline(inputfile,s);    // Skip one line
            getline(inputfile,s);        // Skip one line
            while(getline(inputfile,s) && s.length() > 1){
                len = s.length();
                char *tokenPrt, *ptr = new char [len+1], *newtoken;
                s.copy(ptr,len,0);
                ptr[len] = 0;
                tokenPrt = strtok_s(ptr," ",&newtoken);
                tokenPrt = strtok_s(NULL," ",&newtoken);
                s2 = string(tokenPrt);
                int znaux;
                if ( (znaux = seekzn(s2)) > 0){
                    zn[ncen] = (double) znaux;
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    znaux = atof(tokenPrt);
//                     if (zn[ncen] != znaux){    // Compares charge with atomic number
//                         cout << "zn[" << ncen << "] = " << zn[ncen] << " != znaux = " << znaux << endl;
//                         break;
//                     }
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    x[ncen] = atof(tokenPrt) * unit_transf;
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    y[ncen] = atof(tokenPrt) * unit_transf;
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    z[ncen] = atof(tokenPrt) * unit_transf;
                    zntot += zn[ncen];
                    xc += zn[ncen] * x[ncen];
                    yc += zn[ncen] * y[ncen];
                    zc += zn[ncen] * z[ncen];
                    ncen++;
                    if (ncen >= MXCEN){
                        outimportfile << "Number of centers higher than maximum allowable " << MXCEN << endl ;
                        outimportfile << "Redefine parameter MXCEN and recompile " << endl ;
                        exit(1);
                    }
                }
            }
//             cout << "Geometry no. " << kntgeom << endl;
//             for (int i = 0 ; i < ncen ; i++){
//                 cout << "x = " << x[i] << " y = " << y[i] << " z = " << z[i] << " zn = " << zn[i] << endl;
//             }
            break;
        }
    }

    inputfile.clear();
        inputfile.seekg(0, ios::beg);

//     Reads the basis set from basis file
    facts[0] = .5L * sqrt(PI);
    fact2l1[0] = 1.L;
    for(int i = 1; i < MXLBAS+1 ; i++){ facts[i] = facts[i-1] * .5L * (2*i+1); fact2l1[i] = fact2l1[i-1]*(2*i+1);}
    int *ncontr = new int[MXATDIF];
    double (*primexp)[MXPRIMCENT] = new double [MXATDIF][MXPRIMCENT];
    double (*cfcontr)[MXPRIMCENT] = new double [MXATDIF][MXPRIMCENT];
    int (*nprimit)[MXSHELLAT] = new int [MXATDIF][MXSHELLAT];
    int (*lvec)[MXSHELLAT] = new int [MXATDIF][MXSHELLAT];
    int icontr = 0, kntprim = 0;
    int *lmax = new int[MXATDIF];
    int natdif = 0;
    int kntatdif = -1;
    leecontrol = false;
    while(getline(inputfile,s)){
//         outimportfile << "s = " << s << endl;
        if (strstr(s.c_str(),"Basis") && strstr(s.c_str(),"(spherical)")){
            leecontrol = true;
            break;
        }
    }
    if (!leecontrol){
        outimportfile << "No spherical ao basis set found" << endl;
        exit(1);
    }
    leecontrol = false;
    bool lstart = false;
    bool lblank = false;
    kntprim = 0;
    while(getline(inputfile,s)){
//         outimportfile << " len = " << s.length() << " s = " << s << endl;
        if (strstr(s.c_str(),"Summary") && strstr(s.c_str(),"(spherical)") ){
            leecontrol = true;
            break;
        }
        if (!(s.find(string("("))==string::npos) ) {
            len = s.length();
            char *tokenPrt, *ptr = new char [len+1], *newtoken;
            s.copy(ptr,len,0);
            ptr[len] = 0;
            tokenPrt = strtok_s(ptr," ",&newtoken);
            s2 = string(tokenPrt);
//             cout << "aqui hay un atomo: s2 = " << s2 ;
            kntatdif++;
            lmax[kntatdif] = 0;
            znatdif[kntatdif] = seekzn(s2);
//             cout << " no at = " << znatdif[kntatdif] << endl;
            delete [] ptr;
            ncontr[kntatdif] = 0;
            icontr = 0;
            kntprim = 0;
            nprimit[kntatdif][icontr] = 0;
            lstart = true;
            lblank = false;
        }
        if (lstart && s.length() > 1){
            len = s.length();
            char *tokenPrt, *ptr = new char [len+1], *newtoken;
            s.copy(ptr,len,0);
            ptr[len] = 0;
            tokenPrt = strtok_s(ptr," ",&newtoken);
            s2 = string(tokenPrt);
//             cout << "isdigit(" << s2 << ") = " << isdigit(s2[0]) << endl;
            if (isdigit(s2[0])){
                tokenPrt = strtok_s(NULL," ",&newtoken);
                s2 = string(tokenPrt);
                lvec[kntatdif][icontr] = seekl(s2);
//                 cout << "no cuantico l = " << s2 << " " << lvec[kntatdif][icontr] << endl;
                if (lvec[kntatdif][icontr] > lmax[kntatdif])
                    lmax[kntatdif] = lvec[kntatdif][icontr];
                tokenPrt = strtok_s(NULL," ",&newtoken);
                primexp[kntatdif][kntprim] = atof(tokenPrt);
                tokenPrt = strtok_s(NULL," ",&newtoken);
//                 cfcontr[kntatdif][kntprim] = atof(tokenPrt);
                cfcontr[kntatdif][kntprim] = atof(tokenPrt) *
                    sqrt(sqrt(4.L*pow(2.L*primexp[kntatdif][kntprim],2*lvec[kntatdif][icontr]+3))
                        / facts[lvec[kntatdif][icontr]]);
                kntprim++;
                nprimit[kntatdif][icontr]++;
            }
            lblank = false;
        }
        if (lstart && s.length() == 0 && !lblank){
            icontr++;
            ncontr[kntatdif] = icontr;
            nprimit[kntatdif][icontr] = 0;
            lblank = true;
        }
    }
    natdif = kntatdif + 1;
    for (int i = 0 ; i < natdif ; i++){
//         cout << "Num contraidas en atdif " << i << " = " << ncontr[i] << endl ;
        int knt = 0;
        for (int j = 0 ; j < ncontr[i] ; j++){
//             cout << "lvec en contraida no = " << j << " = " << lvec[i][j] << endl ;
//             cout << "Num primitivas  =" << nprimit[i][j] << endl ;
            for (int k = 0 ; k < nprimit[i][j] ; k++){
//                 cout << "exp = " << primexp[i][knt] << " cf = " << cfcontr[i][knt] <<endl ;
                knt++;
            }
        }
    }

    if (!leecontrol){
        outimportfile << "Input file ended reading basis set" << endl;
        exit(1);
    }


    while(getline(inputfile,s)){
        if (strstr(s.c_str(),"Output coordinates")){
//             Scans for center coordinates and atom symbols
            ncen = 0;
            zntot = 0.;
            kntgeom++;
            double unit_transf = 1.;
            if (strstr(s.c_str(),"angstrom")){
                unit_transf = 1.8897259886;
            }
            else if (strstr(s.c_str(),"nanometer")){
                unit_transf = 18.897259886;
            }
            else if (strstr(s.c_str(),"picometer")){
                unit_transf = 0.018897259886;
            }
            getline(inputfile,s);    // Skip one line
            getline(inputfile,s);    // Skip one line
            getline(inputfile,s);        // Skip one line
            while(getline(inputfile,s) && s.length() > 1){
                len = s.length();
                char *tokenPrt, *ptr = new char [len+1], *newtoken;
                s.copy(ptr,len,0);
                ptr[len] = 0;
                tokenPrt = strtok_s(ptr," ",&newtoken);
                tokenPrt = strtok_s(NULL," ",&newtoken);
                s2 = string(tokenPrt);
                int znaux;
                if ( (znaux = seekzn(s2)) > 0){
                    zn[ncen] = (double) znaux;
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    znaux = atof(tokenPrt);
                    if (zn[ncen] != znaux){    // Compares charge with atomic number
//                         cout << "zn[" << ncen << "] = " << zn[ncen] << " != znaux = " << znaux << endl;
                        break;
                    }
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    x[ncen] = atof(tokenPrt) * unit_transf;
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    y[ncen] = atof(tokenPrt) * unit_transf;
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    z[ncen] = atof(tokenPrt) * unit_transf;
                    zntot += zn[ncen];
                    xc += zn[ncen] * x[ncen];
                    yc += zn[ncen] * y[ncen];
                    zc += zn[ncen] * z[ncen];
                    ncen++;
                    if (ncen >= MXCEN){
                        outimportfile << "Number of centers higher than maximum allowable " << MXCEN << endl ;
                        outimportfile << "Redefine parameter MXCEN and recompile " << endl ;
                        exit(1);
                    }
                }
            }
//             cout << "Geometry no. " << kntgeom << endl;
//             for (int i = 0 ; i < ncen ; i++){
//                 cout << "x = " << x[i] << " y = " << y[i] << " z = " << z[i] << " zn = " << zn[i] << endl;
//             }
        }
        else if(strstr(s.c_str(),"output vectors")){
            len = s.length();
            char *tokenPrt, *ptr = new char [len+1], *newtoken;
            s.copy(ptr,len,0);
            ptr[len] = 0;
            tokenPrt = strtok_s(ptr," ",&newtoken);
            tokenPrt = strtok_s(NULL," ",&newtoken);
            tokenPrt = strtok_s(NULL," ",&newtoken);
            tokenPrt = strtok_s(NULL," ",&newtoken);
            movecs = string(tokenPrt);
            string movecsfile;
            if (movecs.find(inpath) != std::string::npos) {
                movecsfile = movecs;
            }
            else{
                movecsfile = inpath+movecs;
            }
            bool existmovecs = exists_test3(movecsfile);
            if (!existmovecs){
                outimportfile << "File " << movecsfile << " with molecular orbitals does not exist" << endl;
                outimportfile << "Operation is aborted" << endl;
                exit(1);
            }
//             std::size_t found = movecs.find_last_of("/\\");
//             if (found != string::npos) movecs = movecs.substr(found+1);
//             outimportfile << "output vectors = " << movecs << endl;
        }
    }
//     cout << "Number of geometries read = " << kntgeom << " Keeps the coordinates of the last one. "<< endl ;

//    redefines the coordinates with respect to the center of positive charges and writes the to file .ggbs
    xc = xc / zntot; yc = yc / zntot; zc = zc / zntot;
//    xc = 0.; yc = 0.; zc = 0.;  // Keeps original coordinates without translation of origin to center of positive charges
    ggbsfile << ncen << endl ;
    ggbsfile.setf(ios::scientific,ios::floatfield);
    for(int i=0; i<ncen;i++){
        izn = int(zn[i]+0.5L);
        ggbsfile << setprecision(15) << x[i]-xc << " " << y[i]-yc << " " << z[i]-zc << " " ;
        ggbsfile << izn << endl ;
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
    stringstream ss;
    ss << nbasis;
    string mov2asc;
    string movecsfile;
    if (movecs.find(inpath) != std::string::npos) {
        movecsfile = movecs;
    }
    else{
        movecsfile = inpath+movecs;
    }

//    outimportfile << "antes de mov2asc: inpath = " << inpath << " outpath = " << outpath << " newprojectname = " << newprojectname << " ss = " << ss.str() << endl;
    mov2asc = "$NWCHEM_TOP/contrib/mov2asc/mov2asc " + ss.str() + " " + movecsfile + " " + outpath + newprojectname + "_txt";
//   outimportfile << "mov2asc = " << mov2asc << endl;
    int indmovasc = system(mov2asc.c_str());
    if (indmovasc != 0){
        outimportfile << "\nError trying to run \"" <<  mov2asc << "\". Error code = " << indmovasc;
        char* nwchemtop = std::getenv("NWCHEM_TOP");
        outimportfile << "\n\nCheck that variable \"NWCHEM_TOP = ";
        if (nwchemtop) outimportfile << nwchemtop;
        outimportfile << "\" is correctly set.\n\nTrying an alternative\n";
        mov2asc = "mov2asc " + ss.str() + " " + movecsfile + " " + outpath + newprojectname + "_txt";
        indmovasc = system(mov2asc.c_str());
        if (indmovasc != 0){
            char* pPath = std::getenv("PATH");
            outimportfile << "\nError trying to run \"" <<  mov2asc << "\". Error code = " << indmovasc
                    << "\n\nCheck that mov2asc utility is installed in your system."
                    << " (It usually resides in $NWCHEM_TOP/contrib/mov2asc/).\n"
                    << "Check that it has been compiled and the executable is available. Otherwise, run make in that directory.\n"
                    << "\nIf mov2asc is installed, set NWCHEM_TOP in the working dir (that from which DAMQT was launched) "
                    << "to the current home directory of nwchem, export it and restart DAMQT.\n"
                    << "\nIf this does not work, add the directory where mov2asc resides to variable PATH in the working dir and restart DAMQT."
                    << "\nCurrent PATH = " << pPath
                    << "\n\nAlternatively, create a symbolic link or copy the mov2asc executable to the working dir.\n";
            outimportfile << endl;
            exit(1);
        }
        else{
            outimportfile << mov2asc << " succesfully run\n\n";
        }
    }

    ofstream lvecfile;
    s = outpath + newprojectname + ".lvec";
    lvecfile.open(s.c_str(),ios::out);
    lvecfile.setf(ios::dec);
    for (int i = 0 ; i < ncontrtot ; i++){
        lvecfile << lvectot[i] << " " ;
        if ((i+1)%50 == 0) lvecfile << endl;
    }
    if (ncontrtot%50 != 0) lvecfile << endl;

//     Reads the molecular orbitals from movecs_txt file

    s = outpath + newprojectname + "_txt";
//     outimportfile << "s = " << s << endl;
    ifstream inputMO;
    inputMO.open(s.c_str(),ios::in);
    if (!inputMO) {
        cerr << "Unable to open file " << s <<  endl ;
        outimportfile << "Unable to open file " << s <<  endl ;
        exit(1);
    }
    bool lMO = false;
    while(getline(inputMO,s)){
        if(strstr(s.c_str(),"ao basis")){
            lMO = true;
            break;
        }
    };    // skips header records
    if (!lMO){
        cerr << "No MO found " << endl ;
        outimportfile << "No MO found in file " << s <<  endl ;
        exit(1);
    }
    getline(inputMO,s);    // Here reads the number of MO sets (1 for RHF, 2 for UHF)
    int nshells = atoi(s.c_str());
//     cout << "num capas = " << nshells << endl;
    if (nshells == 2){
        lclosedsh = false;
        luhf = true;
    }
    getline(inputMO,s);    // Skips the number of MO
    getline(inputMO,s);    // Reads the number of MO in each shell
    len = s.length();
    char *tokenPrt, *ptr = new char [len+1], *newtoken;
    s.copy(ptr,len,0);
    ptr[len] = 0;
    tokenPrt = strtok_s(ptr," ",&newtoken);
    norb = atoi(tokenPrt);
    nalphaorb = norb;
//     cout << " num of alpha MO  = " << nalphaorb << "  nbasis = " << nbasis << endl;
    nbetaorb = 0;
    if (nshells == 2){
        tokenPrt = strtok_s(NULL," ",&newtoken);
        nbetaorb = atoi(tokenPrt);
//         cout << "num of beta MO  = " << nbetaorb << endl;
    }

    OM = new double[nbasis*nbasis];

//     Read occupancies of alpha orbitals
    fisrtalphaocc = 1;
    lastalphaocc = 0;
    double occ[norb];
    int kntocc = 0;
    for (int i = 0 ; i < norb/3 ; i++){
        getline(inputMO,s);
        len = s.length();
        char *tokenPrt, *ptr = new char [len+1], *newtoken;
        s.copy(ptr,len,0);
        ptr[len] = 0;
        tokenPrt = strtok_s(ptr," ",&newtoken);
        occ[kntocc++] = atof(tokenPrt);
        if (occ[kntocc-1] > 0.) lastalphaocc++;
        tokenPrt = strtok_s(NULL," ",&newtoken);
        occ[kntocc++] = atof(tokenPrt);
        if (occ[kntocc-1] > 0.) lastalphaocc++;
        tokenPrt = strtok_s(NULL," ",&newtoken);
        occ[kntocc++] = atof(tokenPrt);
        if (occ[kntocc-1] > 0.) lastalphaocc++;
    }
    if (norb%3 != 0){
        getline(inputMO,s);
        len = s.length();
        char *tokenPrt, *ptr = new char [len+1], *newtoken;
        s.copy(ptr,len,0);
        ptr[len] = 0;
        tokenPrt = strtok_s(ptr," ",&newtoken);
        occ[kntocc++] = atof(tokenPrt);
        if (occ[kntocc-1] > 0.) lastalphaocc++;
        if (norb%3 > 1){
            tokenPrt = strtok_s(NULL," ",&newtoken);
            occ[kntocc++] = atof(tokenPrt);
            if (occ[kntocc-1] > 0.) lastalphaocc++;
        }
    }
//     cout << "Number of alpha occupied orbitals = " << lastalphaocc << endl;
//     cout << "Occupancies = " ;
//     for (int i = 0 ; i < kntocc ; i++){
//         cout << " " << occ[i] ;
//     }
//     cout << endl;
//     Skips orbital energies of alpha orbitals
    for (int i = 0 ; i < norb/3 ; i++){
        getline(inputMO,s);
    }
    if (norb%3 != 0){
        getline(inputMO,s);
    }

//     Read coefficients of alpha MOs
// cout << "OM = " << endl;
    for (int n = 0 ; n < norb ; n++){
        int j = 0;
        for (int i = 0 ; i < norb/3 ; i++){
            getline(inputMO,s);    // Reads the number of MO in each shell
            len = s.length();
            char *tokenPrt, *ptr = new char [len+1], *newtoken;
            s.copy(ptr,len,0);
            ptr[len] = 0;
            tokenPrt = strtok_s(ptr," ",&newtoken);
            OM[(j++)*nbasis+n] = atof(tokenPrt);
// cout << " " << OM[(j-1)*nbasis+n];
            tokenPrt = strtok_s(NULL," ",&newtoken);
            OM[(j++)*nbasis+n] = atof(tokenPrt);
// cout << " " << OM[(j-1)*nbasis+n];
            tokenPrt = strtok_s(NULL," ",&newtoken);
            OM[(j++)*nbasis+n] = atof(tokenPrt);
// cout << " " << OM[(j-1)*nbasis+n];
        }
        if (norb%3 != 0){
            getline(inputMO,s);    // Reads the number of MO in each shell
            len = s.length();
            char *tokenPrt, *ptr = new char [len+1], *newtoken;
            s.copy(ptr,len,0);
            ptr[len] = 0;
            tokenPrt = strtok_s(ptr," ",&newtoken);
            OM[(j++)*nbasis+n] = atof(tokenPrt);
// cout << " " << OM[(j-1)*nbasis+n];
            if (norb%3 > 1){
                tokenPrt = strtok_s(NULL," ",&newtoken);
                OM[(j++)*nbasis+n] = atof(tokenPrt);
// cout << " " << OM[(j-1)*nbasis+norb];
            }
        }
    }
//     cout << endl << "OM = " << endl;
//     for(int j=0 ; j < nbasis ; j++){
//         for (int i=0 ; i < ncontrtot ; i++){
//             cout << " " << OM[i*nbasis+j];
//         }
//         cout << endl;
//     }


//    Reorders the MO from NWChem to DAM (only P functions must be reordered)
    double px, py, pz, zlmv[2*MXLBAS+1];
    for(int j=0 ; j < nbasis ; j++){
        int knt=0;
        for (int i=0 ; i < ncontrtot ; i++){
            if(lvectot[i] == 1){
                py=OM[(knt+1)*nbasis+j];
                pz=OM[(knt+2)*nbasis+j];
                px=OM[knt*nbasis+j];
                OM[knt*nbasis+j]=py;
                OM[(knt+1)*nbasis+j]=pz;
                OM[(knt+2)*nbasis+j]=px;
            }
            if(lvectot[i] == 2){
                OM[(knt+3)*nbasis+j]=-OM[(knt+3)*nbasis+j]; // Changes sign in case of d+1
            }
            if(lvectot[i] == 3){
                OM[(knt+4)*nbasis+j]=-OM[(knt+4)*nbasis+j]; // Changes sign in case of f+1
            }
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
    aorbfile << nbasis << "  " << lastalphaocc-fisrtalphaocc+1 << "  " << norb << endl;
    aorbfile << setprecision(15) ;
     aorbfile.setf(ios::scientific,ios::floatfield);
    for(int j=0 ; j < nbasis ; j++){
        int knt = 0;
        aorbfile << endl ;
        for (int i=0 ; i < ncontrtot ; i++){
            for (int m = -lvectot[i] ; m <=  lvectot[i] ; m++){
                aorbfile << OM[(knt++)*nbasis+j] << " ";
                if (knt%5 == 0)
                    aorbfile << endl;
            }
        }
        aorbfile << endl;
    }
    aorbfile.close();
// //    Builds the density matrix for the closed shell case and the alpha component for the open shell case
    dmat = new double[nbasis*nbasis];
    double sum;
    for (int i=0 ; i < nbasis ; i++){
        for (int j=0 ; j < nbasis ; j++){
            sum = 0;
            for (int k = fisrtalphaocc-1; k < lastalphaocc; k++){
                sum += occ[k] * OM[i*nbasis+k] * OM[j*nbasis+k];
            }
            dmat[i*nbasis+j] = sum;
        }
    }

//     Reads the beta molecular orbitals component for uhf case
    if (luhf){
//         Read occupancies of beta orbitals
        norb = nbetaorb;
        firstbetaocc = 1;
        lastbetaocc = 0;
        int kntocc = 0;
        for (int i = 0 ; i < norb/3 ; i++){
            getline(inputMO,s);
            len = s.length();
            char *tokenPrt, *ptr = new char [len+1], *newtoken;
            s.copy(ptr,len,0);
            ptr[len] = 0;
            tokenPrt = strtok_s(ptr," ",&newtoken);
            occ[kntocc++] = atof(tokenPrt);
            if (occ[kntocc-1] > 0.) lastbetaocc++;
            tokenPrt = strtok_s(NULL," ",&newtoken);
            occ[kntocc++] = atof(tokenPrt);
            if (occ[kntocc-1] > 0.) lastbetaocc++;
            tokenPrt = strtok_s(NULL," ",&newtoken);
            occ[kntocc++] = atof(tokenPrt);
            if (occ[kntocc-1] > 0.) lastbetaocc++;
        }
        if (norb%3 != 0){
            getline(inputMO,s);
            len = s.length();
            char *tokenPrt, *ptr = new char [len+1], *newtoken;
            s.copy(ptr,len,0);
            ptr[len] = 0;
            tokenPrt = strtok_s(ptr," ",&newtoken);
            occ[kntocc++] = atof(tokenPrt);
            if (occ[kntocc-1] > 0.) lastbetaocc++;
            if (norb%3 > 1){
                tokenPrt = strtok_s(NULL," ",&newtoken);
                occ[kntocc++] = atof(tokenPrt);
                if (occ[kntocc-1] > 0.) lastbetaocc++;
            }
        }
//         cout << "Number of beta occupied orbitals = " << lastbetaocc << endl;
//         cout << "Occupancies = " ;
//         for (int i = 0 ; i < kntocc ; i++){
//             cout << " " << occ[i] ;
//         }
//         cout << endl;
//     Skips orbital energies of beta orbitals
        for (int i = 0 ; i < norb/3 ; i++){
            getline(inputMO,s);
        }
        if (norb%3 != 0){
            getline(inputMO,s);
        }

//     Read coefficients of beta MOs
// cout << "OM = " << endl;
        for (int n = 0 ; n < norb ; n++){
            int j = 0;
            for (int i = 0 ; i < norb/3 ; i++){
                getline(inputMO,s);    // Reads the number of MO in each shell
                len = s.length();
                char *tokenPrt, *ptr = new char [len+1], *newtoken;
                s.copy(ptr,len,0);
                ptr[len] = 0;
                tokenPrt = strtok_s(ptr," ",&newtoken);
                OM[(j++)*nbasis+n] = atof(tokenPrt);
// cout << " " << OM[(j-1)*nbasis+n];
                tokenPrt = strtok_s(NULL," ",&newtoken);
                OM[(j++)*nbasis+n] = atof(tokenPrt);
// cout << " " << OM[(j-1)*nbasis+n];
                tokenPrt = strtok_s(NULL," ",&newtoken);
                OM[(j++)*nbasis+n] = atof(tokenPrt);
// cout << " " << OM[(j-1)*nbasis+n];
            }
            if (norb%3 != 0){
                getline(inputMO,s);    // Reads the number of MO in each shell
                len = s.length();
                char *tokenPrt, *ptr = new char [len+1], *newtoken;
                s.copy(ptr,len,0);
                ptr[len] = 0;
                tokenPrt = strtok_s(ptr," ",&newtoken);
                OM[(j++)*nbasis+n] = atof(tokenPrt);
// cout << " " << OM[(j-1)*nbasis+n];
                if (norb%3 > 1){
                    tokenPrt = strtok_s(NULL," ",&newtoken);
                    OM[(j++)*nbasis+n] = atof(tokenPrt);
// cout << " " << OM[(j-1)*nbasis+norb];
                }
            }
        }
//     cout << endl << "OM = " << endl;
//     for(int j=0 ; j < nbasis ; j++){
//         for (int i=0 ; i < ncontrtot ; i++){
//             cout << " " << OM[i*nbasis+j];
//         }
//         cout << endl;
//     }


//    Reorders the MO from NWChem to DAM (only P functions must be reordered)
        for(int j=0 ; j < nbasis ; j++){
            int knt=0;
            for (int i=0 ; i < ncontrtot ; i++){
                if(lvectot[i] == 1){
                    py=OM[(knt+1)*nbasis+j];
                    pz=OM[(knt+2)*nbasis+j];
                    px=OM[knt*nbasis+j];
                    OM[knt*nbasis+j]=py;
                    OM[(knt+1)*nbasis+j]=pz;
                    OM[(knt+2)*nbasis+j]=px;
                }
                if(lvectot[i] == 2){
                    OM[(knt+3)*nbasis+j]=-OM[(knt+3)*nbasis+j]; // Changes sign in case of d+1
                }
                if(lvectot[i] == 3){
                    OM[(knt+4)*nbasis+j]=-OM[(knt+4)*nbasis+j]; // Changes sign in case of f+1
                }
                knt += 2*lvectot[i]+1;
            }
        }
        s = outpath + newprojectname + ".GAorbb";
        aorbfile.open(s.c_str(),ios::out);
        if (!aorbfile) {
            cerr << "Unable to open file " << s << endl ;
            outimportfile << "Unable to open file " << s << endl ;
            exit(1);
        }
        aorbfile << nbasis << "  " << lastalphaocc-fisrtalphaocc+1 << "  " << norb << endl;
        aorbfile << setprecision(15) ;
        aorbfile.setf(ios::scientific,ios::floatfield);
        for(int j=0 ; j < nbasis ; j++){
            int knt = 0;
            aorbfile << endl ;
            for (int i=0 ; i < ncontrtot ; i++){
                for (int m = -lvectot[i] ; m <=  lvectot[i] ; m++){
                    aorbfile << OM[(knt++)*nbasis+j] << " ";
                    if (knt%5 == 0)
                        aorbfile << endl;
                }
            }
            aorbfile << endl;
        }
        aorbfile.close();
    // //    Builds the density matrix for the closed shell case and the alpha component for the open shell case
        double sum;
        for (int i=0 ; i < nbasis ; i++){
            for (int j=0 ; j < nbasis ; j++){
                sum = 0;
                for (int k = fisrtalphaocc-1; k < lastalphaocc; k++){
                    sum += occ[k] * OM[i*nbasis+k] * OM[j*nbasis+k];
                }
                dmat[i*nbasis+j] += sum;    // Adds the contribution to density matrix of beta orbitals
            }
        }
    }    // End of UHF case


/*   Writes the density matrix to file denfile (.den)   */

    denfile << " " << nbasis << endl;

    denfile << setprecision(15) ;
    denfile.setf(ios::scientific,ios::floatfield);
    int klin=1;
    for(int i=0; i<nbasis;i++){
        for(int j=0; j<=i; j++, klin++){
            denfile << dmat[i*nbasis+j] << " ";
            if (klin == 5){
                denfile << endl ;
                klin = 0;
            }
        }
    }
    denfile.unsetf(ios::scientific);
    clock_t endTime = clock();
    clock_t clockTicksTaken = endTime - startTime;
    double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;
// 
// /*   Prints out some statistics  */
    outimportfile << "NWChem_interface output" << endl;
    outimportfile << "==========================" << endl<< endl;
    outimportfile << "inpath = " << inpath << endl;
    outimportfile << "outpath = " << outpath << endl << endl;

    outimportfile << "PROJECT NAME: " <<  newprojectname << endl << endl;
    outimportfile << "Number of centers: " << ncen << endl;
    outimportfile << "Number of geometries read = " << kntgeom << " Keeps the coordinates of the last one. "<< endl ;
    if (lclosedsh){
        outimportfile << "Number of electrons: " << 2*(lastalphaocc-fisrtalphaocc+1) << endl;
    }
    if (luhf){
        outimportfile << "Number of alpha electrons: " << lastalphaocc-fisrtalphaocc+1 << endl;
        outimportfile << "Number of beta electrons: " << lastbetaocc-firstbetaocc+1 << endl;
    }
    outimportfile << "Number of contractions: " << ncontrtot << endl;
    outimportfile << "Number of basis functions: " << nbasis << endl;
    if (lclosedsh){
        outimportfile << "Total number of MO: " << norb << endl ;
        outimportfile << "Number of occupied MO: " << lastalphaocc-fisrtalphaocc+1 << endl ;
        outimportfile << "Last occupied MO: " << lastalphaocc << endl ;
    }
    else if (luhf){
        outimportfile << "Total number of alpha MO: " << nalphaorb << endl ;
        outimportfile << "Number of occupied alpha MO: " << lastalphaocc-fisrtalphaocc+1 << endl ;
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

// Function editinputMO
string editinputMO(string s){
    size_t found = 0;
    string s2 = "E";
    while (!((found=s.find(string("D"),found+1))==string::npos)){
        s.replace(found,1,string("E"));
        found += 4;
        s.insert(found," ");
    }
    return s;
}
// End of function editinputMO

