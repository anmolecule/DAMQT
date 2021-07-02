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
#include <algorithm>
#include <ctime>
#include <string>
#include <sstream>
#include <vector>
#include <regex>
using namespace std;

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

// Function prototypes
int IndicesToLower(int u,int v);
int seekzn(string s);
vector<int> LowerToIndices(int w);
vector<int> sph_l_and_offset(string s);
/*int seekl(string s);
int seekiatdif(int natdif, int n);*/  
vector<string> split (const string &s, char delim);

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


int main(int argc,char *argv[]) {
    bool lbeta = false, lbetaMO = false, lden = false, lmo = false;
    int izn, knt, kntshell, len, nbasis, ncen, natomindex, nden, new_i, new_j, new_k, nexponents, nmo, norbtypes, npqn ;
    string inpath, newprojectname, outpath, Mopacfile, projectname, s, s2;
    double *aMO, *bMO, *den, *exponents, *x, *y, *z, *zn;
    int *atomindex, *lvec, *pqn, *shftindex, *nshells;
    string *orbtypes;
    
    clock_t startTime = clock();
//  Default values for inpath and outpath
    inpath = "." + slash;
    outpath = "." + slash;
    
    if (argc<2) {
        cout << "Please, supply a project name: ";
        cin >> projectname ;
        
    } else{
        projectname = argv[1]; // First argument (if available): project name (name of files .basis, ... without extension) 
        if (argc > 2) {
            inpath = argv[2]; // Second argument (if available): full path to Mopac .aux file
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
    
    s = outpath + newprojectname + "-MOPAC_aux_interface.out";
    ofstream outimportfile(s.c_str(),ios::out);
    if (!outimportfile){
        cerr << "In Mopac_aux_interface: unable to open file " << s << endl ;
        exit(1);
    }
    Mopacfile = inpath + projectname;
    s = Mopacfile + ".aux";
    
    ifstream inputfile;
    inputfile.open(s.c_str(),ios::in);
    if (!inputfile) {
        cerr << "In Mopac_aux_interface: unable to open file " << s <<  endl ;
        outimportfile << "In Mopac_aux_interface: unable to open file " << s <<  endl ;
        exit(1);
    }

    s = outpath + newprojectname + ".sgbsden";
    ofstream sgbsdenfile(s.c_str(),ios::out);
    if (!sgbsdenfile) {
        cerr << "In Mopac_aux_interface: unable to open file " << s << endl ;
        outimportfile << "In Mopac_aux_interface: unable to open file " << s << endl ;
        exit(1);
    }
    
//     Reads the number of centers
    while(getline(inputfile,s)){
        transform(s.begin(), s.end(),s.begin(),::toupper);
        if (strstr(s.c_str(),"ATOM_EL")){
            std::size_t pos1 = s.find("[");
            std::size_t pos2 = s.find("]");
            s2 = s.substr(pos1+1,pos2-pos1-1);
//             cout << "s2 = " << s2 << std::flush << endl;
            ncen = atoi(s2.c_str());
//             cout << "ncen = " << ncen << std::flush << endl;
            break;
        }
    }
    
//     Reads the symbols of centers and stores atomic numbers
    zn = new double [ncen];
    knt = 0;
    while(knt < ncen && getline(inputfile,s)){
        vector<string> v = split (s, ' ');
        for (auto i : v){
            if (knt >= ncen)
                break;
            if (i == "") continue;
//             cout << "knt = " << knt << i ;
            zn[knt] = (double) seekzn(i);
//             cout << " zn[" << i << "] = " << zn[knt] << endl;
            knt++;
        }
//        cout << endl;
    }
        
//     Reads the number of indices relating basis functions to centers
    while(getline(inputfile,s)){
        transform(s.begin(), s.end(),s.begin(),::toupper);
        if (strstr(s.c_str(),"AO_ATOMINDEX")){
            std::size_t pos1 = s.find("[");
            std::size_t pos2 = s.find("]");
            s2 = s.substr(pos1+1,pos2-pos1-1);
//             cout << "s2 = " << s2 << endl;
            natomindex = atoi(s2.c_str());
//             cout << "natomindex = " << natomindex << endl;
            break;
        }
    }
//     Reads the indices relating basis functions to centers
    atomindex = new int [natomindex];
    knt = 0;
    while(knt < natomindex && getline(inputfile,s)){
        vector<string> v = split (s, ' ');
        for (auto i : v){
            if (knt >= natomindex)
                break;
            if (i == "") continue;
//             cout << "knt = " << knt << i ;
            atomindex[knt] = atoi(i.c_str());
//             cout << " atomindex[" << i << "] = " << atomindex[knt] << endl;
            knt++;
        }
//         cout << endl;
    }
            
//     Reads the number of orbital types
    while(getline(inputfile,s)){
        transform(s.begin(), s.end(),s.begin(),::toupper);
        if (strstr(s.c_str(),"ATOM_SYMTYPE")){
            std::size_t pos1 = s.find("[");
            std::size_t pos2 = s.find("]");
            s2 = s.substr(pos1+1,pos2-pos1-1);
//             cout << "s2 = " << s2 << std::flush << endl;
            norbtypes = atoi(s2.c_str());
//             cout << "norbtypes = " << norbtypes << std::flush << endl;
            break;
        }
    }
//     Reads the orbital types
    orbtypes = new string [norbtypes];
    shftindex = new int [norbtypes];
    nshells = new int [ncen];
    lvec = new int [norbtypes];
    knt = 0;
    kntshell = -1;
    while(knt < norbtypes && getline(inputfile,s)){
        vector<string> v = split (s, ' ');
        for (auto i : v){
            if (knt >= norbtypes)
                break;
//            cout << "knt = " << knt << " i = " << i << std::flush << endl;
            if (i == "") continue;
            orbtypes[knt] = i;
            vector<int> w = sph_l_and_offset(i);
            lvec[knt] = w[0];
            shftindex[knt] = w[1];
            if (lvec[knt] == 0){
                kntshell++;
                nshells[kntshell] = 1;
            }
            else if(lvec[knt] == 1){
                nshells[kntshell] = 2;
            }
            else if(lvec[knt] == 2){
                nshells[kntshell] = 3;
            }

//             cout << " orbtypes[" << knt << "] = " << orbtypes[knt] << " l = " << lvec[knt] << " shift = " << shftindex[knt] << std::flush << endl;
            knt++;
        }
//         cout << std::flush << endl;
    }
//    cout << "norbtypes = " << norbtypes << endl;

    if (!(kntshell == ncen-1)){
        outimportfile << "Error: number of nshells elements = " << kntshell+1 << " different from number of centers = "
            << ncen << " Stop." << endl;
        exit(1);
    }

//     Reads the shell exponents
    while(getline(inputfile,s)){
//cout << "s = " << s << std::flush << endl;
        transform(s.begin(), s.end(),s.begin(),::toupper);
        if (strstr(s.c_str(),"AO_ZETA")){
            std::size_t pos1 = s.find("[");
            std::size_t pos2 = s.find("]");
            s2 = s.substr(pos1+1,pos2-pos1-1);
//             cout << "s2 = " << s2 << std::flush << endl;
            nexponents = atoi(s2.c_str());
//             cout << "nexponents = " << nexponents << std::flush << endl;
            break;
        }
    }
//     Reads the orbital types
    exponents = new double [nexponents];
    knt = 0;
    while(knt < nexponents && getline(inputfile,s)){
        vector<string> v = split (s, ' ');
        for (auto i : v){
            if (knt >= nexponents)
                break;
            if (i == "") continue;
            exponents[knt] = atof(i.c_str());
//             cout << " exponents[" << knt << "] = " << exponents[knt] << std::flush << endl;
            knt++;
        }
//         cout << endl;
    }
                    
//     Reads the principal quantum numbers
    while(getline(inputfile,s)){
        transform(s.begin(), s.end(),s.begin(),::toupper);
        if (strstr(s.c_str(),"ATOM_PQN")){
            std::size_t pos1 = s.find("[");
            std::size_t pos2 = s.find("]");
            s2 = s.substr(pos1+1,pos2-pos1-1);
//             cout << "s2 = " << s2 << std::flush << endl;
            npqn = atoi(s2.c_str());
//             cout << "npqn = " << npqn << std::flush << endl;
            break;
        }
    }
//     Reads the orbital types
    pqn = new int [npqn];
    knt = 0;
    while(knt < npqn && getline(inputfile,s)){
        vector<string> v = split (s, ' ');
        for (auto i : v){
            if (knt >= npqn)
                break;
            if (i == "") continue;
            pqn[knt] = atoi(i.c_str());
//             cout << " pqn[" << knt << "] = " << pqn[knt] << std::flush << endl;
            knt++;
        }
//         cout << endl;
    }
                        
//     Geometry
    while(getline(inputfile,s)){
        transform(s.begin(), s.end(),s.begin(),::toupper);
        if (strstr(s.c_str(),"ATOM_X_OPT:ANGSTROMS")){
            break;
        }
    }
//     Reads the atomic coordinates in angstrom
    x = new double [ncen];
    y = new double [ncen];
    z = new double [ncen];
    double xc=0.L, yc=0.L, zc=0.L, zntot=0.;
    knt = 0;
//     cout << "Geometry" << endl;
    while(knt < ncen && getline(inputfile,s)){
        vector<string> v = split (s, ' ');
        int iaux = 0;
        for (int i = 0 ; i < v.size() ; i++){
            if (v[i] == "") continue;
            if (iaux == 0){
                x[knt] = atof(v[i++].c_str()) * 1.889725989; // conversion to Bohr
                xc += x[knt] * zn[knt];
                iaux++;
                continue;
            }
            if (iaux == 1){
                y[knt] = atof(v[i++].c_str()) * 1.889725989; // conversion to Bohr
                yc += y[knt] * zn[knt];
                iaux++;
                continue;
            }
            if (iaux == 2){
                z[knt] = atof(v[i].c_str()) * 1.889725989; // conversion to Bohr
                zc += z[knt] * zn[knt];
                zntot += zn[knt];
                iaux++;
                continue;
            }
        }
//         cout << " x[" << knt << "] = " << x[knt] << " y[" << knt << "] = " << y[knt] << " z[" << knt << "] = " << z[knt] << endl;
        knt++;
    }

//     Molecular orbitals

    while(getline(inputfile,s)){
        transform(s.begin(), s.end(),s.begin(),::toupper);
        if (strstr(s.c_str(),"EIGENVECTORS") || strstr(s.c_str(),"LMO_VECTORS")){
            lmo = true;
            if (strstr(s.c_str(),"ALPHA")){
                lbetaMO = true;
            }
            std::size_t pos1 = s.find("[");
            std::size_t pos2 = s.find("]");
            s2 = s.substr(pos1+1,pos2-pos1-1);
//             cout << "s2 = " << s2 << std::flush << endl;
            nmo = atoi(s2.c_str());
//             cout << "nmo = " << nmo << std::flush << endl;
            break;
        }
        else if (strstr(s.c_str(),"DENSITY_MATRIX")){
            lmo = false;
            lden = true;
            break;
        }
    }

    if (lmo){
        aMO = new double [nmo];
        for (int i = 0 ; i < nmo ; i++){
            aMO[i] = 0.;
        }
        knt = 0;
//     cout << "Alpha orbitals " << endl;
        int naux = round(sqrt(nmo));
        while(knt < nmo && getline(inputfile,s)){
//         cout << "s = " << s << endl;
            vector<string> v = split (s, ' ');
            for (auto k : v){
                if (knt >= nmo || k == "#")
                     break;
                if (k == "") continue;
//             cout << "k = " << k << std::flush << endl;
                double val = atof(k.c_str());
//             cout << "val = " << val << std::flush << endl;
                int oldi = knt / naux;
                int oldj = knt % naux;
//             cout << "oldi = " << oldi << " oldj = " << oldj << std::flush << endl;
                new_i = oldi + shftindex[oldi];
                new_j = oldj + shftindex[oldj];
                new_k = new_i * naux + new_j;
                aMO[new_k] = val;
                knt++;
            }
        }

        if (lbetaMO){
            lbetaMO = false;
            while(getline(inputfile,s)){
                transform(s.begin(), s.end(),s.begin(),::toupper);
                if (strstr(s.c_str(),"BETA_EIGENVECTORS") || strstr(s.c_str(),"BETA_LMO_VECTORS")){
                    lbetaMO = true;
                    std::size_t pos1 = s.find("[");
                    std::size_t pos2 = s.find("]");
                    s2 = s.substr(pos1+1,pos2-pos1-1);
//             cout << "s2 = " << s2 << std::flush << endl;
                    nmo = atoi(s2.c_str());
//             cout << "nmo = " << nmo << std::flush << endl;
                    break;
                }
                else if (strstr(s.c_str(),"DENSITY_MATRIX")){
                    lbetaMO = false;
                    lden = true;
                    break;
                }
            }

            if (lbetaMO){
                bMO = new double [nmo];
                for (int i = 0 ; i < nmo ; i++){
                    bMO[i] = 0.;
                }
                knt = 0;
//     cout << "Beta orbitals " << endl;
                int naux = round(sqrt(nmo));
                while(knt < nmo && getline(inputfile,s)){
//         cout << "s = " << s << endl;
                    vector<string> v = split (s, ' ');
                    for (auto k : v){
                        if (knt >= nmo || k == "#") break;
                        if (k == "") continue;
//             cout << "k = " << k << endl;
                        double val = atof(k.c_str());
//             cout << "val = " << val << endl;
                        int oldi = knt / naux;
                        int oldj = knt % naux;
//             cout << "oldi = " << oldij[0] << " oldj = " << oldij[1] << endl;
                        new_i = oldi + shftindex[oldi];
                        new_j = oldj + shftindex[oldj];
                        new_k = new_i * naux + new_j;
                        bMO[new_k] = val;
                        knt++;
                    }
                }
            }
        }
    }
                            
//     Density matrix
    if (!lden){
        while(getline(inputfile,s)){
            transform(s.begin(), s.end(),s.begin(),::toupper);
            if (strstr(s.c_str(),"DENSITY_MATRIX")){
                lden = true;
                break;
            }
        }
    }
    if (!lden){
        cerr << "No density matrix found. Stop " << endl ;
        exit(1);
    }
    if (strstr(s.c_str(),"ALPHA")){
        lbeta = true;
    }
    std::size_t pos1 = s.find("[");
    std::size_t pos2 = s.find("]");
    s2 = s.substr(pos1+1,pos2-pos1-1);
//             cout << "s2 = " << s2 << std::flush << endl;
    nden = atoi(s2.c_str());
//             cout << "nden = " << nden << std::flush << endl;
////     Reads the density matrix (lower triangle)
    den = new double [nden];
    for (int i = 0 ; i < nden ; i++){
        den[i] = 0.;
    }
    knt = 0;
//     cout << "Density matrix (lower triangle)" << endl;
    while(knt < nden && getline(inputfile,s)){
//         cout << "s = " << s << std::flush << endl;
        vector<string> v = split (s, ' ');
        for (auto k : v){
            if (knt >= nden || k == "#") break;
            if (k == "") continue;
//             cout << "k = " << k << endl;
            double val = atof(k.c_str());
//             cout << "val = " << val << endl;
            vector<int> oldij = LowerToIndices(knt);
//             cout << "oldi = " << oldij[0] << " oldj = " << oldij[1] << endl;
            new_i = oldij[0] + shftindex[oldij[0]];
            new_j = oldij[1] + shftindex[oldij[1]];
            new_k = IndicesToLower(new_i, new_j);
            den[new_k] = val;
            knt++;
        }
    }
//     for (int i = 0 ; i < nden ; i+=8){
//         for (int j = 0 ; j < 8 ; j++){
//             if (i+j >= nden)
//                 break;
//             cout << den[i+j] << " ";
//         }
//         cout << endl;
//     }

//    Searches for the beta density matrix
    if (lbeta){
        while(getline(inputfile,s)){
            transform(s.begin(), s.end(),s.begin(),::toupper);
            if (strstr(s.c_str(),"BETA_DENSITY_MATRIX")){
                std::size_t pos1 = s.find("[");
                std::size_t pos2 = s.find("]");
                s2 = s.substr(pos1+1,pos2-pos1-1);
//             cout << "s2 = " << s2 << endl;
                nden = atoi(s2.c_str());
//             cout << "nden = " << nden << endl;
                break;
            }
        }
        knt = 0;
    //     cout << "Beta density matrix (lower triangle)" << endl;
        while(knt < nden && getline(inputfile,s)){
    //         cout << "s = " << s << endl;
            vector<string> v = split (s, ' ');
            for (auto k : v){
                if (knt >= nden || k == "#") break;
                if (k == "") continue;
    //             cout << "k = " << k << endl;
                double val = atof(k.c_str());
    //             cout << "val = " << val << endl;
                vector<int> oldij = LowerToIndices(knt);
    //             cout << "oldi = " << oldij[0] << " oldj = " << oldij[1] << endl;
                new_i = oldij[0] + shftindex[oldij[0]];
                new_j = oldij[1] + shftindex[oldij[1]];
                new_k = IndicesToLower(new_i, new_j);
                den[new_k] += val;
                knt++;
            }
        }
    }
//     Writes sgbsden file
    sgbsdenfile << ncen << endl ;
    sgbsdenfile.setf(ios::scientific,ios::floatfield);
//     Transforms coordinates to set the origin in the center of positive charges and writes to file sgbsden
    xc /= zntot ; yc /= zntot ; zc /= zntot ;
    for(int i=0; i<ncen;i++){
        izn = int(zn[i]+0.5L);
        sgbsdenfile << setprecision(15) << x[i]-xc << " " << y[i]-yc << " " << z[i]-zc << " " ;
        sgbsdenfile << izn << " " << nshells[i] << endl ;
    }
//     Writes basis set to file sgbsden
    sgbsdenfile << pqn[0] << " " << lvec[0] << " ";
        sgbsdenfile << setprecision(15) << exponents[0] << endl ;
    for (int i = 1 ; i < npqn ; i++){
        if ((pqn[i] == pqn[i-1]) && (lvec[i] == lvec[i-1]) && (abs(exponents[i]-exponents[i-1]) < 1.e-5)
            && (atomindex[i] == atomindex[i-1]) )
                continue;
        sgbsdenfile << pqn[i] << " " << lvec[i] << " ";
        sgbsdenfile << setprecision(15) << exponents[i] << endl ;
    }
//     Writes the density matrix to file sgbsden
    sgbsdenfile << setprecision(15);
    for (int i = 0 ; i < nden ; i+=8){
        for (int j = 0 ; j < 8 ; j++){
            if (i+j >= nden)
                break;
            sgbsdenfile << den[i+j] << " ";
        }
        sgbsdenfile << endl;
    }
    sgbsdenfile << endl;
    sgbsdenfile.close();
//  Writes files with molecular orbitals
    if (lmo){
        int naux = round(sqrt(nmo));
        s = outpath + newprojectname + ".SLorba";
        ofstream slorbafile(s.c_str(),ios::out);
        if (!slorbafile) {
            cerr << "In Mopac_aux_interface: unable to open file " << s << endl ;
            outimportfile << "In Mopac_aux_interface: unable to open file " << s << endl ;
            exit(1);
        }
        slorbafile << naux << " " << naux << " " << naux << " " << endl;
        slorbafile << setprecision(15) ;
        slorbafile.setf(ios::scientific,ios::floatfield);
        for (int i = 0 ; i < nmo ; i+=8){
            for (int j = 0 ; j < 8 ; j++){
                if (i+j >= nmo)
                    break;
                slorbafile << aMO[i+j] << " ";
            }
            slorbafile << endl;
        }
        slorbafile << endl;
        slorbafile.close();

        if (lbetaMO){
            s = outpath + newprojectname + ".SLorbb";
            ofstream slorbbfile(s.c_str(),ios::out);
            if (!slorbbfile) {
                cerr << "In Mopac_aux_interface: unable to open file " << s << endl ;
                outimportfile << "In Mopac_aux_interface: unable to open file " << s << endl ;
                exit(1);
            }
            slorbbfile << naux << " " << naux << " " << naux << " " << endl;
            slorbbfile << setprecision(15) ;
            slorbbfile.setf(ios::scientific,ios::floatfield);
            for (int i = 0 ; i < nmo ; i+=8){
                for (int j = 0 ; j < 8 ; j++){
                    if (i+j >= nmo)
                        break;
                    slorbbfile << bMO[i+j] << " ";
                }
                slorbbfile << endl;
            }
            slorbbfile << endl;
            slorbbfile.close();
        }
    }
//    Writes summary to outimportfile file

    clock_t endTime = clock();
    clock_t clockTicksTaken = endTime - startTime;
    double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;
    outimportfile << "Mopac_interface output " << endl;
    outimportfile << "====================== " << endl << endl;
    outimportfile << "inpath = " << inpath << endl;
    outimportfile << "outpath = " << outpath << endl;
    outimportfile << "projectname = " << newprojectname << endl;
    outimportfile << "mopacfile = " << Mopacfile << endl;
    outimportfile << "sgbsdenfile = " << outpath + newprojectname + ".sgbsden" << endl;
    outimportfile << "Time (in secs) = " << timeInSeconds << endl;
    outimportfile.close();
    inputfile.close();
}

// Function IndicesToLower
int IndicesToLower(int u,int v){ 
    int i = u+1;
    int j = v+1; 
    if (i >= j)
        return int((i*(i-1)/2.) + j)-1;
    else
        return int((j*(j-1)/2.) + i)-1;
}
// End of function IndicesToLower

//-----------------------------------------------------------------------------------------------------------------------------   

// Function LowerToIndices
vector<int> LowerToIndices(int k){
    vector<int> result;
    int i = int(floor((1+sqrt(8.*k+1))/2.)) - 1;
    int j = int(floor(k - (i*(i+1))/2.));
    result.push_back (i);
    result.push_back (j);
    return result;
}
// End of function LowerToIndices
    
//-----------------------------------------------------------------------------------------------------------------------------
// 
// Function split: splits a string taking a separator characters and return result in a vector
vector<string> split (const string &s, char delim) {
    vector<string> result;
    stringstream ss (s);
    string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }

    return result;
}
// End of function split
//-----------------------------------------------------------------------------------------------------------------------------

// Function sph_l_and_offset: returns the L quantum number and the index offset of a given orbital type
vector<int>  sph_l_and_offset(string s){
    vector<int> result;
    if (s.find("S") !=std::string::npos){
        result.push_back (0);
        result.push_back (0);
        return  result;
    }else if (s.find("PX") !=std::string::npos){
        result.push_back (1);
        result.push_back (2);
        return  result;
    }else if (s.find("PY") !=std::string::npos){
        result.push_back (1);
        result.push_back (-1);
        return  result;
    }else if (s.find("PZ") !=std::string::npos){
        result.push_back (1);
        result.push_back (-1);
        return  result;
    }else if (s.find("X2") !=std::string::npos){
        result.push_back (2);
        result.push_back (4);
        return  result;
    }else if (s.find("XZ") !=std::string::npos){
        result.push_back (2);
        result.push_back (2);
        return  result;
    }else if (s.find("Z2") !=std::string::npos){
        result.push_back (2);
        result.push_back (0);
        return  result;
    }else if (s.find("YZ") !=std::string::npos){
        result.push_back (2);
        result.push_back (-2);
        return  result;
    }else if (s.find("XY") !=std::string::npos){
        result.push_back (2);
        result.push_back (-4);
        return  result;
    }else {
        cerr << "Orbital type = " << s <<  "not allowed";
        exit(1);
    }
}
// End of function sph_l_and_offset
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
// 
// //-----------------------------------------------------------------------------------------------------------------------------
// 
// // Function seekl
// int seekl(string s){
//     int l = 0;
//     transform(s.begin(), s.end(),s.begin(),::tolower);
//     for (int i=0 ; i < 9 ; i++){
//         if (s.compare(0,1,string(lvals[i])) == 0){ // only compares the first character of s (s may contains trailing blanks)
//             l = i;
//             break;
//         }
//     }
//     return l;
// }
// // End of function seekl
// 
// //-----------------------------------------------------------------------------------------------------------------------------
// 
// // Function seekiatdif
// int seekiatdif(int natdif, int n){
//     int j = 0;
//     for (int i = 0 ; i < natdif ; i++){
//         if (zn[n] == znatdif[i]){
//             j = i;
//             break;
//         }
//     }
//     return j;
// }
// // End of function seekiatdif
                        
