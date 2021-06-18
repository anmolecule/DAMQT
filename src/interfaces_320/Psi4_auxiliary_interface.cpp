//  Copyright 2013-2020, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
//  Guillermo Ramirez, Anmol Kumar, Sachin D. Yeole, Shridhar R. Gadre
// 
//  This file is part of DAM320.
// 
//  DAM2017 is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
// 
//  DAM2017 is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with DAM320.  If not, see <http://www.gnu.org/licenses/>.
//
//------------------------------------------------------------------------
//
//  Version of April 2020
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
#include <sstream>
#include <vector>
#include <regex>
using namespace std;

const double angsTobohr = 1.8897259886L;
const double pi = 3.14159265358979L;
const double alpha = 1.e-5;

struct atombasis{
    string symbol;
    double zn;
    double atomcharge;
    vector<int> lv;
    vector<double> expv;
    vector<double> rnorm;       // Radial normalization of the expansion function
    vector<double> rnormdummy;   // Radial normalization of the dummy function
    vector<double> charge;      // Charge or multipole of the expansion function
};

#define abs(n) ((n>0) ? n : -n)
#define minimum(n,m) ((n>m) ? m : n)

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
int seekzn(string s);
vector<string> split (const string &s, char delim);

#if defined(_WIN32) || defined(WIN32) 
    static const std::string slash = "\\";    
#else
    static const std::string slash = "/";
#endif

int main(int argc,char *argv[]) {
    int izn, knt, lmax, natomdif, ncen, nfun;
    string inpath, psi4file, newprojectname, outpath, projectname, s, s2;
    double *x, *y, *z, *zn;
    double fact[30], facts[30];
    double chargealpha, rnormalpha;
    
    clock_t startTime = clock();
    
    fact[0] = 1.;
    facts[0] = 0.5 * sqrt(pi);
    for (int i = 1 ; i < 30 ; i++){
        fact[i] = i*fact[i-1];
        facts[i] = 0.5 * (2.*i+1.) * facts[i-1];
    }
    rnormalpha = sqrt(2.*pow((2.*alpha),1.5) / facts[0]);
    chargealpha = 1. / (2.*pow((2.*alpha),1.5) / facts[0]);
    
//  Default values for inpath and outpath
    inpath = "." + slash;
    outpath = "." + slash;
    
    if (argc<2) {
        cout << "Please, supply a project name: ";
        cin >> projectname ;
        
    } else{
        projectname = argv[1]; // First argument (if available): project name 
        if (argc > 2) {
            inpath = argv[2]; // Second argument (if available): full path to data files
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
    
    s = outpath + newprojectname + "-Psi4_auxiliary_interface.out";
    ofstream outimportfile(s.c_str(),ios::out);
    if (!outimportfile){
        cerr << "In Psi4_auxiliary_interface: unable to open file " << s << endl ;
        exit(1);
    }
    psi4file = inpath + projectname;
    
    s = psi4file + ".psixyz.xyz";    
    ifstream xyzfile;
    xyzfile.open(s.c_str(),ios::in);
    if (!xyzfile) {
        cerr << "In Psi4_auxiliary_interface: unable to open file " << s <<  endl ;
        outimportfile << "In Psi4_auxiliary_interface: unable to open file " << s <<  endl ;
        exit(1);
    }
    
    s = outpath + newprojectname + ".ggbs";
    ofstream ggbsfile(s.c_str(),ios::out);
    if (!ggbsfile) {
        cerr << "In Psi4_auxiliary_interface: unable to open file " << s << endl ;
        outimportfile << "In Psi4_auxiliary_interface: unable to open file " << s << endl ;
        exit(1);
    }
    
//     Reads geometry from .psixyz.xyz file and writes to .ggbs. Transforms distances to Bohr
    getline(xyzfile,s);
    ncen = atoi(s.c_str());
    getline(xyzfile,s); // Reads comment line to dumy
    x = new double [ncen];
    y = new double [ncen];
    z = new double [ncen];
    zn = new double [ncen];
    double xc=0.L, yc=0.L, zc=0.L, zntot=0.;
    
    knt = 0;
    string *auxstring = new string [4];
    natomdif = 0;
    while(knt < ncen && getline(xyzfile,s)){
        vector<string> v = split (s, ' ');
        int knt2 = 0;
        for (auto i : v){
            if (knt2 >= 4)
                break;
            if (i == "") continue;
            auxstring[knt2] = i; 
            knt2++;
        }
        if (knt2 < 4){
            cerr << "In Psi4_auxiliary_interface: error reading xyz file, last line read: " << s << endl ;
            exit(1);
        }
        zn[knt] = (double) seekzn(auxstring[0]);
        x[knt] = atof(auxstring[1].c_str()) * angsTobohr;
        xc += x[knt] * zn[knt];
        y[knt] = atof(auxstring[2].c_str()) * angsTobohr;
        yc += y[knt] * zn[knt];
        z[knt] = atof(auxstring[3].c_str()) * angsTobohr;
        zc += z[knt] * zn[knt];
        bool existatom = false;
        for (int j = 0 ; j < knt ; j++){
            if (zn[knt]==zn[j]){
                existatom = true;
                break;
            }
        }
        if (!existatom){ 
            natomdif++;
        }
        zntot += zn[knt];
        knt++;
    }
    
    if (!(knt == ncen)){
        outimportfile << "Error reading xyz file wrong number of centers in xyz file " << endl;
        exit(1);
    }
    xyzfile.close();
   
    ggbsfile << ncen << endl ;
    ggbsfile.setf(ios::scientific,ios::floatfield);
    xc /= zntot ; yc /= zntot ; zc /= zntot ;
    for(int i = 0; i < ncen; i++){
        izn = int(zn[i]+0.5L);
        ggbsfile << setprecision(15) << x[i]-xc << " " << y[i]-yc << " " << z[i]-zc << " " ;
        ggbsfile << izn << endl ;
    }
    
//     Reads basis set from psigbs file and writes it to .ggbs including a dummy S function for each atom.    
    
    s = psi4file + ".psigbs";
    ifstream basisfile;
    basisfile.open(s.c_str(),ios::in);
    if (!basisfile) {
        cerr << "In Psi4_auxiliary_interface: unable to open file " << s <<  endl ;
        outimportfile << "In Psi4_auxiliary_interface: unable to open file " << s <<  endl ;
        exit(1);
    }
    
    while(getline(basisfile,s)){
        if (!(s.find(string("****"))==string::npos))
            break;
    }
    lmax = 0;
    vector<atombasis> atom;
    while(getline(basisfile,s)){
        atombasis atombas;
        vector<string> v = split (s, ' ');
        for (auto i : v){
            if (i == "") continue;
            atombas.symbol = i;
            break;
        }
        atombas.zn = (double) seekzn(atombas.symbol);
        while(getline(basisfile,s) && (s.find(string("****"))==string::npos)){
            v = split (s, ' ');
            string ftype = v[0];
            getline(basisfile,s);
            v = split (s, ' ');
            double charge, dummyexpon, expon, rnorm, rnormdummy;
            for (auto i : v){
                if (i == "") continue;
                replace( i.begin(), i.end(), 'D', 'e');
                expon = atof(i.c_str());
                dummyexpon = expon - alpha;
                break;
            }
            if (ftype == "S"){
//                 cout << "funcion de tipo S" << endl;
                atombas.lv.push_back(0);
                atombas.expv.push_back(dummyexpon);
                rnormdummy = sqrt(2.*pow((2.*dummyexpon),1.5) / facts[0]);      // Radial normalization of the dummy function
                rnorm = sqrt(2.*pow((2.*expon),1.5) / facts[0]);          // Radial normalization of the expansion function
                charge = 2.*pi / (pow((expon),1.5) / facts[0]);                   // Charge or multipole of the expansion function
                atombas.rnormdummy.push_back(rnormdummy);
                atombas.rnorm.push_back(rnorm);
                atombas.charge.push_back(charge);
            }
            else if (ftype == "SP"){
//                 cout << "funcion de tipo SP" << endl;
                atombas.lv.push_back(0);
                atombas.expv.push_back(dummyexpon);
                rnormdummy = sqrt(2.*pow((2.*dummyexpon),1.5) / facts[0]);      // Radial normalization of the dummy function
                rnorm = sqrt(2.*pow((2.*expon),1.5) / facts[0]);          // Radial normalization of the expansion function
                charge = 2.*pi / (pow((expon),1.5) / facts[0]);                    // Charge or multipole of the expansion function
                atombas.rnormdummy.push_back(rnormdummy);
                atombas.rnorm.push_back(rnorm);
                atombas.charge.push_back(charge);
                atombas.lv.push_back(1);
                atombas.expv.push_back(dummyexpon);
                rnormdummy = sqrt(2.*pow((2.*dummyexpon),2.5) / facts[1]);      // Radial normalization of the dummy function
                rnorm = sqrt(2.*pow((2.*expon),2.5) / facts[1]);          // Radial normalization of the expansion function
                charge = 2.*pi / (pow((expon),2.5) / facts[1]);                    // Charge or multipole of the expansion function
                atombas.rnormdummy.push_back(rnormdummy);
                atombas.rnorm.push_back(rnorm);
                atombas.charge.push_back(charge);
                lmax = max(lmax,1);
            }
            else if (ftype == "P"){
//                 cout << "funcion de tipo P" << endl;
                atombas.lv.push_back(1);
                atombas.expv.push_back(dummyexpon);
                rnormdummy = sqrt(2.*pow((2.*dummyexpon),2.5) / facts[1]);      // Radial normalization of the dummy function
                rnorm = sqrt(2.*pow((2.*expon),2.5) / facts[1]);          // Radial normalization of the expansion function
                charge = 2.*pi / (pow((expon),2.5) / facts[1]);                    // Charge or multipole of the expansion function
                atombas.rnormdummy.push_back(rnormdummy);
                atombas.rnorm.push_back(rnorm);
                atombas.charge.push_back(charge);
                lmax = max(lmax,1);
            }
            else if (ftype == "D"){
//                 cout << "funcion de tipo D" << endl;
                atombas.lv.push_back(2);
                atombas.expv.push_back(dummyexpon);
                rnormdummy = sqrt(2.*pow((2.*dummyexpon),3.5) / facts[2]);      // Radial normalization of the dummy function 
                rnorm = sqrt(2.*pow((2.*expon),3.5) / facts[2]);          // Radial normalization of the expansion function
                charge = 2.*pi / (pow((expon),3.5) / facts[2]);                    // Charge or multipole of the expansion function
                atombas.rnormdummy.push_back(rnormdummy);
                atombas.rnorm.push_back(rnorm);
                atombas.charge.push_back(charge);
                lmax = max(lmax,2);
            }
            else if (ftype == "F"){
//                 cout << "funcion de tipo F" << endl;
                atombas.lv.push_back(3);
                atombas.expv.push_back(dummyexpon);
                rnormdummy = sqrt(2.*pow((2.*dummyexpon),4.5) / facts[3]);      // Radial normalization of the dummy function 
                rnorm = sqrt(2.*pow((2.*expon),4.5) / facts[3]);          // Radial normalization of the expansion function
                charge = 2.*pi / (pow((expon),4.5) / facts[3]);                    // Charge or multipole of the expansion function
                atombas.rnormdummy.push_back(rnormdummy);
                atombas.rnorm.push_back(rnorm);
                atombas.charge.push_back(charge);
                lmax = max(lmax,3);
            }
            else{
                cerr << "Function type " << ftype << " not allowable" << endl;
                exit(1);
            }
        }     
        atom.push_back(atombas);
        
    }
//  Angular normalization factors:
//     angnorm(l*(l+1)/2+m+1) = sqrt( (2*l+1) * fact(l-m) / (2 * pi * (1 + delta(m,0)) * fact(l+m)) )
    vector<double> angnorm;     
    for (int l = 0 ; l < lmax+1 ; l++){
        for (int m = -l ; m < 0 ; m++){
            double anorm = sqrt((2.*l+1.) * fact[l-abs(m)] / (2.*pi * fact[l+abs(m)]));
            angnorm.push_back(anorm);
        }
        angnorm.push_back(sqrt((2.*l+1.) / (4.*pi)));
        for (int m = 1 ; m <= l ; m++){
            double anorm = sqrt((2.*l+1.) * fact[l-abs(m)] / (2.*pi * fact[l+abs(m)]));
            angnorm.push_back(anorm);
        }
    }
    
    vector<atombasis> totalbasis;
    for (int j = 0 ; j < ncen ; j++){
        bool found = false;
        for (int k = 0 ; k < atom.size() ; k++){
            if (abs(zn[j]-atom[k].zn) < 1.e-3){
                found = true;
                ggbsfile << endl << atom[k].lv.size()+1 << endl ;
//                 The first function is an S dummy function with exponent alpha
                ggbsfile << "1 0" << endl ;
                ggbsfile << alpha << endl ;
                ggbsfile << "1." << endl ;
                atombasis auxbasis;
                auxbasis.symbol = atom[k].symbol;
                auxbasis.zn = atom[k].zn;
                auxbasis.lv.push_back(0);
                auxbasis.expv.push_back(alpha);
                auxbasis.rnormdummy.push_back(rnormalpha);
                auxbasis.rnorm.push_back(1. / rnormalpha);
                auxbasis.charge.push_back(chargealpha);
                
                nfun++;
                for (int n = 0 ; n < atom[k].lv.size() ; n++){
                    ggbsfile << "1 " << atom[k].lv[n] << endl ;
                    ggbsfile << atom[k].expv[n] << endl ;
                    ggbsfile << "1." << endl ;
                    nfun += 2*atom[k].lv[n]+1;
                    auxbasis.lv.push_back(atom[k].lv[n]);
                    auxbasis.expv.push_back(atom[k].expv[n]);
                    auxbasis.rnormdummy.push_back(atom[k].rnormdummy[n]);
                    auxbasis.rnorm.push_back(atom[k].rnorm[n]);
                    auxbasis.charge.push_back(atom[k].charge[n]);
                }
                totalbasis.push_back(auxbasis);
                break;
            }
        }
        if (!found){
            cerr << "cannot find a basis set for atom with Znuc = " << zn[j] << endl;
            exit(1);
        }
    }
    ggbsfile.close();
    
    s = psi4file + ".psiauxden";  
    ifstream denpsi4file;
    denpsi4file.open(s.c_str(),ios::in);
    if (!denpsi4file) {
        cerr << "In Psi4_auxiliary_interface: unable to open file " << s <<  endl ;
        outimportfile << "In Psi4_auxiliary_interface: unable to open file " << s <<  endl ;
        exit(1);
    }
    
    s = outpath + newprojectname + ".den";
    ofstream denfile(s.c_str(),ios::out);
    if (!denfile) {
        cerr << "In Psi4_auxiliary_interface: unable to open file " << s << endl ;
        outimportfile << "In Psi4_auxiliary_interface: unable to open file " << s << endl ;
        exit(1);
    }
    
    getline(denpsi4file,s); 
    int kntatom = -1;
    int kntrow = 0;
    vector<double> dmat;
    vector<int> nfbas;
    double vaux[2*lmax+1];
    double normaux[2*lmax+1];
    double totalcharge = 0.;
    while(getline(denpsi4file,s)){
        if (!(s.find(string("number of functions"))==string::npos)){
            vector<string> v = split (s, ':');
            nfbas.push_back(atoi(v[1].c_str())+1);
            kntatom++;
            int kntbas = 0;
            int jk = 0;
            for (int j = 0 ; j < kntatom ; j++){
                for (int k = 0 ; k < nfbas[j] ; k++){
                    dmat.push_back(0.);
                    jk++;
                }
            }
            dmat.push_back(0.);
            kntbas++;
            kntrow++;
            double chargeatom = 0.;
            while(getline(denpsi4file,s) && (s.find(string("Atom"))==string::npos)){
                v = split (s, ' ');
                int l = atoi(v[0].c_str());
                double daux = atof(v[1].c_str());
                double nraux = atof(v[2].c_str());
                vaux[l] = daux;
                normaux[l] = nraux;
                for (int j = 0 ; j < 2*l ; j++){
                    getline(denpsi4file,s);
                    v = split (s, ' ');
                    daux = atof(v[1].c_str());
                    nraux = atof(v[2].c_str());
                    if (j%2 == 0){
                        vaux[l+(j/2)+1]=daux;
                        normaux[l+(j/2)+1]=nraux;
                    }
                    else{
                        vaux[l-(j/2)-1]=daux;
                        normaux[l-(j/2)-1]=nraux;
                    }
                }
                for (int j = 0 ; j < 2*l+1 ; j++){
                    for (int k = 0 ; k < jk ; k++){
                        dmat.push_back(0.);
                    }
                    dmat.push_back( 0.5 * vaux[j] * totalbasis[kntatom].rnorm[kntbas]
                                    / ( totalbasis[kntatom].rnormdummy[kntbas] * rnormalpha * angnorm[0]) );
                    if (l == 0){ 
                        chargeatom += vaux[j] * totalbasis[kntatom].charge[kntbas] * totalbasis[kntatom].rnorm[kntbas] * angnorm[l*l+j];
                        totalcharge += vaux[j] * totalbasis[kntatom].charge[kntbas] * totalbasis[kntatom].rnorm[kntbas] * angnorm[l*l+j];
                    }
                    for (int k = 0 ; k < kntrow-jk ; k++){
                        dmat.push_back(0.);
                    }
                    kntrow++;
                }
                kntbas++;
            }
            totalbasis[kntatom].atomcharge = -chargeatom;
        }
    } 
    int nden = dmat.size();
    denfile << nfun << " ";
    denfile << setprecision(15);
    for (int i = 0 ; i < nden ; i+=8){
        for (int j = 0 ; j < 8 ; j++){
            if (i+j >= nden)
                break;
            denfile << dmat[i+j] << " ";
        }
        denfile << endl;
    }
    denfile.close();
    clock_t endTime = clock();
    clock_t clockTicksTaken = endTime - startTime;
    double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;
    outimportfile << "Psi4_auxiliary_interface output" << endl;
    outimportfile << "===============================" << endl<< endl;
    outimportfile << "inpath = " << inpath << endl;
    outimportfile << "outpath = " << outpath << endl << endl;

    outimportfile << "PROJECT NAME: " <<  newprojectname << endl << endl;
    outimportfile << "Number of centers: " << ncen << endl;
    outimportfile << "Number of different atoms = " << natomdif << endl << endl;

    outimportfile << "Basis set for expansion of atomic densities (l, exponent) " << endl;
    for (int i = 0 ; i < atom.size() ; i++){
        outimportfile << endl << "Atom: " << atom[i].symbol << endl;
        for (int j = 0 ; j < atom[i].lv.size() ; j++){
            outimportfile << atom[i].lv[j] << " " << atom[i].expv[j] << endl;
        }
    }
    outimportfile << endl << "Exponent for dummy (S) basis functions = " << alpha << endl;
    outimportfile << endl << "Atomic electron charges " << endl;
    outimportfile << "-----------------------" << endl;
    for (int i = 0 ; i < totalbasis.size() ; i++){
        outimportfile << "Atom no " << i << "(" << totalbasis[i].symbol << ")"<< " charge = "
                      << totalbasis[i].atomcharge << endl;
    }
    outimportfile << endl << "Total molecular electron charge = " << -totalcharge << endl;
    outimportfile << endl << "Time (in secs) = " << timeInSeconds << endl;
}

    
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
