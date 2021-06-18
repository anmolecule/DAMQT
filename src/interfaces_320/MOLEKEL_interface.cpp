//  Copyright 2013-2017, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
//  Guillermo Ramirez, Anmol Kumar, Sachin D. Yeole, Shridhar R. Gadre
// 
//  This file is part of DAM2017.
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
//  along with DAM2017.  If not, see <http://www.gnu.org/licenses/>.
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
#include <string>
#include <stdlib.h>
#include <math.h>
using namespace std;

#define abs(n) ((n>0) ? n : -n)
#define minimum(n,m) ((n>m) ? m : n)

string limpia(string);
int min(int,int);

const int MXCEN = 10000;			// Maximum number of centers
const int MXSHELL = 20000;		// Maximum total number of contractions
const int MXSHELLAT = 50;		// Maximum number of contractions per atom
const int MXPRIMCENT = 200;   	// Maximum number of primitives per center
const double PI = 3.141592653589793;
int i,j,k,ii,jj,ki,kj,i1,j11,klin;
int len, lastocca, lastoccb, ncen, nbasis, kntshell = 0, izn;
double zn[MXCEN], x[MXCEN], y[MXCEN], z[MXCEN], primexp[MXPRIMCENT], cfcontr[MXPRIMCENT];
double *OMa, *OMb, *occa, *occb;
bool lzn = false, lcoord = false, lbasis = false, lcoefalfa = false, loccalfa = false, lcoefbeta = false, loccbeta = false;


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

int main(int argc,char *argv[])
{
	string MOLEKELfiles, inpath, newprojectname, outpath, projectname, s;
	if (argc<2) {
		cout << "Filename (with or without the extension  .mkl): ";
		cin >> projectname ;
		if (projectname.size() > 3 && projectname.substr(projectname.size()-4,4) == ".mkl")
			projectname = projectname.substr(0,projectname.size()-4);
		inpath = "." + slash;
		outpath = "." + slash;
	}
	else{
		projectname = argv[1]; // First argument (if available): project name (name of file .mkl without extension) 
		if (projectname.size() > 3 && projectname.substr(projectname.size()-4,4) == ".mkl")
			projectname = projectname.substr(0,projectname.size()-4);
		if (argc > 2) {
			inpath = argv[2]; // Second argument (if available): full path to MOLEKEL .mkl files
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
	
	s = outpath + newprojectname + "-MOLEKEL_interface.out";
	ofstream outimportfile(s.c_str(),ios::out);
	if (!outimportfile){
		cerr << "Unable to open file " << s << endl ;
		exit(1);
	}
	MOLEKELfiles = inpath + projectname;
	s = MOLEKELfiles + ".mkl";
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
	
	int *lvec = new int[MXSHELLAT];
	int *lvectot = new int[MXSHELL];
	int *iprimit = new int[MXSHELLAT+1];
	int ncontr, nprimit, ncontrtot = 0;
	while(getline(inputfile,s)){
		len = s.length();
/*    seek for atomic numbers and coordinates of the centers */
		if (!lcoord){
			if( !(s.find("COORD")==string::npos) ) {
				lcoord = true;
				ncen = 0;
				double xc=0., yc=0., zc=0., zntot=0.;
				while(getline(inputfile,s)&&!lzn){
					if ( !(s.find("END")==string::npos) ) 
						lzn = true;
					else {
						if ( !(s.find("$$")==string::npos) ) // If there are several sets of coordinates, restarts storage index to keep the last set at the end
							ncen = 0;
						len = s.length();
						char *tokenPrt, *ptr = new char [len+1], *newtoken;
						s.copy(ptr,len,0);
						ptr[len] = 0;
						tokenPrt = strtok_s(ptr," ",&newtoken);
						zn[ncen] = atof(tokenPrt);
						tokenPrt = strtok_s(NULL," ",&newtoken);
						x[ncen] = atof(tokenPrt)/0.529177249;
						tokenPrt = strtok_s(NULL," ",&newtoken);
						y[ncen] = atof(tokenPrt)/0.529177249;
						tokenPrt = strtok_s(NULL," ",&newtoken);
						z[ncen] = atof(tokenPrt)/0.529177249;
						zntot += zn[ncen];
						xc += zn[ncen] * x[ncen];
						yc += zn[ncen] * y[ncen];
						zc += zn[ncen] * z[ncen];
						ncen++;
						delete [] ptr;
					}
				}
//				redefines the coordinates with respect to the center of positive charges and writes the to file .ggbs
                xc = xc / zntot; yc = yc / zntot; zc = zc / zntot;
//				xc = 0.; yc = 0.; zc = 0.;  // Keeps original coordinates without translation of origin to center of positive charges
				ggbsfile << ncen << endl ;
				ggbsfile.setf(ios::scientific,ios::floatfield);
				for(i=0; i<ncen;i++){
					izn = int(round(zn[i]));
					ggbsfile << setprecision(10) << x[i]-xc << " " << y[i]-yc << " " << z[i]-zc << " " ;
					ggbsfile << izn << endl ;
				}
			}
		}
/*    seek for the basis set */
		if (!lbasis){
			if( !(s.find("$BASIS")==string::npos) ) {
				lbasis = true; nbasis = 0;
				int ncontr = 0, kprim = 0, kcntbas = 1; 
				while(getline(inputfile,s)&&(s.find("$END")==string::npos)){
					len = s.length();
					if (!(s.find("$$")==string::npos) ) {
						iprimit[ncontr] = kprim;
						ncontrtot += ncontr;
						ggbsfile.unsetf(ios::scientific);
						ggbsfile << endl << ncontr << endl;
						for(i=0; i<ncontr;i++){
							nprimit = iprimit[i+1] - iprimit[i];
							ggbsfile << nprimit << " " << lvec[i] << endl;
							ggbsfile.setf(ios::scientific,ios::floatfield);
							ggbsfile << setprecision(10);
							for (j=0; j<nprimit;j++){ 
								ggbsfile << primexp[j+iprimit[i]] << " ";
								if ((j+1)%5==0) ggbsfile << endl;
							}
							if (nprimit%5 != 0)	ggbsfile << endl;
							for (j=0; j<nprimit;j++){ 
								ggbsfile << cfcontr[j+iprimit[i]] << " ";
								if ((j+1)%5==0) ggbsfile << endl;
							}
							if (nprimit%5 != 0)	ggbsfile << endl;
							ggbsfile.unsetf(ios::scientific);
						}
						ncontr = 0; kprim = 1; kcntbas++;
					}
					else {
						if ( !(s.find_first_of("SPDFGHI")==string::npos) ){
							iprimit[ncontr] = kprim;
 							char *tokenPrt, *ptr = new char [len+1], *newtoken;
							s.copy(ptr,len,0);
							tokenPrt = strtok_s(ptr," ",&newtoken);
							tokenPrt = strtok_s(NULL," ",&newtoken);
							char car = *tokenPrt;
							switch(car) {
								case 'S': lvec[ncontr] = 0; break;
								case 'P': lvec[ncontr] = 1; break;
								case 'D': lvec[ncontr] = 2; break;
								case 'F': lvec[ncontr] = 3; break;
								case 'G': lvec[ncontr] = 4; break;
								case 'H': lvec[ncontr] = 5; break;
								case 'I': lvec[ncontr] = 6; break;
							}
							lvectot[kntshell++] = lvec[ncontr];
						     nbasis = nbasis + 2*lvec[ncontr++]+1;
							delete [] ptr;
						}
						else if (len != 0) {
							char *tokenPrt, *ptr = new char [len+1], *newtoken;
							s.copy(ptr,len,0);
							tokenPrt = strtok_s(ptr," ",&newtoken);
							primexp[kprim] = atof(tokenPrt);
							tokenPrt = strtok_s(NULL," ",&newtoken);	
							cfcontr[kprim++] = atof(tokenPrt);
							delete [] ptr;
						}	
					}
				}
				iprimit[ncontr] = kprim;
				ncontrtot += ncontr;
				ggbsfile.unsetf(ios::scientific);
				ggbsfile << endl << ncontr << endl;
				for(i=0; i<ncontr;i++){
					nprimit = iprimit[i+1] - iprimit[i];
					ggbsfile << nprimit << " " << lvec[i] << endl;
					ggbsfile.setf(ios::scientific,ios::floatfield);
					ggbsfile << setprecision(10);
					for (j=0; j<nprimit;j++){ 
						ggbsfile << primexp[j+iprimit[i]] << " ";
						if ((j+1)%5==0) ggbsfile << endl;
					}
					if (nprimit%5 != 0)	ggbsfile << endl;
					for (j=0; j<nprimit;j++){ 
						ggbsfile << cfcontr[j+iprimit[i]] << " ";
						if ((j+1)%5==0) ggbsfile << endl;
					}
					if (nprimit%5 != 0)	ggbsfile << endl;
					ggbsfile.unsetf(ios::scientific);
				}
				ggbsfile << endl;
			}
		}
/*    seek for the molecular orbitals */
		if (!lcoefalfa){
			if( !(s.find("$COEFF_ALPHA")==string::npos) ) {
				lcoefalfa = true;
				//	Allocates array for alpha molecular orbitals, OMa 
				OMa = new double[nbasis*nbasis];
				int ncic = nbasis / 5, lastcic = nbasis%5, iorb = 0;
				for (i=1; i<= ncic; i++){
					getline(inputfile,s); getline(inputfile,s);
					for (j=0; j<nbasis;j++){
						getline(inputfile,s);
						len = s.length();
						char *tokenPrt, *ptr = new char [len+1], *newtoken;
						s.copy(ptr,len,0);
						tokenPrt = strtok_s(ptr," ",&newtoken);	
						OMa[j*nbasis+iorb] = atof(tokenPrt);
						for (k=1;k<=4;k++){
							tokenPrt = strtok_s(NULL," ",&newtoken);	
							OMa[j*nbasis+iorb+k] = atof(tokenPrt);
						}
					}
					iorb = iorb+5;
				}
				if (lastcic != 0){
					getline(inputfile,s); getline(inputfile,s);
					for (j=0; j<nbasis;j++){
						getline(inputfile,s);
						len = s.length();
						char *tokenPrt, *ptr = new char [len+1], *newtoken;
						s.copy(ptr,len,0);
						tokenPrt = strtok_s(ptr," ",&newtoken);	
						OMa[j*nbasis+iorb] = atof(tokenPrt);
						for (k=1;k<=lastcic-1;k++){
							tokenPrt = strtok_s(NULL," ",&newtoken);	
							OMa[j*nbasis+iorb+k] = atof(tokenPrt);
						}
					}	
				}
// 				Reorders alpha MO from MOLEKEL to DAM 
				double px, py, pz, zlmv[25]; 
				for(j=0; j<nbasis;j++){
					int knt=0;
					for (i=0; knt < nbasis;i++){
						if(lvectot[i] == 1){
							py=OMa[(knt+1)*nbasis+j];
							pz=OMa[(knt+2)*nbasis+j];
							px=OMa[knt*nbasis+j];
							OMa[knt*nbasis+j]=py;
							OMa[(knt+1)*nbasis+j]=pz;
							OMa[(knt+2)*nbasis+j]=px;
						}
						else {
							for (k=0; k<= 2*lvectot[i]; k++){
								zlmv[k] = OMa[(knt+k)*nbasis+j];
							}
							ii = 0;
							for (k=2*lvectot[i]; k>=0; k -= 2){
								OMa[(knt+ii++)*nbasis+j] = zlmv[k];
							}
							for (k=1; k<2*lvectot[i]; k += 2){
								OMa[(knt+ii++)*nbasis+j] = zlmv[k];
							}
						}
						knt += 2*lvectot[i]+1;
					}
				}
			}
		}
/*    seek for the occupancies of alpha molecular orbitals */
		if (!loccalfa){
			if( !(s.find("$OCC_ALPHA")==string::npos) ) {
				loccalfa = true;
				occa = new double[nbasis];
				int ncic = nbasis / 5, lastcic = nbasis%5, iorb=0;
				lastocca = 0;
				for (i=1; i<= ncic; i++){
					getline(inputfile,s);
					
					len = s.length();
					char *tokenPrt, *ptr = new char [len+1], *newtoken;
					s.copy(ptr,len,0);
					tokenPrt = strtok_s(ptr," ",&newtoken);	
					occa[iorb] = atof(tokenPrt);
					if (occa[iorb] > 0) lastocca = iorb;
					iorb++;
					for (k=1;k<=4;k++){
						tokenPrt = strtok_s(NULL," ",&newtoken);	
						occa[iorb] = atof(tokenPrt);
						if (occa[iorb] > 0) lastocca = iorb;
						iorb++;
					}
				}
				if (lastcic != 0){
					getline(inputfile,s);
					len = s.length();
					char *tokenPrt, *ptr = new char [len+1], *newtoken;
					s.copy(ptr,len,0);
					tokenPrt = strtok_s(ptr," ",&newtoken);	
					occa[iorb] = atof(tokenPrt);
					if (occa[iorb] > 0) lastocca = iorb;
					iorb++;
					for (k=1;k<=lastcic-1;k++){
						tokenPrt = strtok_s(NULL," ",&newtoken);	
						occa[iorb] = atof(tokenPrt);
						if (occa[iorb] > 0) lastocca = iorb;
						iorb++;
					}	
				}
                string aorbfname(outpath + newprojectname + ".GAorba");
				ofstream aorbfile(aorbfname.c_str(),ios::out);
				if (!aorbfile) {
					cerr << "Unable to open file " << aorbfname << endl ;
					outimportfile << "Unable to open file " << aorbfname << endl ;
					exit(1);
				}
				aorbfile << nbasis << "  " << iorb  << "  " << iorb << endl;
				aorbfile << setprecision(15) ;
				aorbfile.setf(ios::scientific,ios::floatfield);
				for(i=0 ; i < iorb ; i++){
					aorbfile << endl ;
					int knt2 = 0;
					for (j=0 ; j < nbasis ; j++){
						aorbfile << OMa[j*nbasis+i] << " ";
						knt2++;
						if (knt2%5 == 0)
							aorbfile << endl;
					}
					if (knt2%5 != 0) aorbfile << endl;
				}
				aorbfile.close();
			}
		}
		if (!lcoefbeta){
			if( !(s.find("$COEFF_BETA")==string::npos) ) {
				lcoefbeta = true;
				//	Allocates array for beta molecular orbitals, OMb
				OMb = new double[nbasis*nbasis];
				int ncic = nbasis / 5, lastcic = nbasis%5, iorb = 0;
				for (i=1; i<= ncic; i++){
					getline(inputfile,s); getline(inputfile,s);
					for (j=0; j<nbasis;j++){
						getline(inputfile,s);
						len = s.length();
						char *tokenPrt, *ptr = new char [len+1], *newtoken;
						s.copy(ptr,len,0);
						tokenPrt = strtok_s(ptr," ",&newtoken);	
						OMb[j*nbasis+iorb] = atof(tokenPrt);
						for (k=1;k<=4;k++){
							tokenPrt = strtok_s(NULL," ",&newtoken);	
							OMb[j*nbasis+iorb+k] = atof(tokenPrt);
						}
					}
					iorb = iorb+5;
				}
				if (lastcic != 0){
					getline(inputfile,s); getline(inputfile,s);
					for (j=0; j<nbasis;j++){
						getline(inputfile,s);
						len = s.length();
						char *tokenPrt, *ptr = new char [len+1], *newtoken;
						s.copy(ptr,len,0);
						tokenPrt = strtok_s(ptr," ",&newtoken);	
						OMb[j*nbasis+iorb] = atof(tokenPrt);
						for (k=1;k<=lastcic-1;k++){
							tokenPrt = strtok_s(NULL," ",&newtoken);	
							OMb[j*nbasis+iorb+k] = atof(tokenPrt);
						}
					}	
				}
// 					Reorders beta MO from MOLEKEL to DAM 
				double px, py, pz, zlmv[25]; 
				for(j=0; j<nbasis;j++){
					int knt=0;
					for (i=0 ; knt < nbasis ; i++){
						if(lvectot[i] == 1){
							py=OMb[(knt+1)*nbasis+j];
							pz=OMb[(knt+2)*nbasis+j];
							px=OMb[knt*nbasis+j];
							OMb[knt*nbasis+j]=py;
							OMb[(knt+1)*nbasis+j]=pz;
							OMb[(knt+2)*nbasis+j]=px;
						}
						else {
							for (k=0; k<= 2*lvectot[i]; k++){
								zlmv[k] = OMb[(knt+k)*nbasis+j];
							}
							ii = 0;
							for (k=2*lvectot[i]; k>=0; k -= 2){
								OMb[(knt+ii++)*nbasis+j] = zlmv[k];
							}
							for (k=1; k<2*lvectot[i]; k += 2){
								OMb[(knt+ii++)*nbasis+j] = zlmv[k];
							}
						}
						knt += 2*lvectot[i]+1;
					}
				}
			}
		}
/*    seek for the occupancies of beta molecular orbitals */
		if (!loccbeta){
			if( !(s.find("$OCC_BETA")==string::npos) ) {
				loccbeta = true;
				occb = new double[nbasis];
				int ncic = nbasis / 5, lastcic = nbasis%5, iorb=0;
				lastoccb = 0;
				for (i=1; i<= ncic; i++){
					getline(inputfile,s);
					len = s.length();
					char *tokenPrt, *ptr = new char [len+1], *newtoken;
					s.copy(ptr,len,0);
					tokenPrt = strtok_s(ptr," ",&newtoken);	
					occb[iorb] = atof(tokenPrt);
					if (occb[iorb] > 0) lastoccb = iorb;
					iorb++;
					for (k=1;k<=4;k++){
						tokenPrt = strtok_s(NULL," ",&newtoken);	
						occb[iorb] = atof(tokenPrt);
						if (occb[iorb] > 0) lastoccb = iorb;
						iorb++;
					}
				}
				if (lastcic != 0){
					getline(inputfile,s);
					len = s.length();
					char *tokenPrt, *ptr = new char [len+1], *newtoken;
					s.copy(ptr,len,0);
					tokenPrt = strtok_s(ptr," ",&newtoken);	
					occb[iorb] = atof(tokenPrt);
					if (occb[iorb] > 0) lastoccb = iorb;
					iorb++;
					for (k=1;k<=lastcic-1;k++){
						tokenPrt = strtok_s(NULL," ",&newtoken);	
						occb[iorb] = atof(tokenPrt);
						if (occb[iorb] > 0) lastoccb = iorb;
						iorb++;
					}	
				}
                string borbfname(outpath + newprojectname + ".GAorbb");
				ofstream borbfile(borbfname.c_str(),ios::out);
				if (!borbfile) {
					cerr << "Unable to open file " << borbfname << endl ;
					outimportfile << "Unable to open file " << borbfname << endl ;
					exit(1);
				}
				borbfile << nbasis << "  " << iorb << "  " << iorb << endl;
				borbfile << setprecision(15) ;
				borbfile.setf(ios::scientific,ios::floatfield);
				for(i=0 ; i < iorb ; i++){
					int knt2 = 0;
					borbfile << endl ;
					for (j=0 ; j < nbasis ; j++){
						borbfile << OMb[j*nbasis+i] << " ";
						knt2++;
						if (knt2%5 == 0)
							borbfile << endl;
					}
					if (knt2%5 != 0) borbfile << endl;
				}
				borbfile.close();
			}
		}
	}
/*	Builds the density matrix from MO and occupancies */
	double *dmat = new double[nbasis*nbasis];
	double sum;
	for (i=0; i<nbasis; i++){
		for (j=0; j<nbasis; j++){
			sum = 0;
			for (k=0; k<= lastocca; k++){
				sum += OMa[i*nbasis+k] * OMa[j*nbasis+k] * occa[k];
			}
			if (loccbeta){
				for (k=0; k<= lastoccb; k++){
					sum += OMb[i*nbasis+k] * OMb[j*nbasis+k] * occb[k];
				}
			}
			dmat[i*nbasis+j] = sum;
		}
	}

/*   Writes the density matrix to file denfile (.den)   */

	denfile << " " << nbasis << endl;

	denfile << setprecision(10) ;
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

/*   Prints out some statistics  */
	outimportfile << "MOLEKEL_interface output" << endl;
	outimportfile << "=====================" << endl<< endl;
	outimportfile << "inpath = " << inpath << endl;
	outimportfile << "outpath = " << outpath << endl << endl;
    outimportfile << "PROJECT NAME: " <<  newprojectname << endl << endl;
	outimportfile << "Number of centers: " << ncen << endl;
	outimportfile << "Number of contractions: " << ncontrtot << endl;
	outimportfile << "Number of basis functions: " << nbasis << endl;
	if (!loccbeta){
		outimportfile << "Last occupied MO: " << lastocca+1 << endl ;
		outimportfile << "Occupancies: " ;}
	else{
		outimportfile << "Last occupied Alpha MO: " << lastocca+1 << endl ;
		outimportfile << "Occupancies of  Alpha MO: " ;
	}
	for (i=0; i<= lastocca; i++){
		outimportfile << occa[i] << " ";
		if (((i+1)%50)==0) outimportfile << endl ;
	}
	if (loccbeta){
		outimportfile << "Last occupied Beta MO: " << lastoccb+1 << endl ;
		outimportfile << "Occupancies of  Beta MO: " ;
		for (i=0; i<= lastoccb; i++){
			outimportfile << occb[i] << " ";
			if (((i+1)%50)==0) outimportfile << endl ;
		}
	}
	outimportfile << endl;
}
