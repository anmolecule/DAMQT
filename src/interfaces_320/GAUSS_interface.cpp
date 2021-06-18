//  Copyright 2013-2019, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
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
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <ctime>
#include <string>
using namespace std;

#define abs(n) ((n>0) ? n : -n)
#define minimum(n,m) ((n>m) ? m : n)

string limpia(string);
int min(int,int);


#if __cplusplus <= 199711L
    #define nullpointer NULL
#else
//  C++11 compliant compiler
    #define nullpointer nullptr
#endif


const int MXCEN = 10000;
const int MXSHELL = 20000;
const int MXBASIS = 45000;
const int MXPRIMCENT = 20;   // Maximum number of primitives per center
const int MXPRIMTOT = 30000;  // Maximum total number of primitives
const double PI = 3.141592653589793;
int i,j,k,ii,jj,ki,kj,kl,i1,j11,klin;
int lineknt = 0;
int konta, kontadns, lmax = 0, len, spos;
int ncen, nshell, nshelltot, nprimtot, ndens, nbasis, nindep, nbasaux, nsp=0, numMO;
int nalfa, nbeta, ncfalfa, ncfbeta;



bool lallocateden = false, lncen = false, lzn = false, lcoord = false, lshell = false, lprimit = false
   , lcontr = false, lshellat = false, lprimexp = false, lcfcontr = false, lcfP = true
   , lendl = false , ldens = false, lbasis = false, lindep = false, lRHF = true, lROHF = false
   , lDTotSCF = false , lDSpSCF = false , lDTotCC = false , lDSpCC = false, lDTotCI = false , lDSpCI = false, lDTotMP2 = false , lDSpMP2 = false
   , lalfae = false, lbetae = false, lcfalfa = false, lcfbeta = false, lzdo = false;
double facts[50];

#if ( defined(_WIN32) || defined(WIN32) ) && !defined(__MINGW32__)
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
        return nullpointer;
    last = str + strlen(str);
    if( (*save = res = strtok(str, delim)) )
    {
        *save += strlen(res);
        if( *save < last )
            (*save)++;
        else
            *save = nullpointer;
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




int main(int argc,char *argv[])
{
    double *dens=nullpointer;
    double *dmat=nullpointer;
    double *primexp=nullpointer, *cfcontr=nullpointer, *cfcontrP=nullpointer, *zn=nullpointer;
//    double cfcontrP[5000];
    int *lvec=nullpointer, *lvectot=nullpointer, *ncontr=nullpointer, *ncontrtot=nullpointer, *nprimit=nullpointer;
//    int ncontr[MXCEN], ncontrtot[MXCEN];
	int mden = 0;
    string fchkfile, outpath, newprojectname, inpath, s, ss;	// 0=checks file, 1=SCF electron density matrix, 2=SCF spin  density matrix,
											// 3=CI electron density matrix, 4= CI spin  density matrix
											// 5=MP2 electron density matrix, 6= MP2 spin  density matrix
                                            // 7=CC electron density matrix, 8= CC spin  density matrix
	clock_t startTime = clock();
	if (argc<2) {
		cout << "Filename (with or without the extension  .fchk): ";
		cin >> ss;
		mden=1;
		inpath = "." + slash;
		outpath = "." + slash;
	}
	else{ 
		if (argc==2){
		   if (strncmp(argv[1],"?",1)==0){
				cout << "Usage: the following commands are valid:" << endl ;
				cout << "   GAUSS_interface.exe " << endl ;
				cout << "   GAUSS_interface.exe projectname (mden) " << endl ;
				cout << "   GAUSS_interface.exe projectname mden inputdir (outputdir)" << endl ;
                cout << "   GAUSS_interface.exe fchkfile mden inputdir outputdir projectname" << endl ;
				cout << "Available density options (mden): " << endl ;
				cout << "   GAUSS_interface  fimename.fchk   : Total SCF density matrix " << endl ;
				cout << "   GAUSS_interface  fimename.fchk 1 : Total SCF density matrix " << endl ;
				cout << "   GAUSS_interface  fimename.fchk 2 : SCF spin density density matrix " << endl ;
				cout << "   GAUSS_interface  fimename.fchk 3 : Total CI density matrix " << endl ;
				cout << "   GAUSS_interface  fimename.fchk 4 : CI spin density matrix " << endl ;
                cout << "   GAUSS_interface  fimename.fchk 5 : Total CC density matrix " << endl ;
                cout << "   GAUSS_interface  fimename.fchk 6 : CC spin density matrix " << endl ;
				exit(0);
			}
			else{
				ss=argv[1];
				mden=1;
			}
			inpath = "." + slash;
			outpath = "." + slash;
		} 
		else {
			ss=argv[1];
			mden=atoi(argv[2]);
			if (argc >= 4) {
				inpath = argv[3];
				if(inpath.find_last_of(slash) != (inpath.length()-1) ) inpath += slash;
			}
			else inpath = "." + slash;
			if (argc >= 5) {
				outpath = argv[4];
				if(outpath.find_last_of(slash) != (outpath.length()-1) ) outpath += slash;
				if (argc == 6){
                    newprojectname = argv[5];
				}
			}
			else outpath = "." + slash;
		}
	}
	string ext(ss.end()-5,ss.end());
	if(ext==".fchk"){
		fchkfile.assign(ss.begin(),ss.end()-5);
	}
	else{
		fchkfile.assign(ss.begin(),ss.end());
	}
	if (argc < 6)
        newprojectname = fchkfile;

//cout << "inpath = " << inpath << endl;
	string s1 = inpath + fchkfile + ".fchk";
	ifstream inputfile(s1.c_str(),ios::in);
	if (!inputfile) {
		cerr << "Unable to open file " << s1 << endl ;
		exit(1);
	}
    string s4(outpath + newprojectname + "-GAUSS_interface.out");
	ofstream outinterface(s4.c_str(),ios::out);
	if (!outinterface){
		cerr << "Unable to open file " << s4 << endl ;
		exit(1);
	}
    if (mden==0) { //searches CI, MP2 and CC density and Spin matrices
		outinterface << "existSCF" << endl;
		while(getline(inputfile,s)){
		   lineknt++;
			if( !(s.find("SCF")==string::npos) && !(s.find("Density")==string::npos) &&
				  !(s.find("Spin")==string::npos) ) {
                //prints to an ancillary file that SCF spin matrix exists
				outinterface << "existSpSCF" << endl;
			}	
			if( !(s.find("CI")==string::npos) && !(s.find("Density")==string::npos) &&
				  !(s.find("Total")==string::npos) ) {
				//prints to an ancillary file that CI matrix exists
				outinterface << "existCI" << endl;
			}
			if( !(s.find("CI")==string::npos) && !(s.find("Density")==string::npos) &&
				  !(s.find("Spin")==string::npos) ) {
                //prints to an ancillary file that CI spin matrix exists
				outinterface << "existSpCI" << endl;
			}
			if( !(s.find("MP2")==string::npos) && !(s.find("Density")==string::npos) &&
				  !(s.find("Total")==string::npos) ) {
                //prints to an ancillary file that MP2 matrix exists
				outinterface << "existMP2" << endl;
			}
			if( !(s.find("MP2")==string::npos) && !(s.find("Density")==string::npos) &&
				  !(s.find("Spin")==string::npos) ) {
                //prints to an ancillary file that MP2 spin matrix exists
				outinterface << "existSpMP2" << endl;
			}
            if( !(s.find("CC")==string::npos) && !(s.find("Density")==string::npos) &&
                  !(s.find("Total")==string::npos) ) {
                //prints to an ancillary file that CC matrix exists
                outinterface << "existCC" << endl;
            }
            if( !(s.find("CC")==string::npos) && !(s.find("Density")==string::npos) &&
                  !(s.find("Spin")==string::npos) ) {
                //prints to an ancillary file that CC spin matrix exists
                outinterface << "existSpCC" << endl;
            }
		}
		inputfile.clear(); // eof bit must be removed, it was set to read up to the end of file 
		inputfile.seekg(0);
		exit(0);
	}else{
//outinterface << "mden = " << mden << endl;
		switch (mden) {
			case 1:	//use SCF electron density
				lDTotSCF=true;
				break ;
			case 2:  //use SCF spin density
				lDSpSCF=true; 
				break ;
			case 3:  //use CI electron density
				lDTotCI=true;
				break ;
			case 4:  //use CI spin density
				lDSpCI=true;
				break ;
			case 5:  //use MP2 electron density
				lDTotMP2=true;
				break ;
			case 6:  //use MP2 spin density
				lDSpMP2=true;
				break ;
            case 7:  //use CC electron density
                lDTotCC=true;
                break ;
            case 8:  //use CC spin density
                lDSpCC=true;
                break ;
		}
	}
	
    string ggbsfname(outpath + newprojectname + ".ggbs");
	ofstream ggbsfile(ggbsfname.c_str(),ios::out);
	if (!ggbsfile) {
		cerr << "Unable to open file " << ggbsfname << endl ;
		outinterface << "Unable to open file " << ggbsfname << endl ;
		exit(1);
	}
    string denfname(outpath + newprojectname + ".den");
	ofstream denfile(denfname.c_str(),ios::out);
	if (!denfile) {
		cerr << "Unable to open file " << denfname << endl ;
		outinterface << "Unable to open file " << denfname << endl ;
		exit(1);
	}
    string aorbfname(outpath + newprojectname + ".GAorba");
	ofstream aorbfile(aorbfname.c_str(),ios::out);
	if (!aorbfile) {
		cerr << "Unable to open file " << aorbfname << endl ;
		outinterface << "Unable to open file " << aorbfname << endl ;
		exit(1);
	}

	while(getline(inputfile,s)){
		lineknt++;
		len = s.length();

/*    Tests whether is a ZDO calculation */
        if (!lzdo && !(s.find("ZDO")==string::npos)) {
//cout << "aqui lzdo: " <<  std::flush << endl;
            lzdo = true;
            string ss(outpath + "zdo");
            ofstream outzdo(ss.c_str(),ios::out);
            outzdo.close();
            string ss2(outpath + "valence");
            ofstream outzdo2(ss2.c_str(),ios::out);
            outzdo2.close();
//cout << "sale de lzdo: " << std::flush << endl;
        }

/*    searches and processes the line with the number of centers  */

        if (!lncen){
            if( !(s.find("Number of atoms")==string::npos) ) {
//cout << "aqui lncen: " <<  std::flush << endl;
                char *tokenPrt=nullpointer, *ptr = new char [len+1], *tokenPrt2=nullpointer, *newtoken;
                lncen = true;
                s.copy(ptr,len,0);
                ptr[len] = 0;
                tokenPrt = strtok_s(ptr," ",&newtoken);
                while (tokenPrt != nullpointer){
                    tokenPrt2 = tokenPrt;
                    tokenPrt = strtok_s(nullpointer," ",&newtoken);
//cout << "en lncen tokenPrt = " << tokenPrt2 << std::flush << endl;
                }
                ncen = atoi(tokenPrt2);
                ggbsfile << tokenPrt2 << endl;
                delete [] ptr;
//cout << "sale de lncen nalfa: " << ncen << std::flush << endl;
            }
        }           /*  End of if(!lncen)  */
        if(ncen > MXCEN){
            outinterface << endl << "Number of centers: " << ncen << " greater than allowed: " \
                << MXCEN << " increase the value of parameter MXCEN and recompile. Stop" << endl;
            exit(1);
        }


/*    Tests whether is a ROHF calculation */
        if (!lROHF && !(s.find("ROHF")==string::npos) && (s.find("IROHF")==string::npos)) {lRHF = false; lROHF = true;}

/*    searches the number of alpha electrons */
        if (!lalfae){
            if( !(s.find("Number of alpha electrons")==string::npos) ) {
//cout << "aqui lalfae: " <<  std::flush << endl;
                lalfae = true;
                char *tokenPrt=nullpointer, *ptr = new char [len+1], *tokenPrt2=nullpointer, *newtoken;
                lncen = true;
                s.copy(ptr,len,0);
                ptr[len] = 0;
                tokenPrt = strtok_s(ptr," ",&newtoken);
                while (tokenPrt != nullpointer){
                    tokenPrt2 = tokenPrt;
                    tokenPrt = strtok_s(nullpointer," ",&newtoken);
                }
                nalfa = atoi(tokenPrt2);
                delete [] ptr;
//cout << "sale de lalfae nalfa: " << nalfa << std::flush << endl;
            }
        }           /*  End of if(!lalfae)  */


/*    searches the number of beta electrons */

        if (!lbetae){
            if( !(s.find("Number of beta electrons")==string::npos) ) {
//cout << "aqui lbetae: " <<  std::flush << endl;
                lbetae = true;
                char *tokenPrt=nullpointer, *ptr = new char [len+1], *tokenPrt2=nullpointer, *newtoken;
                lncen = true;
                s.copy(ptr,len,0);
                ptr[len] = 0;
                tokenPrt = strtok_s(ptr," ",&newtoken);
                while (tokenPrt != nullpointer){
                    tokenPrt2 = tokenPrt;
                    tokenPrt = strtok_s(nullpointer," ",&newtoken);
                }
                nbeta = atoi(tokenPrt2);
                delete [] ptr;
                if (nalfa != nbeta) lRHF = false ;
//cout << "sale de lbetae nbeta: " << nbeta << std::flush << endl;
            }
        }           /*  End of if(!lbetae)  */


/*    searches and processes the line with the number of basis functions   */

        if (!lbasis){
            if( !(s.find("Number of basis functions")==string::npos) ) {
//cout << "aqui lbasis: " <<  std::flush << endl;
                char *tokenPrt=nullpointer, *ptr = new char [len+1], *tokenPrt2=nullpointer, *newtoken;
                lbasis = true;
                s.copy(ptr,len,0);
                ptr[len] = 0;
                tokenPrt = strtok_s(ptr," ",&newtoken);
                while (tokenPrt != nullpointer){
                    tokenPrt2 = tokenPrt;
                    tokenPrt = strtok_s(nullpointer," ",&newtoken);
                }
                nbasis = atoi(tokenPrt2);
                nindep = nbasis;
                if(nbasis > MXBASIS){
                    outinterface << endl << "Number of basis functions: " << nbasis << " greater than allowed: " \
                         << MXBASIS << " increase the value of parameter MXBASIS and recompile. Stop" << endl;
                    exit(1);
                }
                delete [] ptr;
//cout << "sale de lbasis nbasis: " << nbasis << std::flush << endl;
            }
        }           /*  End of if(!lbasis)  */

/*    searches and processes the line with the number of independent functions  (O.M.)  */

        if (!lindep){
            if( !(s.find("Number of independent functions")==string::npos) ) {
//cout << "aqui lindep: " <<  std::flush << endl;
                char *tokenPrt=nullpointer, *ptr = new char [len+1], *tokenPrt2=nullpointer, *newtoken;
                lindep = true;
                s.copy(ptr,len,0);
                ptr[len] = 0;
                tokenPrt = strtok_s(ptr," ",&newtoken);
                while (tokenPrt != nullpointer){
                    tokenPrt2 = tokenPrt;
                    tokenPrt = strtok_s(nullpointer," ",&newtoken);
                }
                nindep = atoi(tokenPrt2);
                delete [] ptr;
//cout << "sale de lindep nindep: " << nindep << std::flush << endl;
            }
        }           /*  End of if(!lindep)  */


/*    searches and processes the line with the number of atomic numbers   */

        if (!lzn){
            if( !(s.find("Atomic numbers")==string::npos) ) {
//cout << "aqui lzn ncen : " << ncen <<  std::flush << endl;
                konta = 0;
                zn = new double[ncen];
                while(getline(inputfile,s)&&!lzn){
                    lineknt++;
                    len = s.length();
                    char *tokenPrt3=nullpointer, *ptr3 = new char [len+1], *newtoken;
                    s.copy(ptr3,len,0);
                    ptr3[len] = 0;
                    tokenPrt3 = strtok_s(ptr3," ",&newtoken);
                    while (tokenPrt3 != nullpointer){
                        konta++;
                        if (konta>=ncen) lzn=true;
                        zn[konta-1] = atof(tokenPrt3) ;
                        tokenPrt3 = strtok_s(nullpointer," ",&newtoken);
                    }
                    delete [] ptr3;
                }
//cout << "sale de lzn: " <<  std::flush << endl;
            }
        }           /*  End of if(!lzn)  */


/*    searches and processes the line with the coordinates of centers  */

        if (!lcoord){
            if( !(s.find("Current cartesian")==string::npos) ) {
//cout << "aqui lcoord: " <<  std::flush << endl;
                konta = 0;
                int kontc = 1;
                int kontb = 1;
                double (*xyz) = new double[3*ncen];
                int kntxyz = 0;
                while(getline(inputfile,s)&&!lcoord){
                    lineknt++;
                    len = s.length();
                    char *tokenPrt3=nullpointer, *ptr3 = new char [len+1], *newtoken;
                    s.copy(ptr3,len,0);
                    ptr3[len] = 0;
                    tokenPrt3 = strtok_s(ptr3," \r\n",&newtoken);
                    while (tokenPrt3 != nullpointer){
                        konta++;
                        if (konta>=(3*ncen)) lcoord=true;
                        xyz[kntxyz++] = atof(tokenPrt3);
//						ggbsfile << tokenPrt3 << " ";
//						if(kontb++==3){
//							kontb = 1;
//							ggbsfile << " " << zn[kontc-1] << endl ;
//							kontc = kontc + 1;
//						}
                        tokenPrt3 = strtok_s(nullpointer," \r\n",&newtoken);
                    }
                    delete [] ptr3;
                }
                double xC = 0., yC = 0., zC = 0., qtot = 0.;
                for (int i = 0; i < kntxyz ; i += 3){
                    xC += xyz[i] * zn[i/3];
                    yC += xyz[i+1] * zn[i/3];
                    zC += xyz[i+2] * zn[i/3];
                    qtot += zn[i/3];
                }
                xC /= qtot;
                yC /= qtot;
                zC /= qtot;
                ggbsfile.setf(ios::scientific,ios::floatfield);
                ggbsfile << setprecision(10) ;
                for (int i = 0; i < kntxyz ; i += 3){
                    ggbsfile << xyz[i] - xC << " " << xyz[i+1] - yC << " "
                             << xyz[i+2] - zC << " " << zn[i/3] << endl;
                }
                ggbsfile << endl;
                ggbsfile.unsetf(ios::scientific);
                delete [] xyz;
//cout << "sale de lcoord: " << std::flush << endl;
            }
        }           /*  End of if(!lcoord)  */


/*    searches and processes the line with the "l" quantum numbers of shells  */

        if (!lshell){
            if( !(s.find("Shell types")==string::npos) ) {
//cout << "aqui lshell: " << std::flush << endl;
                nbasaux = 0;
                len = s.length();
                char *tokenPrt=nullpointer, *ptr = new char [len+1], *tokenPrt2=nullpointer, *newtoken;
                lncen = true;
                s.copy(ptr,len,0);
                ptr[len] = 0;
                tokenPrt = strtok_s(ptr," ",&newtoken);
                while (tokenPrt != nullpointer){
                    tokenPrt2 = tokenPrt;
                    tokenPrt = strtok_s(nullpointer," ",&newtoken);
                }
                nshell = atoi(tokenPrt2);
                nshelltot = nshell;
                if(nshelltot > MXSHELL){
                    outinterface << endl << "Number of shells: " << nshelltot << " greater than allowed: " \
                        << MXSHELL << " increase the value of parameter MXSHELL and recompile. Stop" << endl;
                    exit(1);
                }
                delete [] ptr;
                konta = 0;
                lvec = new int [nshell];
                lvectot = new int [2*nshell];
                int *pntr = lvectot;
//                pntr = lvectot;
                while(getline(inputfile,s)&&!lshell){
                    lineknt++;
                    len = s.length();
                    char *tokenPrt3=nullpointer, *ptr3 = new char [len+1], *newtoken;
                    s.copy(ptr3,len,0);
                    ptr3[len] = 0;
                    tokenPrt3 = strtok_s(ptr3," ",&newtoken);
                    while (tokenPrt3 != nullpointer){
                        konta++;
                        if (konta>=nshell) lshell=true;
                        lvec[konta-1]=atoi(tokenPrt3);
                        if (abs(lvec[konta-1]) > lmax) lmax = abs(lvec[konta-1]);
                        nbasaux = nbasaux + (2*abs(lvec[konta-1])+1);
                        if(lvec[konta-1]==-1){
                            nsp++; lcfP = false; nshelltot++; *pntr++ = 0;
                        }
                        *pntr++ = abs(lvec[konta-1]);
                        tokenPrt3 = strtok_s(nullpointer," ",&newtoken);
                    }
                    delete [] ptr3;
                }
                if(nshelltot > MXSHELL){
                    outinterface << endl << "Number of shells: " << nshelltot << " greater than allowed: " \
                        << MXSHELL << " increase the value of parameter MXSHELL and recompile. Stop" << endl;
                    exit(1);
                }
//cout << "sale de lbasis nshell: " << nshell << std::flush << endl;
            }
        }           /*  End of if(!lshell)  */


/*    searches and processes the line with the number of primitives for each contraction  */

        if (!lprimit){
            if( !(s.find("Number of primitives")==string::npos) ) {
//cout << "aqui lprimit nshell :" << nshell << std::flush << endl;
                konta = 0;
                nprimtot = 0;
                nprimit = new int [nshell];
                while(getline(inputfile,s)&&!lprimit){
                    lineknt++;
                    len = s.length();
                    char *tokenPrt3=nullpointer, *ptr3 = new char [len+1], *newtoken;
                    s.copy(ptr3,len,0);
                    ptr3[len] = 0;
                    tokenPrt3 = strtok_s(ptr3," ",&newtoken);
                    while (tokenPrt3 != nullpointer){
                        konta++;
                        if (konta>=nshell) lprimit=true;
                        nprimit[konta-1] = atoi(tokenPrt3);
                        nprimtot = nprimtot + nprimit[konta-1];
                        tokenPrt3 = strtok_s(nullpointer," ",&newtoken);
                        if(nprimit[konta-1]>=MXPRIMCENT){
                            outinterface << "Number of primitives exceeded in center " << konta
                                << ": " << nprimit[konta-1] << " exceeds the maximum allowed: " << MXPRIMCENT << endl
                                << "Increase the parameter MXPRIMCENT in file GAUSS.cpp and recompile " << endl;
                            exit(1);
                        }
                    }
                    delete [] ptr3;
                }
                if(nprimtot>=MXPRIMTOT){
                    outinterface << "Total number of primitives: " << nprimtot ;
                    outinterface << "   exceeds the maximum allowed: " << MXPRIMTOT << endl
                        << "Increase the parameter MXPRIMTOT in file GAUSS.cpp and recompile " << endl;
                    exit(1);
                }
//cout << "sale de lprimit: " <<  std::flush << endl;
            }
        }           /*  End of if(!lprimit)  */


/*    searches and processes the lines which associate each contraction with a center  */

        if (!lshellat){
            if( !(s.find("Shell to atom")==string::npos) ) {
//cout << "aqui lshellat ncen :" << ncen << std::flush << endl;
                konta = 0;
                int iatom = 1, katom;
                ncontr = new int[ncen];
                ncontrtot = new int[ncen];
                for(i=0 ; i< ncen ; i++){ncontr[i] = 0; ncontrtot[i] = 0;}
                while(getline(inputfile,s)&&!lshellat){
                    lineknt++;
                    len = s.length();
                    char *tokenPrt3=nullpointer, *ptr3 = new char [len+1], *newtoken;
                    s.copy(ptr3,len,0);
                    ptr3[len] = 0;
                    tokenPrt3 = strtok_s(ptr3," ",&newtoken);
                    while (tokenPrt3 != nullpointer){
                        konta++;
                        if (konta>=nshell) lshellat=true;
                        katom=atoi(tokenPrt3);
                        if(katom != iatom){
                            if(katom == iatom+1){
                                iatom=iatom+1;}
                            else{ outinterface << endl << "Error: katom-iatom = " << katom-iatom
                                << " Stop "; exit(1);}
                        }
                        ncontr[iatom-1]=ncontr[iatom-1]+1;
                        ncontrtot[iatom-1]=ncontrtot[iatom-1]+1;
                        if(lvec[konta-1]==-1) ncontrtot[iatom-1]=ncontrtot[iatom-1]+1;
                        tokenPrt3 = strtok_s(nullpointer," ",&newtoken);
                    }
                    delete [] ptr3;
                }
//cout << "sale de lshellat: " <<  std::flush << endl;
            }
        }           /*  End of if(!lshellat)  */


/*    searches and processes the line with the primitive exponents  */

        if (!lprimexp){
            if( !(s.find("Primitive exponents")==string::npos) ) {
//cout << "aqui lprimexp nprimtot: " << nprimtot << std::flush << endl;
                konta = 0;
                primexp = new double[nprimtot];
                while(getline(inputfile,s)&&!lprimexp){
                    lineknt++;
                    len = s.length();
                    char *tokenPrt3=nullpointer, *ptr3 = new char [len+1], *newtoken;
                    s.copy(ptr3,len,0);
                    ptr3[len] = 0;
                    tokenPrt3 = strtok_s(ptr3," ",&newtoken);
                    while (tokenPrt3 != nullpointer){
                        konta++;
                        if (konta>=nprimtot) lprimexp=true;
                        primexp[konta-1] = atof(tokenPrt3);
                        tokenPrt3 = strtok_s(nullpointer," ",&newtoken);
                    }
                    delete [] ptr3;
                }
//cout << "sale de lprimexp: " <<  std::flush << endl;
            }
        }           /*  End of if(!lprimexp)  */


/*    searches and processes the line with the contraction coefficients  */

        if (!lcfcontr){
            if( !(s.find("Contraction coefficients")==string::npos) ) {
//cout << "aqui lcfcontr nprimtot: " << nprimtot <<  std::flush << endl;
                konta = 0;
                cfcontr = new double[nprimtot];
                while(getline(inputfile,s)&&!lcfcontr){
                    lineknt++;
                    len = s.length();
                    char *tokenPrt3=nullpointer, *ptr3 = new char [len+1], *newtoken;
                    s.copy(ptr3,len,0);
                    ptr3[len] = 0;
                    tokenPrt3 = strtok_s(ptr3," ",&newtoken);
                    while (tokenPrt3 != nullpointer){
                        konta++;
                        if (konta>=nprimtot) lcfcontr=true;
                        double cfaux;
                        cfaux = atof(tokenPrt3);
                        cfcontr[konta-1] = cfaux;
                        tokenPrt3 = strtok_s(nullpointer," ",&newtoken);
                    }
                    delete [] ptr3;
                }
//cout << "sale de lcfcontr: " <<  std::flush << endl;
            }
        }           /*  End of if(!lcfcontr)  */


/*    special case of basis sets with SP shells  */

        if (!lcfP){
            if( !(s.find("P(S=P) Contraction coefficients")==string::npos) ) {
//cout << "aqui lcfP nprimtot : " << nprimtot << std::flush << endl;
                konta = 0;
                cfcontrP = new double[5000];
                while(getline(inputfile,s)&&!lcfP){
                    lineknt++;
                    len = s.length();
                    char *tokenPrt3=nullpointer, *ptr3 = new char [len+1], *newtoken;
                    s.copy(ptr3,len,0);
                    ptr3[len] = 0;
                    tokenPrt3 = strtok_s(ptr3," ",&newtoken);
                    while (tokenPrt3 != nullpointer){
                        konta++;
                        if (konta>=nprimtot) lcfP=true;
                        double cfaux;
                        cfaux = atof(tokenPrt3);
                        cfcontrP[konta-1] = cfaux;
                        tokenPrt3 = strtok_s(nullpointer," ",&newtoken);
                    }
                    delete [] ptr3;
                }
//cout << "sale de lcfP: " <<  std::flush << endl;
            }
        }           /*  End of if(!lcfP)  */

/*    If nalfa != nbeta and it is a ROHF calculation the density matrix is built from the molecular orbitals.
To enable this, stores the MO in dmat temporarily  */
        if (!lcfalfa){
            if( !(s.find("Alpha MO coefficients")==string::npos) && !(s.find("N=")==string::npos) ) {
//cout << "aqui lcfalfa: " <<  std::flush << endl;
/*    Reads the coefficients of the alpha orbitals and loads the elements of the density matrix to cover the case of ROHF calculation */
                len = s.length();
                char *tokenPrt=nullpointer, *ptr = new char [len+1], *tokenPrt2=nullpointer, *newtoken;
                s.copy(ptr,len,0);
                ptr[len] = 0;
                tokenPrt = strtok_s(ptr," ",&newtoken);
                while (tokenPrt != nullpointer){
                    tokenPrt2 = tokenPrt;
                    tokenPrt = strtok_s(nullpointer," ",&newtoken);
                }
                ncfalfa = atoi(tokenPrt2);
                delete [] ptr;
                if (ncfalfa != nbasis*nindep){
                    outinterface << "Error en Alpha MO coefficients " << endl;
                    outinterface << "ncfalfa = " << ncfalfa << " != nbasis*nindep = " << nbasis*nindep << endl ;
                    exit(1);
                }
                double *aorb = new double[nbasis * nindep];
                double *vaux = new double[2*lmax+1];
//	Allocates arrays dens and dmat to be loaded into the heap instead of into the stack. This enables to allocate big arrays
//	without having problems of stack overflow. Currently, stack size is 8Mb, wheras heap size is much higher.
                if (lallocateden){
                    delete [] dens;
                    delete [] dmat;
                }
                if (nbasis <= 0){
                    outinterface << "Error allocating density matrix when reading Alpha MO coefficients ";
                    outinterface << "Number of basis functions = " << nbasis << endl ;
                    exit(1);
                }
                dens = new double[nbasis*(nbasis+1)/2];
                dmat = new double[nbasis*nbasis];
                lallocateden = true;
                kontadns = 0;
                i = 0; j = 0;
                while(getline(inputfile,s) && !lcfalfa){
                    lineknt++;
                    len = s.length();
                    char *tokenPrt3=nullpointer, *ptr3 = new char [len+1], *newtoken;
                    s.copy(ptr3,len,0);
                    ptr3[len] = 0;
                    tokenPrt3 = strtok_s(ptr3," ",&newtoken);
                    while (tokenPrt3 != nullpointer){
                        aorb[j+i*nbasis]=atof(tokenPrt3);
                        dmat[(j++)*nbasis+i]=atof(tokenPrt3);
                        if(j == nbasis){
                            j = 0; i++;
                            if(i == nindep){ lcfalfa = true;}
                        }
                        tokenPrt3 = strtok_s(nullpointer," ",&newtoken);
                    }
                    delete [] ptr3;
                }
                if(i != nindep ){
                    outinterface << "Error in Alpha MO coefficients ";
                    outinterface << "Number of MO read = " << i << endl ;
                    exit(1);
                }
                else{
/*   	Sorts the GAUSSIAN-ordered alpha molecular orbitals: M = (0,+1,-1,+2,-2,...)  to canonical order: (...,-2,-1,0,1,2,...) */
                    for(i=0 ; i < nindep ; i++){
                        kj = 0;
                        for(j=0 ; j<nshelltot ; j++){
                            switch(lvectot[j]){
                                case 0:
                                    vaux[0] = aorb[kj+i*nbasis];
                                    break;
                                case 1:
                                    vaux[0] = aorb[kj+1+i*nbasis] ;
                                    vaux[1] = aorb[kj+2+i*nbasis] ;
                                    vaux[2] = aorb[kj+i*nbasis]   ;
                                    break;
                                default:
                                    j11 = 0;
                                    for(jj = 2*lvectot[j] ; jj >= 0 ; jj -= 2, j11++){
                                        vaux[j11] = aorb[kj+jj+i*nbasis];
                                    }
                                    for(jj = 1 ; jj <= 2*lvectot[j]-1 ; jj += 2, j11++){
                                        vaux[j11] = aorb[kj+jj+i*nbasis];
                                    }
                                    break;
                            }
                            for (int row = 0; row < 2*lvectot[j]+1 ; row++){
                                aorb[kj+row+i*nbasis] = vaux[row];
                            }
                            kj += 2*lvectot[j]+1;
                        }
                    }
                    aorbfile << nbasis << "  " << nalfa  << "  " << nindep << endl;
                    aorbfile << setprecision(15) ;
                    aorbfile.setf(ios::scientific,ios::floatfield);
                    int knt = 0;
                    for(i=0 ; i < nindep ; i++){
                        aorbfile << endl ;
                        int knt2 = 0;
                        for (j=0 ; j < nbasis ; j++){
                            aorbfile << aorb[knt++] << " ";
                            knt2++;
                            if (knt2%5 == 0)
                                aorbfile << endl;
                        }
                        if (knt2%5 != 0) aorbfile << endl;
                    }
                    aorbfile.close();
                    if (lDTotSCF && !lRHF){
                        double sum; int r, s, knt=0 , kntprt=0;
                        for(r = 0; r < nbasis ; r++){
                            for(s = 0; s <= r ; s++, knt++){
                                for(i = 0, sum = 0. ; i < nalfa ; i++){
                                    sum += dmat[r*nbasis+i] * dmat[s*nbasis+i];
                                }
                                dens[knt] = sum;
                                kntprt++;
                            }
                        }
                    }
                }
                delete [] aorb;
                delete [] vaux;
//cout << "sale de lcfalfa: " <<  std::flush << endl;
            }
        }    /* End of if(lcfalfa)  */

/*   If there are not beta orbital coefficients and it is a ROHF calculation:   */
/*   it will build the density matrix from alpha molecular orbitals 			*/
/*   If there are beta orbitals, reads density matrix from file .fchk			*/

        if( !lcfbeta){
            if(!(s.find("Total")==string::npos) &&
                !(s.find("SCF")==string::npos) &&
                !(s.find("Density")==string::npos) &&
                !(s.find("N=")==string::npos) ) {

                if (!lRHF){
                    lROHF = true ;
                    ldens = true ;
                }
                lcfbeta = true ;
            }
            else{
                if (!(s.find("Beta MO coefficients")==string::npos) && !(s.find("N=")==string::npos) ) {
/*    Reads the coefficients of the beta orbitals  */
                    len = s.length();
                    char *tokenPrt=nullpointer, *ptr = new char [len+1], *tokenPrt2=nullpointer, *newtoken;
                    s.copy(ptr,len,0);
                    ptr[len] = 0;
                    tokenPrt = strtok_s(ptr," ",&newtoken);
                    while (tokenPrt != nullpointer){
                        tokenPrt2 = tokenPrt;
                        tokenPrt = strtok_s(nullpointer," ",&newtoken);
                    }
                    ncfbeta = atoi(tokenPrt2);
                    delete [] ptr;
                    if (ncfbeta != nbasis*nindep){
                        outinterface << "Error en Beta MO coefficients " << endl;
                        outinterface << "ncfbeta = " << ncfbeta << " != nbasis*nindep = " << nbasis*nindep << endl ;
                        exit(1);
                    }
                    double *borb = new double[nbasis * nindep];
                    double *vaux = new double[2*lmax+1];

                    kontadns = 0;
                    i = 0; j = 0;
                    while(getline(inputfile,s) && !lcfbeta){
                        lineknt++;
                        len = s.length();
                        char *tokenPrt3=nullpointer, *ptr3 = new char [len+1], *newtoken;
                        s.copy(ptr3,len,0);
                        ptr3[len] = 0;
                        tokenPrt3 = strtok_s(ptr3," ",&newtoken);
                        while (tokenPrt3 != nullpointer){
                            borb[j+i*nbasis]=atof(tokenPrt3);
                            j++;
                            if(j == nbasis){
                                j = 0; i++;
                                if(i == nindep){ lcfbeta = true;};
                            }
                            tokenPrt3 = strtok_s(nullpointer," ",&newtoken);
                        }
                        delete [] ptr3;
                    }
                    if(i != nindep ){
                        outinterface << "Error in Beta MO coefficients ";
                        outinterface << "Number of MO read = " << i << endl ;
                        exit(1);
                    }
                    else{
/*   	Sorts the GAUSSIAN-ordered beta molecular orbitals: M = (0,+1,-1,+2,-2,...)  to canonical order: (...,-2,-1,0,1,2,...) */
                        for(i=0 ; i < nindep ; i++){
                            kj = 0;
                            for(j=0 ; j<nshelltot ; j++){
                                switch(lvectot[j]){
                                    case 0:
                                        vaux[0] = borb[kj+i*nbasis];
                                        break;
                                    case 1:
                                        vaux[0] = borb[kj+1+i*nbasis] ;
                                        vaux[1] = borb[kj+2+i*nbasis] ;
                                        vaux[2] = borb[kj+i*nbasis]   ;
                                        break;
                                    default:
                                        j11 = 0;
                                        for(jj = 2*lvectot[j] ; jj >= 0 ; jj -= 2, j11++){
                                            vaux[j11] = borb[kj+jj+i*nbasis];
                                        }
                                        for(jj = 1 ; jj <= 2*lvectot[j]-1 ; jj += 2, j11++){
                                            vaux[j11] = borb[kj+jj+i*nbasis];
                                        }
                                        break;
                                }
                                for (int row = 0; row < 2*lvectot[j]+1 ; row++){
                                    borb[kj+row+i*nbasis] = vaux[row];
                                }
                                kj += 2*lvectot[j]+1;
                            }
                        }
                        string borbfname(outpath + newprojectname + ".GAorbb");
                        ofstream borbfile(borbfname.c_str(),ios::out);
                        if (!borbfile) {
                            cerr << "Unable to open file " << borbfname << endl ;
                            outinterface << "Unable to open file " << borbfname << endl ;
                            exit(1);
                        }
                        borbfile << nbasis << "  " << nbeta << "  " << nindep << endl;
                        borbfile << setprecision(15) ;
                        borbfile.setf(ios::scientific,ios::floatfield);
                        int knt = 0;
                        for(i=0 ; i < nindep ; i++){
                            int knt2 = 0;
                            borbfile << endl ;
                            for (j=0 ; j < nbasis ; j++){
                                borbfile << borb[knt++] << " ";
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
        }


/*    If nalfa == nbeta Total RHF density matrix   */
        if(lDTotSCF &&  !(s.find("Total")==string::npos) &&
            !(s.find("SCF")==string::npos) &&
            !(s.find("Density")==string::npos) &&
            !(s.find("N=")==string::npos) ) {
//outinterface << "aqui lDTotSCF: " <<  std::flush << endl;
            len = s.length();
            char *tokenPrt=nullpointer, *ptr = new char [len+1], *tokenPrt2=nullpointer, *newtoken;
            lncen = true;
            s.copy(ptr,len,0);
            ptr[len] = 0;
            tokenPrt = strtok_s(ptr," ",&newtoken);
            while (tokenPrt != nullpointer){
                tokenPrt2 = tokenPrt;
                tokenPrt = strtok_s(nullpointer," ",&newtoken);
            }
            ndens = atoi(tokenPrt2);

            if (lallocateden){
                delete [] dens;
                delete [] dmat;
            }
            if (ndens != nbasis*(nbasis+1)/2){
                outinterface << "Error allocating RHF density matrix ";
                outinterface << "Wrong number of elements, nbasis = " << nbasis << " ndens = " << ndens << endl ;
                exit(1);
            }
            dens = new double[ndens];
            dmat = new double[nbasis*nbasis];
            lallocateden = true;

            delete [] ptr;
            kontadns = 0;
            while(getline(inputfile,s)&&!ldens){
                lineknt++;
                len = s.length();
                char *tokenPrt3=nullpointer, *ptr3 = new char [len+1], *newtoken;
                s.copy(ptr3,len,0);
                ptr3[len] = 0;
                tokenPrt3 = strtok_s(ptr3," ",&newtoken);
                while (tokenPrt3 != nullpointer){
                    kontadns++;
                    if (kontadns>=ndens) ldens=true;
                    dens[kontadns-1]=atof(tokenPrt3);
                    tokenPrt3 = strtok_s(nullpointer," ",&newtoken);
                }
                delete [] ptr3;
            }
//cout << "sale de lDTotSCF: " <<  std::flush << endl;
        }

/*    Searches SCF spin density matrix      */
        if( lDSpSCF  && !(s.find("Spin")==string::npos)
            && !(s.find("SCF")==string::npos)
            && !(s.find("Density")==string::npos)
            &&  !(s.find("N=")==string::npos) ) {
            ldens = false ;
            len = s.length();
            char *tokenPrt=nullpointer, *ptr = new char [len+1], *tokenPrt2=nullpointer, *newtoken;
            lncen = true;
            s.copy(ptr,len,0);
            ptr[len] = 0;
            tokenPrt = strtok_s(ptr," ",&newtoken);
            while (tokenPrt != nullpointer){
                tokenPrt2 = tokenPrt;
                tokenPrt = strtok_s(nullpointer," ",&newtoken);
            }
            ndens = atoi(tokenPrt2);

            if (lallocateden){
                delete [] dens;
                delete [] dmat;
            }
            if (ndens != nbasis*(nbasis+1)/2){
                outinterface << "Error allocating RHF density matrix ";
                outinterface << "Wrong number of elements, nbasis = " << nbasis << " ndens = " << ndens << endl ;
                exit(1);
            }
            dens = new double[ndens];
            dmat = new double[nbasis*nbasis];
            lallocateden = true;

            delete [] ptr;
            kontadns = 0;
            while(getline(inputfile,s)&&!ldens){
                lineknt++;
                len = s.length();
                char *tokenPrt3=nullpointer, *ptr3 = new char [len+1], *newtoken;
                s.copy(ptr3,len,0);
                ptr3[len] = 0;
                tokenPrt3 = strtok_s(ptr3," ",&newtoken);
                while (tokenPrt3 != nullpointer){
                    kontadns++;
                    if (kontadns>=ndens) ldens=true;
                    dens[kontadns-1]=atof(tokenPrt3);
                    tokenPrt3 = strtok_s(nullpointer," ",&newtoken);
                }
                delete [] ptr3;
            }
        }
/*    Searches CI density matrix      */
        if(lDTotCI &&  !(s.find("Total")==string::npos)
            &&  !(s.find("CI")==string::npos)
            && !(s.find("Density")==string::npos)
            && !(s.find("N=")==string::npos) ) {
            ldens = false ;
            len = s.length();
            char *tokenPrt=nullpointer, *ptr = new char [len+1], *tokenPrt2=nullpointer, *newtoken;
            lncen = true;
            s.copy(ptr,len,0);
            ptr[len] = 0;
            tokenPrt = strtok_s(ptr," ",&newtoken);
            while (tokenPrt != nullpointer){
                tokenPrt2 = tokenPrt;
                tokenPrt = strtok_s(nullpointer," ",&newtoken);
            }
            ndens = atoi(tokenPrt2);

            if (lallocateden){
                delete [] dens;
                delete [] dmat;
            }
            if (ndens != nbasis*(nbasis+1)/2){
                outinterface << "Error allocating CI density matrix ";
                outinterface << "Wrong number of elements, nbasis = " << nbasis << " ndens = " << ndens << endl ;
                exit(1);
            }
            dens = new double[ndens];
            dmat = new double[nbasis*nbasis];
            lallocateden = true;

            delete [] ptr;
            kontadns = 0;
            while(getline(inputfile,s)&&!ldens){
                lineknt++;
                len = s.length();
                char *tokenPrt3=nullpointer, *ptr3 = new char [len+1], *newtoken;
                s.copy(ptr3,len,0);
                ptr3[len] = 0;
                tokenPrt3 = strtok_s(ptr3," ",&newtoken);
                while (tokenPrt3 != nullpointer){
                    kontadns++;
                    if (kontadns>=ndens) ldens=true;
                    dens[kontadns-1]=atof(tokenPrt3);
                    tokenPrt3 = strtok_s(nullpointer," ",&newtoken);
                }
                delete [] ptr3;
            }
        }
/*    Searches CI spin density matrix        */
        if(lDSpCI && !(s.find("Spin")==string::npos)
            &&  !(s.find("CI")==string::npos)
            && !(s.find("Density")==string::npos)
            && !(s.find("N=")==string::npos) ) {
            ldens = false ;
            len = s.length();
            char *tokenPrt=nullpointer, *ptr = new char [len+1], *tokenPrt2=nullpointer, *newtoken;
            lncen = true;
            s.copy(ptr,len,0);
            ptr[len] = 0;
            tokenPrt = strtok_s(ptr," ",&newtoken);
            while (tokenPrt != nullpointer){
                tokenPrt2 = tokenPrt;
                tokenPrt = strtok_s(nullpointer," ",&newtoken);
            }
            ndens = atoi(tokenPrt2);

            if (lallocateden){
                delete [] dens;
                delete [] dmat;
            }
            if (ndens != nbasis*(nbasis+1)/2){
                outinterface << "Error allocating CI density matrix ";
                outinterface << "Wrong number of elements, nbasis = " << nbasis << " ndens = " << ndens << endl ;
                exit(1);
            }
            dens = new double[ndens];
            dmat = new double[nbasis*nbasis];
            lallocateden = true;

            delete [] ptr;
            kontadns = 0;
            while(getline(inputfile,s)&&!ldens){
                lineknt++;
                len = s.length();
                char *tokenPrt3=nullpointer, *ptr3 = new char [len+1], *newtoken;
                s.copy(ptr3,len,0);
                ptr3[len] = 0;
                tokenPrt3 = strtok_s(ptr3," ",&newtoken);
                while (tokenPrt3 != nullpointer){
                    kontadns++;
                    if (kontadns>=ndens) ldens=true;
                    dens[kontadns-1]=atof(tokenPrt3);
                    tokenPrt3 = strtok_s(nullpointer," ",&newtoken);
                }
                delete [] ptr3;
            }
        }
/*    Searches MP2 density matrix      */
        if(lDTotMP2 &&  !(s.find("Total")==string::npos)
            &&  !(s.find("MP2")==string::npos)
            && !(s.find("Density")==string::npos)
            && !(s.find("N=")==string::npos) ) {
            ldens = false ;
            len = s.length();
            char *tokenPrt=nullpointer, *ptr = new char [len+1], *tokenPrt2=nullpointer, *newtoken;
            lncen = true;
            s.copy(ptr,len,0);
            ptr[len] = 0;
            tokenPrt = strtok_s(ptr," ",&newtoken);
            while (tokenPrt != nullpointer){
                tokenPrt2 = tokenPrt;
                tokenPrt = strtok_s(nullpointer," ",&newtoken);
            }
            ndens = atoi(tokenPrt2);

            if (lallocateden){
                delete [] dens;
                delete [] dmat;
            }
            if (ndens != nbasis*(nbasis+1)/2){
                outinterface << "Error allocating MP2 density matrix ";
                outinterface << "Wrong number of elements, nbasis = " << nbasis << " ndens = " << ndens << endl ;
                exit(1);
            }
            dens = new double[ndens];
            dmat = new double[nbasis*nbasis];
            lallocateden = true;

            delete [] ptr;
            kontadns = 0;
            while(getline(inputfile,s)&&!ldens){
                lineknt++;
                len = s.length();
                char *tokenPrt3=nullpointer, *ptr3 = new char [len+1], *newtoken;
                s.copy(ptr3,len,0);
                ptr3[len] = 0;
                tokenPrt3 = strtok_s(ptr3," ",&newtoken);
                while (tokenPrt3 != nullpointer){
                    kontadns++;
                    if (kontadns>=ndens) ldens=true;
                    dens[kontadns-1]=atof(tokenPrt3);
                    tokenPrt3 = strtok_s(nullpointer," ",&newtoken);
                }
                delete [] ptr3;
            }
        }
/*    Searches MP2 spin density matrix        */
        if(lDSpMP2 && !(s.find("Spin")==string::npos)
            &&  !(s.find("MP2")==string::npos)
            && !(s.find("Density")==string::npos)
            && !(s.find("N=")==string::npos) ) {
            ldens = false ;
            len = s.length();
            char *tokenPrt=nullpointer, *ptr = new char [len+1], *tokenPrt2=nullpointer, *newtoken;
            lncen = true;
            s.copy(ptr,len,0);
            ptr[len] = 0;
            tokenPrt = strtok_s(ptr," ",&newtoken);
            while (tokenPrt != nullpointer){
                tokenPrt2 = tokenPrt;
                tokenPrt = strtok_s(nullpointer," ",&newtoken);
            }
            ndens = atoi(tokenPrt2);

            if (lallocateden){
                delete [] dens;
                delete [] dmat;
            }
            if (ndens != nbasis*(nbasis+1)/2){
                outinterface << "Error allocating MP2 density matrix ";
                outinterface << "Wrong number of elements, nbasis = " << nbasis << " ndens = " << ndens << endl ;
                exit(1);
            }
            dens = new double[ndens];
            dmat = new double[nbasis*nbasis];
            lallocateden = true;

            delete [] ptr;
            kontadns = 0;
            while(getline(inputfile,s)&&!ldens){
                lineknt++;
                len = s.length();
                char *tokenPrt3=nullpointer, *ptr3 = new char [len+1], *newtoken;
                s.copy(ptr3,len,0);
                ptr3[len] = 0;
                tokenPrt3 = strtok_s(ptr3," ",&newtoken);
                while (tokenPrt3 != nullpointer){
                    kontadns++;
                    if (kontadns>=ndens) ldens=true;
                    dens[kontadns-1]=atof(tokenPrt3);
                    tokenPrt3 = strtok_s(nullpointer," ",&newtoken);
                }
                delete [] ptr3;
            }
        }
        /*    Searches CC density matrix      */
        if(lDTotCC &&  !(s.find("Total")==string::npos)
            &&  !(s.find("CC")==string::npos)
            && !(s.find("Density")==string::npos)
            && !(s.find("N=")==string::npos) ) {
//outinterface << "aqui lDTotCC: " <<  std::flush << endl;
            ldens = false ;
            len = s.length();
            char *tokenPrt=nullpointer, *ptr = new char [len+1], *tokenPrt2=nullpointer, *newtoken;
            lncen = true;
            s.copy(ptr,len,0);
            ptr[len] = 0;
            tokenPrt = strtok_s(ptr," ",&newtoken);
            while (tokenPrt != nullpointer){
                tokenPrt2 = tokenPrt;
                tokenPrt = strtok_s(nullpointer," ",&newtoken);
            }
            ndens = atoi(tokenPrt2);

            if (lallocateden){
                delete [] dens;
                delete [] dmat;
            }
            if (ndens != nbasis*(nbasis+1)/2){
                outinterface << "Error allocating CC density matrix ";
                outinterface << "Wrong number of elements, nbasis = " << nbasis << " ndens = " << ndens << endl ;
                exit(1);
            }
            dens = new double[ndens];
            dmat = new double[nbasis*nbasis];
            lallocateden = true;

            delete [] ptr;
            kontadns = 0;
            while(getline(inputfile,s)&&!ldens){
                lineknt++;
                len = s.length();
                char *tokenPrt3=nullpointer, *ptr3 = new char [len+1], *newtoken;
                s.copy(ptr3,len,0);
                ptr3[len] = 0;
                tokenPrt3 = strtok_s(ptr3," ",&newtoken);
                while (tokenPrt3 != nullpointer){
                    kontadns++;
                    if (kontadns>=ndens) ldens=true;
                    dens[kontadns-1]=atof(tokenPrt3);
                    tokenPrt3 = strtok_s(nullpointer," ",&newtoken);
                }
                delete [] ptr3;
            }
        }
/*    Searches CC spin density matrix        */
        if(lDSpCC && !(s.find("Spin")==string::npos)
            &&  !(s.find("CC")==string::npos)
            && !(s.find("Density")==string::npos)
            && !(s.find("N=")==string::npos) ) {
            ldens = false ;
            len = s.length();
            char *tokenPrt=nullpointer, *ptr = new char [len+1], *tokenPrt2=nullpointer, *newtoken;
            lncen = true;
            s.copy(ptr,len,0);
            ptr[len] = 0;
            tokenPrt = strtok_s(ptr," ",&newtoken);
            while (tokenPrt != nullpointer){
                tokenPrt2 = tokenPrt;
                tokenPrt = strtok_s(nullpointer," ",&newtoken);
            }
            ndens = atoi(tokenPrt2);

            if (lallocateden){
                delete [] dens;
                delete [] dmat;
            }
            if (ndens != nbasis*(nbasis+1)/2){
                outinterface << "Error allocating CC density matrix ";
                outinterface << "Wrong number of elements, nbasis = " << nbasis << " ndens = " << ndens << endl ;
                exit(1);
            }
            dens = new double[ndens];
            dmat = new double[nbasis*nbasis];
            lallocateden = true;

            delete [] ptr;
            kontadns = 0;
            while(getline(inputfile,s)&&!ldens){
                lineknt++;
                len = s.length();
                char *tokenPrt3=nullpointer, *ptr3 = new char [len+1], *newtoken;
                s.copy(ptr3,len,0);
                ptr3[len] = 0;
                tokenPrt3 = strtok_s(ptr3," ",&newtoken);
                while (tokenPrt3 != nullpointer){
                    kontadns++;
                    if (kontadns>=ndens) ldens=true;
                    dens[kontadns-1]=atof(tokenPrt3);
                    tokenPrt3 = strtok_s(nullpointer," ",&newtoken);
                }
                delete [] ptr3;
            }
        }
    }              /*  End of main loop (while) */

    if (!ldens){
        outinterface << "Insufficient number of density matrix elements" << endl;
        outinterface << "Number of elements read: " << kontadns << endl ;
        outinterface << "It should be: " << ndens << endl ;
        exit(1);
    }

    if((nbasaux+nsp) != nbasis){
        outinterface << endl << "nbasaux = " << nbasaux << "   nsp = " << nsp ;
        outinterface << endl << "nbasis = " << nbasis << " != nbasaux+nsp = " << nbasaux+nsp << endl ;
        outinterface << "Check whether the calculation has been done in Cartesians " << endl;
        outinterface << "G-DAM is only prepared for spherical functions " << endl;
        outinterface << "Try to run Gaussian with options 5D 7F in the basis set tag " << endl;
        exit(1);
    }

/* In  ROHF  builds the density matrix from the coefficients 	*/
/* of  alpha  MO temporarily stored in  dmat					*/

    if (lDTotSCF &&  lROHF){
        double sum; int r, s, i, knt=0 , kntprt=0;
        for(r = 0; r < nbasis ; r++){
            for(s = 0; s <= r ; s++, knt++){
                for(i = 0, sum = 0. ; i < nbeta ; i++){
                    sum += 2. * dmat[r*nbasis+i] * dmat[s*nbasis+i];
                }
                for(i = nbeta ; i < nalfa ; i++){
                    sum += dmat[r*nbasis+i] * dmat[s*nbasis+i];
                }
                dens[knt] = sum;
                kntprt++;
            }
        }
    }

/*   Writes the basis set to file .ggbs  */

    double* normprim=nullpointer;
    facts[0] = .5 * sqrt(PI);
    for(i = 1; i < 50 ; i++){ facts[i] = facts[i-1] * .5 * (2*i+1);}
    konta = 0;
    int kontb = 0;
    for(i=1; i<=ncen;i++){
        ggbsfile << endl << ncontrtot[i-1] << endl;
        for(j=1; j<=ncontr[i-1]; j++){
            konta++;
            int kloop = 0;
            if(lvec[konta-1]==-1) kloop = 1;
            lvec[konta-1]=abs(lvec[konta-1]);
            for(kl = 0; kl <= kloop; kl++){
                int laux;
                if (kl == (kloop-1)) laux = 0;
                else laux = lvec[konta-1];
                ggbsfile << nprimit[konta-1] << " " << laux << endl;
                normprim=new double[nprimit[konta-1]];
                for(k=0; k < nprimit[konta-1] ; k++){
                    normprim[k] = sqrt(2. * pow(2*primexp[k+kontb],laux+1.5) / facts[laux]) ;
                }
                ggbsfile.setf(ios::scientific,ios::floatfield);
                ggbsfile << setprecision(10) ;
                for(k=0,klin = 1; k < nprimit[konta-1] ; k++,klin++){
                    lendl = true;
                    ggbsfile << setprecision(10) << primexp[k+kontb] << " " ;
                    if (klin == 5){
                        ggbsfile << endl ;
                        klin = 0;
                        lendl = false;
                    }
                }
                if(lendl) ggbsfile << endl ;
                for(k=0, klin = 1; k < nprimit[konta-1] ; k++,klin++){
                    lendl = true;
                    if (kl == 0) ggbsfile << setprecision(10) << cfcontr[k+kontb] * normprim[k] << " " ;
                    else ggbsfile << setprecision(10) << cfcontrP[k+kontb] * normprim[k] << " " ;
                    if (klin == 5){
                        ggbsfile << endl ;
                        klin = 0;
                        lendl = false;
                    }
                }
                delete [] normprim;
                if(lendl) ggbsfile << endl ;
                ggbsfile.unsetf(ios::scientific);
            }
            kontb = kontb + nprimit[konta-1];
        }
    }

/*   Loads the array with density matrix into an auxiliary matrix */

    int ijmax, ijmin;
    ggbsfile << setprecision(10) ;
    ggbsfile.setf(ios::scientific,ios::floatfield);
    for(i=0; i<nbasis;i++){
        for(j=0; j<nbasis; j++){
            if( i > j){ ijmax = i; ijmin = j;}
            else{ ijmax = j; ijmin = i;}
            dmat[i*nbasis+j] = dens[((ijmax+1)*ijmax)/2+ijmin];
        }
    }
    ggbsfile << endl;

/*   Sorts the GAUSSIAN-ordered density matrix: M = (0,+1,-1,+2,-2,...)
     to canonical order: (...,-2,-1,0,1,2,...) */

    ki = 0;
    double *mataux= new double[(2*lmax+1)*(2*lmax+1)];
    for(i=0 ; i<nshelltot ; i++){
        kj = 0;
        for(j=0 ; j<nshelltot ; j++){
            i1 = 0;
            switch(lvectot[i]){
                case 0:
                switch(lvectot[j]){
                    case 0:
                        mataux[0] = dmat[ki*nbasis+kj];
                        break;
                    case 1:
                        mataux[0] = dmat[ki*nbasis+(kj+1)] ;
                        mataux[1] = dmat[ki*nbasis+(kj+2)] ;
                        mataux[2] = dmat[ki*nbasis+kj]   ;
                        break;
                    default:
                        j11 = 0;
                        for(jj = 2*lvectot[j] ; jj >= 0 ; jj -= 2, j11++){
                            mataux[j11] = dmat[ki*nbasis+(kj+jj)];
                        }
                        for(jj = 1 ; jj <= 2*lvectot[j]-1 ; jj += 2, j11++){
                            mataux[j11] = dmat[ki*nbasis+(kj+jj)];
                        }
                    }
                    break;
                case 1:
                switch(lvectot[j]){
                    case 0:
                        mataux[0] = dmat[(ki+1)*nbasis+kj] ;
                        mataux[(2*lvectot[j]+1)]   = dmat[(ki+2)*nbasis+kj] ;
                        mataux[2*(2*lvectot[j]+1)] = dmat[ki*nbasis+kj]   ;
                        break;
                    case 1:
                        mataux[0] = dmat[(ki+1)*nbasis+(kj+1)]  ;
                        mataux[1] = dmat[(ki+1)*nbasis+(kj+2)]  ;
                        mataux[2] = dmat[(ki+1)*nbasis+kj]    ;
                        mataux[(2*lvectot[j]+1)]     = dmat[(ki+2)*nbasis+(kj+1)]  ;
                        mataux[(2*lvectot[j]+1)+1]   = dmat[(ki+2)*nbasis+(kj+2)]  ;
                        mataux[(2*lvectot[j]+1)+2]   = dmat[(ki+2)*nbasis+kj]    ;
                        mataux[2*(2*lvectot[j]+1)]   = dmat[ki*nbasis+(kj+1)]    ;
                        mataux[2*(2*lvectot[j]+1)+1] = dmat[ki*nbasis+(kj+2)]    ;
                        mataux[2*(2*lvectot[j]+1)+2] = dmat[ki*nbasis+kj]      ;
                        break;
                    default:
                        j11 = 0;
                        for(jj = 2*lvectot[j] ; jj >= 0 ; jj -= 2, j11++){
                            mataux[j11] = dmat[(ki+1)*nbasis+(kj+jj)];
                            mataux[(2*lvectot[j]+1)+j11]   = dmat[(ki+2)*nbasis+(kj+jj)];
                            mataux[2*(2*lvectot[j]+1)+j11] = dmat[ki*nbasis+(kj+jj)]  ;
                        }
                        for(jj = 1 ; jj <= 2*lvectot[j]-1 ; jj += 2, j11++){
                            mataux[j11] = dmat[(ki+1)*nbasis+(kj+jj)];
                            mataux[(2*lvectot[j]+1)+j11]   = dmat[(ki+2)*nbasis+(kj+jj)];
                            mataux[2*(2*lvectot[j]+1)+j11] = dmat[ki*nbasis+(kj+jj)]  ;
                        }
                    }
                    break;
                default:
                    switch(lvectot[j]){
                        case 0:
                        for(ii = 2*lvectot[i] ; ii >= 0 ; ii -= 2, i1++){
                            mataux[i1*(2*lvectot[j]+1)] = dmat[(ki+ii)*nbasis+kj];
                        }
                        for(ii = 1; ii <= 2*lvectot[i]-1 ; ii += 2, i1++){
                            mataux[i1*(2*lvectot[j]+1)] = dmat[(ki+ii)*nbasis+kj];
                        }
                        break;
                    case 1:
                        for(ii = 2*lvectot[i] ; ii >= 0 ; ii -= 2, i1++){
                            mataux[i1*(2*lvectot[j]+1)]   = dmat[(ki+ii)*nbasis+(kj+1)];
                            mataux[i1*(2*lvectot[j]+1)+1] = dmat[(ki+ii)*nbasis+(kj+2)];
                            mataux[i1*(2*lvectot[j]+1)+2] = dmat[(ki+ii)*nbasis+kj]  ;
                        }
                        for(ii = 1; ii <= 2*lvectot[i]-1 ; ii += 2, i1++){
                            mataux[i1*(2*lvectot[j]+1)]   = dmat[(ki+ii)*nbasis+(kj+1)];
                            mataux[i1*(2*lvectot[j]+1)+1] = dmat[(ki+ii)*nbasis+(kj+2)];
                            mataux[i1*(2*lvectot[j]+1)+2] = dmat[(ki+ii)*nbasis+kj]  ;
                        }
                        break;
                    default:
                        for(ii = 2*lvectot[i] ; ii >= 0 ; ii -= 2, i1++){
                            j11 = 0;
                            for(jj = 2*lvectot[j] ; jj >= 0 ; jj -= 2, j11++){
                                mataux[i1*(2*lvectot[j]+1)+j11] = dmat[(ki+ii)*nbasis+(kj+jj)];
                            }
                            for(jj = 1 ; jj <= 2*lvectot[j]-1 ; jj += 2, j11++){
                                mataux[i1*(2*lvectot[j]+1)+j11] = dmat[(ki+ii)*nbasis+(kj+jj)];
                            }
                        }
                        for(ii = 1; ii <= 2*lvectot[i]-1 ; ii += 2, i1++){
                            j11 = 0;
                            for(jj = 2*lvectot[j] ; jj >= 0 ; jj -= 2, j11++){
                               mataux[i1*(2*lvectot[j]+1)+j11] = dmat[(ki+ii)*nbasis+(kj+jj)];
                            }
                            for(jj = 1 ; jj <= 2*lvectot[j]-1 ; jj += 2, j11++){
                                mataux[i1*(2*lvectot[j]+1)+j11] = dmat[(ki+ii)*nbasis+(kj+jj)];
                            }
                        }
                    }
                    break;
                }
                for(ii=0 ; ii <= 2*lvectot[i] ; ii++){
                    for(jj=0 ; jj <= 2*lvectot[j] ; jj++){
                        dmat[(ki+ii)*nbasis+(kj+jj)] = mataux[ii*(2*lvectot[j]+1)+jj];
                    }
                }
                kj = kj + 2*lvectot[j]+1;
            }
        ki = ki + 2*lvectot[i]+1;
    }
    delete [] mataux;

/*   Writes the density matrix to file .den   */
    denfile << nbasis << endl;
    denfile << setprecision(10) ;
    denfile.setf(ios::scientific,ios::floatfield);
    for(i=0; i<nbasis;i++){
        for(j=0,klin=1; j<=i; j++, klin++){
            denfile << dmat[i+j*nbasis] << " ";
            if (klin == 5){
                denfile << endl ;
                klin = 0;
            }
        }
        if (klin != 1) denfile << endl ;
    }
    denfile.unsetf(ios::scientific);

/*   Prints out some statistics  */
    outinterface << "GAUSS_interface output" << endl;
    outinterface << "======================" << endl<< endl;
    outinterface << "input path = " << inpath << "\noutput path = " << outpath  << endl << endl;
    outinterface << "PROJECT NAME: " <<  newprojectname << endl << endl;
    outinterface << "No Alpha Electrons: " << nalfa << endl ;
    outinterface << "No Beta Electrons: " << nbeta << endl ;
    if (lROHF) outinterface << "ROHF" << endl ;
    outinterface << "Num. of centers = " << ncen << endl;
    outinterface << "Num. of primitives = " << nprimtot << endl;
    outinterface << "Num. of contractions = " << nshelltot << endl;
    outinterface << "Num. of basis functions = " << nbasis << endl;
    if (!lcfbeta){
        outinterface << "Num. of Molecular Orbitals = " << nindep << endl;
    }
    else
    {
        outinterface << "Num. of alpha Molecular Orbitals = " << nindep << endl;
        outinterface << "Num. of beta Molecular Orbitals = " << nindep << endl;
    }
    outinterface << endl ;
    switch (mden) {
        case 1:
            outinterface << "Total SCF density matrix in file " << denfname << endl ;
            break ;
        case 2:
            outinterface << "SCF spin density matrix in file " << denfname << endl ;
            break ;
        case 3:
            outinterface << "Total CI density matrix in file " << denfname << endl ;
            break ;
        case 4:
            outinterface << "CI spin density matrix in file " << denfname << endl ;
            break ;
        case 5:
            outinterface << "Total MP2 density matrix in file " << denfname << endl ;
            break ;
        case 6:
            outinterface << "MP2 spin density matrix in file " << denfname << endl ;
            break ;
    }
    outinterface.close();
	
	return(0);
}
