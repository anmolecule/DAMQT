#ifndef GLOBALINFO_H
#define GLOBALINFO_H
#include <QString>

extern bool iswindows;  // true if Windows operative system, false otherwise
extern bool mpi;        // true if MPI available
extern int maxnumprocessors; // Highest number of processors
extern QString mpicommand;
extern QString mpiflags;

#endif // GLOBALINFO_H
