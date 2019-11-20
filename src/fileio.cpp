/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Oct 23, 2019

  Standalone cpp Poisson solver

*/

//****************** fileio.cpp **************************

#include "fileio.h"

double GetDoubleParam(string fname, string parnam, double defval, int verbose)
{
  string parval = ReadPar(fname, parnam);
  if (parval.compare("") == 0)
    {
      if (verbose > 1)
	{
	  printf("Parameter %s not found in %s. Using default: %f.\n", parnam.c_str(), fname.c_str(), defval);
	}
      return defval;
    }
  else
    {
      return atof(parval.c_str());
    }
}

int GetIntParam(string fname, string parnam, int defval, int verbose)
{
  string parval = ReadPar(fname, parnam);
  if (parval.compare("") == 0)
    {
      if (verbose > 1)
	{
	  printf("Parameter %s not found in %s. Using default: %d.\n", parnam.c_str(), fname.c_str(), defval);
	}
      return defval;
    }
  else
    {
      return atoi(parval.c_str());
    }
}

double* GetDoubleList(string fname, string parnam, int numvals, double* defval, int verbose)
{
  string parval = ReadPar(fname, parnam);
  size_t found;
  string first;
  int i;
  if (parval.compare("") == 0)
    {
      if (verbose > 1)
	{
	  printf("Parameter %s not found in %s. Using default.\n", parnam.c_str(), fname.c_str());
	}
      return defval;
    }
  else
    {
      for (i=0; i<numvals-1; i++)
	{
	    found = parval.find_first_of(" \t");
	    first = parval.substr(0,found);
	    defval[i] = atof(first.c_str());
	    parval = parval.substr(found+1,parval.size());
	}
      defval[numvals-1] = atof(parval.c_str());
      return defval;
    }
}

int* GetIntList(string fname, string parnam, int numvals, int* defval, int verbose)
{
  string parval = ReadPar(fname, parnam);
  size_t found;
  string first;
  int i;
  if (parval.compare("") == 0)
    {
      if (verbose > 1)
	{
	  printf("Parameter %s not found in %s. Using default.\n", parnam.c_str(), fname.c_str());
	}
      return defval;
    }
  else
    {
      for (i=0; i<numvals-1; i++)
	{
	    found = parval.find_first_of(" \t");
	    first = parval.substr(0,found);
	    defval[i] = atoi(first.c_str());
	    parval = parval.substr(found+1,parval.size());
	}
      defval[numvals-1] = atoi(parval.c_str());
      return defval;
    }
}

string  GetStringParam(string fname, string parnam, string defval, int verbose)
{
  string parval = ReadPar(fname, parnam);
  if (parval.compare("") == 0)
    {
      if (verbose > 1)
	{
	  printf("Parameter %s not found in %s. Using default: %s.\n", parnam.c_str(), fname.c_str(), defval.c_str());
	}
      return defval;
    }
  else
    {
      return parval;
    }
}

string ReadPar(string fname, string parnam)
{
  string line;
  ifstream myfile(fname.c_str());
  string first, second;
  size_t found;
    if (myfile.is_open())
      {
	while (! myfile.eof() )
	  {
	    getline (myfile,line);
	    // First strip out comments
	    found = line.find_first_of("#");
	    if (found != string::npos && found < line.size()) line.erase(found, line.size());
	    // Now split the line at the equals sign (if it exists)
	    found = line.find("=");
	    if (found == string::npos) continue;
	    first = line.substr(0,found);
	    second = line.substr(found+1,line.size());
	    // Now strip leading spaces or tabs
	    found = first.find_first_not_of(" \t");
	    if (found > 0) first.erase(0,found);
	    // Now strip trailing spaces or tabs
	    found = first.find_first_of(" \t");
	    if (found >0 && found < first.size()) first.erase(found, first.size());
	    if (first.compare(parnam) == 0) // if it matches the one we're searching for
	      {
		// Now strip leading spaces or tabs
		found = second.find_first_not_of(" \t");
		if (found > 0) second.erase(0,found);
		return second;
	      }
	  }
	//close the stream:
	myfile.close();
      }
    else printf("Unable to open file");
    second = "";
  return second;
}

