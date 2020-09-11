/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Oct 23, 2019

  Standalone cpp Poisson solver

*/

//****************** multigrid.cpp **************************

#include "multigrid.h"

MultiGrid::MultiGrid(string inname) //Constructor
{
  // This reads in the data from the poisson.cfg
  // file, sets the initial conditions, and solves
  // Poisson's equation using SOR and Multi-Grid methods
  double setup_time, solution_time, efield_time, trace_time;
  time_t time1, time2;
  string StepNum;
  string underscore = "_";
  int n, m, m_init = 0;
  time1 = time(NULL);
  //Set the random number seed

  // First we read in the configuration information
  ReadConfigurationFile(inname);
  printf("Finished Reading config file\n");
  unsigned int seed = (unsigned int)Seed;
  printf("Seed = %d\n",seed);
  srand48(seed);

  // Then, we build the multigrid arrays and set the initial conditions
  phi = new Array3D*[nsteps+1];
  rho = new Array3D*[nsteps+1];
  elec = new Array3D*[nsteps+1];
  hole = new Array3D*[nsteps+1];
  eps = new Array3D*[nsteps+1];
  BCType = new Array2DInt*[nsteps+1];
  QFe = new Array2D*[nsteps+1];
  QFh = new Array2D*[nsteps+1];
  Ckmin = new Array2DInt*[nsteps+1];
  Vkmin = new Array2DInt*[nsteps+1];
  E = new Array3D*[3];

  BuildArrays(phi, rho, elec, hole, eps, E, BCType, QFe, QFh, Ckmin, Vkmin);
  // Save coordinate grid along each axis.
  SaveGrid();
  // If requested, the loop below saves the grids for all scales
  if (SaveMultiGrids == 1)
    {
      SaveGridMulti();  
    }
  printf("Finished Building Arrays. \n");
  //CountCharges(rho, elec, hole);  
  Setkmins(rho, phi, eps, Ckmin, Vkmin);
  Channelkmin = rho[0]->ChannelCkmin;
  ChannelStopkmin = rho[0]->ChannelStopCkmin;
  for (n=0; n<nsteps+1; n++)
    {
      SetInitialVoltages(phi[n], eps[n], BCType[n], Vkmin[n], Ckmin[n]);
      SetFixedCharges(rho[n], Ckmin[n]); // Place fixed charges
    }
  Set_QFh(QFh); // Set hole Quasi-Fermi level
  time2 = time(NULL);
  setup_time = difftime(time2, time1);
  printf("Finished Setting ICs. Setup time = %.3f seconds.\n",setup_time);
  fflush(stdout);
  if (NumberofPixelRegions > 0 && ElectronMethod == 1)
    {
      if (BuildQFeLookup == 1)
	{
	  printf("Building QFe Lookup table\n");
	  Adjust_QFe(QFe, elec, Ckmin);
	  VCycle_Inner(phi, rho, elec, hole, eps, BCType, QFe, QFh, Ckmin, Vkmin, 0, nsteps);
	  Adjust_QFe(QFe, elec, Ckmin);        
	  WriteQFeLookup(outputfiledir, outputfilebase, "QFe");
	  StepNum = std::to_string(999);      
	  WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Elec", elec[0]);
	  WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Hole", hole[0]);
	  BuildQFeLookup = 0;
	}
      else
	{
	  ReadQFeLookup(outputfiledir, outputfilebase, "QFe");
	}
    }
  Adjust_QFe(QFe, elec, Ckmin); // Set QFe	  	    
  // Now we run NumSteps cycle, adding NumElec electrons each step and re-solving
  // Poisson's equation at each step.
  // If a Continuation is requested, we read in the existing results and start
  // from where we left off.
  if (Continuation == 1)
    {
      StepNum = std::to_string(LastContinuationStep);
      ReadOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "phi", phi[0]);      
      ReadOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Elec", elec[0]);
      ReadOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Hole", hole[0]);        
      m_init = LastContinuationStep + 1;
    }

  for (m=m_init; m<NumSteps; m++)
    {
      time1 = time(NULL);
      // Now we cycle through the VCycles to solve Poisson's equation
      for (int n=0; n<iterations; n++)
	{
          printf("--- Starting VCycle %d / %d of step %d / %d.\n",
            n + 1, iterations, m + 1, NumSteps);
	  if (ElectronMethod == 1)
	    {
	      Adjust_QFe(QFe, elec, Ckmin); // Set QFe	  	  
	    }
	  if (ElectronMethod == 0 && m > 0)
	    {
	      // In this case, only iterate at the finest grid.
	      VCycle_Inner(phi, rho, elec, hole, eps, BCType, QFe, QFh, Ckmin, Vkmin, m, 0);
	    }
	  else
	    {
	      VCycle_Inner(phi, rho, elec, hole, eps, BCType, QFe, QFh, Ckmin, Vkmin, m, nsteps);
	    }
	  //CountCharges(rho, elec, hole);
	}
      StepNum = std::to_string(m);
      time2 = time(NULL);
      solution_time = difftime(time2, time1);
      printf("Finished solving Poisson's equation. Solution time = %.3f seconds\n",solution_time);
      time1 = time(NULL);
      // Next we calculate the E fields
      Gradient(phi[0], E);
      time2 = time(NULL);
      efield_time = difftime(time2, time1);
      printf("Finished calculating E Fields. E Field time = %.3f seconds\n",efield_time);
      // Now, we write out the potential and charge density results
      if (SaveElec !=0 && m % SaveElec == 0)
	{
	  WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Elec", elec[0]);
	  WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Hole", hole[0]);
	}
      
      if (SaveData !=0 && m % SaveData == 0)
	{
	  WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "phi", phi[0]);
	  WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "rho", rho[0]);
	  if (VerboseLevel > 2)
	    {
	      WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "eps", eps[0]);	  
	      Write2DFile(outputfiledir, outputfilebase+underscore+StepNum, "QFe", QFe[0]);
	      Write2DFile(outputfiledir, outputfilebase+underscore+StepNum, "QFh", QFh[0]);
	      Write2DIntFile(outputfiledir, outputfilebase+underscore+StepNum, "BCType", BCType[0]);
	      Write2DIntFile(outputfiledir, outputfilebase+underscore+StepNum, "Vkmin", Vkmin[0]);
	      Write2DIntFile(outputfiledir, outputfilebase+underscore+StepNum, "Ckmin", Ckmin[0]);
	    }
	      if(LogEField > 0)
	    {
	      
	      WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Ex", E[0]);
	      WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Ey", E[1]);
	      WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Ez", E[2]);
	    }
	}
      time1 = time(NULL);
      // Now we trace the electrons.
      if (NumberofPixelRegions > 0)
	{
	  if (PixelBoundaryTestType == 0)
	    {
	      TraceGrid(m);
	    }
	  if (PixelBoundaryTestType == 1)
	    {
	      TraceSpot(m);
	      WriteCollectedCharge(outputfiledir, outputfilebase+underscore+StepNum, "CC");
	    }
	  if (PixelBoundaryTestType == 2 || PixelBoundaryTestType == 4)
	    {
	      TraceRegion(m);
	    }
	  if (PixelBoundaryTestType == 3)
	    {
	      TraceFringes(m);
	      WriteCollectedCharge(outputfiledir, outputfilebase+underscore+StepNum, "CC");
	    }
	  if (PixelBoundaryTestType == 5)
	    {
	      TraceList(m);
	      WriteCollectedCharge(outputfiledir, outputfilebase+underscore+StepNum, "CC");		    }
	  if (PixelBoundaryTestType == 6)
	    {
	      TraceFe55Cloud(m);
	      WriteCollectedCharge(outputfiledir, outputfilebase+underscore+StepNum, "CC");		    }
	  // Calculate pixel areas after tracing electrons, if requested.
	  if (PixelAreas >= 0 && (m % PixelAreas) == 0)
	    {
	      CalculatePixelAreas(m);
	    }
	  time2 = time(NULL);
	  trace_time = difftime(time2, time1);
	  printf("Finished tracing electrons. Trace time = %.3f seconds\n",trace_time);
	}
    }
  return;
}

MultiGrid::~MultiGrid() //Destructor                                                                                            
{
  int n, m;
  for (n=0; n<nsteps+1; n++)
    {
      delete phi[n];
      delete rho[n];
      delete elec[n];
      delete hole[n];
      delete eps[n];            
      delete BCType[n];
      delete QFe[n];
      delete QFh[n];                  
      delete Ckmin[n];
      delete Vkmin[n];
    }
  for (n=0; n<3; n++)
    {
      delete E[n];
    }
  delete[] E;
  delete[] phi;
  delete[] rho;
  delete[] elec;
  delete[] hole;
  delete[] eps;
  delete[] BCType;
  delete[] QFe;
  delete[] QFh;
  delete[] Ckmin;
  delete[] Vkmin;
  delete[] SimulationRegionLowerLeft;
  if (ChannelProfile > 0)
    {
      delete[] ChannelDose;
      delete[] ChannelSigma;
      delete[] ChannelPeak;      
    }
  if (ChannelStopProfile > 0)
    {
      delete[] ChannelStopDose;
      delete[] ChannelStopSigma;
      delete[] ChannelStopPeak;      
    }
  delete[] PixelBoundaryLowerLeft;
  delete[] PixelBoundaryUpperRight;
  delete[] PixelBoundaryStepSize;
  if (ElectronMethod == 1)
    {
      delete[] QFeLookup;
    }
  for (n=0; n<NumberofFixedRegions; n++)
    {
      delete FixedRegionLowerLeft[n];
      delete FixedRegionUpperRight[n];      
    }
  delete[] FixedRegionLowerLeft;
  delete[] FixedRegionUpperRight;      
  delete[] FixedRegionVoltage;
  delete[] FixedRegionDoping;
  delete[] FixedRegionOxide;
  delete[] FixedRegionQFe;
  delete[] FixedRegionQFh;
  delete[] FixedRegionBCType;
  for (n=0; n<NumberofPixelRegions; n++)
    {
      delete PixelRegionLowerLeft[n];
      delete PixelRegionUpperRight[n];
      for (m=0; m<NumberofFilledWells[n]; m++)
	{
	  delete FilledPixelCoords[n][m];
	}
      delete[] FilledPixelCoords[n];
      delete[] CollectedCharge[n];
      delete[] ElectronCount[n];
      if (ElectronMethod == 1)
	{
	  delete[] PixelQFe[n];
	}
    }
  delete[] NumberofFilledWells;
  delete[] PixelRegionLowerLeft;
  delete[] PixelRegionUpperRight;
  delete[] FilledPixelCoords;
  delete[] CollectedCharge;
  delete[] ElectronCount;
  if (ElectronMethod == 1)
    {
      delete[] PixelQFe;
    }
  return;
}

void MultiGrid::SaveGrid() {
    printf("Saving coordinate grids.\n");
    string grid_name = outputfiledir + "/grid_";

    // The phi, rho, Ex, Ey, Ez arrays all use the same grid.
    Array3D *A = phi[0];

    string xgrid_name = grid_name + "x.dat";
    ofstream xgrid_out(xgrid_name.c_str());
    xgrid_out << setw(8) << "index" << setw(16) << "lower edge" << setw(16) << "grid center" << setw(16) << "upper edge" << endl;
    for(int i = 0; i < A->nx; ++i) {
      xgrid_out << setw(8) << i << setw(16) << A->x[i] - 0.5 * A->dx << setw(16) << A->x[i] << setw(16) << A->x[i] + 0.5 * A->dx << endl;
    }
    xgrid_out.close();

    string ygrid_name = grid_name + "y.dat";
    ofstream ygrid_out(ygrid_name.c_str());
    ygrid_out << setw(8) << "index" << setw(16) << "lower edge" << setw(16) << "grid center" << setw(16) << "upper edge" << endl;
    for(int i = 0; i < A->ny; ++i) {
      ygrid_out << setw(8) << i << setw(16) << A->y[i] - 0.5 * A->dy << setw(16) << A->y[i] << setw(16) << A->y[i] + 0.5 * A->dy << endl;
    }
    ygrid_out.close();

    string zgrid_name = grid_name + "z.dat";
    ofstream zgrid_out(zgrid_name.c_str());
    zgrid_out << setw(8) << "index" << setw(16) << "lower edge" << setw(16) << "grid center" << setw(16) << "upper edge" << endl;
    for(int i = 0; i < A->nz; ++i) {
      zgrid_out << setw(8) << i << setw(16) << A->zmz[i] << setw(16) << A->z[i] << setw(16) << A->zpz[i] << endl;
    }
    zgrid_out.close();

}

void MultiGrid::SaveGridMulti() {
    printf("Saving coordinate grids at all scales.\n");
    string grid_name = outputfiledir + "/multigrid_";
    Array3D *A;
    for (int ii = 0; ii < nsteps+1; ii++)
      {
    
	// The phi, rho, Ex, Ey, Ez arrays all use the same grid.
	A = phi[ii];
	
	string xgrid_name = grid_name + std::to_string(ii) + "_x.dat";
	ofstream xgrid_out(xgrid_name.c_str());
	xgrid_out << setw(8) << "index" << setw(16) << "lower edge" << setw(16) << "grid center" << setw(16) << "upper edge" << endl;
	for(int i = 0; i < A->nx; ++i) {
	  xgrid_out << setw(8) << i << setw(16) << A->x[i] - 0.5 * A->dx << setw(16) << A->x[i] << setw(16) << A->x[i] + 0.5 * A->dx << endl;
	}
	xgrid_out.close();
	
	string ygrid_name = grid_name + std::to_string(ii) + "_y.dat";
	ofstream ygrid_out(ygrid_name.c_str());
	ygrid_out << setw(8) << "index" << setw(16) << "lower edge" << setw(16) << "grid center" << setw(16) << "upper edge" << endl;
	for(int i = 0; i < A->ny; ++i) {
	  ygrid_out << setw(8) << i << setw(16) << A->y[i] - 0.5 * A->dy << setw(16) << A->y[i] << setw(16) << A->y[i] + 0.5 * A->dy << endl;
	}
	ygrid_out.close();
	
	string zgrid_name = grid_name + std::to_string(ii) + "_z.dat";
	ofstream zgrid_out(zgrid_name.c_str());
	zgrid_out << setw(8) << "index" << setw(16) << "lower edge" << setw(16) << "grid center" << setw(16) << "upper edge" << endl;
	for(int i = 0; i < A->nz; ++i) {
	  zgrid_out << setw(8) << i << setw(16) << A->zmz[i] << setw(16) << A->z[i] << setw(16) << A->zpz[i] << endl;
	}
	zgrid_out.close();
      }
}

void MultiGrid::ReadPhotonList(string outputfiledir, string name)
{
  // This reads the PhotonList
  int i = 0, success = 0;
  double x, y, dxdz, dydz, lambda;
  string line, filename = outputfiledir+'/'+name;
  printf("filename = %s\n", filename.c_str());
  fflush(stdout);

  ifstream pl_in(filename.c_str());
  if (pl_in.is_open())
    {
      getline (pl_in,line); //Skip first line, which is labels
      while (!pl_in.eof())
	{
	  getline (pl_in,line);
	  i += 1;
	}
      NumPhotons = i-1;
      pl_in.clear();
      pl_in.seekg (0, ios::beg); //Reset back to line 1
      getline (pl_in,line);//Skip first line, which is labels
      PhotonListx = new double[NumPhotons];
      PhotonListy = new double[NumPhotons];
      PhotonListdxdz = new double[NumPhotons];
      PhotonListdydz = new double[NumPhotons];
      PhotonListlambda = new double[NumPhotons];  
      
      while (!pl_in.eof())
	{
	  getline (pl_in,line);
	  istringstream iss(line);
	  //printf("Line %d = %s\n",i,line.c_str());
	  if (!(iss >> i >> x >> y >> dxdz >> dydz >> lambda))
	    {
	      success = 0;
	      break; // error
	    }
	  if (i < 0 || i > NumPhotons - 1)
	    {
	      success = 0;
	      break; // error
	    }
	  PhotonListx[i] = x;
	  PhotonListy[i] = y;
	  PhotonListdxdz[i] = dxdz;
	  PhotonListdydz[i] = dydz;
	  PhotonListlambda[i] = lambda;      
	}
      if (i == NumPhotons - 1) success = 1;
      pl_in.close();
    }
  if (success == 1)
    {
      printf("File %s successfully read. Contains %d photons\n", filename.c_str(), NumPhotons);
      return;
    }
  else
    {
      printf("Problem reading file %s. Quitting\n", filename.c_str());
      exit(0);
    }
}


void MultiGrid::WriteQFeLookup(string outputfiledir, string filenamebase, string name)
{
  // This writes the QFeLookup file
  int i, n;
  double qf, dqf;
  dqf = (QFemax - QFemin) / ((double)NQFe - 1.0);
  string underscore = "_", slash = "/", filename, step;
  filename = outputfiledir+slash+filenamebase+underscore+name+".dat";
  ofstream qf_out(filename.c_str());
  qf_out << setw(16) << "QFe" ;
  for (n=0; n<nsteps+1; n++)
    {
      qf_out << setw(16) << "Ne("+ std::to_string(n)+")";
    }
  qf_out << endl;	 
  for (i=0; i<NQFe; i++)
    {
      qf = QFemin + dqf * (double)i;
      qf_out << setw(16) << qf ;
      for (n=0; n<nsteps+1; n++)
	{
	  qf_out << setw(16) << QFeLookup[n * NQFe + i];
	}
        qf_out << endl;	 
    }
  qf_out.close();
  printf("File %s successfully written\n", filename.c_str());
}

void MultiGrid::ReadQFeLookup(string outputfiledir, string filenamebase, string name)
{
  // This reads the QFeLookup file
  int n, i = 0, success = 0;
  double qfin, dqf, qfcalc, ne;
  dqf = (QFemax - QFemin) / ((double)NQFe - 1.0);  
  string line, underscore = "_", slash = "/", filename;
  filename = outputfiledir+slash+filenamebase+underscore+name+".dat";
  ifstream qf_in(filename.c_str());
  if (qf_in.is_open())
    {
      getline (qf_in,line); // First line is labels      
      while (!qf_in.eof())
	{
	  getline (qf_in,line);
	  istringstream iss(line);
	  if (VerboseLevel > 2) printf("In ReadQFeLookup. Line %d = %s\n",i,line.c_str());
	  iss >> qfin;
	  qfcalc = QFemin + dqf * (double)i;
	  if (VerboseLevel > 2) printf("In ReadQFeLookup. Line %d. qfin = %f, qfcalc = %f\n",i,qfin, qfcalc);
	  if (fabs(qfcalc - qfin) > 1.0E-4) break; // error	  
	  for (n=0; n<nsteps+1; n++)
	    {
	      iss >> ne;
	      QFeLookup[n * NQFe + i] = ne;
	    }
	  i += 1;
	}
      if (i == NQFe) success = 1;
    }
  qf_in.close();
  if (success == 1)
    {
      printf("File %s successfully read\n", filename.c_str());
      return;
    }
  else
    {
      printf("Problem reading file %s. Quitting\n", filename.c_str());
      exit(0);
    }
}


void MultiGrid::WriteCollectedCharge(string outputfiledir, string filenamebase, string name)
{
  // This writes the CollectedCharge file
  int m, q, PixX, PixY;
  string underscore = "_", slash = "/", filename;
  filename = outputfiledir+slash+filenamebase+underscore+name+".dat";
  ofstream cc_out(filename.c_str());
  cc_out << setw(8) << "PixX" << setw(16) << "PixY" << setw(16) << "electrons" << endl;
  for (m=0; m<NumberofPixelRegions; m++)
    {
      for (q=0; q<NumberofFilledWells[m]; q++)
	{
	  // First figure out the pixel coordinates
	  PixX = (int)floor((FilledPixelCoords[m][q][0] - PixelRegionLowerLeft[m][0]) / PixelSizeX);
	  PixY = (int)floor((FilledPixelCoords[m][q][1] - PixelRegionLowerLeft[m][1]) / PixelSizeY);
	  cc_out << setw(8) << PixX << setw(16) << PixY << setw(16) << CollectedCharge[m][q] << endl;
	}
    }
  cc_out.close();
  printf("File %s successfully written\n", filename.c_str());
}


void MultiGrid::ReadConfigurationFile(string inname)
{
  VerboseLevel = GetIntParam(inname, "VerboseLevel", 1, 2); // 0 - minimal output, 1 - normal, 2 - more verbose.
  outputfilebase  = GetStringParam(inname,"outputfilebase", "Test", VerboseLevel); //Output filename base
  outputfiledir  = GetStringParam(inname,"outputfiledir", "data", VerboseLevel); //Output filename directory
  int i, j, k, n, jj;
  ScaleFactor =  GetIntParam(inname, "ScaleFactor", 1, VerboseLevel);     // Power of 2 that sets the grid size
  // ScaleFactor = 1 means grid size is 5/6 micron, 128 grids in the z-direction
  nsteps = 2 + (int)(log2(ScaleFactor));  
  printf("There will be %d total multi-grids\n",(nsteps+1));  
  // nsteps is the number of reduction steps in the Vcycle_Inner Cycle.
  qfh = GetDoubleParam(inname, "qfh", -100.0, VerboseLevel);
  qfh1 = GetDoubleParam(inname, "qfh1", -100.0, VerboseLevel);
  UseDoubleQFh =  GetIntParam(inname, "UseDoubleQFh", 0, VerboseLevel);
  DoubleQFhZmax = GetDoubleParam(inname, "DoubleQFhZmax", 0.0, VerboseLevel);  
  // Poisson solver constants
  w = GetDoubleParam(inname, "w", 1.9, VerboseLevel);			// Successive Over-Relaxation factor
  ncycle = GetIntParam(inname, "ncycle", 100, VerboseLevel);		// Number of SOR cycles at each resolution
  iterations =  GetIntParam(inname, "iterations", 1, VerboseLevel);	// Number of VCycles
  Seed =  GetIntParam(inname, "Seed", 77, VerboseLevel);	// Seed
  NZExp = GetDoubleParam(inname, "NZExp", 10.0, VerboseLevel);        // Non-linear z axis exponent
  
  // Overall setup
  NumSteps = GetIntParam(inname, "NumSteps", 100, VerboseLevel);
  SaveData =  GetIntParam(inname, "SaveData", 1, VerboseLevel);     // 0 - Save only Pts, N save phi,rho,E every Nth step
  SaveElec =  GetIntParam(inname, "SaveElec", 1, VerboseLevel);     // 0 - Save only Pts, N save Elec every Nth step
  SaveMultiGrids =  GetIntParam(inname, "SaveMultiGrids", 0, VerboseLevel);     // 1 - Save all of the grids at all scales  
  PixelSizeX = GetDoubleParam(inname, "PixelSizeX", -1.0, VerboseLevel);    // Pixel size in microns
  PixelSizeY = GetDoubleParam(inname, "PixelSizeY", -1.0, VerboseLevel);    // Pixel size in microns    
  GridsPerPixelX = GetIntParam(inname, "GridsPerPixelX", 16, VerboseLevel); // Grids per pixel at ScaleFactor = 1
  GridsPerPixelX = GridsPerPixelX * ScaleFactor;
  GridSpacingX = PixelSizeX / (double)GridsPerPixelX;
  GridsPerPixelY = GetIntParam(inname, "GridsPerPixelY", 16, VerboseLevel); // Grids per pixel at ScaleFactor = 1
  GridsPerPixelY = GridsPerPixelY * ScaleFactor;
  GridSpacingY = PixelSizeY / (double)GridsPerPixelY;

  if (PixelSizeX < 0.0)
    {
      // This is for backward compatibility with square only pixels
      PixelSizeX = GetDoubleParam(inname, "PixelSize", 10.0, VerboseLevel);    // Pixel size in microns
      PixelSizeY = GetDoubleParam(inname, "PixelSize", 10.0, VerboseLevel);    // Pixel size in microns      
      GridsPerPixelX = GetIntParam(inname, "GridsPerPixel", 12, VerboseLevel); // Grids per pixel at ScaleFactor = 1
      GridsPerPixelX = GridsPerPixelX * ScaleFactor;
      GridSpacingX = PixelSizeX / (double)GridsPerPixelX;
      GridsPerPixelY = GetIntParam(inname, "GridsPerPixel", 12, VerboseLevel); // Grids per pixel at ScaleFactor = 1
      GridsPerPixelY = GridsPerPixelY * ScaleFactor;
      GridSpacingY = PixelSizeY / (double)GridsPerPixelY;

    }
  SensorThickness = GetDoubleParam(inname, "SensorThickness", 100.0, VerboseLevel);    // Sensor Thickness in microns  

  Nx = GetIntParam(inname, "Nx", 160, VerboseLevel);                // Number of grids in x at ScaleFactor = 1
  Nx = Nx * ScaleFactor;
  Ny = GetIntParam(inname, "Ny", 160, VerboseLevel);                // Number of grids in y at ScaleFactor = 1
  Ny = Ny * ScaleFactor;
  Nz = GetIntParam(inname, "Nz", 160, VerboseLevel);                // Number of grids in z at ScaleFactor = 1
  Nz = Nz * ScaleFactor;
  Nzelec = GetIntParam(inname, "Nzelec", 32, VerboseLevel);                // Number of grids in z in electron and hole grids at ScaleFactor = 1
  Nzelec = Nzelec * ScaleFactor;
  XBCType = GetIntParam(inname, "XBCType", 1, VerboseLevel);        // 0 - Free BC, 1 - Periodic BC
  YBCType = GetIntParam(inname, "YBCType", 1, VerboseLevel);        // 0 - Free BC, 1 - Periodic BC
  SimulationRegionLowerLeft = new double[2];
  SimulationRegionLowerLeft[0] = 0.0;  SimulationRegionLowerLeft[1] = 0.0;
  SimulationRegionLowerLeft = GetDoubleList(inname, "SimulationRegionLowerLeft", 2, SimulationRegionLowerLeft, VerboseLevel);
  // Voltages and Charges
  Vbb = GetDoubleParam(inname, "Vbb", -50.0, VerboseLevel);		        // Back bias
  Vparallel_lo = GetDoubleParam(inname, "Vparallel_lo", -8.0, VerboseLevel);	// Parallel Low Voltage
  Vparallel_hi = GetDoubleParam(inname, "Vparallel_hi", 4.0, VerboseLevel);	// Parallel High Voltage
  GateOxide = GetDoubleParam(inname, "GateOxide", 0.15, VerboseLevel);
  FieldOxide = GetDoubleParam(inname, "FieldOxide", 0.40, VerboseLevel);
  FieldOxideTaper = GetDoubleParam(inname, "FieldOxideTaper", 0.50, VerboseLevel);  
  ChannelStopWidth = GetDoubleParam(inname, "ChannelStopWidth", 1.0, VerboseLevel);
  BackgroundDoping = GetDoubleParam(inname, "BackgroundDoping", -1.0E12, VerboseLevel);
  ChannelSurfaceCharge = GetDoubleParam(inname, "ChannelSurfaceCharge", 0.0, VerboseLevel);
  ChannelStopSurfaceCharge = GetDoubleParam(inname, "ChannelStopSurfaceCharge", 0.0, VerboseLevel);
  ChannelStopSideDiff = GetDoubleParam(inname, "ChannelStopSideDiff", FieldOxideTaper, VerboseLevel);  
  NumPhases = GetIntParam(inname, "NumPhases", 3, VerboseLevel);
  CollectingPhases = GetIntParam(inname, "CollectingPhases", 1, VerboseLevel);
  GateGap = GetDoubleParam(inname, "GateGap", 0.0, VerboseLevel);    
  if(CollectingPhases > NumPhases)
    {
      printf("Cannot have more collecting phases than phases!  Quitting\n");
      exit(0);
    }
  string profilenum, regionnum;
  ChannelProfile = GetIntParam(inname, "ChannelProfile", 0, VerboseLevel);
  if (ChannelProfile == 0)
    {
      ChannelDoping = GetDoubleParam(inname, "ChannelDoping", 5.0E11, VerboseLevel);
      ChannelDepth = GetDoubleParam(inname, "ChannelDepth", 1.0, VerboseLevel);
    }
  else
    {
      ChannelDose = new double[ChannelProfile];
      ChannelSigma = new double[ChannelProfile];
      ChannelPeak = new double[ChannelProfile];            
      for (i=0; i<ChannelProfile; i++)
	{
	  profilenum = std::to_string(i);
	  ChannelDose[i] = GetDoubleParam(inname, "ChannelDose_"+profilenum, 5.0E11, VerboseLevel);
	  ChannelSigma[i] = GetDoubleParam(inname, "ChannelSigma_"+profilenum, 0.5, VerboseLevel);
	  ChannelPeak[i] = GetDoubleParam(inname, "ChannelPeak_"+profilenum, 0.0, VerboseLevel);	  
	}
    }
  ChannelStopProfile = GetIntParam(inname, "ChannelStopProfile", 0, VerboseLevel);
  if (ChannelStopProfile == 0)
    {
      ChannelStopDoping = GetDoubleParam(inname, "ChannelStopDoping", 5.0E11, VerboseLevel);
      ChannelStopDepth = GetDoubleParam(inname, "ChannelStopDepth", 1.0, VerboseLevel);
    }
  else
    {
      ChannelStopDose = new double[ChannelStopProfile];
      ChannelStopSigma = new double[ChannelStopProfile];
      ChannelStopPeak = new double[ChannelStopProfile];            
      for (i=0; i<ChannelStopProfile; i++)
	{
	  profilenum = std::to_string(i);
	  ChannelStopDose[i] = GetDoubleParam(inname, "ChannelStopDose_"+profilenum, 5.0E11, VerboseLevel);
	  ChannelStopSigma[i] = GetDoubleParam(inname, "ChannelStopSigma_"+profilenum, 0.5, VerboseLevel);
	  ChannelStopPeak[i] = GetDoubleParam(inname, "ChannelStopPeak_"+profilenum, 0.0, VerboseLevel);	  
	}
    }
  // Channel Stop "Dots", if used
  ChannelStopDotHeight = GetDoubleParam(inname, "ChannelStopDotHeight", 0.0, VerboseLevel);    
  if (ChannelStopDotHeight > 0.001)
    {
      ChannelStopDotCenter = GetDoubleParam(inname, "ChannelStopDotCenter", 5.0, VerboseLevel);  
      ChannelStopDotSurfaceCharge = GetDoubleParam(inname, "ChannelStopDotSurfaceCharge", 0.0, VerboseLevel);
      ChannelStopDotProfile = GetIntParam(inname, "ChannelStopDotProfile", 0, VerboseLevel);
      if (ChannelStopDotProfile == 0)
	{
	  ChannelStopDotDoping = GetDoubleParam(inname, "ChannelStopDotDoping", 5.0E11, VerboseLevel);
	  ChannelStopDotDepth = GetDoubleParam(inname, "ChannelStopDotDepth", 1.0, VerboseLevel);
	}
      else
	{
	  ChannelStopDotDose = new double[ChannelStopDotProfile];
	  ChannelStopDotSigma = new double[ChannelStopDotProfile];
	  ChannelStopDotPeak = new double[ChannelStopDotProfile];            
	  for (i=0; i<ChannelStopDotProfile; i++)
	    {
	      profilenum = std::to_string(i);
	      ChannelStopDotDose[i] = GetDoubleParam(inname, "ChannelStopDotDose_"+profilenum, 5.0E11, VerboseLevel);
	      ChannelStopDotSigma[i] = GetDoubleParam(inname, "ChannelStopDotSigma_"+profilenum, 0.5, VerboseLevel);
	      ChannelStopDotPeak[i] = GetDoubleParam(inname, "ChannelStopDotPeak_"+profilenum, 0.0, VerboseLevel);	  
	    }
	}
    }
  AddTreeRings = GetIntParam(inname, "AddTreeRings", 0, VerboseLevel);
  if (AddTreeRings == 1)
    {
      TreeRingAngle = GetDoubleParam(inname, "TreeRingAngle", 0.0, VerboseLevel);	  
      TreeRingAmplitude = GetDoubleParam(inname, "TreeRingAmplitude", 0.0, VerboseLevel);	  
      TreeRingPeriod = GetDoubleParam(inname, "TreeRingPeriod", 0.0, VerboseLevel);	  
    }
  else
    {
      TreeRingAngle = 0.0;
      TreeRingAmplitude = 0.0;
      TreeRingPeriod = 100.0;
    }
  
  // Continuations
  Continuation = GetIntParam(inname, "Continuation", 0, VerboseLevel);
  LastContinuationStep = GetIntParam(inname, "LastContinuationStep", 0, VerboseLevel);    

  // Temperature and diffusion parameters 
  CCDTemperature = GetDoubleParam(inname, "CCDTemperature", 173.0, VerboseLevel);
  Ni = NS * pow(CCDTemperature / 300.0, 1.5) * exp( - QE * EG / (2.0 * KBOLTZMANN * CCDTemperature)) / pow(MICRON_PER_CM, 3.0);   // Intrinsic carrier concentration per micron^3
  ktq = .026 * CCDTemperature / 300.0;//KBOLTZMANN * CCDTemperature / QE;
  printf("Intrinsic carrier concentration Ni = %g in code units, kT/q = %f Volts\n",Ni,ktq);
  DiffMultiplier = GetDoubleParam(inname, "DiffMultiplier", 2.30, VerboseLevel);
  TopAbsorptionProb = GetDoubleParam(inname, "TopAbsorptionProb", 0.0, VerboseLevel);
  SaturationModel = GetIntParam(inname, "SaturationModel", 0, VerboseLevel);
  NumDiffSteps = GetIntParam(inname, "NumDiffSteps", 1, VerboseLevel);
  EquilibrateSteps = GetIntParam(inname, "EquilibrateSteps", 100, VerboseLevel);
  BottomSteps = GetIntParam(inname, "BottomSteps", 1000, VerboseLevel);
  NumVertices = GetIntParam(inname,"NumVertices",2, VerboseLevel);
  CalculateZ0 = GetIntParam(inname,"CalculateZ0",0, VerboseLevel);
  ElectronMethod = GetIntParam(inname,"ElectronMethod",0, VerboseLevel);    
  ElectronZ0Fill = GetDoubleParam(inname,"ElectronZ0Fill",100.0, VerboseLevel);
  ElectronZ0Area = GetDoubleParam(inname,"ElectronZ0Area",100.0, VerboseLevel);
  LogEField = GetIntParam(inname, "LogEField", 0, VerboseLevel);
  LogPixelPaths = GetIntParam(inname, "LogPixelPaths", 0, VerboseLevel);
  PixelAreas = GetIntParam(inname, "PixelAreas", 0, VerboseLevel);
  

  // Pixel Regions
  NumberofPixelRegions = GetIntParam(inname, "NumberofPixelRegions", 0, VerboseLevel);
  if (NumberofPixelRegions > 0)
    {
      bool PixelBoundaryInsideAnyPixelRegion = false, PixelBoundaryInsideThisPixelRegion;
      string fillednum;
      PixelRegionLowerLeft = new double*[NumberofPixelRegions];
      PixelRegionUpperRight = new double*[NumberofPixelRegions];
      NumberofFilledWells = new int[NumberofPixelRegions];
      FilledPixelCoords = new double**[NumberofPixelRegions];
      CollectedCharge = new int*[NumberofPixelRegions];
      ElectronCount = new double*[NumberofPixelRegions];
      if (ElectronMethod == 1)
	{
	  BuildQFeLookup = GetIntParam(inname, "BuildQFeLookup", 0, VerboseLevel);
	  NQFe = GetIntParam(inname, "NQFe", 81, VerboseLevel);
	  QFemin = GetDoubleParam(inname, "QFemin", 5.0, VerboseLevel);
	  QFemax = GetDoubleParam(inname, "QFemax", 10.0, VerboseLevel);
	  QFeLookup = new double[NQFe * (nsteps + 1)];
	  PixelQFe = new double*[NumberofPixelRegions];
	}
      for (i=0; i<NumberofPixelRegions; i++)
	{
	  PixelRegionLowerLeft[i] = new double[2];
	  PixelRegionUpperRight[i] = new double[2];
	  for (j=0; j<2; j++)
	    {
	      PixelRegionLowerLeft[i][j] = 0.0;
	      PixelRegionUpperRight[i][j] = 100.0;
	    }
	}
      
      // Pixel Boundary Tests
      PixelBoundaryTestType = GetIntParam(inname, "PixelBoundaryTestType", 0, VerboseLevel);
      PixelBoundaryLowerLeft = new double(2);
      PixelBoundaryUpperRight = new double(2);
      PixelBoundaryStepSize = new double(2);
      for (j=0; j<2; j++)
	{
	  PixelBoundaryLowerLeft[j] = 0.0;
	  PixelBoundaryUpperRight[j] = 100.0;
	  PixelBoundaryStepSize[j] = 1.0;
	}
      PixelBoundaryLowerLeft = GetDoubleList(inname, "PixelBoundaryLowerLeft", 2, PixelBoundaryLowerLeft, VerboseLevel);
      PixelBoundaryUpperRight = GetDoubleList(inname, "PixelBoundaryUpperRight", 2, PixelBoundaryUpperRight, VerboseLevel);
      PixelBoundaryNx = GetIntParam(inname, "PixelBoundaryNx", 9, VerboseLevel);
      PixelBoundaryNy = GetIntParam(inname, "PixelBoundaryNy", 9, VerboseLevel);
      
      if (PixelBoundaryTestType == 0)
	{
	  for (j=0; j<2; j++)
	    {
	      PixelBoundaryStepSize[j] = 1.0;
	    }
	  PixelBoundaryStepSize = GetDoubleList(inname, "PixelBoundaryStepSize", 2, PixelBoundaryStepSize, VerboseLevel);
	}
      if (PixelBoundaryTestType == 1)
	{
	  NumElec = GetIntParam(inname, "NumElec", 1000, VerboseLevel);
	  Sigmax = GetDoubleParam(inname, "Sigmax", 1.0, VerboseLevel);
	  Sigmay = GetDoubleParam(inname, "Sigmay", 1.0, VerboseLevel);
	  Xoffset = GetDoubleParam(inname, "Xoffset", 0.0, VerboseLevel);
	  Yoffset = GetDoubleParam(inname, "Yoffset", 0.0, VerboseLevel);
	}
      if((PixelBoundaryTestType == 2) || (PixelBoundaryTestType == 4))
	{
	  NumElec = GetIntParam(inname, "NumElec", 1000, VerboseLevel);
	}
      if (PixelBoundaryTestType == 3)
	{
	  NumElec = GetIntParam(inname, "NumElec", 1000, VerboseLevel);	  
	  FringeAngle = GetDoubleParam(inname, "FringeAngle", 0.0, VerboseLevel);	  
	  FringePeriod = GetDoubleParam(inname, "FringePeriod", 0.0, VerboseLevel);	  
	  Xoffset = GetDoubleParam(inname, "Xoffset", 0.0, VerboseLevel);
	  Yoffset = GetDoubleParam(inname, "Yoffset", 0.0, VerboseLevel);
	}
      if (PixelBoundaryTestType == 6)
	{
	  NumElec = GetIntParam(inname, "NumElec", 1620, VerboseLevel);
	  Fe55CloudRadius = GetDoubleParam(inname, "Fe55CloudRadius", 0.2, VerboseLevel);	  
	  Fe55ElectronMult = GetDoubleParam(inname, "Fe55ElectronMult", 1.0, VerboseLevel);
	  Fe55HoleMult = GetDoubleParam(inname, "Fe55HoleMult", 1.0, VerboseLevel);	  	  
	}
      
      for (i=0; i<NumberofPixelRegions; i++)
	{
	  int NewNumberofFilledWells;
	  PixelBoundaryInsideThisPixelRegion = false;
	  regionnum = std::to_string(i);
	  PixelRegionLowerLeft[i] = GetDoubleList(inname, "PixelRegionLowerLeft_"+regionnum, 2, PixelRegionLowerLeft[i], VerboseLevel);
	  PixelRegionUpperRight[i] = GetDoubleList(inname, "PixelRegionUpperRight_"+regionnum, 2, PixelRegionUpperRight[i], VerboseLevel);
	  NumberofFilledWells[i] = GetIntParam(inname, "NumberofFilledWells_"+regionnum, 0, VerboseLevel);
	  NewNumberofFilledWells = NumberofFilledWells[i];
	  // Will reduce number of filled wells if some are inside PixelBoundary
	  if ((PixelBoundaryLowerLeft[0] >= PixelRegionLowerLeft[i][0]) && (PixelBoundaryLowerLeft[1] >= PixelRegionLowerLeft[i][1]) && (PixelBoundaryUpperRight[0] <= PixelRegionUpperRight[i][0]) && (PixelBoundaryUpperRight[1] <= PixelRegionUpperRight[i][1]))
	    {
	      PixelBoundaryInsideAnyPixelRegion = true;
	      PixelBoundaryInsideThisPixelRegion = true;
	      NumberofFilledWells[i] += PixelBoundaryNx * PixelBoundaryNy;
	      // PixelBoundary added to number of wells filled in .cfg file.
	      // Use first PixelBoundaryNx * PixelBoundaryNy locations for PixelBoundary tests,
	      // then any cells filled in .cfg file which are outside pixel boundary
	      NewNumberofFilledWells = NumberofFilledWells[i];
	      // Will reduce number of filled wells if some are inside PixelBoundary
	      if (VerboseLevel >2) printf("In ReadConfig. NewNumberofFilledWells = %d\n",NewNumberofFilledWells);
	    }
	  CollectedCharge[i] = new int[NumberofFilledWells[i]];
	  FilledPixelCoords[i] = new double*[NumberofFilledWells[i]];
	  ElectronCount[i] = new double[NumberofFilledWells[i] * (nsteps + 1)];
	  if (ElectronMethod == 1)
	    {
	      PixelQFe[i] = new double[NumberofFilledWells[i] * (nsteps + 1)];
	    }
	  for (j=0; j<NumberofFilledWells[i]; j++)
	    {
	      int PixX, PixY;
	      PixX = j % PixelBoundaryNx;
	      PixY = (j - PixX) / PixelBoundaryNx;
	      for (n=0; n<nsteps+1; n++)
		{
		  ElectronCount[i][n*NumberofFilledWells[i]+j] = 0.0;
		  if (ElectronMethod == 1)
		    {
		      PixelQFe[i][n*NumberofFilledWells[i]+j] = 100.0;
		    }
		}
	      FilledPixelCoords[i][j] = new double[2];
	      for (k=0; k<2; k++)
		{
		  FilledPixelCoords[i][j][k] = 0.0;
		}
	      if (!PixelBoundaryInsideThisPixelRegion)
		{
		  fillednum = std::to_string(j);
		  CollectedCharge[i][j] = GetIntParam(inname,"CollectedCharge_"+regionnum+"_"+fillednum,0, VerboseLevel);
		  FilledPixelCoords[i][j] = GetDoubleList(inname, "FilledPixelCoords_"+regionnum+"_"+fillednum, 2, FilledPixelCoords[i][j], VerboseLevel);
		}
	      else
		{
		  if (j < PixelBoundaryNx * PixelBoundaryNy)
		    {
		      CollectedCharge[i][j] = 0.0;
		      FilledPixelCoords[i][j][0] = PixelBoundaryLowerLeft[0] + ((double)PixX + 0.5) * PixelSizeX;  
		      FilledPixelCoords[i][j][1] = PixelBoundaryLowerLeft[1] + ((double)PixY + 0.5) * PixelSizeY;   		    
		    }
		  else
		    {
		      fillednum = std::to_string(j -  PixelBoundaryNx * PixelBoundaryNy);
		      CollectedCharge[i][j] = GetIntParam(inname,"CollectedCharge_"+regionnum+"_"+fillednum,0, VerboseLevel);
		      FilledPixelCoords[i][j] = GetDoubleList(inname, "FilledPixelCoords_"+regionnum+"_"+fillednum, 2, FilledPixelCoords[i][j], VerboseLevel);
		      // Check if this is inside the PixelRegion
		      if ((FilledPixelCoords[i][j][0] >= PixelBoundaryLowerLeft[0]) && (FilledPixelCoords[i][j][1] >= PixelBoundaryLowerLeft[1]) && (FilledPixelCoords[i][j][0] <= PixelBoundaryUpperRight[0]) && (FilledPixelCoords[i][j][1] <= PixelBoundaryUpperRight[1]))
			{
			  // This pixel is inside the pixel boundary.  Delete it from the list and move the charge to the appropriate pixel
			  NewNumberofFilledWells -= 1;
			  PixX = (int)floor((FilledPixelCoords[i][j][0] - PixelBoundaryLowerLeft[0]) / PixelSizeX);
			  PixY = (int)floor((FilledPixelCoords[i][j][1] - PixelBoundaryLowerLeft[1]) / PixelSizeY);
			  jj = PixX + PixelBoundaryNx * PixY;
			  CollectedCharge[i][jj] = CollectedCharge[i][j];
			  CollectedCharge[i][j] = 0.0;
			  if (VerboseLevel >2) printf("In ReadConfig. PixX = %d, PixY = %d, jj = %d, NewNumberofFilledWells = %d\n",PixX, PixY, jj, NewNumberofFilledWells);
			}
		    }
		}
	    }
	  NumberofFilledWells[i] = NewNumberofFilledWells;
	}
      if(!PixelBoundaryInsideAnyPixelRegion)
	{
	  printf("Pixel Boundary not Inside any Pixel Region!!!  This may cause errors!!!\n");
	  //exit(0);
	}
    }
  if (PixelBoundaryTestType == 5)
    {
      NumElec = GetIntParam(inname, "NumElec", 1000, VerboseLevel);      
      PhotonList  = GetStringParam(inname,"PhotonList", "None", VerboseLevel); //Photon List file
      if (PhotonList == "None")
	{
	  printf("No PhotonList specified.  Quitting.");
	}
      else
	{
	  printf("Reading Photon List from file %s\n", PhotonList.c_str());	  
	  ReadPhotonList(outputfiledir, PhotonList);
	}
    }

  // Fixed Voltage Regions
  NumberofFixedRegions = GetIntParam(inname, "NumberofFixedRegions", 0, VerboseLevel);
  FixedRegionLowerLeft = new double*[NumberofFixedRegions];
  FixedRegionUpperRight = new double*[NumberofFixedRegions];
  FixedRegionVoltage = new double[NumberofFixedRegions];
  FixedRegionQFe = new double[NumberofFixedRegions];
  FixedRegionQFh = new double[NumberofFixedRegions];    
  FixedRegionDoping = new int[NumberofFixedRegions];
  FixedRegionOxide = new int[NumberofFixedRegions];
  FixedRegionBCType = new int[NumberofFixedRegions];

  for (i=0; i<NumberofFixedRegions; i++)
    {
      FixedRegionLowerLeft[i] = new double[2];
      FixedRegionUpperRight[i] = new double[2];
      for (j=0; j<2; j++)
	{
	  FixedRegionLowerLeft[i][j] = 0.0;
	  FixedRegionUpperRight[i][j] = 100.0;
	}
    }
  for (i=0; i<NumberofFixedRegions; i++)
    {
      regionnum = std::to_string(i);
      FixedRegionLowerLeft[i] = GetDoubleList(inname, "FixedRegionLowerLeft_"+regionnum, 2, FixedRegionLowerLeft[i], VerboseLevel);
      FixedRegionUpperRight[i] = GetDoubleList(inname, "FixedRegionUpperRight_"+regionnum, 2, FixedRegionUpperRight[i], VerboseLevel);
      FixedRegionVoltage[i] = GetDoubleParam(inname, "FixedRegionVoltage_"+regionnum,0.0, VerboseLevel);
      FixedRegionQFe[i] = GetDoubleParam(inname, "FixedRegionQFe_"+regionnum,100.0, VerboseLevel);
      FixedRegionQFh[i] = GetDoubleParam(inname, "FixedRegionQFh_"+regionnum,-100.0, VerboseLevel);            
      FixedRegionDoping[i] = GetIntParam(inname, "FixedRegionDoping_"+regionnum,0, VerboseLevel);
      FixedRegionOxide[i] = GetIntParam(inname, "FixedRegionOxide_"+regionnum,0, VerboseLevel);
      FixedRegionBCType[i] = GetIntParam(inname, "FixedRegionBCType_"+regionnum,0, VerboseLevel);
    }

  if (CalculateZ0 == 1)
    {
      // Filter band configuration.
      FilterBand = GetStringParam(inname, "FilterBand", "none", VerboseLevel);
      if(FilterBand == "u") {
        FilterIndex = 0;
      }
      else if(FilterBand == "g") {
        FilterIndex = 1;
      }
      else if(FilterBand == "r") {
        FilterIndex = 2;
      }
      else if(FilterBand == "i") {
        FilterIndex = 3;
      }
      else if(FilterBand == "z") {
        FilterIndex = 4;
      }
      else if(FilterBand == "y") {
        FilterIndex = 5;
      }
      else {
        printf("No filter response will be used.\n");
        FilterIndex = -1;
      }
      FilterFile = GetStringParam(inname, "FilterFile", "notebooks/depth_pdf.dat", VerboseLevel);
      if(FilterIndex >= 0) {
        ifstream filter_input(FilterFile.c_str());
        string header;
        getline(filter_input, header); // Skip header line
        for(int i = 0; i < n_filter_cdf; i++) {
	  for(int j = 0; j < n_band; j++) {
	    assert(filter_input >> filter_cdf[j * n_filter_cdf + i]);
	  }
        }
        filter_input.close();
      }
    }

  return;
}

void MultiGrid::BuildArrays(Array3D** phi, Array3D** rho, Array3D** elec, Array3D** hole, Array3D** eps, Array3D** E, Array2DInt** BCType, Array2D** QFe, Array2D** QFh, Array2DInt** Ckmin, Array2DInt** Vkmin)
{
  // Builds the multigrid arrays
  int nx, ny, nz, nzelec, nze, nxx, nyy, nzz, n;
  double xmin, xmax, ymin, ymax, zmin, zmax, zmaxelec, dx, dy, dz;
  nxx = Nx + 1;
  nyy = Ny + 1;
  nzz = Nz + 1;
  if (Nzelec == Nz)
    {
      nzelec = nzz;
    }
  else
    {
      nzelec = Nzelec;
    }
  dx = PixelSizeX / (double)GridsPerPixelX;
  dy = PixelSizeY / (double)GridsPerPixelY;  
  dz = SensorThickness / (double)Nz;
  Xmin = SimulationRegionLowerLeft[0] - dx/2;
  Ymin = SimulationRegionLowerLeft[1] - dx/2;
  Zmin = -dz / 2.0;
  Xmax = dx * (double)nxx + Xmin;
  Ymax = dy * (double)nyy + Ymin;
  Zmax = SensorThickness + dz / 2.0;
  phi[0] = new Array3D(Xmin,Xmax,nxx,Ymin,Ymax,nyy,Zmin,Zmax,nzz,NZExp,SensorThickness);
  rho[0] = new Array3D(Xmin,Xmax,nxx,Ymin,Ymax,nyy,Zmin,Zmax,nzz,NZExp,SensorThickness);
  zmaxelec = rho[0]->Z(rho[0]->zp[nzelec] + rho[0]->dzp / 2.0);
  elec[0] = new Array3D(Xmin,Xmax,nxx,Ymin,Ymax,nyy,Zmin,zmaxelec,nzelec,NZExp,SensorThickness);
  hole[0] = new Array3D(Xmin,Xmax,nxx,Ymin,Ymax,nyy,Zmin,zmaxelec,nzelec,NZExp,SensorThickness);
  eps[0] = new Array3D(Xmin,Xmax,nxx,Ymin,Ymax,nyy,Zmin,zmaxelec,nzelec,NZExp,SensorThickness);
  BCType[0] = new Array2DInt(Xmin,Xmax,nxx,Ymin,Ymax,nyy);
  QFe[0] = new Array2D(Xmin,Xmax,nxx,Ymin,Ymax,nyy);
  QFh[0] = new Array2D(Xmin,Xmax,nxx,Ymin,Ymax,nyy);
  Ckmin[0] = new Array2DInt(Xmin,Xmax,nxx,Ymin,Ymax,nyy);
  Vkmin[0] = new Array2DInt(Xmin,Xmax,nxx,Ymin,Ymax,nyy);    

  for (n=1; n<nsteps+1; n++)
    {
      nx = (phi[0]->nx - 1) / (int)pow(2,n) + 1;
      ny = (phi[0]->ny - 1) / (int)pow(2,n) + 1;
      nz = (phi[0]->nz - 1) / (int)pow(2,n) + 1;
      if (Nzelec == Nz)
	{
	  nze = nz;
	}
      else
	{
	  nze = Nzelec/ (int)pow(2,n);      
	}
      dx = phi[0]->dx * (int)pow(2,n);
      dy = phi[0]->dy * (int)pow(2,n);
      dz = phi[0]->dzp * (int)pow(2,n);
      xmin = phi[0]->xmin + phi[0]->dx / 2.0 - dx / 2.0;
      ymin = phi[0]->ymin + phi[0]->dy / 2.0 - dy / 2.0;
      zmin = phi[0]->zmin + phi[0]->dzp / 2.0 - dz / 2.0;
      xmax = phi[0]->xmax - phi[0]->dx / 2.0 + dx / 2.0;
      ymax = phi[0]->ymax - phi[0]->dy / 2.0 + dy / 2.0;
      zmax = phi[0]->zmax - phi[0]->dzp / 2.0 + dz / 2.0;
      phi[n] = new Array3D(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz,NZExp,SensorThickness);
      rho[n] = new Array3D(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz,NZExp,SensorThickness);
      elec[n] = new Array3D(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nze,NZExp,SensorThickness);
      hole[n] = new Array3D(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nze,NZExp,SensorThickness);
      eps[n] = new Array3D(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nze,NZExp,SensorThickness);      
      BCType[n] = new Array2DInt(xmin,xmax,nx,ymin,ymax,ny);
      QFe[n] = new Array2D(xmin,xmax,nx,ymin,ymax,ny);
      QFh[n] = new Array2D(xmin,xmax,nx,ymin,ymax,ny);            
      Ckmin[n] = new Array2DInt(xmin,xmax,nx,ymin,ymax,ny);
      Vkmin[n] = new Array2DInt(xmin,xmax,nx,ymin,ymax,ny);    
    }

  for (n=0; n<3; n++)
    {
      E[n] = new Array3D(Xmin,Xmax,nxx,Ymin,Ymax,nyy,Zmin,Zmax,nzz,NZExp,SensorThickness);
    }
  return;
}

void MultiGrid::SetInitialVoltages(Array3D* phi, Array3D* eps, Array2DInt* BCType, Array2DInt* Vkmin, Array2DInt* Ckmin)
{
  // This sets up the initial voltages on the boundaries
  int i, j, k, n, m, p, index, index2;
  int PixY;
  double PixYmin, GateWidth, eps_factor_gap;
  double Ylower[NumPhases], Yupper[NumPhases];
  double GateVoltage[NumPhases], GapVoltage[NumPhases];
  // Potential on top

  // These "fudge factors" for epsilon in the gate gap are still experimental
  if (phi->dx > 1.20) eps_factor_gap = 10.0;
  else if (phi->dx > 0.60)  eps_factor_gap = 5.5;
  else if (phi->dx > 0.30)  eps_factor_gap = 2.0;
  else eps_factor_gap = 1.3;  
    
  
  for (i=0; i<phi->nx; i++)
    {
      for (j=0; j<phi->ny; j++)
	{
	  index = i + j * phi->nx + (phi->nz - 1) * phi->nx * phi->ny;
	  phi->data[index] = Vbb;
	}
    }
  // Potentials in Pixel Region
  GateWidth = PixelSizeY / (double)NumPhases;
  for (m=0; m<NumPhases; m++)
    {
      GateVoltage[m] = Vparallel_hi;
      GapVoltage[m] = Vparallel_lo;
      if ((NumPhases + CollectingPhases) % 2 == 0)
	{
	  // Even case - gate boundary coincident with pixel edge
	  Ylower[m] = (double)m * GateWidth - phi->dy / 2.0;	
	  Yupper[m] = Ylower[m] + GateWidth;
	  // Need to shift last gate by 1 grid to keep things centered.
	  if (m == NumPhases - 2)
	    {
	      Yupper[m] += phi->dy;
	    }
	  if (m == NumPhases - 1)
	    {
	      Ylower[m] += phi->dy;
	    }
	}
      else
	{
	  // Odd case - gate boundary straddles pixel edge
	  Ylower[m] = ((double)m - 0.5) * GateWidth;
	  Yupper[m] = Ylower[m] + GateWidth;
	}
      //printf("m = %d, Ylower = %f, Yupper = %f\n", m, Ylower[m], Yupper[m]);
    }
  for (p=0; p<(CollectingPhases-1); p++)
    {
      // Gap voltage is only high if both sides are high.
      GapVoltage[p+2] = Vparallel_hi;
    }
  for (p=0; p<(NumPhases-CollectingPhases); p++)
    {
      // Set the outer phases low, bottom, then top, then bottom, then top, ...
      if (p % 2 == 0)
	{
	  GateVoltage[p/2] = Vparallel_lo;
	}
      else
	{
	  GateVoltage[NumPhases - 1 - (p - 1)/2] = Vparallel_lo;
	}
    }
  
  if (VerboseLevel > 2)
    {
      for (m=0; m<NumPhases; m++)
	{
	  printf("In SetInitialVoltages. m = %d, Ylower = %f, Yupper = %f, GateVoltage = %f, GapVoltage = %f\n",m, Ylower[m],Yupper[m],GateVoltage[m], GapVoltage[m]);
	}
    }
  for (n=0; n<NumberofPixelRegions; n++)
    {
      for (i=0; i<phi->nx; i++)
	{
	  if (phi->x[i] < PixelRegionLowerLeft[n][0] || phi->x[i] > PixelRegionUpperRight[n][0])
	    {
	      continue; // If not in PixelRegion, continue
	    }
	  for (j=0; j<phi->ny; j++)
	    {
	      if (phi->y[j] < PixelRegionLowerLeft[n][1] || phi->y[j] > PixelRegionUpperRight[n][1])
		{
		  continue; // If not in PixelRegion, continue
		}
	      index = i + j * phi->nx;
	      BCType->data[index] = -1; // Used to identify the pixel regions	      
	      PixY = (int)floor((phi->y[j] - PixelRegionLowerLeft[n][1]) / PixelSizeY);
	      PixYmin = PixelRegionLowerLeft[n][1] + (double)PixY * PixelSizeY;
	      for (m=0; m<NumPhases; m++)
		{
 		  if (((phi->y[j] >= PixYmin + Ylower[m] - GateGap / 2.0 && phi->y[j] <= PixYmin + Ylower[m] + GateGap / 2.0) || (Ylower[m] < 0.0 && (phi->y[j] >= PixYmin + Ylower[m] + PixelSizeY - GateGap / 2.0 && phi->y[j] <= PixYmin + Ylower[m] + PixelSizeY + GateGap / 2.0))) && GateGap > 0.001) 
 		    {
		      // In the gate gap
		      Vkmin->data[index] = 1;
		      //BCType->data[index] = 1;  // Free BC here
		      phi->data[index] = GapVoltage[m];
		      for (k=Vkmin->data[index]; k<Ckmin->data[index]; k++)		      
			{
			  index2 = index + k * phi->nx * phi->ny;
			  eps->data[index2] = EPSILON_OX / EPSILON_SI * eps_factor_gap;
			}
 		    }
 		  else if ((phi->y[j] >= PixYmin + Ylower[m] && phi->y[j] <= PixYmin + Yupper[m]) || (Ylower[m] < 0.0 && (phi->y[j] >= PixYmin + Ylower[m] + PixelSizeY && phi->y[j] <= PixYmin + Yupper[m] + PixelSizeY)))
 		    {
		      for (k=0; k<Vkmin->data[index]; k++)
			{
			  index2 = index + k * phi->nx * phi->ny;
			  phi->data[index2] = GateVoltage[m];
			}
 		    }
		}
	      if (VerboseLevel >2 && i == (phi->nx-1)/2) printf("In SetInitialVoltages. j = %d, Vkmin = %d, phi[0] = %f\n",j,Vkmin->data[index],phi->data[index]);		      
 	    }
 	}
    }
  // This sets the potentials and the boundary conditions in the Fixed Voltage regions.
  for (n=0; n<NumberofFixedRegions; n++)
    {
      for (i=0; i<phi->nx; i++)
	{
	  for (j=0; j<phi->ny; j++)
	    {
	      index = i + j * phi->nx;
	      if (phi->x[i] >= FixedRegionLowerLeft[n][0] && phi->x[i] <= FixedRegionUpperRight[n][0] && phi->y[j] >= FixedRegionLowerLeft[n][1] && phi->y[j] <= FixedRegionUpperRight[n][1])
		{
		  if (FixedRegionBCType[n] == 0) //  Fixed potential region
		    {
		      for (k=0; k<Vkmin->data[index]; k++)
			{
			  index2 = index + k * phi->nx * phi->ny;
			  phi->data[index2] = FixedRegionVoltage[n];
			}
		    }
		  else // Free Boundary Condition region
		    {
		      BCType->data[index] = 1;
		    }
		}
	    }
	}
    }
  if(VerboseLevel > 2)
    {
      printf("Finished setting Boundary Potentials, \n");
      fflush(stdout);
    }
  return;
}


void MultiGrid::SetFixedCharges(Array3D* rho, Array2DInt* Ckmin)
{
  // This sets the fixed lattice charges
  int i, j, k, n, index, index2, PixX, PixY;
  double PixXmin, PixYmin, DotYmin, DotYmax, ChargeFactor = (QE*MICRON_PER_M/(EPSILON_0*EPSILON_SI)) / pow(MICRON_PER_CM, 3);
  bool AlreadyDoped = false;  
  double TaperRatio;
  int NLeft, NRight, NTaper, LastCS, LastC;
  // ChargeFactor converts doping in cm^-3 into the code units
  // Set the background charge:
  double tr_xprime, tr_background;
  for (i=0; i<rho->nx; i++)
    {
      for (j=0; j<rho->ny; j++)
	{
	  index = i + j * rho->nx;
	  tr_xprime = rho->x[i] * cos(TreeRingAngle * pi / 180.0) + rho->y[j] * sin(TreeRingAngle * pi / 180.0);
	  tr_background = BackgroundDoping * ChargeFactor * (1.0 + TreeRingAmplitude * cos(tr_xprime / TreeRingPeriod));
	  //if (rho->nz == 161 && rho->y[j] == 80) printf("i = %d, rho = %f\n",i,tr_background);
	  for (k=Ckmin->data[index]; k<rho->nz-1; k++)
	    {
	      index2 = index + k * rho->nx * rho->ny;
	      rho->data[index2] = tr_background;
	    }
	}
    }
  // Charges in Pixel Regions// This uses the same basic logic as Setkmins to find the edges and TaperRatio,
  // But here the taper is set by ChannelStopSideDiff instead of FieldOxideTaper.
  for (n=0; n<NumberofPixelRegions; n++)
    {
      NLeft = 0; NRight = 0;
      for (i=0; i<rho->nx; i++)
	{
	  if (rho->x[i] >= PixelRegionLowerLeft[n][0] && rho->x[i] <= PixelRegionUpperRight[n][0])
	    {
	      PixX = (int)floor((rho->x[i] - PixelRegionLowerLeft[n][0]) / PixelSizeX);
	      PixXmin = PixelRegionLowerLeft[n][0] + (double)PixX * PixelSizeX;
	      if (PixX !=0) continue;
	      if (rho->x[i] <= PixXmin + ChannelStopWidth/2.0 || rho->x[i] >= PixXmin + PixelSizeX - ChannelStopWidth/2.0)
		{
		  // Channel Stop region
		  continue;
		}
	      else if (rho->x[i] <= PixXmin + ChannelStopWidth/2.0 + ChannelStopSideDiff)
		{
		  // Left side taper
		  if (VerboseLevel > 2)
		    {
		    printf("In SetFixedCharges. Setting Tapers. Nz = %d, i = %d, Left\n",rho->nz,i);
		    }
		  NLeft += 1;
		}
	      else if (rho->x[i] >= PixXmin + PixelSizeX - ChannelStopWidth/2.0 - ChannelStopSideDiff)
		{
		  // Right side taper
		  if (VerboseLevel > 2)
		    {
		      printf("In SetFixedCharges. Setting Tapers. Nz = %d, i = %d, Right\n",rho->nz,i);
		    }
		  NRight += 1;
		}
	      else
		{
		  // Channel region
		  continue;
		}
	    }
	}
      if (NLeft == NRight)  NTaper = NLeft;
      else NTaper = max(NLeft, NRight);
      LastCS = 0, LastC = i;	  
      if (VerboseLevel > 2) printf("In SetFixedCharges. Nz = %d, NTaper = %d \n",rho->nz,NTaper);
      for (i=0; i<rho->nx; i++)
	{
	  if (rho->x[i] >= PixelRegionLowerLeft[n][0] && rho->x[i] <= PixelRegionUpperRight[n][0])
	    {
	      PixX = (int)floor((rho->x[i] - PixelRegionLowerLeft[n][0]) / PixelSizeX);
	      PixXmin = PixelRegionLowerLeft[n][0] + (double)PixX * PixelSizeX;
	      if (rho->x[i] <= PixXmin + ChannelStopWidth/2.0 || rho->x[i] >= PixXmin + PixelSizeX - ChannelStopWidth/2.0)
		{
		  // Channel Stop region
		  TaperRatio = 1.0;
		  LastCS = i;
		}
	      else if (rho->x[i] <= PixXmin + ChannelStopWidth/2.0 + ChannelStopSideDiff)
		{
		  // Left side taper		      
		  TaperRatio = 1.0 - (double)(i - LastCS) / (double)(NTaper + 1);
		}
	      else if (rho->x[i] >= PixXmin + PixelSizeX - ChannelStopWidth/2.0 - ChannelStopSideDiff)
		{
		  // Right side taper		      		      
		  TaperRatio = (double)(i - LastC) / (double)(NTaper + 1);
		}
	      else
		{
		  // Channel region
		  TaperRatio = 0.0;
		  LastC = i;		  
		  for (j=0; j<rho->ny; j++)
		    {
		      if (rho->y[j] < PixelRegionLowerLeft[n][1] || rho->y[j] > PixelRegionUpperRight[n][1])
			{
			  continue; // If not in PixelRegion, continue
			}
		      SetCharge(rho, Ckmin, i, j, 1, 1.0);
		      if (VerboseLevel >2 && j == (rho->ny-1)/2) printf("In SetFixedCharges. Nz = %d, i = %d, PixX = %d, Ckmin = %d, Channel\n",rho->nz,i,PixX,Ckmin->data[i+j*rho->nx]);
		    }
		  continue;
		}
	      for (j=0; j<rho->ny; j++)
		{
		  if (rho->y[j] < PixelRegionLowerLeft[n][1] || rho->y[j] > PixelRegionUpperRight[n][1])
		    {
		      continue; // If not in PixelRegion, continue
		    }
		  PixY = (int)floor((rho->y[j] - PixelRegionLowerLeft[n][1]) / PixelSizeY);
		  PixYmin = PixelRegionLowerLeft[n][1] + (double)PixY * PixelSizeY;
		  DotYmin = PixYmin + ChannelStopDotCenter - ChannelStopDotHeight / 2.0;
		  DotYmax = PixYmin + ChannelStopDotCenter + ChannelStopDotHeight / 2.0;
		  // Channel Stop region
		  SetCharge(rho, Ckmin, i, j, 2, TaperRatio);		  
		  if (VerboseLevel >2 && j == (rho->ny-1)/2) printf("in SetFixedCharges. Nz = %d, i = %d, PixX = %d, Ckmin = %d, Channel Stop\n",rho->nz,i,PixX,Ckmin->data[i+j*rho->nx]);
		  if (DotYmin > PixYmin) // "Dot" wholly within one pixel
		    {
		      if  ((rho->y[j] >= DotYmin) && (rho->y[j] <= DotYmax))
			{
			  // Channel Stop Dot region
			  SetCharge(rho, Ckmin, i, j, 3, TaperRatio);		  
			  if (VerboseLevel >2) printf("in SetFixedCharges. Nz = %d, i = %d, y=%.2f, PixY = %d, Ckmin = %d, Channel Stop Dot Center\n",rho->nz,i,rho->y[j],PixY,Ckmin->data[i+j*rho->nx]);
			}
		    }
		  else // "Dot" straddles pixel edge
		    {
		      if  ((rho->y[j] <= DotYmax) || (rho->y[j] >= PixelSizeY + DotYmin))
			{
			  // Channel Stop Dot region
			  SetCharge(rho, Ckmin, i, j, 3, TaperRatio);		  
			  if (VerboseLevel >2) printf("in SetFixedCharges. Nz = %d, i = %d, y=%.2f, PixY = %d, Ckmin = %d, Channel Stop Dot Corner\n",rho->nz,i,rho->y[j],PixY,Ckmin->data[i+j*rho->nx]);
			}
		    }
		}
	    }
	}
    }
  //Charges in Fixed regions
  for (n=0; n<NumberofFixedRegions; n++)
    {
      for (i=0; i<rho->nx; i++)
	{
	  for (j=0; j<rho->ny; j++)
	    {
	      if (rho->x[i] >= FixedRegionLowerLeft[n][0] && rho->x[i] <= FixedRegionUpperRight[n][0] && rho->y[j] >= FixedRegionLowerLeft[n][1] && rho->y[j] <= FixedRegionUpperRight[n][1])
		{
		  index = i + j * rho->nx;
		  AlreadyDoped = false;
		  for (k=Ckmin->data[index]; k<rho->nz-1; k++)
		    {
		      // Check if this region is already doped.
		      // Don't want double doping of fixed regions.
		      index2 = index + k * rho->nx * rho->ny;
		      if (fabs(rho->data[index2]) > 1.1 * fabs(BackgroundDoping * ChargeFactor))
			{
			  AlreadyDoped = true;
			  break;
			}
		    }
		  if (AlreadyDoped) continue;
		  if (FixedRegionDoping[n] == 1) // Channel Doping
		    {
		      SetCharge(rho, Ckmin, i, j, 1, 1.0);		  			  
		    }
		  else if (FixedRegionDoping[n] == 2) // Channel Stop Doping
		    {
		      SetCharge(rho, Ckmin, i, j, 2, 1.0);		  			  
		    }
		}
	    }
	}
    }
  if(VerboseLevel > 2)
    {
      printf("Finished setting fixed charges\n");
    }
  return;
}


double MultiGrid::SOR_Inner(Array3D* phi, Array3D* rho, Array3D* elec, Array3D* hole, Array3D* eps, Array2DInt* BCType, Array2D* QFe, Array2D* QFh, Array2DInt* Ckmin, Array2DInt* Vkmin)
{
  // This is the main loop that calculates the potentials in the array through Successive Over Relaxation.
  // An inner Newton's method loop is used to solve the non-linear Quasi-Fermi level equations.
  // Assumes fixed potentials on the top, mixture of fixed and free BC on bottom, free or periodic BC on the sides.
  double newphi, oldnewphi, tol, exponent = 0.0;
  double omw, w6, hsquared, zhsquaredp, zhsquaredm;
  double apx, amx, apy, amy, apz, amz; // Dielectric constant averages
  double SORChargeFactor =  (QE*MICRON_PER_M/(EPSILON_0*EPSILON_SI));
  double MaxExponent = log(1.0E10 / (Ni * SORChargeFactor)); // This limits mobile carrier density from getting too large
  double MinElec = 1.0E-6; // electron concentrations are not adjusted below this value
  double MinDeltaPhi = ktq * log(MinElec / (Ni * SORChargeFactor)); // This controls where we calculate mobile carriers
  double DeltaPhi, MaxDeltaPhi = 0.1;   // This limits the change in phi per cycle and helps convergence.
  double Emax = 100.0; // This controls how quickly the electron concentrations are adjusted.
  double PreFactor = phi->dx / Emax;  
  double nmx, npx, nmy, npy, nmz, npz, dnmx, dnpx, dnmy, dnpy, dnmz, dnpz, deltan;
  bool InnerLoop;
  double ElecCharge=0.0, HoleCharge=0.0, Term1=0.0, Term2=0.0, CellVolume=0.0, AveIterations = 0.0;
  double NumElec=0.0, NumHoles=0.0, TotalHoles=0.0, TotalElectrons=0.0;
  double localQFh = -100.0;
  int nn = 0, mm = 0, iter_counter, iter_limit = 10000, red_black;
  int i, j, k, kstart, kmax, im, ip, j0, jm, jp, nxy, ind, ind2, indmx, indpx, indmy, indpy, indmz, indpz;
  kmax = min(phi->nz - 1, elec->nz-1);
  nxy = phi->nx * phi->ny;
  hsquared =  phi->dx * phi->dy;
  omw = 1.0 - w;
      for (red_black=0; red_black<2; red_black++)
	{
	  for (i=0; i<phi->nx; i++)
	    {
	      if (XBCType == 0)
		{
		  if (i == 0) {im = 0;} else {im = 1;} // Free BC
		  if (i == phi->nx-1) {ip = 0;} else {ip = 1;} // Free BC
		}
	      else
		{
		  if (i == 0) {im = -phi->nx + 2;} else {im = 1;} // Periodic BC
		  if (i == phi->nx-1) {ip = -phi->nx + 2;} else {ip = 1;} // Periodic BC
		}
	      for (j=0; j<phi->ny; j++)
		{
		  j0 = j * phi->nx;
		  if (YBCType == 0)
		    {
		      if (j == 0) {jm = 0;} else {jm = phi->nx;} // Free BC
		      if (j == phi->ny-1) {jp = 0;} else {jp = phi->nx;} // Free BC
		    }
		  else
		    {
		      if (j == 0) {jm = (-phi->ny + 2) * phi->nx;} else {jm = phi->nx;} // Periodic BC
		      if (j == phi->ny-1) {jp = (-phi->ny + 2) * phi->nx;} else {jp = phi->nx;} // Periodic BC
		    }
		  ind2 = i + j0;
		  for (k=Vkmin->data[ind2]; k<kmax; k++)
		    {
		      InnerLoop = false;
		      if ((i + j + k + red_black) % 2 == 0) continue; // Implements Red-Black alternation		  
		      ind = ind2 + k * nxy;
		      indmx = ind - im;
		      indpx = ind + ip;
		      indmy = ind - jm;
		      indpy = ind + jp;
		      if (k == 0 && BCType->data[ind2] > 0) // Free BC at z = 0
			{
			  indmz = ind;
			}
		      else
			{
			  indmz = ind - nxy;
			}
		      indpz = ind + nxy;
		      if (k < eps->nz - 1)
			{
			  apx = (eps->data[ind] + eps->data[indpx]) / 2.0;
			  amx = (eps->data[ind] + eps->data[indmx]) / 2.0;	      
			  apy = (eps->data[ind] + eps->data[indpy]) / 2.0;
			  amy = (eps->data[ind] + eps->data[indmy]) / 2.0;	      
			  apz = (eps->data[ind] + eps->data[indpz]) / 2.0;
			  amz = (eps->data[ind] + eps->data[indmz]) / 2.0;	      
			}
		      else
			{
			  apx = 1.0; amx = 1.0; apy = 1.0; amy = 1.0; apz = 1.0; amz = 1.0;
			}
		      w6 = w / (apx + amx + apy + amy + apz * phi->zplus[k] + amz * phi->zminus[k]);
		      CellVolume = rho->dx * rho->dy * rho->zw[k];
		      Term1 = omw * phi->data[ind] + w6 * (amx * phi->data[indmx] + apx * phi->data[indpx] + amy * phi->data[indmy] + apy * phi->data[indpy] + amz * phi->zminus[k] * phi->data[indmz] + apz * phi->zplus[k] * phi->data[indpz] + hsquared * rho->data[ind]);
		      newphi = phi->data[ind];
		      if (k>=Ckmin->data[ind2] && QFe->data[ind2] > 1.0E5 && elec->data[ind] > MinElec)
			{
			  // Inner Newton loop for electrons when using ElectronMethod == 2
			  // Move electrons to drive QFe toward a constant value
			  // This conserves electrons.
			  // These two terms are needed to compensate for non-linear z-axis
			  zhsquaredp = 4.0 * phi->dzp / (phi->dzpdz[k] + phi->dzpdz[k+1]) * phi->dzp / (phi->dzpdz[k] + phi->dzpdz[k+1]);
			  zhsquaredm = 4.0 * phi->dzp / (phi->dzpdz[k] + phi->dzpdz[k-1]) * phi->dzp / (phi->dzpdz[k] + phi->dzpdz[k-1]);
			  InnerLoop = true;
			  mm++;
			  // n** is the value at the cell boundary
			  // dn** is the number of electrons to be moved across the boundary
			  nmx = (elec->data[ind] + elec->data[indmx]) / 2.0;
			  dnmx = max(-elec->data[indmx], PreFactor / hsquared * (ktq * (elec->data[ind] - elec->data[indmx]) - nmx * (phi->data[ind] - phi->data[indmx])));
			  if (fabs(dnmx) < MinElec) dnmx = 0.0;
			  npx = (elec->data[ind] + elec->data[indpx]) / 2.0;
			  dnpx = max(-elec->data[indpx], PreFactor / hsquared * (ktq * (elec->data[ind] - elec->data[indpx]) - npx * (phi->data[ind] - phi->data[indpx])));
			  if (fabs(dnpx) < MinElec) dnpx = 0.0;
			  nmy = (elec->data[ind] + elec->data[indmy]) / 2.0;
			  dnmy = max(-elec->data[indmy], PreFactor / hsquared * (ktq * (elec->data[ind] - elec->data[indmy]) - nmy * (phi->data[ind] - phi->data[indmy])));
			  if (fabs(dnmy) < MinElec) dnmy = 0.0;
			  npy = (elec->data[ind] + elec->data[indpy]) / 2.0;
			  dnpy = max(-elec->data[indpy], PreFactor / hsquared * (ktq * (elec->data[ind] - elec->data[indpy]) - npy * (phi->data[ind] - phi->data[indpy])));
			  if (fabs(dnpy) < MinElec) dnpy = 0.0;
			  if (k > Ckmin->data[ind2])
			    {
			      // Don't move electrons below Ckmin
			      nmz = (elec->data[ind] + elec->data[indmz]) / 2.0;
			      dnmz = max(-elec->data[indmz], PreFactor / zhsquaredm * (ktq * (elec->data[ind] - elec->data[indmz]) - nmz * (phi->data[ind] - phi->data[indmz])));
			      if (fabs(dnmz) < MinElec) dnmz = 0.0;
			    }
			  else
			    {
			      dnmz = 0.0;
			    }
			  if (k < elec->nz)
			    {
			      // Don't move electrons above elec->nz-1
			      npz = (elec->data[ind] + elec->data[indpz]) / 2.0;
			      dnpz = max(-elec->data[indpz], PreFactor / zhsquaredp * (ktq * (elec->data[ind] - elec->data[indpz]) - npz * (phi->data[ind] - phi->data[indpz])));
			      if (fabs(dnpz) < MinElec) dnpz = 0.0;
			    }
			  else
			    {
			      dnpz = 0.0;
			    }
			  deltan = (dnmx + dnpx + dnmy + dnpy + dnmz + dnpz);
			  
			    if (VerboseLevel > 2 && i == 20 && j == 20)
			    {
			    printf("In SOR_Inner, ElectronMethod=2. i = %d, j = %d, k = %d, kmax = %d, deltan = %f\n",i,j,k,elec->nz-1,deltan);
			    printf("In SOR_Inner, ElectronMethod=2. dnmx=%f,dnpx=%f,dnmy=%f,dnpy=%f,dnmz=%f,dnpz=%f\n",dnmx,dnpx,dnmy,dnpy,dnmz,dnpz);
			    }
			  if (deltan > elec->data[ind])
			    {
			      // Don't allow the cell to go negative
			      // If it would, ratio down the values until the cell is zero.
			      dnmx *= elec->data[ind] / deltan;
			      dnpx *= elec->data[ind] / deltan;
			      dnmy *= elec->data[ind] / deltan;
			      dnpy *= elec->data[ind] / deltan;
			      dnmz *= elec->data[ind] / deltan;
			      dnpz *= elec->data[ind] / deltan;
			      deltan = (dnmx + dnpx + dnmy + dnpy + dnmz + dnpz);
			    }
			  elec->data[ind] -= deltan;
			  // Now update the surrounding cells
			  elec->data[indmx] += dnmx;
			  elec->data[indpx] += dnpx;
			  elec->data[indmy] += dnmy;
			  elec->data[indpy] += dnpy;
			  elec->data[indmz] += dnmz;
			  elec->data[indpz] += dnpz;		      
			  TotalElectrons += fabs(deltan);
			  // Now calculate the new phi.
			  ElecCharge = - elec->data[ind] * SORChargeFactor / CellVolume;
			  newphi = Term1 + w6 * hsquared * ElecCharge;
			  DeltaPhi = max(-MaxDeltaPhi, min(MaxDeltaPhi, newphi - phi->data[ind]));		  
			  newphi = phi->data[ind] + DeltaPhi;
			  if (isnan(newphi) || isnan(elec->data[ind]))
			    {
			      printf("Nan encountered in SOR_Inner electrons at point i,j,k = (%d,%d,%d)\n",i,j,k);
			      printf("Phi = %f, NumElec = %f, ElecCharge = %f, exponent = %f, Term1 = %f, Term2 = %f\n",phi->data[ind], NumElec, ElecCharge, exponent, Term1, Term2);
			      exit(0);
			    }
			}
		      if (k>=Ckmin->data[ind2] && phi->data[ind]>QFe->data[ind2]+MinDeltaPhi)
			{
			  // Inner Newton loop for electrons in Fixed Regions
			  // and in Pixel regions when using ElectronMethod == 1
			  InnerLoop = true;
			  mm++;
			  iter_counter = 0;
			  tol = 1.0;
			  while (tol > 1.0E-9 && iter_counter < iter_limit)
			    {
			      oldnewphi = newphi;
			      exponent = min(MaxExponent, (newphi - QFe->data[ind2]) / ktq);// Prevents ElecCharge getting too large
			      ElecCharge = - Ni * exp(exponent) * SORChargeFactor;
			      Term2 = hsquared * w6 * ElecCharge / ktq;
			      if (fabs(1.0 - Term2) < 1.0E-18) break; 
			      newphi = newphi - (newphi - (w6 * hsquared * ElecCharge + Term1)) / (1.0 - Term2);
			      tol = fabs(newphi - oldnewphi);		      
			      iter_counter++;
			    }
			  nn += iter_counter;
			  NumElec = -ElecCharge / SORChargeFactor * CellVolume;
			  elec->data[ind] = NumElec;
			  TotalElectrons += NumElec;
			  if (isnan(newphi) || isnan(elec->data[ind]))
			    {
			      printf("Nan encountered in SOR_Inner electrons at point i,j,k = (%d,%d,%d)\n",i,j,k);
			      printf("Phi = %f, NumElec = %f, ElecCharge = %f, exponent = %f, Term1 = %f, Term2 = %f\n",phi->data[ind], NumElec, ElecCharge, exponent, Term1, Term2);
			      exit(0);
			    }
			  if (iter_counter >= iter_limit)
			    {
			      printf("Warning electron inner loop iteration exceeded at point i,j,k = (%d,%d,%d)\n",i,j,k);
			    }
			  DeltaPhi = max(-MaxDeltaPhi, min(MaxDeltaPhi, newphi - phi->data[ind]));		  
			  newphi = phi->data[ind] + DeltaPhi;
			}
		      if (!InnerLoop && elec->data[ind] > MinElec)
			{
			  // Inner Newton loop (not really a loop) when using ElectronMethod == 0
			  InnerLoop = true;
			  mm++;
			  // Now calculate the new phi.
			  TotalElectrons += elec->data[ind];
			  ElecCharge = - elec->data[ind] * SORChargeFactor / CellVolume;
			  newphi = Term1 + w6 * hsquared * ElecCharge;
			  DeltaPhi = max(-MaxDeltaPhi, min(MaxDeltaPhi, newphi - phi->data[ind]));		  
			  newphi = phi->data[ind] + DeltaPhi;
			}			  
		      // Allow use of two values of QFh.  A bit of a hack, but works.
		      if (UseDoubleQFh == 0)
			{
			  localQFh = QFh->data[ind2];
			}
		      else
			{
			  if (phi->x[i] >= PixelRegionLowerLeft[0][0] && phi->x[i] <= PixelRegionUpperRight[0][0] && phi->y[j] >= PixelRegionLowerLeft[0][1] && phi->y[j] <= PixelRegionUpperRight[0][1] && phi->z[k] <= DoubleQFhZmax)
			    {
			      localQFh = qfh1;
			    }
			  else
			    {
			      localQFh = QFh->data[ind2];
			    }
			}
		      if (k>=Ckmin->data[ind2] && phi->data[ind]<localQFh-MinDeltaPhi)
			{
			  // Inner Newton loop for holes
			  InnerLoop = true;
			  mm++;
			  iter_counter = 0;
			  tol = 1.0;
			  while (tol > 1.0E-9 && iter_counter < iter_limit)
			    {
			      oldnewphi = newphi;
			      exponent = min(MaxExponent,(localQFh - newphi) / ktq);// Prevents HoleCharge getting too large
			      HoleCharge = Ni * exp(exponent) * SORChargeFactor;
			      Term2 = hsquared * w6 * HoleCharge / ktq;
			      if (fabs(1.0 + Term2) < 1.0E-18) break;
			      newphi = newphi - (newphi - (w6 * hsquared * HoleCharge + Term1)) / (1.0 + Term2);	  
			      tol = fabs(newphi - oldnewphi);		      
			      iter_counter++;
			    }
			  nn += iter_counter;
			  NumHoles = HoleCharge / SORChargeFactor * CellVolume;
			  hole->data[ind] = NumHoles;
			  TotalHoles += NumHoles;
			  if (isnan(newphi) || isnan(hole->data[ind]))
			    {
			      printf("Nan encountered in SOR_Inner holes at point i,j,k = (%d,%d,%d)\n",i,j,k);
			      printf("Phi = %f, NumHoles = %f, HoleCharge = %f\n",phi->data[ind], NumHoles, HoleCharge);
			      exit(0);
			    }
			  if (iter_counter >= iter_limit)
			    {
			      printf("Warning hole inner loop iteration exceeded at point i,j,k = (%d,%d,%d)\n",i,j,k);
			      printf("Phi = %f, NumHoles = %f, HoleCharge = %f\n",phi->data[ind], NumHoles, HoleCharge);
			    }
			  DeltaPhi = max(-MaxDeltaPhi, min(MaxDeltaPhi, newphi - phi->data[ind]));		  
			  newphi = phi->data[ind] + DeltaPhi;
			}
		      if (!InnerLoop)
			{
			  // No inner Newton loop if no free carriers.
			  newphi = Term1;
			  if (isnan(newphi))
			    {
			      printf("Nan encountered no free carriers at point i,j,k = (%d,%d,%d)\n",i,j,k);
			      exit(0);
			    }
			  elec->data[ind] = 0.0;			  
			  hole->data[ind] = 0.0;
			}
		      phi->data[ind] = newphi;
		    }
		  kstart = max(kmax, Vkmin->data[ind2]);
		  for (k=kstart; k<phi->nz-1; k++)
		    {
		      // Above the elec and hole arrays, no inner loop is needed
		      if ((i + j + k + red_black) % 2 == 0) continue; // Implements Red-Black alternation		  
		      ind = ind2 + k * nxy;
		      indmx = ind - im;
		      indpx = ind + ip;
		      indmy = ind - jm;
		      indpy = ind + jp;
		      if (k == 0 && BCType->data[ind2] > 0) // Free BC at z = 0
			{
			  indmz = ind;
			}
		      else
			{
			  indmz = ind - nxy;
			}
		      indpz = ind + nxy;
		      w6 = w / (4.0 + phi->zplus[k] + phi->zminus[k]);
		      newphi = omw * phi->data[ind] + w6 * (phi->data[indmx] + phi->data[indpx] + phi->data[indmy] + phi->data[indpy] + phi->zminus[k] * phi->data[indmz] + phi->zplus[k] * phi->data[indpz] + hsquared * rho->data[ind]);
		      if (isnan(newphi))
			{
			  printf("Nan encountered above mobile arrays at point i,j,k = (%d,%d,%d)\n",i,j,k);
			  exit(0);
			}
		      phi->data[ind] = newphi;
		    }
		}
	    }
	}
      if (VerboseLevel > 1)
	{
	  if (ElectronMethod == 2) printf("Finished SOR_Inner, nz = %d, TotalElectronsMoved = %.1f, TotalHoles = %.1f\n",rho->nz,TotalElectrons,TotalHoles);
	  else printf("Finished SOR_Inner, nz = %d, TotalElectrons = %.1f, TotalHoles = %.1f\n",rho->nz,TotalElectrons,TotalHoles);
	}
  if (mm > 0) AveIterations = (double)nn / (double)mm;
  return AveIterations;
}


double MultiGrid::Error_Inner(Array3D* phi, Array3D* rho, Array3D* elec, Array3D* hole, Array3D* eps, Array2DInt* BCType, Array2DInt* Vkmin)
{
  // This calculates the residual error after completing an SOR cycle
  double RhoChargeFactor =  (QE*MICRON_PER_M/(EPSILON_0*EPSILON_SI)) / (rho->dx * rho->dy);
  double newphi, rhosum, error = 0.0;
  double hsquared, w6;
  double apx, amx, apy, amy, apz, amz; // Dielectric constant averages  
  int i, j, k, nxy, ind, ind2, indmx, indpx, indmy, indpy, indmz, indpz;
  nxy = phi->nx * phi->ny;
  hsquared =  phi->dx * phi->dy;
  for (i=1; i<phi->nx-1; i++)
    { 
      for (j=1; j<phi->ny-1; j++)
	{
	  ind2 = i + j * phi->nx;
	  for (k=Vkmin->data[ind2]; k<phi->nz-1; k++)
	    {
	      ind = ind2 + k * nxy;
	      indmx = ind - 1;
	      indpx = ind + 1;
	      indmy = ind - phi->nx;
	      indpy = ind + phi->nx;
	      if (k == 0 && BCType->data[ind2] > 0) // Free BC at z = 0
		{
		  indmz = ind;
		}
	      else
		{
		  indmz = ind - nxy;
		}
	      indpz = ind + nxy;
	      if (k < elec->nz)
		{
		  rhosum = rho->data[ind] + (hole->data[ind] - elec->data[ind]) * RhoChargeFactor / rho->zw[k];
		}
	      else
		{
		  rhosum = rho->data[ind];
		}
	      if (k < eps->nz - 1)
		{
		  apx = (eps->data[ind] + eps->data[indpx]) / 2.0;
		  amx = (eps->data[ind] + eps->data[indmx]) / 2.0;	      
		  apy = (eps->data[ind] + eps->data[indpy]) / 2.0;
		  amy = (eps->data[ind] + eps->data[indmy]) / 2.0;	      
		  apz = (eps->data[ind] + eps->data[indpz]) / 2.0;
		  amz = (eps->data[ind] + eps->data[indmz]) / 2.0;	      
		}
	      else
		{
		  apx = 1.0; amx = 1.0; apy = 1.0; amy = 1.0; apz = 1.0; amz = 1.0;
		}
	      w6 = 1.0 / (apx + amx + apy + amy + apz * phi->zplus[k] + amz * phi->zminus[k]);
	      newphi = w6 * (amx * phi->data[indmx] + apx * phi->data[indpx] + amy * phi->data[indmy] + apy * phi->data[indpy] + amz * phi->zminus[k] * phi->data[indmz] + apz * phi->zplus[k] * phi->data[indpz] + hsquared * rhosum);
	      error = max(error, fabs(phi->data[ind] - newphi));
	    }
	}
    }
  if(VerboseLevel > 2) printf("In Error_Inner. Grid Size = %d, error = %.3g\n",phi->nx, error);
  return error;
}

void MultiGrid::Prolongate(Array3D* phi, Array3D* newphi, Array3D* elec, Array3D* newelec, Array2DInt* newBCType, Array2DInt* newVkmin, Array2DInt* newCkmin)
{
  // This propagates a solution at one grid scale up to the next finer scale.
  // Assumes fixed potentials on the top, mixture of fixed and free BC on bottom, free or periodic BC on the sides.
  int i, j, k, im, ip, jm, jp, km=0, kp=0, nxy, newnxy;
  int newindex, newindex2=0;
  int indpmm, indpmp, indppm, indppp;
  int indmmm, indmmp, indmpm, indmpp;
  double zrp=0.0, zrm=0.0; // z-axis height ratios
  double newelecsum;
  if (VerboseLevel >2) printf("In Prolongate. Starting.kmax = %d, newkmax = %d\n",elec->nz, newelec->nz);
  nxy = phi->nx * phi->ny;
  newnxy = newphi->nx * newphi->ny;
  for (i=0; i<newphi->nx; i++)
    {
      im = max(0.0,rint((float)i / 2.0 - 0.1));
      ip = min((double)newphi->nx-1,rint((float)i / 2.0 + 0.1));
      for (j=0; j<newphi->ny; j++)
	{
	  jm = max(0.0,rint((float)j / 2.0 - 0.1));
	  jp = min((double)newphi->ny-1,rint((float)j / 2.0 + 0.1));
	  newindex = i + j * newphi->nx;
	  if (newBCType->data[i + j * newphi->nx] == 1) // Free BC at z = 0
	    {
	      indmmm = im + jm * phi->nx;
	      indmmp = im + jm * phi->nx;
	      indmpm = im + jp * phi->nx;
	      indmpp = im + jp * phi->nx;
	      indpmm = ip + jm * phi->nx;
	      indpmp = ip + jm * phi->nx;
	      indppm = ip + jp * phi->nx;
	      indppp = ip + jp * phi->nx;
	      newphi->data[newindex] = (phi->data[indmmm] + phi->data[indmmp] + phi->data[indmpm] + phi->data[indmpp] + phi->data[indpmm] + phi->data[indpmp] + phi->data[indppm] + phi->data[indppp]) / 8.0;
	    }
	  for (k=newVkmin->data[newindex]; k<newphi->nz-1; k++)
	    {
	      km = rint((float)k / 2.0 - 0.1);
	      kp = rint((float)k / 2.0 + 0.1);
	      newindex2 = newindex + k * newnxy;
	      indmmm = im + jm * phi->nx + km * nxy;
	      indmmp = im + jm * phi->nx + kp * nxy;
	      indmpm = im + jp * phi->nx + km * nxy;
	      indmpp = im + jp * phi->nx + kp * nxy;
	      indpmm = ip + jm * phi->nx + km * nxy;
	      indpmp = ip + jm * phi->nx + kp * nxy;
	      indppm = ip + jp * phi->nx + km * nxy;
	      indppp = ip + jp * phi->nx + kp * nxy;
	      newphi->data[newindex2] = (phi->data[indmmm] + phi->data[indmmp] + phi->data[indmpm] + phi->data[indmpp] + phi->data[indpmm] + phi->data[indpmp] + phi->data[indppm] + phi->data[indppp]) / 8.0;
	    }
	  if (ElectronMethod == 2)
	    {
	      // Prolongate the electrons up to the next level in this case.
	      for (k=newCkmin->data[newindex]-1; k<newelec->nz-1; k++)
		{
		  km = rint((float)k / 2.0 - 0.1);
		  kp = rint((float)k / 2.0 + 0.1);
		  zrp = (newphi->zpz[k] - newphi->z[k]) / (16.0 * phi->zw[kp]);
		  zrm = (newphi->z[k] - newphi->zmz[k]) / (16.0 * phi->zw[km]);
		  newindex2 = newindex + k * newnxy;
		  indmmm = im + jm * phi->nx + km * nxy;
		  indmmp = im + jm * phi->nx + kp * nxy;
		  indmpm = im + jp * phi->nx + km * nxy;
		  indmpp = im + jp * phi->nx + kp * nxy;
		  indpmm = ip + jm * phi->nx + km * nxy;
		  indpmp = ip + jm * phi->nx + kp * nxy;
		  indppm = ip + jp * phi->nx + km * nxy;
		  indppp = ip + jp * phi->nx + kp * nxy;
		  newelecsum = (elec->data[indmmm] + elec->data[indmpm] + elec->data[indpmm] + elec->data[indppm]) * zrm + (elec->data[indmmp] + elec->data[indmpp] + elec->data[indpmp] + elec->data[indppp]) * zrp;
		  if (k < newCkmin->data[newindex])
		    {
		      newindex2 = newindex + newCkmin->data[newindex] * newnxy;
		      newelec->data[newindex2] += newelecsum;
		    }
		  else
		    {
		      newelec->data[newindex2] += newelecsum;
		    }
		}
	      if (isnan(newelec->data[newindex2]))
		{
		  printf("Nan encountered in electron prolongate! Quitting!\n");
		  printf("i=%d,j=%d,k=%d,im=%d,jm=%d,km=%d,ip=%d,jp=%d,kp=%d,zrm=%f,zrp=%f\n",i,j,k,im,jm,km,ip,jp,kp,zrm,zrp);
		  exit(0);
		  }
	      
	      if (VerboseLevel > 2 && i == 40 && j == 40 && newelec->data[newindex2] > 1.0E-18)
		{
		  printf("In Prolongate. i=%d,j=%d,k=%d,km=%d,kp=%d,zrp=%f,zrm=%f, newelec=%f\n",i,j,k,km,kp,zrp,zrm,newelec->data[newindex2]);
		}
	    }
	}
    }
  return;
}



void MultiGrid::VCycle_Inner(Array3D** phi, Array3D** rho, Array3D** elec, Array3D** hole, Array3D** eps, Array2DInt** BCType, Array2D** QFe, Array2D** QFh, Array2DInt** Ckmin, Array2DInt** Vkmin, int run, int stepstart)
{
  // This controls the running of the SOR iterations atall grid scales.
  // It does a given number of steps (niter) at each scale  as it moves to finer and finer scales.
  int i, j, niter, NumSOR;
  double error = 100.0, AveIterations;
  bool PixelsFilled = false;  // Stores whether the pixels have been filled with electrons in ElectronMethod = 2
  for (i=stepstart; i>-1; i--)
    {
      niter = ncycle * (int)pow(4,i);      
      NumSOR = 0;
      AveIterations = 0;
      if (ElectronMethod == 2)
	{
	  if (PixelsFilled)
	    {
	      //printf("Pixels are already filled, nz = %d, i = %d\n",rho[i]->nz,i);
	      // If the pixels have been filled at a coarser grid, Prolongate the electrons up to
	      // the current grid and find the solution.
	      FillElectronWells(rho[i], elec[i], Ckmin[i], 0.0);// Empty the wells 	  	  
	      Prolongate(phi[i+1], phi[i], elec[i+1], elec[i], BCType[i], Vkmin[i], Ckmin[i]);
	    }
	  else if (GridsPerPixelX / (int)pow(2,i) < 8 || GridsPerPixelY / (int)pow(2,i) < 8)
	    {
	      // If the pixels haven't been filled, we don't want to fill them until the pixels
	      // are at least 8 grid cells wide.  Otherwise, there is not yet a well-defined
	      // collection region, and the electrons will "leak out".
	      //printf("Pixel less than 8 grid cells, nz = %d, i = %d\n",rho[i]->nz,i);	      
	      if (i < stepstart)
		{
		  // After we have done the corsest grid, we need to Prolongate the solution up and re-solve.
		  FillElectronWells(rho[i], elec[i], Ckmin[i], 0.0);// Empty the wells 	  	  		  
		  Prolongate(phi[i+1], phi[i], elec[i+1], elec[i], BCType[i], Vkmin[i], Ckmin[i]);
		}
	      else
		{
		  // No need to Prolongate if it's the coarsest grid
		  FillElectronWells(rho[i], elec[i], Ckmin[i], 0.0);// Empty the wells 	  	  		  
		}
	    }
	  else
	    {
	      //printf("Pixel >= 8 grid cells, but not yet filled, nz = %d, i = %d\n",rho[i]->nz,i);	      	      
	      if (i < stepstart)
		{
		  // This is the case when we are filling the pixels for the first time.
		  // If it's not the coarsest grid, we Prolongate the solution up and re-solve
		  // to generate an approximate solution before filling the wells.
		  Prolongate(phi[i+1], phi[i], elec[i+1], elec[i], BCType[i], Vkmin[i], Ckmin[i]);
		  for (j=0; j<niter; j++)
		    {
		      AveIterations = SOR_Inner(phi[i], rho[i], elec[i], hole[i], eps[i], BCType[i], QFe[i], QFh[i], Ckmin[i], Vkmin[i]);
		      NumSOR += 1;
		    }
		  FillElectronWells(rho[i], elec[i], Ckmin[i], 1.0);// Fill the wells 	  	  		  
		}
	      else
		{
		  // This is the case when we are filling the pixels for the first time.
		  // If it's the coarsest grid, we generate an approximate solution before filling the wells.
		  for (j=0; j<niter; j++)
		    {
		      AveIterations = SOR_Inner(phi[i], rho[i], elec[i], hole[i], eps[i], BCType[i], QFe[i], QFh[i], Ckmin[i], Vkmin[i]);
		      NumSOR += 1;
		    }
		  FillElectronWells(rho[i], elec[i], Ckmin[i], 1.0);// Fill the wells 	  	  		  
		}
	      PixelsFilled = true; // Set this to true so we don't fill them again.
	    }
	}
      else if (i < stepstart)
	{
	  // Only Prolongate if it's not the coarsest grid.
	  Prolongate(phi[i+1], phi[i], elec[i+1], elec[i], BCType[i], Vkmin[i], Ckmin[i]);
	}
      for (j=0; j<niter; j++)
	{
	  // Now here's where we really generate the solution at this grid.
	  AveIterations = SOR_Inner(phi[i], rho[i], elec[i], hole[i], eps[i], BCType[i], QFe[i], QFh[i], Ckmin[i], Vkmin[i]);
	  NumSOR += 1;
	}
      // And here's where we evaluate the residual error.
      error = Error_Inner(phi[i], rho[i], elec[i], hole[i], eps[i], BCType[i], Vkmin[i]);      
      if(VerboseLevel > 1)
	{
	  printf("Completed iterations at resolution %dx%dx%d. Number of steps = %d. Error = %.3g Ave inner iterations = %.2f\n",phi[i]->nx-1,phi[i]->ny-1,phi[i]->nz-1,j,error,AveIterations/(double)NumSOR);
	  CountCharges(rho, elec, hole);      	  
	}
      if (SaveMultiGrids == 1)
	{
	  // If requested, save all of the coarser multi-grid data
	  string istepnum = std::to_string(i);            
	  string StepNum = std::to_string(run);
	  WriteOutputFile(outputfiledir, outputfilebase+"_Multi_"+istepnum+"_"+StepNum, "phi", phi[i]);
	  WriteOutputFile(outputfiledir, outputfilebase+"_Multi_"+istepnum+"_"+StepNum, "rho", rho[i]);
	  WriteOutputFile(outputfiledir, outputfilebase+"_Multi_"+istepnum+"_"+StepNum, "Elec", elec[i]);
	  WriteOutputFile(outputfiledir, outputfilebase+"_Multi_"+istepnum+"_"+StepNum, "Hole", elec[i]);
	  if (VerboseLevel > 2)
	    {
	      WriteOutputFile(outputfiledir, outputfilebase+"_Multi_"+istepnum+"_"+StepNum, "eps", eps[i]);	  
	      Write2DFile(outputfiledir, outputfilebase+"_Multi_"+istepnum+"_"+StepNum, "QFe", QFe[i]);
	      Write2DFile(outputfiledir, outputfilebase+"_Multi_"+istepnum+"_"+StepNum, "QFh", QFh[i]);
	      Write2DIntFile(outputfiledir, outputfilebase+"_Multi_"+istepnum+"_"+StepNum, "BCType", BCType[i]);
	      Write2DIntFile(outputfiledir, outputfilebase+"_Multi_"+istepnum+"_"+StepNum, "Vkmin", Vkmin[i]);
	      Write2DIntFile(outputfiledir, outputfilebase+"_Multi_"+istepnum+"_"+StepNum, "Ckmin", Ckmin[i]);
	    }
	  
	}
    }
  return;
}


void MultiGrid::WriteOutputFile(string outputfiledir, string filenamebase, string name, Array3D* array)
{
  // This writes the data to the HDF files
  string underscore = "_", slash = "/", hdfname, filename;
  int* int_attr_data  = new int[3];
  double* double_attr_data  = new double[3];
  double* flipped_data  = new double[array->nx*array->ny*array->nz];
  int i, j, k, index, flipped_index;

  // There must be a better way to change from x fast to z fast, but this works.
  for (i=0; i<array->nx; i++)
    {
      for (j=0; j<array->ny; j++)
	{
	  for (k=0; k<array->nz; k++)
	    {
	      index = i + j * array->nx + k * array->nx * array->ny;
	      flipped_index = k + j * array->nz + i * array->nz * array->ny;
	      flipped_data[flipped_index] = array->data[index];
	    }
	}
    }
  hdfname = filenamebase+underscore+name;
  filename = outputfiledir+slash+hdfname+".hdf5";
  WriteHDF5File3(filename, hdfname, array->nx, array->ny, array->nz, flipped_data);
  // Now we write the attributes
  int_attr_data[0] = array->nx; int_attr_data[1] = array->ny; int_attr_data[2] = array->nz;
  WriteHDF5IntAttribute(filename, hdfname, "Dimension", 3, int_attr_data);
  double_attr_data[0] = array->xmin; double_attr_data[1] = array->ymin; double_attr_data[2] = array->zmin;
  WriteHDF5DoubleAttribute(filename, hdfname, "Lower_Left", 3, double_attr_data);
  double_attr_data[0] = array->xmax; double_attr_data[1] = array->ymax; double_attr_data[2] = array->zmax;
  WriteHDF5DoubleAttribute(filename, hdfname, "Upper_Right", 3, double_attr_data);

  delete[] flipped_data;
  return;
}

void MultiGrid::Write3DFile(string outputfiledir, string filenamebase, string name, Array3D* phi)
{
  // This writes a 3D .dat output file
  int i, j, k, index;
  string underscore = "_", slash = "/", filename;
  filename = outputfiledir+slash+filenamebase+underscore+name+".dat";
  ofstream file_out(filename.c_str());
  file_out << setw(16) << "X" << setw(16) << "Y" << setw(16) << "Z" << setw(16) << "data" << endl;
  for (k=0; k<phi->nz; k++)
    {
      for (i=0; i<phi->nx; i++)	        
	{
	  for (j=0; j<phi->ny; j++)
	    {
	      index = i + j * phi->nx + k * phi->nx * phi->ny;
	      file_out << setw(16) << phi->x[i] << setw(16) << phi->y[j] << setw(16) << phi->z[k] << setw(16) << phi->data[index] << endl;
	    }
	}
    }
  file_out.close();
  printf("File %s successfully written\n", filename.c_str());
}

void MultiGrid::Write2DFile(string outputfiledir, string filenamebase, string name, Array2D* QFe)
{
  // This writes a 2D .dat output file
  int i, j, index;
  string underscore = "_", slash = "/", filename;
  filename = outputfiledir+slash+filenamebase+underscore+name+".dat";
  ofstream file_out(filename.c_str());
  file_out << setw(16) << "X" << setw(16) << "Y" << setw(16) << "data" << endl;
  for (i=0; i<QFe->nx; i++)	        
    {
      for (j=0; j<QFe->ny; j++)
	{
	  index = i + j * QFe->nx;
	  file_out << setw(16) << QFe->x[i] << setw(16) << QFe->y[j] << setw(16) << QFe->data[index] << endl;
	}
    }
  file_out.close();
  printf("File %s successfully written\n", filename.c_str());
}

void MultiGrid::Write2DIntFile(string outputfiledir, string filenamebase, string name, Array2DInt* Ckmin)
{
  // This writes a 2D integer .dat output file
  int i, j, index;
  string underscore = "_", slash = "/", filename;
  filename = outputfiledir+slash+filenamebase+underscore+name+".dat";
  ofstream file_out(filename.c_str());
  file_out << setw(16) << "X" << setw(16) << "Y" << setw(16) << "data" << endl;
  for (i=0; i<Ckmin->nx; i++)	        
    {
      for (j=0; j<Ckmin->ny; j++)
	{
	  index = i + j * Ckmin->nx;
	  file_out << setw(16) << Ckmin->x[i] << setw(16) << Ckmin->y[j] << setw(16) << Ckmin->data[index] << endl;
	}
    }
  file_out.close();
  printf("File %s successfully written\n", filename.c_str());
}


void MultiGrid::ReadOutputFile(string outputfiledir, string filenamebase, string name, Array3D* array)
{
  // This reads the data from the HDF files
  string underscore = "_", slash = "/", hdfname, filename;
  double* flipped_data  = new double[array->nx*array->ny*array->nz];
  int i, j, k, index, flipped_index;
  hdfname = filenamebase+underscore+name;
  filename = outputfiledir+slash+hdfname+".hdf5";
  ReadHDF5File3(filename, hdfname, array->nz, array->ny, array->nx, flipped_data);
  // There must be a better way to change from x fast to z fast, but this works.
  for (i=0; i<array->nx; i++)
    {
      for (j=0; j<array->ny; j++)
	{
	  for (k=0; k<array->nz; k++)
	    {
	      index = i + j * array->nx + k * array->nx * array->ny;
	      flipped_index = k + j * array->nz + i * array->nz * array->ny;
	      array->data[index] = flipped_data[flipped_index];
	    }
	}
    }
  delete[] flipped_data;
  return;
}

void MultiGrid::Gradient(Array3D* phi, Array3D** E)
{
  // This calculates the E-field as the gradient of the potential
  int i, j, k, nxy, ind;
  nxy = phi->nx * phi->ny;

  // Ex
  for (j=0; j<phi->ny; j++)
    {
      for (k=0; k<phi->nz; k++)
	{
	  ind = j * phi->nx + k * nxy;
	  E[0]->data[ind] = (phi->data[ind+1] - phi->data[ind]) / phi->dx;
	  ind = 1 + j * phi->nx + k * nxy;
	  E[0]->data[ind] = (phi->data[ind+1] - phi->data[ind-1]) / (2.0 * phi->dx);
	  ind = phi->nx-2 + j * phi->nx + k * nxy;
	  E[0]->data[ind] = (phi->data[ind+1] - phi->data[ind-1]) / (2.0 * phi->dx);
	  ind = phi->nx-1 + j * phi->nx + k * nxy;
	  E[0]->data[ind] = (phi->data[ind] - phi->data[ind-1]) / phi->dx;
	  for (i=2; i<phi->nx-2; i++)
	    {
	      ind = i + j * phi->nx + k * nxy;
	      E[0]->data[ind] = (-phi->data[ind+2] + 8.0 * phi->data[ind+1] - 8.0 * phi->data[ind-1] + phi->data[ind-2]) / (12.0 * phi->dx);
	    }
	}
    }

    // Ey
  for (i=0; i<phi->nx; i++)
    {
      for (k=0; k<phi->nz; k++)
	{
	  ind = i + k * nxy;
	  E[1]->data[ind] = (phi->data[ind+phi->nx] - phi->data[ind]) / phi->dy;
	  ind = i + phi->nx + k * nxy;
	  E[1]->data[ind] = (phi->data[ind+phi->nx] - phi->data[ind-phi->nx]) / (2.0 * phi->dy);
	  ind = i + (phi->ny-2) * phi->nx + k * nxy;
	  E[1]->data[ind] = (phi->data[ind+phi->nx] - phi->data[ind-phi->nx]) / (2.0 * phi->dy);
	  ind = i + (phi->ny-1) * phi->nx + k * nxy;
	  E[1]->data[ind] = (phi->data[ind] - phi->data[ind-phi->nx]) / phi->dy;
	  for (j=2; j<phi->ny-2; j++)
	    {
	      ind = i + j * phi->nx + k * nxy;
	      E[1]->data[ind] = (-phi->data[ind+2*phi->nx] + 8.0 * phi->data[ind+phi->nx] - 8.0 * phi->data[ind-phi->nx] + phi->data[ind-2*phi->nx]) / (12.0 * phi->dy);
	    }
	}
    }

  // Ez
  for (i=0; i<phi->nx; i++)
    {
      for (j=0; j<phi->ny; j++)
	{
	  ind = i + j * phi->nx;
	  E[2]->data[ind] = (phi->data[ind+nxy] - phi->data[ind]) / phi->dzp * phi->dzpdz[0];
	  ind = i + j * phi->nx + nxy;
	  E[2]->data[ind] = (phi->data[ind+nxy] - phi->data[ind-nxy]) / (2.0 * phi->dzp)  * phi->dzpdz[1];
	  ind = i + j * phi->nx + (phi->nz-2) * nxy;
	  E[2]->data[ind] = (phi->data[ind+nxy] - phi->data[ind-nxy]) / (2.0 * phi->dzp) * phi->dzpdz[phi->nz-2];
	  ind = i + j * phi->nx + (phi->nz-1) * nxy;
	  E[2]->data[ind] = (phi->data[ind] - phi->data[ind-nxy]) / phi->dzp * phi->dzpdz[phi->nz-1];
	  for (k=2; k<phi->nz-2; k++)
	    {
	      ind = i + j * phi->nx + k * nxy;
	      E[2]->data[ind] = (-phi->data[ind+2*nxy] + 8.0 * phi->data[ind+nxy] - 8.0 * phi->data[ind-nxy] + phi->data[ind-2*nxy]) / (12.0 * phi->dzp) * phi->dzpdz[k];
	    }
	}
    }
  return;
}

void MultiGrid::Trace(double* point, int bottomsteps, bool savecharge, double bottomcharge, ofstream& file)
{
  // This traces an electron down to the bottom, saving path info if requested
  // Diffusion has now been added. This version recalculates mu at each point.
  // And iterates bottomsteps steps after reaching the bottom.
  // If savecharge is true, it finds and stores the self-consistent charge locations
  // If ElectronMethod=0, it will save bottomcharge in each location for bottomsteps steps.
  // If ElectronMethod=1 or 2, it will only determine the pixel and rely on QFe to determine the charge location.
  // id is a unique integer identifier for each electron that we track.
  // phase encodes the tracking phase and is recorded to the Pts file when LogPixelPaths != 0.
  // 0 - initial position.
  // 1 - endpoint of a step through the bulk.
  // 2 - just reached bottom and settling to equilibrium.
  // 3 - equilibrium motion at the bottom, logging charge.
  // 4 - final position.

  int i, j, k, tracesteps = 0, tracestepsmax = 10000;
  bool ReachedBottom = false;
  double mu, E2, Emag, ve, vth, tau, Tscatt;
  double theta, phiangle, zmin, zmax, zbottom;
  zmax = SensorThickness;
  zmin = E[0]->z[Channelkmin] + 2.0 * FieldOxide; // Roughly the top of the collection region
  zbottom = E[0]->zmz[Channelkmin] + 0.01;
  double*  E_interp = new double[3];
  E2 = 0.0;
  for (i=0; i<3; i++)
    {
      E_interp[i] = E[i]->DataInterpolate3D(point[0],point[1],point[2]);
      E2 += E_interp[i] * E_interp[i];
    }
  Emag = max(0.1, sqrt(E2));
  mu = mu_Si(Emag * MICRON_PER_CM, CCDTemperature); // Mobility
  // Thermal Velocity - DiffMultiplier factor greater than 1 as explained in Green, 1990
  vth = sqrt(8.0 * KBOLTZMANN * CCDTemperature / (ME * pi))  * MICRON_PER_M * DiffMultiplier; 
  vth = vth / sqrt((double)NumDiffSteps);
  tau  = ME / QE * mu * METER_PER_CM * METER_PER_CM; // scattering time

  static int id = 0;
  int phase = 0;
  // Log initial position.
  file << setw(8) << id << setw(8) << tracesteps << setw(3) << phase
       << setw(15) << point[0] << setw(15) << point[1] << setw(15) << point[2] << endl;

  phase = 1;
  while (tracesteps < tracestepsmax)
    {
      tracesteps += 1;
      E2 = 0.0;
      for (i=0; i<3; i++)
	{
	  E_interp[i] = E[i]->DataInterpolate3D(point[0],point[1],point[2]);
	  E2 += E_interp[i] * E_interp[i];
	}
      Emag = max(0.1, sqrt(E2));
      mu = mu_Si(Emag * MICRON_PER_CM, CCDTemperature); // Mobility
      ve = mu * MICRON_PER_CM * MICRON_PER_CM; // Drift Velocity Factor (Velocity / E)
      phiangle = 2.0 * pi * drand48();
      theta = acos(-1.0 + 2.0 * drand48());
      Tscatt = -tau * log(1.0 - drand48()) * (double)NumDiffSteps;
      point[0] += (vth * sin(theta) * cos(phiangle) + E_interp[0] * ve) * Tscatt;
      point[1] += (vth * sin(theta) * sin(phiangle) + E_interp[1] * ve) * Tscatt;
      point[2] += (vth * cos(theta) + E_interp[2] * ve) * Tscatt;

      if (point[2] > zmax)
	{
	  // If the electron hits the top surface, it is reflected off,
	  // but has some probability of recombining, defined TopAbsorptionProb.
	  point[2] = 2.0 * zmax - point[2];
	  if (drand48() > (1.0 - TopAbsorptionProb))
	    {
	      break;
	    }
	}
      if (point[2] < zmin && !ReachedBottom)
	{
	  ReachedBottom = true;
	  tracestepsmax = tracesteps + bottomsteps + EquilibrateSteps;
	  phase = 2;
	  // After reaching bottom, iterate bottomsteps+EquilibrateSteps more steps.
	  // The first EquilibrateSteps steps are to let it settle to an
	  // equilibrium location, then we start logging the charge location
	}
      if (ReachedBottom && tracesteps > tracestepsmax - bottomsteps)
	{
	  //  Start logging location after EquilibrateSteps.
	  phase = (tracesteps < tracestepsmax) ? 3 : 4;
	  if (point[2] <= zbottom && SaturationModel == 1)
	    {
	      break; // Electron recombines and is lost if it reaches the gate oxide
	    }
	  point[2] = max(zbottom, point[2]);
	  i = E[0]->XIndex(point[0]);
	  j = E[0]->YIndex(point[1]);
	  k = E[0]->ZIndex(point[2]);
	  if ((hole[0]->data[i + j * hole[0]->nx + k * hole[0]->nx * hole[0]->ny] > 0.1) && SaturationModel == 1)
	    {
	      break; // Electron recombines if it encounters free holes.
	    }
	  if (savecharge && ElectronMethod != 0)
	    {
	      // Find the pixel the charge is in and add 1 electron to it.
	      int PixX = (int)floor((point[0] - PixelBoundaryLowerLeft[0]) / PixelSizeX);
	      int PixY = (int)floor((point[1] - PixelBoundaryLowerLeft[1]) / PixelSizeY);
	      j = PixX + PixelBoundaryNx * PixY;
	      if (j >= 0 && j < PixelBoundaryNx * PixelBoundaryNy)   CollectedCharge[0][j] += 1;
	      phase = 4;
	      break;
	    }
	  if (i > 0 && i < elec[0]->nx-1 && j > 0 && j < elec[0]->ny-1 && k < elec[0]->nz-1 && savecharge && ElectronMethod == 0)
	    {
	      elec[0]->data[i + j * elec[0]->nx + k * elec[0]->nx * elec[0]->ny] += bottomcharge;// Add bottomcharge to this grid cell
	    }
	}
      if(LogPixelPaths == 1) 
      {
        // Log latest position update.
        file << setw(8) << id << setw(8) << tracesteps << setw(3) << phase
	     << setw(15) << point[0] << setw(15) << point[1] << setw(15) << point[2] << endl;
      }
    point[2] = max(zbottom, point[2]);
    }
  delete[] E_interp;
  
  if(LogPixelPaths == 0) 
    {
      // Log final position
      phase = 4;
      file << setw(8) << id << setw(8) << tracesteps << setw(3) << phase
	   << setw(15) << point[0] << setw(15) << point[1] << setw(15) << point[2] << endl;
    }
  id += 1;
  return;
}


double MultiGrid::GetElectronInitialZ()
{
  // This draws an intitial z value from a pdf given a filter and an SED
  if(FilterIndex >= 0 && FilterIndex < n_band && CalculateZ0 == 1) {
        int cdf_index = (int)floor(n_filter_cdf * drand48());
        return SensorThickness - filter_cdf[FilterIndex * n_filter_cdf + cdf_index];
    }
    else {
        return ElectronZ0Fill;
    }
}

void MultiGrid::TraceSpot(int m)
{
  // This builds up a Gaussian spot with given center (Xoffset, Yoffset) and SigmaX and SigmaY
  double x, y, z, rsq, v1, v2, fac, xcenter, ycenter;
  int n;
  double bottomcharge = 1.0 / (double)BottomSteps;
  double* point = new double[3];
  string underscore = "_", slash = "/", name = "Pts";
  string StepNum = std::to_string(m);
  string filename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+name+".dat");
  ofstream file;
  // If LogPixelPaths > 1, then we only log some pixel paths
  int OldLogPixelPaths = LogPixelPaths;
  if (LogPixelPaths != 0)
    {
      if (m % LogPixelPaths == 0)
	{
	  LogPixelPaths = 1;
	}
      else
	{
	  LogPixelPaths = 0;
	}
    }
  file.open(filename.c_str());
  file.setf(ios::fixed);
  file.setf(ios::showpoint);
  file.setf(ios::left);
  file.precision(4);
  // Write header line.
  file << setw(8) << "id" << setw(8) << "step" << setw(3) << "ph"
       << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;

  xcenter = (PixelBoundaryUpperRight[0] + PixelBoundaryLowerLeft[0]) / 2.0 + Xoffset;
  ycenter = (PixelBoundaryUpperRight[1] + PixelBoundaryLowerLeft[1]) / 2.0 + Yoffset;
  for (n=0; n<NumElec; n++)
    {
      //  Use Box-Muller algorithm to generate two Gaussian random numbers
      rsq = 1000.0;
      while (rsq >= 1.0 || rsq == 0.0)
	{
	  v1 = 2.0 * drand48() - 1.0;
	  v2 = 2.0 * drand48() - 1.0;
	  rsq = v1*v1 + v2 *v2;
	}
      fac = sqrt(-2.0 * log(rsq) / rsq);
      x = xcenter + Sigmax * v1 * fac;
      y = ycenter + Sigmay * v2 * fac;
      point[0] = x;
      point[1] = y;
      z = GetElectronInitialZ();
      point[2] = z;
      Trace(point, BottomSteps, true, bottomcharge, file);
    }
  file.close();
  printf("Finished writing grid file - %s\n",filename.c_str());
  fflush(stdout);
  delete[] point;
  LogPixelPaths = OldLogPixelPaths;
  return;
}


void MultiGrid::TraceFringes(int m)
{
  // This traces a random set of starting electron locations in a fringe pattern
  double x, y, z, boxx, boxy, FringeLength, theta;
  double bottomcharge = 1.0 / (double)BottomSteps;
  int n, ntrials;
  bool Reject;
  boxx = PixelBoundaryUpperRight[0] - PixelBoundaryLowerLeft[0];
  boxy = PixelBoundaryUpperRight[1] - PixelBoundaryLowerLeft[1];
  double* point = new double[3];
  string underscore = "_", slash = "/", name = "Pts";
  string StepNum = std::to_string(m);
  string filename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+name+".dat");
  ofstream file;
  file.open(filename.c_str());
  file.setf(ios::fixed);
  file.setf(ios::showpoint);
  file.setf(ios::left);
  file.precision(4);
  // Write header line.
  file << setw(8) << "id" << setw(8) << "step" << setw(3) << "ph"
       << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;
  // If LogPixelPaths > 1, then we only log some pixel paths
  int OldLogPixelPaths = LogPixelPaths;
  if (LogPixelPaths != 0)
    {
      if (m % LogPixelPaths == 0)
	{
	  LogPixelPaths = 1;
	}
      else
	{
	  LogPixelPaths = 0;
	}
    }
  ntrials = 0;
  for (n=0; n<NumElec; n++)
    {
      if (n%1000==0)
	{
	  if (VerboseLevel > 1) printf("In TraceRegion. Finished %d electrons\n",n);
	}
      Reject=true;
      while (Reject) // Rejection sampling
	{
	  ntrials += 1;
	  x = PixelBoundaryLowerLeft[0] + drand48() * boxx + Xoffset;
	  y = PixelBoundaryLowerLeft[1] + drand48() * boxy + Yoffset;
	  theta = atan2(y, x);
	  FringeLength = sqrt(x*x + y*y) * cos(theta - FringeAngle * pi / 180.0);

	  if (1.0 + cos(FringeLength / FringePeriod * 2.0 * pi) > 2.0 * drand48()) Reject = false;
	}

      z = ElectronZ0Fill;
      point[0] = x;
      point[1] = y;
      point[2] = z;
      Trace(point, BottomSteps, true, bottomcharge, file);
    }
  file.close();
  printf("Finished writing grid file - %s. n = %d, ntrials = %d\n",filename.c_str(),n, ntrials);
  fflush(stdout);
  delete[] point;
  LogPixelPaths = OldLogPixelPaths;  
  return;
}

void MultiGrid::TraceFe55Cloud(int m)
{
  // This generates an Fe55 charge cloud and traces all electrons in the cloud
  // down to the final pixels.  Mutual repulsion of the electrons in the cloud is included.
  // Currently only supports ElectronMethod = 2
  
  int i, j, n, nn, phase, bottomphase = 4, bottomcount, topcount, tracesteps = 0, tracestepsmax = 4000;
  double mu, E2, Emag, r2, sqrtr2, ve=0.0, vth=0.0, tau=0.0, Tscatt=0.0;
  double rsq, rcloud, v1, v2, fac, xcenter, ycenter, zcenter;
  double theta, phiangle, zmin, zmax;
  double SORChargeFactor =  (QE*MICRON_PER_M/(EPSILON_0*EPSILON_SI));
  double ElectronRepulsionFactor = SORChargeFactor / (4.0 * pi) * Fe55ElectronMult;
  double HoleRepulsionFactor = SORChargeFactor / (4.0 * pi) * Fe55HoleMult;
  double MinRadiusSquared = 1.0E-6;
  // This sets the minimum radius for calculating the mutual
  // attraction and repulsion.  this prevents these fields from getting too large
  // Currently set at 0.001 micron
  zmax = SensorThickness - 2.0;
  zmin = E[0]->z[Channelkmin];
  double*  E_interp = new double[3];
  double* r = new double[3]; 
  double PenetrationDepth = 30.0; // Exponential penetration depth in microns
  string underscore = "_", slash = "/", name = "Pts";
  string StepNum = std::to_string(m);
  string filename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+name+".dat");
  ofstream file;
  file.open(filename.c_str());
  file.setf(ios::fixed);
  file.setf(ios::showpoint);
  file.setf(ios::left);
  file.precision(4);
  // Write header line.
  file << setw(8) << "id" << setw(8) << "step" << setw(3) << "ph"
       << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;

  name = "HPts";
  string hfilename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+name+".dat");
  ofstream hfile;
  hfile.open(hfilename.c_str());
  hfile.setf(ios::fixed);
  hfile.setf(ios::showpoint);
  hfile.setf(ios::left);
  hfile.precision(4);
  // Write header line.
  hfile << setw(8) << "id" << setw(8) << "step" << setw(3) << "ph"
       << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;
  
  // First we generate the initial charge cloud position
  // Center the cloud somewhere in the center pixel
  xcenter = (PixelBoundaryUpperRight[0] + PixelBoundaryLowerLeft[0]) / 2.0 ;
  ycenter = (PixelBoundaryUpperRight[1] + PixelBoundaryLowerLeft[1]) / 2.0;
  xcenter += PixelSizeX * (0.5 - drand48());
  ycenter += PixelSizeY * (0.5 - drand48());  
  // Draw zcenter from an exponential distribution. Keep it in the silicon.
  zcenter = SensorThickness + PenetrationDepth * log(1.0 - drand48());
  zcenter = max(zmin, min(zmax, zcenter));

  // Now generate the initial charge cloud
  phase = 0;
  double** Ecloud = new double*[NumElec]; // Electrons
  double** Hcloud = new double*[NumElec]; // Holes
  for (n=0; n<NumElec; n++)
    {
      Ecloud[n] = new double[3];
      Hcloud[n] = new double[3];      
      //  Use Box-Muller algorithm to generate two Gaussian random numbers
      // Only using one
      rsq = 1000.0;
      while (rsq >= 1.0 || rsq == 0.0)
	{
	  v1 = 2.0 * drand48() - 1.0;
	  v2 = 2.0 * drand48() - 1.0;
	  rsq = v1*v1 + v2 *v2;
	}
      fac = sqrt(-2.0 * log(rsq) / rsq);
      // Choose a radius from  Gaussian distribution and add a random angle (electrons)
      rcloud = Fe55CloudRadius * v1 * fac;
      phiangle = 2.0 * pi * drand48();
      theta = acos(-1.0 + 2.0 * drand48());
      Ecloud[n][0] = xcenter + rcloud * sin(theta) * cos(phiangle);
      Ecloud[n][1] = ycenter + rcloud * sin(theta) * sin(phiangle);
      Ecloud[n][2] = zcenter + rcloud * cos(theta);            
      // Choose a radius from  Gaussian distribution and add a random angle (holes)
      rcloud = Fe55CloudRadius * v2 * fac;
      phiangle = 2.0 * pi * drand48();
      theta = acos(-1.0 + 2.0 * drand48());
      Hcloud[n][0] = xcenter + rcloud * sin(theta) * cos(phiangle);
      Hcloud[n][1] = ycenter + rcloud * sin(theta) * sin(phiangle);
      Hcloud[n][2] = zcenter + rcloud * cos(theta);            
      // Log initial position.
      file << setw(8) << n << setw(8) << tracesteps << setw(3) << phase
	   << setw(15) << Ecloud[n][0] << setw(15) << Ecloud[n][1] << setw(15) << Ecloud[n][2] << endl;
      // Log initial position.
      hfile << setw(8) << n << setw(8) << tracesteps << setw(3) << phase
	   << setw(15) << Hcloud[n][0] << setw(15) << Hcloud[n][1] << setw(15) << Hcloud[n][2] << endl;
    }
  // Thermal Velocity - DiffMultiplier factor greater than 1 as explained in Green, 1990
  vth = sqrt(8.0 * KBOLTZMANN * CCDTemperature / (ME * pi))  * MICRON_PER_M * DiffMultiplier; 
  vth = vth / sqrt((double)NumDiffSteps);

  // Now trace the cloud down, adding in the mutual repulsion at each step
  // All electrons trace together
  phase = 1;
  bottomcount = 0;
  topcount = 0;
  while (tracesteps < tracestepsmax && bottomcount < NumElec)
    {
      tracesteps += 1;
      for (n=0; n<NumElec; n++)
	{
	  if (Ecloud[n][2] < zmin) continue; // This electron has finished
	  // Get the background field
	  for (i=0; i<3; i++)
	    {
	      E_interp[i] = E[i]->DataInterpolate3D(Ecloud[n][0],Ecloud[n][1],Ecloud[n][2]);
	    }
	  // Now add in the mutual repulsion, one electron at a time
	  for (nn=0; nn<NumElec; nn++)
	    {
	      if (n == nn || Ecloud[nn][2] < zmin) continue;
	      // Skip yourself, and anything which has reached the bottom
	      r2 = 0.0;
	      for (i=0; i<3; i++)
		{
		  r[i] = Ecloud[n][i] - Ecloud[nn][i];
		  r2 += max(MinRadiusSquared,r[i] * r[i]);
		}
	      sqrtr2 = sqrt(r2);
	      for (i=0; i<3; i++)
		{
		  E_interp[i] += r[i] / (r2 * sqrtr2) * ElectronRepulsionFactor;
		}
	    }
	  // Now add in the mutual attraction, one hole at a time
	  for (nn=0; nn<NumElec; nn++)
	    {
	      if (Hcloud[nn][2] > zmax) continue;
	      // Skip any hole which has reached the top side
	      r2 = 0.0;
	      for (i=0; i<3; i++)
		{
		  r[i] = Hcloud[nn][i] - Ecloud[n][i];
		  r2 += max(MinRadiusSquared,r[i] * r[i]);
		}
	      sqrtr2 = sqrt(r2);
	      for (i=0; i<3; i++)
		{
		  E_interp[i] += r[i] / (r2 * sqrtr2) * HoleRepulsionFactor;
		}
	    }
	  // Now calculate the step as a combination of diffusion and drift
	  // Just like in Trace().
	  E2 = 0.0;
	  for (i=0; i<3; i++)
	    {
	      E2 += E_interp[i] * E_interp[i];
	    }
	  Emag = max(0.1, sqrt(E2));
	  mu = mu_Si(Emag * MICRON_PER_CM, CCDTemperature); // Mobility
	  ve = mu * MICRON_PER_CM * MICRON_PER_CM; // Drift Velocity Factor (Velocity / E)
	  tau  = ME / QE * mu * METER_PER_CM * METER_PER_CM; // scattering time
	  Tscatt = -tau * log(1.0 - drand48()) * (double)NumDiffSteps;
	  phiangle = 2.0 * pi * drand48();
	  theta = acos(-1.0 + 2.0 * drand48());
	  Ecloud[n][0] += (vth * sin(theta) * cos(phiangle) + E_interp[0] * ve) * Tscatt;
	  Ecloud[n][1] += (vth * sin(theta) * sin(phiangle) + E_interp[1] * ve) * Tscatt;
	  Ecloud[n][2] += (vth * cos(theta) + E_interp[2] * ve) * Tscatt;
	  if(LogPixelPaths == 1) 
	    {
	      // Log latest position update.
	      file << setw(8) << n << setw(8) << tracesteps << setw(3) << phase
		   << setw(15) << Ecloud[n][0] << setw(15) << Ecloud[n][1] << setw(15) << Ecloud[n][2] << endl;
	    }
	  if (Ecloud[n][2] < zmin)
	    {
	      bottomcount += 1;
	      // Find the pixel the charge is in and add 1 electron to it.
	      int PixX = (int)floor((Ecloud[n][0] - PixelBoundaryLowerLeft[0]) / PixelSizeX);
	      int PixY = (int)floor((Ecloud[n][1] - PixelBoundaryLowerLeft[1]) / PixelSizeY);
	      j = PixX + PixelBoundaryNx * PixY;
	      if (j >= 0 && j < PixelBoundaryNx * PixelBoundaryNy)   CollectedCharge[0][j] += 1;
	      if (VerboseLevel > 1)
		{
		  printf("In TraceFe55Cloud. Electron %d finished in Pixel (%d, %d)\n",n,PixX,PixY);
		}
		// Log the last position update.
	      file << setw(8) << n << setw(8) << tracesteps << setw(3) << bottomphase
		   << setw(15) << Ecloud[n][0] << setw(15) << Ecloud[n][1] << setw(15) << Ecloud[n][2] << endl;
	    }
	} // ends n
      // Now track the holes.
      for (n=0; n<NumElec; n++)
	{
	  if (Hcloud[n][2] > zmax) continue; // This hole has reached the top
	  // Get the background field
	  for (i=0; i<3; i++)
	    {
	      E_interp[i] = E[i]->DataInterpolate3D(Hcloud[n][0],Hcloud[n][1],Hcloud[n][2]);
	    }
	  // Now add in the mutual attraction, one electron at a time
	  for (nn=0; nn<NumElec; nn++)
	    {
	      if (Ecloud[nn][2] < zmin) continue;
	      // Skip anything which has reached the bottom
	      r2 = 0.0;
	      for (i=0; i<3; i++)
		{
		  r[i] = Hcloud[n][i] - Ecloud[nn][i];
		  r2 += max(MinRadiusSquared,r[i] * r[i]);
		}
	      sqrtr2 = sqrt(r2);
	      for (i=0; i<3; i++)
		{
		  E_interp[i] += r[i] / (r2 * sqrtr2) * HoleRepulsionFactor;
		}
	    }
	  // Now add in the mutual repulsion, one hole at a time
	  for (nn=0; nn<NumElec; nn++)
	    {
	      if (n == nn || Hcloud[nn][2] > zmax) continue;	      
	      // Skip yourself and any hole which has reached the top side
	      r2 = 0.0;
	      for (i=0; i<3; i++)
		{
		  r[i] = Hcloud[nn][i] - Hcloud[n][i];
		  r2 += max(MinRadiusSquared,r[i] * r[i]);
		}
	      sqrtr2 = sqrt(r2);
	      for (i=0; i<3; i++)
		{
		  E_interp[i] += r[i] / (r2 * sqrtr2) * HoleRepulsionFactor;
		}
	    }
	  E2 = 0.0;
	  for (i=0; i<3; i++)
	    {
	      E2 += E_interp[i] * E_interp[i];
	    }
	  // Now calculate the step as a combination of diffusion and drift
	  // Just like in Trace().
	  Emag = max(0.1, sqrt(E2));
	  mu = mu_Si(Emag * MICRON_PER_CM, CCDTemperature) / 3.0; // Mobility
	  // For now, use the hack that the hole mobility is 1/3 the electrons
	  // Hole motion is not critical to device operation
	  // Reverse sign on ve because of opposite charge
	  ve = -mu * MICRON_PER_CM * MICRON_PER_CM; // Drift Velocity Factor (Velocity / E)
	  tau  = ME / QE * mu * METER_PER_CM * METER_PER_CM; // scattering time
	  Tscatt = -tau * log(1.0 - drand48()) * (double)NumDiffSteps;
	  phiangle = 2.0 * pi * drand48();
	  theta = acos(-1.0 + 2.0 * drand48());
	  Hcloud[n][0] += (vth * sin(theta) * cos(phiangle) + E_interp[0] * ve) * Tscatt;
	  Hcloud[n][1] += (vth * sin(theta) * sin(phiangle) + E_interp[1] * ve) * Tscatt;
	  Hcloud[n][2] += (vth * cos(theta) + E_interp[2] * ve) * Tscatt;
	  if(LogPixelPaths == 1) 
	    {
	      // Log latest position update.
	      hfile << setw(8) << n << setw(8) << tracesteps << setw(3) << phase
		   << setw(15) << Hcloud[n][0] << setw(15) << Hcloud[n][1] << setw(15) << Hcloud[n][2] << endl;
	    }
	  if (Hcloud[n][2] > zmax)
	    {
	      topcount += 1;
	      if (VerboseLevel > 1)
		{
		  printf("In TraceFe55Cloud. Hole %d has reached the top after %d steps.\n",n,tracesteps);
		}
	      // Log the last position update.
	      hfile << setw(8) << n << setw(8) << tracesteps << setw(3) << bottomphase
		    << setw(15) << Hcloud[n][0] << setw(15) << Hcloud[n][1] << setw(15) << Hcloud[n][2] << endl;
	    }
	} // ends n
    } // ends tracesteps
  file.close();
  hfile.close();  
  printf("Finished writing grid file - %s\n",filename.c_str());
  printf("Finished writing grid file - %s\n",hfilename.c_str());  
  fflush(stdout);
  
  for (n=0; n<NumElec; n++)
    {
      delete[] Ecloud[n];
      delete[] Hcloud[n];      
    }
  delete[] Ecloud;
  delete[] Hcloud;  
  delete[] E_interp;  
  delete[] r;  
  return;
}

void MultiGrid::TraceList(int m)
{
  // This runs a list of photons from a list
  printf("Running TraceList, n = %d\n", NumElec);
  fflush(stdout);

  double x, y, z, zbottom, xcenter, ycenter, abs_length, path_length;
  zbottom = E[0]->zmz[Channelkmin] + 0.01;
  int n, nlist;
  double bottomcharge = 1.0 / (double)BottomSteps;
  double* point = new double[3];
  string underscore = "_", slash = "/", name = "Pts";
  string StepNum = std::to_string(m);
  string filename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+name+".dat");
  ofstream file;
  // If LogPixelPaths > 1, then we only log some pixel paths
  int OldLogPixelPaths = LogPixelPaths;
  if (LogPixelPaths != 0)
    {
      if (m % LogPixelPaths == 0)
	{
	  LogPixelPaths = 1;
	}
      else
	{
	  LogPixelPaths = 0;
	}
    }
  file.open(filename.c_str());
  file.setf(ios::fixed);
  file.setf(ios::showpoint);
  file.setf(ios::left);
  file.precision(4);
  // Write header line.
  file << setw(8) << "id" << setw(8) << "step" << setw(3) << "ph"
       << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;

  xcenter = (PixelBoundaryUpperRight[0] + PixelBoundaryLowerLeft[0]) / 2.0 + Xoffset;
  ycenter = (PixelBoundaryUpperRight[1] + PixelBoundaryLowerLeft[1]) / 2.0 + Yoffset;
  for (n=0; n<NumElec; n++)
    {
      nlist = m * NumElec + n;
      if (nlist > NumPhotons)
	{
	  printf("Reached end of photon list in step %d.\n",m);
	  return;
	}

      abs_length = pow(10.0,(-4.0 + (PhotonListlambda[nlist] - 500.0) / 250.0)) * 1.0E4; //Approximate formula in micron^-1
      path_length = -abs_length * log(1.0 - drand48());
      x = xcenter + PixelSizeX * PhotonListx[nlist]; // in microns
      y = ycenter + PixelSizeY * PhotonListy[nlist]; // in microns      
      x += PhotonListdxdz[nlist] * path_length;
      y += PhotonListdydz[nlist] * path_length;
      z = SensorThickness - path_length / sqrt(1.0 + PhotonListdxdz[nlist]*PhotonListdxdz[nlist] + PhotonListdydz[nlist]*PhotonListdydz[nlist]); // in microns
      if (x < PixelBoundaryLowerLeft[0] || x > PixelBoundaryUpperRight[0] || y < PixelBoundaryLowerLeft[1] || y > PixelBoundaryUpperRight[1] || z > SensorThickness || z < zbottom) continue;
      if (nlist % 1000 == 0)
	{
	  printf("Nlist = %d, (x,y,z) = (%f,%f,%f)\n",nlist,x,y,z);
	}
      point[0] = x;
      point[1] = y;
      point[2] = z;
      Trace(point, BottomSteps, true, bottomcharge, file);
    }
  file.close();
  printf("Finished writing grid file - %s\n",filename.c_str());
  fflush(stdout);
  delete[] point;
  LogPixelPaths = OldLogPixelPaths;
  return;
}

void MultiGrid::TraceGrid(int m)
{
  // This traces a grid of starting electron locations.
  double x, y, z;
  double* point = new double[3];
  string underscore = "_", slash = "/", name = "Pts";
  string StepNum = std::to_string(m);
  string filename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+name+".dat");
  ofstream file;
  file.open(filename.c_str());
  file.setf(ios::fixed);
  file.setf(ios::showpoint);
  file.setf(ios::left);
  file.precision(4);
  // Write header line.
  file << setw(8) << "id" << setw(8) << "step" << setw(3) << "ph"
       << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;

  // If LogPixelPaths > 1, then we only log some pixel paths
  int OldLogPixelPaths = LogPixelPaths;
  if (LogPixelPaths != 0)
    {
      if (m % LogPixelPaths == 0)
	{
	  LogPixelPaths = 1;
	}
      else
	{
	  LogPixelPaths = 0;
	}
    }
  x = PixelBoundaryLowerLeft[0] + PixelBoundaryStepSize[0] / 2.0;
  while (x < PixelBoundaryUpperRight[0])
    {
      y = PixelBoundaryLowerLeft[1] + PixelBoundaryStepSize[1] / 2.0;
      while (y < PixelBoundaryUpperRight[1])
	{
	  point[0] = x;
	  point[1] = y;
	  z = GetElectronInitialZ();
	  point[2] = z;
	  Trace(point, 100, false, 0.0, file);
	  y += PixelBoundaryStepSize[1];
	}
      x += PixelBoundaryStepSize[0];
    }
  file.close();
  printf("Finished writing grid file - %s\n",filename.c_str());
  fflush(stdout);
  delete[] point;
  LogPixelPaths = OldLogPixelPaths;
  return;
}

void MultiGrid::TraceRegion(int m)
{
  // This traces a random set of starting electron locations within the PixelBoundary.
  double x, y, z, boxx, boxy;
  int n;
  boxx = PixelBoundaryUpperRight[0] - PixelBoundaryLowerLeft[0];
  boxy = PixelBoundaryUpperRight[1] - PixelBoundaryLowerLeft[1];
  double* point = new double[3];
  string underscore = "_", slash = "/", name = "Pts";
  string StepNum = std::to_string(m);
  string filename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+name+".dat");
  ofstream file;
  file.open(filename.c_str());
  file.setf(ios::fixed);
  file.setf(ios::showpoint);
  file.setf(ios::left);
  file.precision(4);
  // Write header line.
  file << setw(8) << "id" << setw(8) << "step" << setw(3) << "ph"
       << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;
  // If LogPixelPaths > 1, then we only log some pixel paths
  int OldLogPixelPaths = LogPixelPaths;
  if (LogPixelPaths != 0)
    {
      if (m % LogPixelPaths == 0)
	{
	  LogPixelPaths = 1;
	}
      else
	{
	  LogPixelPaths = 0;
	}
    }
  // Initialize for PixelBoundaryTestType == 4 mode.
  double bottomcharge = 1.0 / (double)BottomSteps;
  double x_center = 0.5 * (PixelBoundaryLowerLeft[0] + PixelBoundaryUpperRight[0]);
  double y_center = 0.5 * (PixelBoundaryLowerLeft[1] + PixelBoundaryUpperRight[1]);

  for (n=0; n<NumElec; n++)
    {
      if (n%1000==0)
	{
	  if (VerboseLevel > 1) printf("In TraceRegion. Finished %d electrons\n",n);
	}
      if(PixelBoundaryTestType == 4) {
          x = x_center + (drand48() - 0.5) * PixelSizeX;
          y = y_center + (drand48() - 0.5) * PixelSizeY;
          z = GetElectronInitialZ();
      }
      else {
	x = PixelBoundaryLowerLeft[0] + drand48() * boxx;
	y = PixelBoundaryLowerLeft[1] + drand48() * boxy;
	z = ElectronZ0Fill;
      }
      point[0] = x;
      point[1] = y;
      point[2] = z;
      if(PixelBoundaryTestType == 4) {
          // Accumulate charge, the same as TraceSpot().
          Trace(point, BottomSteps, true, bottomcharge, file);
      }
      else {
          // Do not accumulate charge, for backwards compatibility.
          Trace(point, 100, false, 0.0, file);
      }
    }
  file.close();
  printf("Finished writing grid file - %s\n",filename.c_str());
  fflush(stdout);
  delete[] point;
  LogPixelPaths = OldLogPixelPaths;  
  return;
}

void MultiGrid::FindEdge(double* point, double theta, ofstream& file)
{
  // This finds the edge of the pixel through binary search given a starting point and a line angle
  int edgesteps, pixx, pixy, newpixx, newpixy;
  double sinth, costh, x, y, z0, deltar, tolerance;
  sinth = sin(theta);
  costh = cos(theta);
  z0 = point[2];
  x = point[0]; y = point[1];
  pixx = (int)floor((point[0] - PixelBoundaryLowerLeft[0]) / PixelSizeX);
  pixy = (int)floor((point[1] - PixelBoundaryLowerLeft[1]) / PixelSizeY);
  deltar = 0.5 * PixelSizeX;
  tolerance = 0.0001;
  edgesteps = 0;
  while (fabs(deltar) > tolerance)
    {
      edgesteps += 1;
      if (edgesteps >200)
	{
	  printf("Too many steps in edge finding, pixx = %d, pixy = %d, theta = %.3f, %d steps, x = %.3f, y = %.3f\n",pixx, pixy, theta, edgesteps,x,y);
	  fflush(stdout);
	  break;
	}
      x += deltar * costh;
      y += deltar * sinth;
      point[0] = x;
      point[1] = y;
      point[2] = z0;
      Trace(point, 100, false, 0.0, file);
      newpixx = (int)floor((point[0] - PixelBoundaryLowerLeft[0]) / PixelSizeX);
      newpixy = (int)floor((point[1] - PixelBoundaryLowerLeft[1]) / PixelSizeY);
      if (VerboseLevel > 2) printf("Finding edge, newpixx = %d, newpixy = %d, theta = %.3f, %d steps, x = %.3f, y = %.3f\n",newpixx, newpixy, theta, edgesteps,point[0],point[1]);
      if (newpixx != pixx || newpixy != pixy)
	{
	  x -= deltar * costh;
	  y -= deltar * sinth;
	  deltar /= 2.0;
	}
    }
  if (VerboseLevel > 2) printf("Found edge, pixx = %d, pixy = %d, theta = %.3f, %d steps, x = %.3f, y = %.3f\n",pixx, pixy, theta, edgesteps,x,y);
  fflush(stdout);
  point[0] = x;
  point[1] = y;
  point[2] = z0;
  return;
}

void MultiGrid::FindCorner(double* point, double* theta, ofstream& file)
{
  // This finds the corner of the pixel through binary search given a starting point and a line angle
  int edgesteps;
  double theta0, delta_theta, r, rm, r0, rp, deltar, tolerance, x0, y0, z0;
  theta0 = *theta;
  delta_theta = 0.10;
  z0 = point[2];
  x0 = point[0]; y0 = point[1];
  deltar = 1.0;
  tolerance = 0.0001;
  edgesteps = 0;
  while (fabs(deltar) > tolerance)
    {
      edgesteps += 1;
      if (edgesteps >200)
	{
	  printf("Too many steps in corner finding, theta = %.3f, %d steps, x0 = %.3f, y0 = %.3f\n",theta0, edgesteps,x0,y0);
	  fflush(stdout);
	  break;
	}
      point[0] = x0; point[1] = y0; point[2] = z0;
      FindEdge(point, theta0, file);
      r = sqrt((point[0] - x0) * (point[0] - x0) + (point[1] - y0) * (point[1] - y0));
      r0 = r;
      point[0] = x0; point[1] = y0; point[2] = z0;
      FindEdge(point, theta0 + delta_theta, file);
      r = sqrt((point[0] - x0) * (point[0] - x0) + (point[1] - y0) * (point[1] - y0));
      rp = r;
      point[0] = x0; point[1] = y0; point[2] = z0;
      FindEdge(point, theta0 - delta_theta, file);
      r = sqrt((point[0] - x0) * (point[0] - x0) + (point[1] - y0) * (point[1] - y0));
      rm = r;
      if (VerboseLevel > 2) printf("In FindCorner. Step = %d, theta = %.5f, r0 = %.4f, rm = %.4f, rp = %.4f\n",edgesteps,theta0,r0,rm,rp);
      if (rp > r0)
	{
	  theta0 = theta0 + delta_theta;
	  deltar = fabs(r0 - rp);
	}
      else if (rm > r0)
	{
	  theta0 = theta0 - delta_theta;
	  deltar = fabs(r0 - rm);
	}
      else if (rp > rm)
	{
	  delta_theta = delta_theta / 2.0;
	  deltar = fabs(r0 - rp);
	}
      else
	{
	  delta_theta = delta_theta / 2.0;
	  deltar = fabs(r0 - rm);
	}
    }
  if (VerboseLevel > 2) printf("Found corner, theta = %.4f, x = %.4f, y = %.4f, after %d steps\n",theta0, point[0], point[1], edgesteps);
  *theta = theta0;
  return ;
}

void MultiGrid::CalculatePixelAreas(int m)
{
  // This finds the pixel vertices and the areas of the pixel grid.
  int OldLogPixelPaths, k, n, pixx, pixy;
  //Turn off Pixel Path logging
  OldLogPixelPaths = LogPixelPaths;
  LogPixelPaths = 0;
  //Turn off Diffusion
  double OldDiffMultiplier;
  OldDiffMultiplier = DiffMultiplier;
  DiffMultiplier = 0.0;
  double x, y, xb, yb, theta, theta0, theta1, dtheta, area;
  double* point = new double[3];
  string dummy = "dummy";
  string ptsfilename = (dummy);
  ofstream ptsfile; //Not needed, but we need to pass something to the Trace subroutine


  // Now calculate the pixel vertices
  Polygon** polyarray = new Polygon*[PixelBoundaryNx * PixelBoundaryNy];
  x = PixelBoundaryLowerLeft[0] + PixelSizeX / 2.0;
  while (x < PixelBoundaryUpperRight[0])
    {
      pixx = (int)floor((x - PixelBoundaryLowerLeft[0]) / PixelSizeX);
      y = PixelBoundaryLowerLeft[1] + PixelSizeY / 2.0;
      while (y < PixelBoundaryUpperRight[1])
	{
	  pixy = (int)floor((y - PixelBoundaryLowerLeft[1]) / PixelSizeY);
	  polyarray[pixx + PixelBoundaryNx * pixy] = new Polygon(4 * NumVertices + 4);
	  //First, find the four corners
	  for (n=1; n<8; n+=2)
	    {
	      theta = (double)n * pi / 4.0;
	      point[0] = x;
	      point[1] = y;
	      point[2] = ElectronZ0Area;
	      FindCorner(point, &theta, ptsfile);
	      Point* two_d_point = new Point(point[0], point[1], theta);
	      polyarray[pixx + PixelBoundaryNx * pixy]->AddPoint(two_d_point);
	    }
	  // Now find NumVertices points along each edge
	  for (n=0; n<4; n++) // Four corners
	    {
	      for (k=0; k<NumVertices; k++)
		{
		  if (n == 3) // Need this to make angles continuous through zero
		    {
		      theta0 = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[n]->theta - 2.0 * pi;
		      theta1 = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[0]->theta;
		    }
		  else
		    {
		      theta0 = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[n]->theta;
		      theta1 = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[(n+1)]->theta;
		    }
		  dtheta = (theta1 - theta0) / ((double)NumVertices + 1.0);
		  theta  = theta0  + ((double)k + 1.0) * dtheta;
		  point[0] = x;
		  point[1] = y;
		  point[2] = ElectronZ0Area;
		  FindEdge(point, theta, ptsfile);
		  if (VerboseLevel > 2) printf("Found edge, pixx = %d, pixy = %d, theta = %.3f, x = %.3f, y = %.3f\n",pixx, pixy, theta, point[0],point[1]);
		  Point* two_d_point = new Point(point[0], point[1], theta);
		  polyarray[pixx + PixelBoundaryNx * pixy]->AddPoint(two_d_point);
		}
	    }
	  if (VerboseLevel > 1) printf("Finished vertex finding for pixel, pixx = %d, pixy = %d\n",pixx, pixy);
	  if (VerboseLevel > 1) printf("Corners at (%.3f, %.3f),(%.3f, %.3f),(%.3f, %.3f),(%.3f, %.3f)\n",
		 polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[0]->x, polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[0]->y,
		 polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[1]->x, polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[1]->y,
		 polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[2]->x, polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[2]->y,
		 polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[3]->x, polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[3]->y);
	  y += PixelSizeY;
	}
      x += PixelSizeX;
    }

  // Now calculate and print out the pixel areas
  string underscore = "_", slash = "/", vertexname = "Vertices", areaname = "Area";
  string StepNum = std::to_string(m);
  string areafilename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+areaname+".dat");
  ofstream areafile;
  areafile.open(areafilename.c_str());
  areafile.setf(ios::fixed);
  areafile.setf(ios::showpoint);
  areafile.setf(ios::left);
  areafile.precision(4);
  areafile  << setw(15) << "Nx" << setw(15) << "Ny" << setw(15) << "Area" << endl;
  for (pixx=0; pixx<PixelBoundaryNx; pixx++)
    {
      for (pixy=0; pixy<PixelBoundaryNy; pixy++)
	{
	  area = polyarray[pixx + PixelBoundaryNx * pixy]->Area();
	  if (VerboseLevel > 1) printf("Found area, pixx = %d, pixy = %d, area = %.3f\n",pixx, pixy, area);
	  areafile  << setw(15) << pixx << setw(15) << pixy << setw(15) << area << endl;
	}
    }
  areafile.close();
  printf("Finished writing grid file - %s\n",areafilename.c_str());

  // Now print out the pixel vertices
  string vertexfilename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+vertexname+".dat");
  ofstream vertexfile;
  vertexfile.open(vertexfilename.c_str());
  vertexfile.setf(ios::fixed);
  vertexfile.setf(ios::showpoint);
  vertexfile.setf(ios::left);
  vertexfile.precision(4);
  vertexfile  << setw(15) << "X0" << setw(15) << "Y0" << setw(15)<< "Theta" << setw(15) << "X" << setw(15) << "Y" << endl;
  x = PixelBoundaryLowerLeft[0] + PixelSizeX / 2.0;
  while (x < PixelBoundaryUpperRight[0])
    {
      pixx = (int)floor((x - PixelBoundaryLowerLeft[0]) / PixelSizeX);
      y = PixelBoundaryLowerLeft[1] + PixelSizeY / 2.0;
      while (y < PixelBoundaryUpperRight[1])
	{
	  pixy = (int)floor((y - PixelBoundaryLowerLeft[1]) / PixelSizeY);
	  for (n=0; n<polyarray[pixx + PixelBoundaryNx * pixy]->npoints; n++)
	    {
	      xb = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[n]->x;
	      yb = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[n]->y;
	      theta = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[n]->theta;
	      vertexfile  << setw(15) << x << setw(15) << y << setw(15)<< theta << setw(15)<< xb << setw(15) << yb << endl;
	    }
	  y += PixelSizeY;
	}
      x += PixelSizeX;
    }
  vertexfile.close();
  ptsfile.close();
  if (VerboseLevel > 1) printf("Finished writing grid file - %s\n",vertexfilename.c_str());
  // Clean up
  delete[] point;
  for (pixx=0; pixx<PixelBoundaryNx; pixx++)
    {
      for (pixy=0; pixy<PixelBoundaryNy; pixy++)
	{
	  delete polyarray[pixx + PixelBoundaryNx * pixy];
	}
    }
  delete[] polyarray;
  //Put Pixel path logging back where it was
  LogPixelPaths = OldLogPixelPaths;
  //Turn Diffusion back on
  DiffMultiplier = OldDiffMultiplier;
  return;
}

double MultiGrid::mu_Si (double E,double T)
{
  // Calculates the electron mobility given E-field and Temperature
  // Shamelessly copied from phosim
  // Jacobini et al. (1977) equation (9)
  double vm=1.53e9 * pow(T,-0.87); // cm/s
  double Ec = 1.01 * pow(T,1.55); // V/cm
  double beta = 2.57e-2 * pow(T,0.66); // index
  return((vm/Ec)/pow(1 + pow(fabs(E)/Ec,beta),1/beta));
}

void MultiGrid::Set_QFh(Array2D** QFh)
{
  //Set hole quasi-Fermi level in the region.
  int i, j, m, n, index;
  for (n=0; n<nsteps+1; n++)
    {
      for (i=0; i<QFh[n]->nx; i++)
	{
	  for (j=0; j<QFh[n]->ny; j++)
	    {
	      index = i + j * QFh[n]->nx;
	      QFh[n]->data[index] = qfh;
	    }
	}
    }
  // Set QFh in the Fixed Voltage regions
  for (m=0; m<NumberofFixedRegions; m++)
    {
      for (n=0; n<nsteps+1; n++)
	{
	  for (i=0; i<QFh[n]->nx; i++)
	    {
	      for (j=0; j<QFh[n]->ny; j++)
		{
		  index = i + j * QFh[n]->nx;
		  if (QFh[n]->x[i] >= FixedRegionLowerLeft[m][0] && QFh[n]->x[i] <= FixedRegionUpperRight[m][0] && QFh[n]->y[j] >= FixedRegionLowerLeft[m][1] && QFh[n]->y[j] <= FixedRegionUpperRight[m][1])
		    if (FixedRegionQFh[m] > -99.9)
		      {
			QFh[n]->data[index] = FixedRegionQFh[m];
		      }
		}
	    }
	}
    }
  fflush(stdout);
  return;
  }

void MultiGrid::Adjust_QFe(Array2D** QFe, Array3D** elec, Array2DInt** Ckmin)
{
  // Set electron quasi-Fermi level in the region to give the desired number of electrons
  int i, j, k, n, m, q, index, index2, imin, imax, jmin, jmax, PixX, PixY;
  double PixXmin, PixXmax, PixYmin, PixYmax, TotalElectrons = 0.0;
  // First clear the array
  for (n=0; n<nsteps+1; n++)
    {
      for (i=0; i<QFe[n]->nx; i++)
	{
	  for (j=0; j<QFe[n]->ny; j++)
	    {
	      index = i + j * QFe[n]->nx;
	      if (ElectronMethod == 2)
		{
		  // Use this large value to recognize use of ElectronMethod 2.
		  QFe[n]->data[index] = 1.0E6;
		}
	      else
		{
		  QFe[n]->data[index] = 100.0;
		}
	    }
	}
    }
  // Set QFe in the Fixed Voltage regions
  for (m=0; m<NumberofFixedRegions; m++)
    {
      for (n=0; n<nsteps+1; n++)
	{
	  for (i=0; i<QFe[n]->nx; i++)
	    {
	      for (j=0; j<QFe[n]->ny; j++)
		{
		  index = i + j * QFe[n]->nx;
		  if (QFe[n]->x[i] >= FixedRegionLowerLeft[m][0] && QFe[n]->x[i] <= FixedRegionUpperRight[m][0] && QFe[n]->y[j] >= FixedRegionLowerLeft[m][1] && QFe[n]->y[j] <= FixedRegionUpperRight[m][1])
		  QFe[n]->data[index] = FixedRegionQFe[m];
		}
	    }
	}
    }
  if (ElectronMethod != 1)
    {
      //Not using QFe to calculate electrons in pixels.
      return;
    }
  // If ElectronMethod == 1, continue with setting QFe
  // Now count the electrons in the pixel regions.
  for (m=0; m<NumberofPixelRegions; m++)
    {
      for (q=0; q<NumberofFilledWells[m]; q++)
	{
	  // First figure out the pixel coordinates
	  PixX = (int)floor((FilledPixelCoords[m][q][0] - PixelRegionLowerLeft[m][0]) / PixelSizeX);
	  PixY = (int)floor((FilledPixelCoords[m][q][1] - PixelRegionLowerLeft[m][1]) / PixelSizeY);
	  PixXmin = PixelRegionLowerLeft[m][0] + (double)PixX * PixelSizeX;
	  PixYmin = PixelRegionLowerLeft[m][1] + (double)PixY * PixelSizeY;
	  PixXmax = PixXmin + PixelSizeX;
	  PixYmax = PixYmin + PixelSizeY;	      
	  for (n=0; n<nsteps+1; n++)
	    {
	      imin = elec[n]->XIndex(PixXmin);
	      imax = elec[n]->XIndex(PixXmax);
	      jmin = elec[n]->YIndex(PixYmin);
	      jmax = elec[n]->YIndex(PixYmax);
	      // Then count how many electrons currently in the pixel
	      ElectronCount[m][n * NumberofFilledWells[m] + q] = 0.0;
	      for (i=imin; i<imax; i++)
		{
		  for (j=jmin; j<jmax; j++)
		    {
		      index = i + j * elec[n]->nx;
		      for (k=Ckmin[n]->data[index]; k<elec[n]->nz; k++)
			{
			  index2 = index + k * elec[n]->nx * elec[n]->ny;
			  ElectronCount[m][n * NumberofFilledWells[m] + q] += elec[n]->data[index2];
			  TotalElectrons +=  elec[n]->data[index2];
			  if (BuildQFeLookup == 1)
			    {
			      elec[n]->data[index2] = 0.0; // Clear the electrons if building the lookup table
			    }
			}
		    }
		}

	      if (VerboseLevel > 2) printf("In Adjust_QFe. n = %d,m=%d,q=%d,PixX = %.1f, Pixy = %.1f, Desired Electrons = %d, ElecCount = %.4f, QFe = %.4f\n",n,m,q,FilledPixelCoords[m][q][0], FilledPixelCoords[m][q][1], CollectedCharge[m][q], ElectronCount[m][n*NumberofFilledWells[m]+q], PixelQFe[m][n*NumberofFilledWells[m]+q]);
	      if (BuildQFeLookup == 1)
		{
		  PixelQFe[m][n * NumberofFilledWells[m] + q] = QFemin + (QFemax - QFemin) / (double)NQFe * (double)q;
		  QFeLookup[n * NumberofFilledWells[m] + q] = ElectronCount[m][n * NumberofFilledWells[m] + q];
		  if (VerboseLevel > 2) printf("In Adjust_QFe. PixX = %.1f, Pixy = %.1f, ElecCount = %.4f, QFe = %.4f\n",FilledPixelCoords[m][q][0], FilledPixelCoords[m][q][1], ElectronCount[m][n*NumberofFilledWells[m]+q], PixelQFe[m][n*NumberofFilledWells[m]+q]);
		}
	      else
		{
		  PixelQFe[m][n * NumberofFilledWells[m] + q] = ElectronQF(CollectedCharge[m][q], n);
		}
	      if (VerboseLevel > 2) printf("In Adjust_QFe. n=%d,m=%d,q=%d,PixX = %.1f, Pixy = %.1f, Desired Electrons = %d, ElecCount = %.4f, QFe = %.4f\n",n,m,q,FilledPixelCoords[m][q][0], FilledPixelCoords[m][q][1], CollectedCharge[m][q], ElectronCount[m][n * NumberofFilledWells[m] + q], PixelQFe[m][n * NumberofFilledWells[m] + q]);  
	      // Now set pixel QFe to calculated value
	      for (i=imin; i<imax; i++)
		{
		  for (j=jmin; j<jmax; j++)
		    {
		      index = i + j * QFe[n]->nx;
		      QFe[n]->data[index] = PixelQFe[m][n * NumberofFilledWells[m] + q];
		    }
		}
	    }
	}
    }
  printf("Finished AdjustQFe, Total Electrons = %.0f\n",TotalElectrons);
  return;
}

double MultiGrid::ElectronQF(int ne, int n)

// Looks up QFe from the QFeLookup table
// With interpolation.
{
  int i = 0;
  double dne, qf, dqf;
  dqf = (QFemax - QFemin) / ((double)NQFe - 1.0);  
  if (ne < 1)
    {
      return 100.0;
    }
  dne = (double)ne;
  if (ne > QFeLookup[n * NQFe])
    {
      // Ne larger than largest value in lookup table
      // Extrapolate from last two values
      qf = QFemin - dqf / (QFeLookup[n * NQFe] - QFeLookup[n * NQFe + 1]) * (dne - QFeLookup[n * NQFe]);
      if (VerboseLevel > 2) printf("In QFe lookup 0, ne = %d, i = %d, qf = %f\n", ne, i, qf); 
    }
  else
    {
      // Find Ne in the table
      for (i=0; i<NQFe-1; i++)
	{
	  if (ne < QFeLookup[n * NQFe + i] && ne > QFeLookup[n * NQFe + i + 1])
	    {
	      break;
	    }
	}
      if (ne > 1000)  
	{
	  // Ne in table, but larger than 1000
	  // Use linear interpolation
	  qf = QFemin + (double)i * dqf - dqf / (QFeLookup[n * NQFe + i] - QFeLookup[n * NQFe + i + 1]) * (dne - QFeLookup[n * NQFe + i]);
	  if (VerboseLevel > 2) printf("In QFe lookup 1, ne = %d, i = %d, qfi = %f, qfi+1 = %f, qf = %f\n", ne, i, QFeLookup[i], QFeLookup[i+1], qf); 
	}
      else
	{
	  // Ne in table, but less than 1000
	  // Use log interpolation
	  // TODO need to add a test for log(0)
	  qf = QFemin + (double)i * dqf - dqf / (log(QFeLookup[n * NQFe + i]) - log(QFeLookup[n * NQFe + i + 1])) * (log(dne) - log(QFeLookup[n * NQFe + i]));
	  if (VerboseLevel > 2) printf("In QFe lookup 2, ne = %d, i = %d, qfi = %f, qfi+1 = %f, qf = %f\n", ne, i, QFeLookup[i], QFeLookup[i+1], qf); 
	}
    }
  return max(0.0,min(100.0,qf));
}



void MultiGrid::Setkmins(Array3D** rho, Array3D** phi, Array3D** eps, Array2DInt** Ckmin, Array2DInt** Vkmin)
{
  // This sets the kmin values, as well as calculating the epsilon values in the oxide regions
  // Vkmin is the lowest z location where the potential is updated.
  // Below Vkmin are the gates, which are at a fixed, unchanging potential.
  // Ckmin is the lowest z location where charges can accumulate.
  // Between Ckmin and Vkmin is the oxide region.
  int i, j, k, m, n, index, index2, PixX, Chkmin, Vokmin, NLeft=0, NRight=0, NTaper=0, LastCS=0, LastC=0;
  int gox_min_k, fox_min_k;
  double gox_error, fox_error;
  double PixXmin, TaperRatio, gox_eff=0.0, fox_eff=0.0, eps_factor_g[nsteps+1], eps_factor_f[nsteps+1], eps_factor, NGoxCells, NFoxCells;
  // First, set to epsilon to nominal values everywhere
  // The silicon dielectric constant is incorporated elsewhere, so epsilon = 1.0 in the silicon,
  // and epsilon = EPSILON_OX / EPSILON_SI in the oxide region.
  for (n=0; n<nsteps+1; n++)
    {
      for (i=0; i<phi[n]->nx; i++)
	{
	  for (j=0; j<phi[n]->ny; j++)
	    {
	      index = i + j * phi[n]->nx;
	      Ckmin[n]->data[index] = 1;
	      Vkmin[n]->data[index] = 1;	  
	      for (k=0; k<eps[n]->nz; k++)
		{
		  index2 = index + k * eps[n]->nx * eps[n]->ny;
		  eps[n]->data[index2] = 1.0;
		}
	    }
	}
    }
  // Set up the basic geometry
  for (n=0; n<nsteps+1; n++)
    {
      rho[n]->ChannelStopVkmin = 1;
      rho[n]->ChannelVkmin = max(1,rho[n]->ZIndex(FieldOxide / 2.0));
      // Find Ckmin which minimizes the Gox and Fox errors
      // Force thicknesses to be greater than requested
      gox_error = 100.0; fox_error = 100.0; gox_min_k = 2; fox_min_k = 2;
      for (k=2; k<eps[n]->nz; k++)
	{
	  gox_eff = rho[n]->zmz[k] - rho[n]->zpz[rho[n]->ChannelVkmin-1];
	  fox_eff = (rho[n]->zmz[k] - rho[n]->zpz[rho[n]->ChannelStopVkmin-1]);
	  if (k > rho[n]->ChannelVkmin && gox_eff > GateOxide && fabs(gox_eff - GateOxide) / GateOxide < gox_error)
	    {
	      gox_error = fabs(gox_eff - GateOxide) / GateOxide;
	      gox_min_k = k;
	    }
	  if (fox_eff > FieldOxide && fabs(fox_eff - FieldOxide) / FieldOxide < fox_error)
	    {
	      fox_error = fabs(fox_eff - FieldOxide) / FieldOxide;
	      fox_min_k = k;
	    }
	}
      rho[n]->ChannelStopCkmin = fox_min_k;
      rho[n]->ChannelCkmin = gox_min_k;
      gox_eff = rho[n]->zmz[rho[n]->ChannelCkmin] - rho[n]->zpz[rho[n]->ChannelVkmin-1];
      fox_eff = (rho[n]->zmz[rho[n]->ChannelStopCkmin] - rho[n]->zpz[rho[n]->ChannelStopVkmin-1]);
      // These are the effective oxide thicknesses, which differ from what is requested
      // due to finite grid resolution.  The below semi-empirical "fudge-factors" are
      // used to tweak the epsilon values to compensate for the finite grid resolution of the
      // oxide thicknesses.  This speeds convergence considerably
      // and makes sure that the effective capacitance matches what it should be.
      // The first term compensates for any error in the effective thickness, and
      // the second term (small) compensates for the finite width of the surface charge layer.
      NGoxCells = (double)(rho[n]->ChannelCkmin - rho[n]->ChannelVkmin);
      NFoxCells = (double)(rho[n]->ChannelStopCkmin - rho[n]->ChannelStopVkmin);
      eps_factor_g[n] = 1.0 + 3.8 * (gox_eff - GateOxide) / (GateOxide * NGoxCells) + 0.09 * rho[n]->zw[rho[n]->ChannelCkmin] / GateOxide;
      eps_factor_f[n] = 1.0 + 3.8 * (fox_eff - FieldOxide) / (FieldOxide * NFoxCells) + 0.09 * rho[n]->zw[rho[n]->ChannelCkmin] / FieldOxide;
      // Make sure these correction factors don't get too wild
      eps_factor_g[n] = max(0.25, min(100.0, eps_factor_g[n]));
      eps_factor_f[n] = max(0.50, min(2.0, eps_factor_f[n]));
      if (VerboseLevel > 2) printf("In Setkmins. n = %d, gox_eff = %f, fox_eff = %f, NGoxCells = %f, NFoxCells = %f, eps_factor_g = %f, eps_factor_f = %f\n",n,gox_eff, fox_eff, NGoxCells, NFoxCells, eps_factor_g[n], eps_factor_f[n]);
      if (VerboseLevel > 2) printf("In Setkmins. n = %d, CSVkmin = %d, CSCkmin = %d CVkmin = %d, CCkmin = %d\n",n,rho[n]->ChannelStopVkmin,rho[n]->ChannelStopCkmin,rho[n]->ChannelVkmin,rho[n]->ChannelCkmin);
      phi[n]->ChannelStopVkmin = rho[n]->ChannelStopVkmin;
      phi[n]->ChannelStopCkmin = rho[n]->ChannelStopCkmin;
      phi[n]->ChannelVkmin = rho[n]->ChannelVkmin;
      phi[n]->ChannelCkmin = rho[n]->ChannelCkmin;  
      if (VerboseLevel > 2) printf("In Setkmins. n = %d, CSVkmin = %d, CSCkmin = %d CVkmin = %d, CCkmin = %d\n",n,rho[n]->ChannelStopVkmin,rho[n]->ChannelStopCkmin,rho[n]->ChannelVkmin,rho[n]->ChannelCkmin);
    }
  for (m=0; m<NumberofPixelRegions; m++)
    {
      for (n=0; n<nsteps+1; n++)
	{
	  // Count the width in grid of the taper regions
	  NLeft = 0; NRight = 0;
	  for (i=0; i<rho[n]->nx; i++)
	    {
	      if (rho[n]->x[i] >= PixelRegionLowerLeft[m][0] && rho[n]->x[i] <= PixelRegionUpperRight[m][0])
		{
		  PixX = (int)floor((rho[n]->x[i] - PixelRegionLowerLeft[m][0]) / PixelSizeX);
		  PixXmin = PixelRegionLowerLeft[m][0] + (double)PixX * PixelSizeX;
		  if (PixX !=0) continue;
		  if (rho[n]->x[i] <= PixXmin + ChannelStopWidth/2.0 || rho[n]->x[i] >= PixXmin + PixelSizeX - ChannelStopWidth/2.0)
		    {
		      // Channel Stop region
		      continue;
		    }
		  else if (rho[n]->x[i] <= PixXmin + ChannelStopWidth/2.0 + FieldOxideTaper)
		    {
		      // Left side taper
		      if (VerboseLevel > 2) printf("In Setkmins. Setting Tapers. n = %d, i = %d, Left\n",n,i);
		      NLeft += 1;
		    }
		  else if (rho[n]->x[i] >= PixXmin + PixelSizeX - ChannelStopWidth/2.0 - FieldOxideTaper)
		    {
		      // Right side taper
		      if (VerboseLevel > 2) printf("In Setkmins. Setting Tapers. n = %d, i = %d, Right\n",n,i);
		      NRight += 1;
		    }
		  else
		    {
		      // Channel region
		      continue;
		    }
		}
	    }
	  if (NLeft == NRight)  NTaper = NLeft;
	  else NTaper = max(NLeft, NRight);
	  if (VerboseLevel > 2) printf("In Setkmins. n = %d, NTaper = %d \n",n,NTaper);
	  for (i=0; i<rho[n]->nx; i++)
	    {
	      if (rho[n]->x[i] >= PixelRegionLowerLeft[m][0] && rho[n]->x[i] <= PixelRegionUpperRight[m][0])
		{
		  PixX = (int)floor((rho[n]->x[i] - PixelRegionLowerLeft[m][0]) / PixelSizeX);
		  PixXmin = PixelRegionLowerLeft[m][0] + (double)PixX * PixelSizeX;
		  if (rho[n]->x[i] <= PixXmin + ChannelStopWidth/2.0 || rho[n]->x[i] >= PixXmin + PixelSizeX - ChannelStopWidth/2.0)
		    {
		      // Channel Stop region
		      TaperRatio = 1.0;
		      LastCS = i;
		    }
		  else if (rho[n]->x[i] <= PixXmin + ChannelStopWidth/2.0 + FieldOxideTaper)
		    {
		      // Left side taper		      
		      TaperRatio = 1.0 - (double)(i - LastCS) / (double)(NTaper + 1);
		    }
		  else if (rho[n]->x[i] >= PixXmin + PixelSizeX - ChannelStopWidth/2.0 - FieldOxideTaper)
		    {
		      // Right side taper		      		      
		      TaperRatio = (double)(i - LastC) / (double)(NTaper + 1);
		    }
		  else
		    {
		      // Channel region
		      TaperRatio = 0.0;
		      LastC = i;
		    }
		  Chkmin = rho[n]->ChannelCkmin + (int)round((TaperRatio * (double)(rho[n]->ChannelStopCkmin - rho[n]->ChannelCkmin)));
		  Vokmin = rho[n]->ChannelVkmin + (int)round((TaperRatio * (double)(rho[n]->ChannelStopVkmin - rho[n]->ChannelVkmin)));
		  if (VerboseLevel > 2) printf("In Setkmins. n = %d, i = %d, x=%f, TaperRatio = %f, LastC=%d, LastCS=%d, Ckmin = %d, Vkmin = %d\n",n,i,rho[n]->x[i],TaperRatio,LastC, LastCS,Chkmin, Vokmin);
		  // Now fill the epsilon, Ckmin, and Vkmin arrays
		  eps_factor =  eps_factor_g[n] + TaperRatio * (eps_factor_f[n] - eps_factor_g[n]);
		  if (VerboseLevel > 2) printf("In Setkmins. n = %d, i = %d, x = %f, TaperRatio = %f, eps_factor = %f\n",n,i,rho[n]->x[i],TaperRatio,eps_factor);

		  for (j=0; j<eps[n]->ny; j++)
		    {
		      index = i + j * rho[n]->nx;
		      if (rho[n]->y[j] >= PixelRegionLowerLeft[m][1] && rho[n]->y[j] <= PixelRegionUpperRight[m][1])
			{
			  Ckmin[n]->data[index] = Chkmin;
			  Vkmin[n]->data[index] = Vokmin;	  
			  for (k=Vkmin[n]->data[index]; k<Ckmin[n]->data[index]; k++)
			    {
			      index2 = index + k * eps[n]->nx * eps[n]->ny;
			      eps[n]->data[index2] = EPSILON_OX / EPSILON_SI * eps_factor;
			      if (VerboseLevel > 2 && j==48 && (i>41 && i<55)) printf("In Setkmins. n=%d, i=%d, j=%d, k=%d, Ckmin = %d, Vkmin = %d, eps = %f\n",n,i,j,k,Ckmin[n]->data[index], Vkmin[n]->data[index], eps[n]->data[index2]);
			    }
			}
		    }
		}
	    }
	}
    }
  // Now fill Ckmin and Vkmin in the fixed regions
  for (m=0; m<NumberofFixedRegions; m++)
    {
      for (n=0; n<nsteps+1; n++)
	{
	  for (i=0; i<rho[n]->nx; i++)
	    {
	      for (j=0; j<rho[n]->ny; j++)
		{
		  index = i + j * rho[n]->nx;
		  if (rho[n]->x[i] >= FixedRegionLowerLeft[m][0] && rho[n]->x[i] <= FixedRegionUpperRight[m][0] && rho[n]->y[j] >= FixedRegionLowerLeft[m][1] && rho[n]->y[j] <= FixedRegionUpperRight[m][1])
		    {
		      if (FixedRegionOxide[m] == 0) // No oxide
			{
			  Ckmin[n]->data[index] = rho[n]->ChannelCkmin;	  
			  Vkmin[n]->data[index] = rho[n]->ChannelCkmin;	  
			  for (k=0; k<eps[n]->nz; k++)
			    {
			      index2 = index + k * eps[n]->nx * eps[n]->ny;
			      eps[n]->data[index2] = 1.0;
			    }
			}
		      else if (FixedRegionOxide[m] == 1) // Channel region - thin oxide
			{
			  Ckmin[n]->data[index] = rho[n]->ChannelCkmin;	  
			  Vkmin[n]->data[index] = rho[n]->ChannelVkmin;	  
			  for (k=0; k<eps[n]->nz; k++)
			    {
			      index2 = index + k * eps[n]->nx * eps[n]->ny;
			      if (k < Ckmin[n]->data[index]) eps[n]->data[index2] = EPSILON_OX / EPSILON_SI * eps_factor_g[n];
			      else  eps[n]->data[index2] = 1.0;
			    }
			}
		      else if (FixedRegionOxide[m] == 2) // Channel Stop region - thick gate oxide
			{
			  Ckmin[n]->data[index] = rho[n]->ChannelStopCkmin;	  
			  Vkmin[n]->data[index] = rho[n]->ChannelStopVkmin;	  
			  for (k=0; k<eps[n]->nz; k++)
			    {
			      index2 = index + k * eps[n]->nx * eps[n]->ny;
			      if (k < Ckmin[n]->data[index]) eps[n]->data[index2] = EPSILON_OX / EPSILON_SI * eps_factor_f[n];
			      else  eps[n]->data[index2] = 1.0;
			    }
			}
		      
		      if (FixedRegionBCType[m] > 0) // Free BC at z=0
			{
			  Vkmin[n]->data[index] = 0;	  
			}
		    }
		}
	    }
	}
    }
  // Check that there are no negative values
  for (n=0; n<nsteps+1; n++)
    {
      for (i=0; i<rho[n]->nx; i++)
	{
	  for (j=0; j<rho[n]->ny; j++)
	    {
	      index = i + j * rho[n]->nx;
	      if (Ckmin[n]->data[index] < 0 || Vkmin[n]->data[index] < 0)
		{
		  printf("Negative kmin value at n=%d, (i,j)=(%d,%d),Ckmin = %d, Vkmin = %d\n",n,i,j,Ckmin[n]->data[index],Vkmin[n]->data[index]);
		  exit(0);
		}
	    }
	}
    }
  return;
}

void MultiGrid::CountCharges(Array3D** rho, Array3D** elec, Array3D** hole)
{
  // Count total Charge.  Used for debug only
  int i, j, k, n, index;
  double CellVolume;
  double TotalElectrons, TotalHoles, PosCharge, NegCharge, ElectronCharge, HoleCharge;
  int NumCells;
  for (n=0; n<nsteps+1; n++)
    {
      NumCells = 0;
      double RhoChargeFactor =  (QE*MICRON_PER_M/(EPSILON_0*EPSILON_SI)) / (rho[n]->dx * rho[n]->dy);
      TotalElectrons = 0.0; TotalHoles = 0.0; PosCharge = 0.0; NegCharge = 0.0; ElectronCharge = 0.0; HoleCharge = 0.0;
      for (i=0; i<rho[n]->nx; i++)
	{
	  for (j=0; j<rho[n]->ny; j++)
	    {
	      for (k=0; k<rho[n]->nz; k++)
		{
		  CellVolume = rho[n]->dx * rho[n]->dy * rho[n]->zw[k];
		  index = i + j * rho[n]->nx + k * rho[n]->nx * rho[n]->ny;
		  if (rho[n]->data[index] > 0.0)
		    {
		      PosCharge += rho[n]->data[index] * CellVolume;
		    }
		  else
		    {
		      NegCharge += rho[n]->data[index] * CellVolume;
		    }
		    if  (k < elec[n]->nz)
		    {
		      TotalHoles += hole[n]->data[index];
		      TotalElectrons += elec[n]->data[index];		      
		      HoleCharge += hole[n]->data[index] * RhoChargeFactor / rho[n]->zw[k] * CellVolume;
		      ElectronCharge += -elec[n]->data[index] * RhoChargeFactor / rho[n]->zw[k] * CellVolume;
		      if (elec[n]->data[index] > 1.0E-18) NumCells += 1;
		    }
		}
	    }
	}
      if (VerboseLevel > 1) printf("n = %d, PosFixedCharge = %.2f, NegFixedCharge = %.2f, HoleCharge = %.2f, ElectronCharge = %.2f, TotalElectrons = %.1f, TotalHoles = %.1f, NumCells = %d\n",n,PosCharge, NegCharge, HoleCharge, ElectronCharge, TotalElectrons, TotalHoles,NumCells);
    }
  return;
}

void MultiGrid::SetCharge(Array3D* rho, Array2DInt* Ckmin, int i, int j, int region, double TaperRatio)
{
  // This sets the charge in the region at location i,j
  // if region=1 => Channel, region=2 => ChannelStop
  int m, k, index, index2, Profile, kmax;
  index = i + j * rho->nx;
  double Depth, ChargeDepth, ChargeFactor =  (QE*MICRON_PER_M/(EPSILON_0*EPSILON_SI)) / pow(MICRON_PER_CM, 3);
  // ChargeFactor converts doping in cm^-3 into the appropriate units
  double Charge=0.0, Sigma=0.0, Dose=0.0, Peak=0.0, Sum = 0.0;
  double Zmin, Zmax, Total;
  double Qss=0.0;
  
  if (region == 1)  // Channel region
    {
      Profile = ChannelProfile;
      Depth = ChannelDepth;
      Charge = ChannelDoping * MICRON_PER_CM  * ChargeFactor;
    }
  else if (region == 2)  // Channel stop region
    {
      Profile = ChannelStopProfile;
      Depth = ChannelStopDepth;
      Charge = ChannelStopDoping * MICRON_PER_CM  * ChargeFactor;
    }
  else if (region == 3) // Channel stop "dot" region
    {
      Profile = ChannelStopDotProfile;
      Depth = ChannelStopDotDepth;
      Charge = ChannelStopDotDoping * MICRON_PER_CM  * ChargeFactor;
    }
  else
    {
      return;
    }
  
  if (Profile == 0) // Square Well
    {
      kmax = rho->ZIndex(rho->z[Ckmin->data[index]] + Depth);
      ChargeDepth = rho->zpz[kmax] - rho->zmz[Ckmin->data[index]];
      for (k=Ckmin->data[index]; k<kmax+1; k++)
	{
	  index2 = index + k * rho->nx * rho->ny;
	  rho->data[index2] += Charge / ChargeDepth;
	}
    }
  else // N Gaussians
    {
      for (m=0; m<Profile; m++)
	{
	  if (region == 1)  // Channel region  
	    {
	      Sigma = ChannelSigma[m];
	      Peak = ChannelPeak[m];	      
	      Dose = ChannelDose[m];
	      Qss = ChannelSurfaceCharge;
	    }
	  else if (region == 2)  // Channel stop region
	    {
	      Sigma = ChannelStopSigma[m];
	      Peak = ChannelStopPeak[m];	      
	      Dose = ChannelStopDose[m];
	      Qss = ChannelStopSurfaceCharge;
	    }
	  else if (region == 3) // Channel stop "dot" region
	    {
	      Sigma = ChannelStopDotSigma[m];
	      Peak = ChannelStopDotPeak[m];	      
	      Dose = ChannelStopDotDose[m];
	      Qss = ChannelStopDotSurfaceCharge;
	    }
	  if (VerboseLevel > 2 && i == rho->nx/2 && j == rho->ny/2) Sum = 0.0;
	  kmax = rho->ZIndex(rho->z[Ckmin->data[index]] + 4.0 * Sigma + Peak);
	  Charge =  Dose * MICRON_PER_CM  * ChargeFactor * TaperRatio;
	  Zmin = rho->zmz[Ckmin->data[index]];
	  Zmax = rho->zpz[kmax];
	  Total = erf((Zmax - Zmin - Peak) / (sqrt(2.0) * Sigma)) + erf(Peak / (sqrt(2.0) * Sigma));
	  if (m == 0)
	    {
	      // Add in surface charge
	      k = Ckmin->data[index];
	      index2 = index + k * rho->nx * rho->ny;
	      rho->data[index2] += Qss * MICRON_PER_CM * ChargeFactor / rho->zw[k];
	    }
	  for (k=Ckmin->data[index]; k<kmax+1; k++)
	    {
	      index2 = index + k * rho->nx * rho->ny;
	      rho->data[index2] += Charge / Total / rho->zw[k] * (erf((rho->zpz[k] - Zmin - Peak) / (sqrt(2.0) * Sigma)) - erf((rho->zmz[k] - Zmin - Peak) / (sqrt(2.0) * Sigma)));
	      if (VerboseLevel > 2 && i == rho->nx/2 && j == rho->ny/2) Sum += rho->data[index2] * rho->zw[k];
	      if (isnan(rho->data[index2]))
		{
		  printf("Nan encountered in SetCharge at (%d,%d,%d), Nz = %d\n",i,j,k,rho->nz);
		  exit(0);
		}
	    }
	  if (VerboseLevel > 2 && i == rho->nx/2 && j == rho->ny/2) printf("In SetCharge. Nz = %d, m = %d, Total = %f, Sum = %f\n",rho->nz, m, Total, Sum);	      
	}
    }
  return;
}

void MultiGrid::FillElectronWells(Array3D* rho, Array3D* elec, Array2DInt* Ckmin, double ChargeMult)
{
  // This fills a pixel with charge from the ColletedCharge array.
  // It assumes the charge is uniformly distributed in the approximate collection region
  // ChargeMult allows you to adjust the charge put in.
  int i, j, k, n, q, index, index2;
  int CollectedChargeimin, CollectedChargeimax, CollectedChargeXWidth;
  int CollectedChargejmin, CollectedChargejmax, CollectedChargeYWidth;
  int CollectedChargekmin, CollectedChargekmax, CollectedChargeZWidth;

  int PixX, PixY, FilledPixX, FilledPixY;
  double PixXmin, PixYmin, CollectCharge=0.0;
  double TotalElectrons = 0.0;
  // First, zero out the electron count.
  for (n=0; n<NumberofPixelRegions; n++)
    {
      for (i=0; i<rho->nx; i++)
	{
	  if (rho->x[i] < PixelRegionLowerLeft[n][0] || rho->x[i] > PixelRegionUpperRight[n][0])
	    {
	      continue; // If not in PixelRegion, continue
	    }
	  for (j=0; j<rho->ny; j++)
	    {
	      if (rho->y[j] < PixelRegionLowerLeft[n][1] || rho->y[j] > PixelRegionUpperRight[n][1])
		{
		  continue; // If not in PixelRegion, continue
		}
	      index = i + j * elec->nx;	      
	      // Zero out the electrons in the pixel regions
	      for (k=0; k<elec->nz; k++)
		{
		  index2 = index + k * elec->nx * elec->ny;
		  elec->data[index2] =0.0;
		}
	    }
	}
    }
  // Now fill each well with the right amount of charge
  for (n=0; n<NumberofPixelRegions; n++)
    {
      for (q=0; q<NumberofFilledWells[n]; q++)
	{
	  if (CollectedCharge[n][q] < 1) continue;
	  FilledPixX = (int)((FilledPixelCoords[n][q][0] - PixelRegionLowerLeft[n][0]) / PixelSizeX);
	  FilledPixY = (int)((FilledPixelCoords[n][q][1] - PixelRegionLowerLeft[n][1]) / PixelSizeY);
	  for (i=0; i<elec->nx; i++)
	    {
	      if (elec->x[i] < PixelRegionLowerLeft[n][0] || elec->x[i] > PixelRegionUpperRight[n][0])
		{
		  continue; // If not in PixelRegion, continue
		}
	      for (j=0; j<elec->ny; j++)
		{
		  if (elec->y[j] < PixelRegionLowerLeft[n][1] || elec->y[j] > PixelRegionUpperRight[n][1])
		    { 
		      continue; // If not in PixelRegion, continue
		    }
		  index = i + j * elec->nx;
		  PixX = (int)((elec->x[i] - PixelRegionLowerLeft[n][0]) / PixelSizeX);
		  PixY = (int)((elec->y[j] - PixelRegionLowerLeft[n][1]) / PixelSizeY);
		  PixXmin = PixelRegionLowerLeft[n][0] + (double)PixX * PixelSizeX;
		  PixYmin = PixelRegionLowerLeft[n][1] + (double)PixY * PixelSizeY;
		  // First set the charges
		  if (!(FilledPixX == PixX && FilledPixY == PixY))
		    {
		      continue;
		    } 
		  CollectedChargeimin = rho->XIndex(PixXmin + ChannelStopWidth/2.0 + max(FieldOxideTaper, 1.2 * ChannelStopSideDiff));
		  CollectedChargeimax = rho->XIndex(PixXmin + PixelSizeX - ChannelStopWidth/2.0 - max(FieldOxideTaper, 1.2 * ChannelStopSideDiff));
		  CollectedChargeXWidth = CollectedChargeimax - CollectedChargeimin + 1;
		  CollectedChargejmin = rho->YIndex(PixYmin + PixelSizeY * (1.0 - (double)CollectingPhases / (double)NumPhases) / 2.0);
		  CollectedChargejmax = rho->YIndex(PixYmin + PixelSizeY * (1.0 + (double)CollectingPhases / (double)NumPhases) / 2.0);
		  CollectedChargeYWidth = CollectedChargejmax - CollectedChargejmin + 1;
		  CollectedChargekmin = Ckmin->data[index] + 1;
		  CollectedChargekmax = Ckmin->data[index] + 1;
		  CollectedChargeZWidth = CollectedChargekmax - CollectedChargekmin + 1;
		  CollectCharge = (double)CollectedCharge[n][q] * ChargeMult / ((double)CollectedChargeXWidth * (double)CollectedChargeYWidth * (double)CollectedChargeZWidth);
		  
		  if (VerboseLevel >2 && ChargeMult > 0.1)
		    {
		      printf("In FillElectronWells, n = %d, q = %d, i = %d, j = %d\n",n,q,i,j);
		      printf("In FillElectronWells, PixX = %d, PixY = %d, FilledPixX = %d, FilledPixY = %d, Charge = %d\n",PixX, PixY, FilledPixX, FilledPixY, CollectedCharge[n][q]);		  
		      printf("In FillElectronWells, Ckmin = %d, CollectedChargekmin = %d, CollectedChargekmax = %d CollectedChargeZWidth = %d grid cells\n",Ckmin->data[index], CollectedChargekmin, CollectedChargekmax, CollectedChargeZWidth);  
		      printf("In FillElectronWells, CollectedChargeimin = %d, CollectedChargeimax = %d CollectedChargeXWidth = %d grid cells\n",CollectedChargeimin, CollectedChargeimax, CollectedChargeXWidth);		  
		      printf("In FillElectronWells, CollectedChargejmin = %d, CollectedChargejmax = %d CollectedChargeYWidth = %d grid cells\n",CollectedChargejmin, CollectedChargejmax, CollectedChargeYWidth);
		      }
		  if (i >= CollectedChargeimin && i <= CollectedChargeimax && j >= CollectedChargejmin && j <= CollectedChargejmax)
		    {
		      for (k=CollectedChargekmin; k<CollectedChargekmax + 1; k++)
			{
			  index2 = index + k * elec->nx * elec->ny;
			  elec->data[index2] = CollectCharge;
			  TotalElectrons += CollectCharge;
			}
		    }
		}
	    }
	}
    }
  if (VerboseLevel > 1 && ChargeMult > 1.0E-6 )
    {
      printf("Finished Filling Electron Wells, %.1f total electrons placed.\n", TotalElectrons);
    }
      fflush(stdout);
  return;
}
