/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Oct 23, 2019

  Standalone cpp Poisson solver

*/

//****************** multigrid.h **************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

using namespace std;

#include <assert.h>
#include <globals.h>
#include <fileio.h>
#include <hdf5write.h>
#include <hdf5read.h>
#include <array3d.h>
#include <array2d.h>
#include <array2dint.h>
#include <polygon.h>

class MultiGrid
{
 public:

  // Grid set-up
  Array3D** phi;      // Phi arrays
  Array3D** rho;      // Rho arrays
  Array3D** E;        // Electric field
  Array3D** elec;      // Number of stored electrons
  Array3D** hole;      // Number of mobile holes
  Array2DInt** BCType;      // BCType - 0->fixed potential; 1->Enormal = 0
  Array2D** QFe;         // Electron Quasi-Fermi level
  Array2D** QFh;         // Hole Quasi-Fermi level  
  Array2DInt** Ckmin;
  Array2DInt** Vkmin;
  MultiGrid(string);
  ~MultiGrid();
  
  double w;		// Successive Over-Relaxation factor
  int ncycle;		// Number of SOR cycles at each resolution
  int iterations;	// Number of VCycles
  int Seed;             // Pseudo-random number generator seed
  double NZExp;		// Non-linear z axis exponent
  int ElectronMethod;             // Electron method (not really used here)  
  
  int ScaleFactor;       // Power of 2 that sets the grid size
  double PixelSizeX;      // Pixel size in microns
  double PixelSizeY;      // Pixel size in microns  
  double SensorThickness; // Thickness of the sensor in microns
  int GridsPerPixelX;     // Number of grids per pixel at ScaleFactor = 1
  int GridsPerPixelY;     // Number of grids per pixel at ScaleFactor = 1  
  double GridSpacingX;    // Grid size in microns
  double GridSpacingY;    // Grid size in microns  
  int Nx;                // Number of grids in x at ScaleFactor = 1
  int Ny;                // Number of grids in y at ScaleFactor = 1
  int Nz;                // Number of grids in z at ScaleFactor = 1
  int Nzelec;            // Number of grids in z in electron arrays at ScaleFactor = 1
  double Xmin;
  double Xmax;
  double Ymin;
  double Ymax;
  double Zmin;
  double Zmax;
  double* SimulationRegionLowerLeft;
  int nsteps;            // Number of steps in Multigrid
  int XBCType;            // Free or periodic BC in X-direction
  int YBCType;            // Free or periodic BC in Y-direction

  // Voltages and Charges
  double qfe;           // Quasi-Fermi level for electrons
  double qfh;           // Quasi-Fermi level for holes
  double Vbb;           // Backside voltage
  double logNi;            // Intrinsic carrier concentration at operating temperature
  double ktq;           // kT/q
  double Vcontact;    // Nominal voltage at contact

  int Contactkmin;              // Bottom of contact region doping
  int Contactkmax;              // Top of contact region doping
  int ContactProfile;           // 0 = Square well, N = N Gaussians

  double ContactWidth;		// Contact region width
  double ContactHeight;		// Contact region height
  double ContactDoping;		// Contact doping
  double ContactDepth;		// Contact depth in microns
  double BackgroundDoping; 	// Background doping
  double TopSurfaceDoping;      // Doping of top surface
  double TopDopingThickness;    // Thickness of top doping layer
  double* ContactDose;		// Contact doping
  double* ContactSigma;		// Contact depth in microns
  double* ContactPeak;           // Depth of peak of contact implant below silicon surface in microns
  double ContactCapacitance;  // Used for calculating deltaV due to stored electrons.
  double BottomOxide;             // Bottom oxide thickness in microns

  // Pixel Regions
  int NumberofPixelRegions;	  	  // 
  double** PixelRegionLowerLeft;	  //
  double** PixelRegionUpperRight;	  //
  int* NumberofContactDeltaVs;		  //
  double*** DeltaVPixelCoords;            // (x,y) coords of pixel center
  double** DeltaV;                 // Voltage deviation in contact region
  int NewNumberofContactDeltaVs;		  //
  double*** NewDeltaVPixelCoords;            // (x,y) coords of pixel center
  double** NewDeltaV;                 // Voltage deviation in contact region
  int** CollectedCharge;		  // Collected charge in e-
  // Constant Voltage Regions
  int NumberofFixedRegions;
  double** FixedRegionLowerLeft;
  double** FixedRegionUpperRight;
  double* FixedRegionVoltage;
  double* FixedRegionQFe;
  double* FixedRegionQFh;    
  int* FixedRegionDoping;
  int* FixedRegionOxide;
  int* FixedRegionBCType;

  // Pixel Boundary Tests
  int PixelBoundaryTestType;
  int LogEField;
  int LogPixelPaths;
  int PixelAreas;
  double* PixelBoundaryLowerLeft;
  double* PixelBoundaryUpperRight;
  double* PixelBoundaryStepSize;
  double Fe55CloudRadius;
  double Fe55ElectronMult;
  double Fe55HoleMult;

  int PixelBoundaryNx;
  int PixelBoundaryNy;
  int NumVertices;
  int CalculateZ0;
  double Lambda;
  double ElectronZ0Fill;
  double ElectronZ0Area;  
  double CCDTemperature;
  double DiffMultiplier;
  double TopAbsorptionProb;
  double RecombinationLifetime; 
  int NumDiffSteps;
  int EquilibrateSteps;
  int BottomSteps;
  int SaturationModel;
  int NumElec;
  int NumSteps;
  double Sigmax;
  double Sigmay;
  double Xoffset;
  double Yoffset;

  // File I/O
  string outputfilebase; // Output filename base
  string outputfiledir; // Output filename directory
  string PhotonList; // Photon list filename
  int NumPhotons;  //Number of photons in list
  double* PhotonListx;
  double* PhotonListy;
  double* PhotonListdxdz;
  double* PhotonListdydz;
  double* PhotonListlambda;  

  int VerboseLevel;
  int SaveData;
  int SaveElec;
  int SaveMultiGrids;
  
  // Continuation info
  int Continuation;
  int LastContinuationStep;  

  
  // Methods
  void ReadConfigurationFile(string);
  void ReadPhotonList(string, string);
  void BuildArrays(Array3D**, Array3D**, Array3D**, Array3D**, Array3D**, Array2DInt**, Array2D**, Array2D**, Array2DInt**, Array2DInt**);
  void SaveGrid();
  void SaveGridMulti();  
  void SetInitialVoltages(Array3D*, Array2DInt*, Array2D*);
  void SetFixedAndMobileCharges(Array3D*, Array3D*, Array3D*, Array2DInt*);
  double SOR_Inner(Array3D*, Array3D*, Array3D*, Array3D*, Array2DInt*, Array2D*, Array2D*, Array2DInt*, Array2DInt*);
  double Error_Inner(Array3D*, Array3D*, Array3D*, Array3D*, Array2DInt*, Array2DInt*);
  void Prolongate(Array3D*, Array3D*, Array2DInt*, Array2DInt*);
  void VCycle_Inner(Array3D**, Array3D**, Array3D**, Array3D**, Array2DInt**, Array2D**, Array2D**, Array2DInt**, Array2DInt**, int, int);
  void WriteOutputFile(string, string, string, Array3D*);
  void ReadOutputFile(string, string, string, Array3D*);
  void Gradient(Array3D*, Array3D**);
  double GetElectronInitialZ();
  void Trace(double*, int, bool, double, ofstream&);
  void TraceSpot(int);
  void TraceList(int);        
  void TraceGrid(int);
  void TraceRegion(int);
  //void TraceFe55Cloud(int);  
  void FindEdge(double*, double, ofstream&);
  void FindCorner(double*, double*, ofstream&);
  void CalculatePixelAreas(int);
  double mu_Si (double, double);
  void Set_QFh(Array2D**);
  //void Set_QFe(Array2D**);
  void UpdateDeltaVs();
  void WriteCollectedCharge(string, string, string);  
  void ReadQFeLookup(string, string, string);
  void Setkmins(Array3D**, Array2DInt**, Array2DInt**);
  void CountCharges(Array3D**, Array3D**, Array3D**);
  void Write3DFile(string, string, string, Array3D*);
  void Write2DFile(string, string, string, Array2D*);
  void Write2DIntFile(string, string, string, Array2DInt*);
  void SetCharge(Array3D*, Array3D*, Array3D*, Array2DInt*, int, int, int);
};
