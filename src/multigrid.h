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
  Array3D** eps;      // Dielectric constant array    
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
  double qfh;           // Quasi-Fermi level in bulk
  double qfh1;           // Quasi-Fermi level near bottom
  int UseDoubleQFh;     // 0 - only 1 qfh, 1 - two values of qfh
  double DoubleQFhZmax; // In PixelRegion, below this value, qfh1 applies
  double Ni;            // Intrinsic carrier concentration at operating temperature
  double ktq;           // kT/q
  double Vbb;		// Back bias
  double Vparallel_lo;	// Parallel Low Voltage
  double Vparallel_hi;	// Parallel High Voltage

  int Channelkmin;              // Bottom of channel region doping
  int Channelkmax;              // Top of channel region doping
  int ChannelStopkmin;          // Bottom of channel stop region doping
  int ChannelStopkmax;          // Top of channel stop region doping
  int ChannelProfile;           // 0 = Square well, N = N Gaussians
  int ChannelStopProfile;       // 0 = Square well, N = N Gaussian
  int ChannelStopDotProfile;    // 0 = Square well, N = N Gaussian

  double ChannelStopDoping;	// Channel Stop doping
  double ChannelStopDepth;     	// Channel stop depth in microns
  double ChannelStopSideDiff; 
  double ChannelDoping;		// Channel doping
  double ChannelDepth;		// Channel depth in microns
  double BackgroundDoping; 	// Background doping
  double* ChannelDose;		// Channel doping
  double* ChannelSigma;		// Channel depth in microns
  double* ChannelPeak;           // Depth of peak of channel implant below silicon surface in microns
  double* ChannelStopDose;	// Channel Stop doping
  double* ChannelStopSigma;     	// Channel stop depth in microns
  double* ChannelStopPeak;       // Depth of peak of channel stop implant below silicon surface in microns
  double ChannelStopWidth;     	// Channel stop width in microns
  double ChannelStopSurfaceCharge;  // Channel stop surface charge in cm^-2
  double ChannelSurfaceCharge;  // Channel surface charge in cm^-2

  double ChannelStopDotCenter;  // Center of channel stop "Dots" above pixel bottom in microns
  double ChannelStopDotHeight;  // Height of channel stop "Dots" in microns
  double ChannelStopDotDoping;	// Square profile doping in cm^-2
  double ChannelStopDotDepth;   // Square profile depth in microns
  double* ChannelStopDotDose;   // Doping in cm^-2
  double* ChannelStopDotPeak;   // Location of peak below silicon surface in microns
  double* ChannelStopDotSigma;  // Sigma in microns
  double ChannelStopDotSurfaceCharge;  // Surface charge density in cm^-2

  double GateOxide;             // Gate oxide thickness in microns
  double FieldOxide;            // Field oxide thickness in microns
  double FieldOxideTaper;       // Field oxide taper width in microns
  int NTaper0;                  // Field oxide taper in grid cells at finest grid
  int NumPhases;             // Number of phases
  int CollectingPhases;      // Number of collecting phases
  double GateGap;            // Gap between gates in microns
  
  // Tree Ring Parameters
  int AddTreeRings;              // 0 - No tree rings; 1 - Add tree rings
  double TreeRingAngle;          // Rotation angle in  degrees
  double TreeRingAmplitude;      // Amplitude in multiples of background charge
  double TreeRingPeriod;         // Period of oscillations in microns

  // Pixel Regions
  int NumberofPixelRegions;	  	  // 
  double** PixelRegionLowerLeft;	  //
  double** PixelRegionUpperRight;	  //
  int* NumberofFilledWells;		  //
  int** CollectedCharge;		  // Collected charge in e-
  double*** FilledPixelCoords;            // (x,y) coords of pixel center
  double** ElectronCount;                 // Number of electrons in each filled pixel
  double** PixelQFe;                      // QFe in each filled pixel

  // QFe Look-up table and electron method
  int ElectronMethod;        // 0 - Leave electrons where they land from tracking
			     // 1 - Set QFe (QFe is always used in Fixed Regions)
			     // 2 - Electron conservation and constant QFe
  int BuildQFeLookup;        // 0 - No QFe lookup; 1 - Build QFeLookup
  int NQFe;                  // Number of entries in QFELookup table             
  double QFemin;             // Minimum QFe in table
  double QFemax;             // Maximum QFe in table
  double* QFeLookup;         // QFeLookup table


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
  double ElectronZ0Fill;
  double ElectronZ0Area;  
  double CCDTemperature;
  double DiffMultiplier;
  double TopAbsorptionProb;
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

  // Fringe Parameters
  double FringeAngle;          // Rotation angle in  degrees
  double FringePeriod;         // Period of oscillations in microns

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

  // Wavelength info
  string FilterBand;    // One of "u", "g", "r", "i", "z", "y".
  int FilterIndex;      // 0=u, 1=g, 2=r, 3=i, 4=z, 5=y
  string FilterFile;    // location of tabulated per-band depth probabilities
  static const int n_band = 6, n_filter_cdf = 5000;
  double filter_cdf[n_band * n_filter_cdf];
  
  // Methods
  void ReadConfigurationFile(string);
  void ReadPhotonList(string, string);
  void BuildArrays(Array3D**, Array3D**, Array3D**, Array3D**, Array3D**, Array3D**, Array2DInt**, Array2D**, Array2D**, Array2DInt**, Array2DInt**);
  void SaveGrid();
  void SaveGridMulti();  
  void SetInitialVoltages(Array3D*, Array3D*, Array2DInt*, Array2DInt*, Array2DInt*);
  void SetFixedCharges(Array3D*, Array2DInt*);
  void FillElectronWells(Array3D*, Array3D*, Array2DInt*, double);  
  double SOR_Inner(Array3D*, Array3D*, Array3D*, Array3D*, Array3D*, Array2DInt*, Array2D*, Array2D*, Array2DInt*, Array2DInt*);
  double Error_Inner(Array3D*, Array3D*, Array3D*, Array3D*, Array3D*, Array2DInt*, Array2DInt*);
  void Prolongate(Array3D*, Array3D*, Array3D*, Array3D*, Array2DInt*, Array2DInt*, Array2DInt*);
  void VCycle_Inner(Array3D**, Array3D**, Array3D**, Array3D**, Array3D**, Array2DInt**, Array2D**, Array2D**, Array2DInt**, Array2DInt**, int, int);
  void WriteOutputFile(string, string, string, Array3D*);
  void ReadOutputFile(string, string, string, Array3D*);
  void Gradient(Array3D*, Array3D**);
  double GetElectronInitialZ();
  void Trace(double*, int, bool, double, ofstream&);
  void TraceSpot(int);
  void TraceList(int);        
  void TraceFringes(int);  
  void TraceGrid(int);
  void TraceRegion(int);
  void TraceFe55Cloud(int);  
  void FindEdge(double*, double, ofstream&);
  void FindCorner(double*, double*, ofstream&);
  void CalculatePixelAreas(int);
  double mu_Si (double, double);
  void Set_QFh(Array2D**);
  void Adjust_QFe(Array2D**, Array3D**, Array2DInt**);
  double ElectronQF(int, int);
  void WriteQFeLookup(string, string, string);
  void WriteCollectedCharge(string, string, string);  
  void ReadQFeLookup(string, string, string);
  void Setkmins(Array3D**, Array3D**, Array3D**, Array2DInt**, Array2DInt**);
  void CountCharges(Array3D**, Array3D**, Array3D**);
  void Write3DFile(string, string, string, Array3D*);
  void Write2DFile(string, string, string, Array2D*);
  void Write2DIntFile(string, string, string, Array2DInt*);
  void SetCharge(Array3D*, Array2DInt*, int, int, int, double);
};
