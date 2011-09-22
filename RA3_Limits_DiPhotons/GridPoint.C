#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

class GridPoint {

 public:

  GridPoint() { Init(); }
  ~GridPoint() {}

  void Init();
  void Print();

  int mGluino;
  int mSquark;
  int mNeutralino;
  int nGen;
  int nAcc;
  int nuisanceModel;
  double accValue;
  double accErrorStat;
  double accErrorSyst;
  double accErrorPDF;     // PDF errors on XS
  double accErrorTot;
  double kfactor;         // not used for 2011 because NLO XS is taken from Prospino directly
  double xsecValue;       // NLO XS from Prospino in pb
  double xsecPDFError;    // PDF error
  double xsecRSErrorNeg;  // renormalization scale error
  double xsecRSErrorPos;  // renormalization scale error
  double limit;           // observed limit
  double explimit;        // expected limit
  double explimit_1L;     // expected limit -1 sigma
  double explimit_1H;     // expected limit +1 sigma
  double explimit_2L;     // expected limit -2 sigma
  double explimit_2H;     // expected limit +2 sigma
};


void GridPoint::Init() {

  mGluino = 0;
  mSquark = 0;
  mNeutralino = 0;
  nGen = 0;
  nAcc = 0;
  nuisanceModel = -1;
  accValue = 0;
  accErrorStat = 0;
  accErrorSyst = 0;
  accErrorPDF = 0;
  accErrorTot = 0;
  kfactor = 0;
  xsecValue = 0;
  xsecPDFError = 0;
  xsecRSErrorNeg = 0;
  xsecRSErrorPos = 0;
  limit = 0;
  explimit = 0;
  explimit_1L = 0;
  explimit_1H = 0;
  explimit_2L = 0;
  explimit_2H = 0;

}


void GridPoint::Print() {

  std::cout << "----------------------------------" << std::endl;
  std::cout << "mGluino : " << mGluino << std::endl;
  std::cout << "mSquark : " << mSquark << std::endl;
  std::cout << "mNeutralino : " << mNeutralino << std::endl;
  std::cout << "nGen : " << nGen << std::endl;
  std::cout << "nAcc : " << nAcc << std::endl;
  std::cout << "nuisanceModel : " << nuisanceModel << std::endl;
  std::cout << "accValue : " << accValue << std::endl;
  std::cout << "accErrorStat : " << accErrorStat << std::endl;
  std::cout << "accErrorSyst : " << accErrorSyst << std::endl;
  std::cout << "accErrorPDF : " << accErrorPDF << std::endl;
  std::cout << "accErrorTot : " << accErrorTot << std::endl;
  std::cout << "kfactor : " << kfactor << std::endl;
  std::cout << "xsecValue : " << xsecValue << std::endl;
  std::cout << "xsecPDFError : " << xsecPDFError << std::endl;
  std::cout << "xsecRSErrorNeg : " << xsecRSErrorNeg << std::endl;
  std::cout << "xsecRSErrorPos : " << xsecRSErrorPos << std::endl;
  std::cout << "limit : " << limit << std::endl;
  std::cout << "explimit : " << explimit << std::endl;
  std::cout << "explimit_1L : " << explimit_1L << std::endl;
  std::cout << "explimit_1H : " << explimit_1H << std::endl;
  std::cout << "explimit_2L : " << explimit_2L << std::endl;
  std::cout << "explimit_2H : " << explimit_2H << std::endl;
}

#pragma link C++ class  GridPoint+;


void CreateGrid(std::vector<GridPoint>& grid, TString datafile, TString metCut, TString jetCut)
{
  grid.clear();
  ifstream fin;
  fin.open(datafile.Data());
  while(1){
    int mG, mS, mChi, nGen, nAcc;
    double dummy;
    //     mG mS mChi nGen Met50_nojet Met50_1jet Met100_nojet Met100_1jet Met70_nojet Met70_1jet
    if(metCut.Contains("met50") && jetCut.Contains("nojet"))
      fin >> mG >> mS >> mChi >> nGen >> nAcc >> dummy >> dummy >> dummy >> dummy >> dummy;
    else if(metCut.Contains("met50") && jetCut.Contains("1jet"))
      fin >> mG >> mS >> mChi >> nGen >> dummy >> nAcc >> dummy >> dummy >> dummy >> dummy;
    else if(metCut.Contains("met100") && jetCut.Contains("nojet"))
      fin >> mG >> mS >> mChi >> nGen >> dummy >> dummy >> nAcc >> dummy >> dummy >> dummy;
    else if(metCut.Contains("met100") && jetCut.Contains("1jet"))
      fin >> mG >> mS >> mChi >> nGen >> dummy >> dummy >> dummy >> nAcc >> dummy >> dummy;
    else if(metCut.Contains("met70") && jetCut.Contains("nojet"))
      fin >> mG >> mS >> mChi >> nGen >> dummy >> dummy >> dummy >> dummy >> nAcc >> dummy;
    else if(metCut.Contains("met70") && jetCut.Contains("1jet"))
      fin >> mG >> mS >> mChi >> nGen >> dummy >> dummy >> dummy >> dummy >> dummy >> nAcc;

    if(!fin.good()) break;
    GridPoint p;
    p.mGluino        = mG;
    p.mSquark        = mS;
    p.mNeutralino    = mChi;
    p.nGen           = nGen;
    p.nAcc           = nAcc;
    grid.push_back(p);
  }
  fin.close();
  std::cout << " Loaded " << grid.size() << " grid points " << std::endl;
}


void GetXSection(std::vector<GridPoint>& grid, TString datafile)
{
  ifstream fin;
  // RS errors are absolute.
  fin.open(datafile.Data());
  while(1){
    int mG, mS, mN;
    double xsec, rsp, rsm, dummy;
    std::string dummyS;
    // printf("%6.0f %6.0f %6.0f LO: %9.3e + %9.3e - %9.3e NLO: %9.3e + %9.3e - %9.3e\n", $ms, $mg, $mchi, $loxsec{$ms}{$mg}{$mchi},$lohi,$lolo,$nloxsec{$ms}{$mg}{$mchi},$nlohi,$nlolo);
    fin >> mS >> mG >> mN >> dummyS >> dummy >> dummyS >> dummy >> dummyS >> dummy >> dummyS >> xsec >> dummyS >> rsp >> dummyS >> rsm;
    if(!fin.good()) break;
    for(unsigned int ig=0; ig<grid.size(); ig++){
      if( (mG ==grid[ig].mGluino) && (mS ==grid[ig].mSquark) ){
	grid[ig].xsecValue = std::abs(xsec);
	grid[ig].xsecRSErrorPos = std::abs(rsp/xsec);  // convert to relative error
	grid[ig].xsecRSErrorNeg = std::abs(rsm/xsec);  // convert to relative error
      }
    }
  }
  fin.close();
}


void GetPDFErrorsXSection(std::vector<GridPoint>& grid, TString datafile)
{
  ifstream fin;
  fin.open(datafile.Data());
  while(1){
    int mG, mS;
    double xs;
    fin >> mG >> mS >> xs;
    if(!fin.good()) break;
    for(unsigned int ig=0; ig<grid.size(); ig++){
      if( (mG == grid[ig].mGluino) && (mS == grid[ig].mSquark) ) {
	// errors are in %. convert to relative error
	grid[ig].xsecPDFError = 0.01*xs;
      }
    } 
  }
  fin.close();
}


void GetPDFErrorsAcceptance(std::vector<GridPoint>& grid, TString datafile)
{
  ifstream fin;
  fin.open(datafile.Data());
  while(1){
    int mG, mS;
    double exs;
    fin >> mG >> mS >> exs;
    if(!fin.good()) break;
    for(unsigned int ig=0; ig<grid.size(); ig++){
      if( (mG == grid[ig].mGluino) && (mS == grid[ig].mSquark) ) {
	// errors are in %. convert to relative error
	grid[ig].accErrorPDF = 0.01*exs;
      }
    } 
  }
  fin.close();
}


void CalculateAcceptanceCorrectedJuly2011(std::vector<GridPoint>& grid)
{

  // Introduce all known corrections to MC acceptance that are flat
  // across the grid

  // value for single photon Data/MC correction
  // from Rachel's study
  // 0.953 +/- 0.033, for >=1 jet
  // 0.945 +/- 0.068, before approval
  double effCorrValue = 0.945;

  // relative error for single photon Data/MC correction
  // 3.9% error effCorrValue is a sum of the following errors (see https://indico.cern.ch/getFile.py/access?contribId=4&resId=0&materialId=slides&confId=149009)
  // 1.4% : stat
  // 0.5% : e/g difference
  // 2.4% : PU
  // 1.8% : signal fit
  // 1.7% : S/B PDF assumption
  // 5% : +1jet difference in Data/MC
  double effCorrError = 0.072;

  // relative Jet Energy Scale error on acceptance is generous 2%
  double JESSyst = 0.02;

  for(unsigned int ig=0; ig<grid.size(); ig++){
    double nG = grid[ig].nGen;    
    double nA = grid[ig].nAcc;
    if(nG>0.5 && nA>-0.5)
      {
	double acc = nA/nG;
	double accEStat = sqrt(nG*acc*(1.0-acc))/nG;
	double corrVal = effCorrValue*effCorrValue;
	acc *= corrVal;
	accEStat *= corrVal;

	double accESyst = 0;
	accESyst += 4.*effCorrError*effCorrError;
	accESyst += JESSyst*JESSyst;
	accESyst += grid[ig].accErrorPDF*grid[ig].accErrorPDF;
	accESyst = sqrt(accESyst);
	accESyst *= acc;

	double accETot = sqrt( accEStat*accEStat + accESyst*accESyst);

	grid[ig].accValue = acc;
	grid[ig].accErrorStat = accEStat;
	grid[ig].accErrorSyst = accESyst;
	grid[ig].accErrorTot = accETot;

      }
  }	
}
