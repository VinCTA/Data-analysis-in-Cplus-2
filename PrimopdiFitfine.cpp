/*
c++ -o ReadData1AndFit ReadData1AndFit.cpp `root-config --glibs --cflags`
*/
//ll

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "TApplication.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h" // matrice di covarianza 
#include "TMatrixDSym.h" // matrice di covarianza 

using namespace std;

double myGauss (double* x, double* par)
{
  return par[0] * exp(-0.5*((x[0]-par[1])*(x[0]-par[1])/(par[2]*par[2])));
}

double myParabola (double* x, double* par)
{
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

double sum (double* x, double* par)
{
  double val = myGauss(x, par) + myParabola(x, &par[3]);
  return val;
}

bool readData(char* fileName, vector<double>& dataList, double& xMin, double& xMax)
{
  double num;
  bool first = true;
  std::fstream in;
  in.open(fileName,std::ios::in);
  if (in.good() == false) return false;
  
  while (true)
    {
      in >> num;
      if (in.eof() == true) break;
      dataList.push_back(num);
      if (first == true)
	{
	  xMin   = num;
	  xMax   = num;
	  first = false;
	}
      else
	{
	  if (num < xMin) xMin = num;
	  if (num > xMax) xMax = num;
	}
    }
  
  return true;
}
 
//ha Fatto la funzione e ha scritto un po per legger i dati


//##### chi quadro,numero di gradi di libertà(NDF) ,pvalue  calcolati da noi 
void computeChi2(TH1D* myHisto, TF1* myFun, double& Chi2, double& NDF, double& pvalue)
{
  double result = 0;
  NDF = 0.;
  for (unsigned int i = 1; i <= myHisto->GetNbinsX(); i++)
    {
      if (myHisto->GetBinContent(i) != 0.)
	{
	  result += pow((myHisto->GetBinContent(i) - myFun->Eval(myHisto->GetBinCenter(i))),2.) / myHisto->GetBinContent(i);
	  NDF++; // Tiene conto dei soli bin NON vuoti
	}
    }

  Chi2   = result;
  NDF   -= myFun->GetNpar();
  pvalue = TMath::Prob(Chi2,NDF);
}

int main (int argc, char** argv)
{
  gStyle->SetOptFit(1112); // Istruzione per far comparire i risultati del fit nella legenda dell’istogramma

  if (argc < 2)
    {
      cout << "Insufficient number of parameters: ./ReadData1AndFit fileName.txt" << endl;
      return 1;
    }

  
  // #####################
  // # Reading data file #
  // #####################
  double xMin, xMax;
  vector<double> dataList;
  if (readData(argv[1],dataList,xMin,xMax) == false)
    {
      cout << "Error reading data file: " << argv[1] << endl;
      return 1;
    }

  cout << "The file " << argv[1] << " contains " << dataList.size() << " data" << endl;
  cout << "Minimum: " << xMin << endl;
  cout << "Maximum: " << xMax << endl;
   
   // ho dichiarato un Tapplication 
  
  TApplication* myApp = new TApplication("myApp", NULL, NULL);

  
  // ###########################
  // # Preparing the histogram #
  // ###########################
  int nBins = 100;
  int nPar  = 6; // parametri della funzione somma 
  
  TH1D* h1 = new TH1D("h1","Data distribution", nBins, xMin, xMax);   // definisco l istogramma e lo chiamo h1
  h1->SetFillColor(kAzure+6); // colore isto 
  h1->GetXaxis()->SetTitle("Variable x"); // nome asse x
  h1->GetYaxis()->SetTitle("Counts"); // nome asse y 
  for (unsigned int i = 0; i < dataList.size(); i++) h1->Fill(dataList[i]); // ciclo for che riempe l istogramma 

  TCanvas* myCanv1 = new TCanvas("myCanv1","myCanv1",0,0,700,500); /// disegna istogramma

  
  // ##############################
  // # Preparing the fit function #
  // ##############################
  TF1* myFun = new TF1("myFun", sum, xMin, xMax, nPar); // il secondo parametro deve essere sempre la pdf 
  myFun->SetParameter(0, 10000);        // Gaussian maximum value, range di valore massimo che può assumere la gaussiana
  myFun->SetParameter(1, 5.);           // Gaussian mean value, range di valore che puo assumere la media
  myFun->SetParameter(2, h1->GetRMS()); // Gaussian standard deviation
  myFun->SetParameter(3, 1000.);        // Intercept of the parabola
   // ho inizializzato i parametri 
  myFun->SetParName(0,"Ampl");
  myFun->SetParName(1,"Mean");
  myFun->SetParName(2,"Sigma");

  myFun->SetParName(3,"Zero"); 
  myFun->SetParName(4,"First");
  myFun->SetParName(5,"Second");
// ho dato i nomi ai parametri

  // ###########
  // # Fitting #
  // ###########
  myCanv1->cd(); // ci si sposta nel canvas
  h1->Draw(); // disegna l istp
  h1->Fit("myFun"); // fitta i dati attraverso myfun con l ausilio di sum

  
  // #############
  // # Computing #
  // #############
  // NDF è il numero di gradi di libertà
  cout << "\nReduced Chi-2: " << myFun->GetChisquare()/myFun->GetNDF() << endl; 
  /* ho usato (get chi-square) per calcolre il chi quadro normale e poi l ho diviso 
  per il numero di gradi di libertà,per calcolare quello ridotto,in prati scrivo  a schermo
il chi-2 ridotto  */
  cout << "p-value: " << myFun->GetProb() << endl ;
  for (unsigned int i = 0; i < myFun->GetNpar(); i++)
    cout << "Parameter-" << i << " = " << myFun->GetParameter(i) << " +/- " << myFun->GetParError(i) << endl;

/*Vengono usati i metodi GetParameter e
GetParError della classe TF1 per accedere
ai valori dei parametri ed alle loro incertezze*/

  // Compute number of entries for the signal (i.e. Gaussian)
  TF1* mySignal = new TF1("mySignal", myGauss, xMin, xMax, 3);
  mySignal->SetLineColor(kBlack);
  mySignal->SetParameter(0,myFun->GetParameter(0));// numero associato al parametro e quanto vale il parametro,il parametro è 1 su radice di 2mpi*sigma
  mySignal->SetParameter(1,myFun->GetParameter(1));// numero associato al parametro e quanto vale il parametro,il parametro è mu
  mySignal->SetParameter(2,myFun->GetParameter(2));//numero associato al parametro e quanto vale il parametro,il param è sigma
 // in pratica con le funzioni getParametrer mi calcola lui(Fit) i parametri cioè in questo caso mu etc.
 cout << "\nNumber of entries for the signal (i.e. Gaussian): ";
  cout << nBins / (xMax - xMin) * mySignal->Integral(xMin, xMax);
  cout << endl;
  //dichiaro la funzione di fit mySignal che fitta attraverso la mia Gaussiana (myGauss) e lo stesso per myBackgroun i.e. Parabola 
  
  // Compute number of entries for the background (i.e. Parabola)
  TF1* myBackground = new TF1("myBackground", myParabola, xMin, xMax, 3);
  myBackground->SetLineColor(kBlue);
  myBackground->SetParameter(0,myFun->GetParameter(3));
  myBackground->SetParameter(1,myFun->GetParameter(4));
  myBackground->SetParameter(2,myFun->GetParameter(5));
  cout << "Number of entries for the background (i.e. Parabola): ";
  cout << nBins / (xMax - xMin) * myBackground->Integral(xMin, xMax);
  cout << endl;

/*quindi signal coincide col gaussiana e background con polinomio di 2° */


  double Chi2, NDF, pvalue;
  computeChi2(h1,myFun,Chi2,NDF,pvalue);// sta dicendo calcola il chi 2 dell istogramma h1 che fitta con la funzione myFun
  cout << "\nMy Chi-2: " << Chi2 << endl;
  cout << "My NDF: " << NDF << endl;
  cout << "My reduced Chi-2: " << Chi2/NDF << endl;
  cout << "My p-value: " << pvalue << endl;
// ha scritto i valori a schermo 
  myBackground->Draw("same"); // disegna background sullo stesso (same) canvas

  myCanv1->Modified();
  myCanv1->Update();
// refresh del canvas 

  // ###########################
  // # Print covariance matrix #
  // ###########################
  TFitResultPtr r = h1->Fit("myFun", "S");
  r->Print("V");
  TMatrixDSym cov = r->GetCovarianceMatrix();

myApp->Run(); 

  return 0;
  
}   