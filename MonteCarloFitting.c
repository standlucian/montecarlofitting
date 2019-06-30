#include "TGraphErrors.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TMath.h"
#include "iostream"

void fitLinearb1()
{
    //Declarations
    Int_t i, number_of_datapoints, int_value;
    double a, b, c;
    Int_t *histogram = new Int_t[100000];
    ifstream input;
    ofstream output;
    TCanvas *myc = new TCanvas("myc", "Fitting nuclear lifetime data");
    myc->SetGrid();
    TRandom *par1 = new TRandom();
    
    //Initialize output histogram
    for (i=0; i<100000; i++)
        histogram[i]=0;
    
    //Declare number of datapoints and their vectors
    number_of_datapoints=6;
    Double_t *x = new Double_t[number_of_datapoints];
    Double_t *y = new Double_t[number_of_datapoints];
    Double_t *e = new Double_t[number_of_datapoints];
    
    //Read data from the file
    input.open("data.txt");
    for (i=0; i<number_of_datapoints; i++) {
        input>>a >> b >> c;
        x[i]=a;
        y[i]=b;
        e[i]=c;
    }
    
    //Declare the formula
    TFormula *bateman1 = new TFormula("bateman1","[0]*[1]*exp(-(x)/[2])+[0]*(exp(-(x)/[3])-exp(-(x)/[2]))*[2]/([3]-[2])+[0]*exp(-(x)/[3])");
    TF1 *b1 = new TF1("b1","bateman1",0,100);
    
    //Declare the graph and draw it
    TGraphErrors *Graph = new TGraphErrors(number_of_datapoints, x, y, 0, e);
    Graph->Draw("a*");
    
    //Make sure output is closed
    if(output.is_open())
        output.close();
    
    //Fit the graph with the defined function a million times
    for(int j=0;j<=1000000;j++)
    {
        //Set the parameters before each fit. Set parameters are varried in the fit, Fixed parameters can't be changed by the fit.
        b1->SetParameter(0,0.001);
        b1->FixParameter(1,par1->Gaus(0.558,0.106));
        b1->SetParameter(2,6.313);
        b1->FixParameter(3,par1->Gaus(12.514,1.388));
        
        //Make the fit
        Graph->Fit("b1","MQ");
        
        //Check if the Chi Square is acceptable. If it is, seed 100 points according to a Gaussian and add them to histogram.
        if(b1->GetChisquare()<3)
            for(int k=0;k<100;k++)
            {
                int_value= (Int_t) (par1->Gaus(b1->GetParameter(2),b1->GetParError(2))*1000);
                if(int_value>0)
                    if(int_value<100000)
                        histogram[int_value]++;
            }
        
        //Every one hundred thousand fits, write an update on the screen and write the histogram to disk.
        if(j%100000==0)
        {
            if(! output.is_open())
                output.open ("output.txt");
            std::cout<<j/100000<<std::endl;
            for(i=0;i<100000;i++)
                output<<histogram[i]<<std::endl;
        }
    }
    
    //Fitting the Gaussian to the results
    //Declare number of points in histogram
    Int_t number_of_histpoints=100000;
    Double_t *z = new Double_t[number_of_histpoints];
    Double_t *t = new Double_t[number_of_histpoints];
    
    //Initialize vectors for the fitting
    for (i=0; i<number_of_histpoints; i++) {
        z[i]=i;
        t[i]=histogram[i];
    }
    
    //Define the formula
    TFormula *gauss = new TFormula("gauss","[0]*exp(-0.5*(x-[1])^2/[2]^2)");
    TF1 *g1 = new TF1("g1","gauss",0,100);
    
    //Draw the function
    TGraph *GraphGauss = new TGraph(number_of_histpoints, z, t);
    GraphGauss->Draw("a*");
    
    //Set initial values for the parameters
    g1->SetParameter(0,10000);
    g1->SetParameter(1,10000);
    g1->SetParameter(2,3000);
    
    //Fit and output center and sigma
    GraphGauss->Fit("g1","MQ");
    printf(" %lf\n", g1->GetParameter(1));
    printf(" %lf\n", g1->GetParameter(2));
}
