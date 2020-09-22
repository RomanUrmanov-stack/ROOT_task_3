#include <iostream>
#include <math.h>
#include "TF1.h"

using namespace std;

Double_t U = -0.5, parameter = 0;

//this is template for our wave function
Double_t wave_function_template(Double_t *x, Double_t *a)
{   
    //fist task
    return a[0] * TMath::Pi()*exp(-(x[0]*x[0])/(a[0]*a[0]))/sqrt(2);

    //second task
    // return exp(-(abs(x[0]))/(a[0]))/(x[0]/a[0] + a[0]/x[0]);
}

//we'll need to consider two integrands
//<psi|H|psi> = int_(-inf)^inf {psi(x)*H(psi(x))dx} = int_-inf^inf{-h**2/(2m)*psi(x)*d^2(psi(x))/(dx)^2 dx} + int_-10^10 {U * psi**2 dx}
Double_t dot_product_first_integrand(Double_t *x, Double_t *a)
{
    Double_t A = a[0]; //characteristic length
    Double_t h = 1, m = 1; //we'll take these constants as 1
    
    TF1 wave_function("wave_function", wave_function_template, -10.*A, 10.*A, 1); //initialization of wave function
    wave_function.SetParameter(0, A);

    wave_function.SetNormalized(true); //wave function normalization

    return -h*h/(2*m) * wave_function.Eval(x[0]) * wave_function.Derivative2(x[0]); 
}

Double_t dot_product_second_integrand(Double_t *x, Double_t *a)
{
    Double_t A = a[0];
    
    TF1 wave_function("wave_function", wave_function_template, -10.*A, 10.*A, 1);
    wave_function.SetParameter(0, A);

    wave_function.SetNormalized(true);

    return U * wave_function.Eval(x[0]) * wave_function.Eval(x[0]);
}

Double_t dot_product(Double_t *x, Double_t *a)
{
    TF1 first_term("first_term", dot_product_first_integrand, -x[0], x[0], 1);
    first_term.SetParameter(0, x[0]);

    TF1 second_term("second_term", dot_product_second_integrand, -10., 10., 1);
    second_term.SetParameter(0, x[0]);

    return first_term.Integral(-x[0], x[0]) + second_term.Integral(-10, 10);
}

Double_t average_x_squared_integrand(Double_t *x, Double_t *a)
{
    Double_t A = a[0]; 

    TF1 wave_function("wave_function", wave_function_template, -10.*A, 10.*A, 1);
    wave_function.SetParameter(0, parameter);

    return wave_function.Eval(x[0]) * x[0] * x[0] * wave_function.Eval(x[0]);
}

Double_t average_x_squared(Double_t *x, Double_t *a)
{
    TF1 integrand("integrand", average_x_squared_integrand, -x[0], x[0], 1);
    integrand.SetParameter(0, x[0]);

    return integrand.Integral(-x[0], x[0]);
}

Int_t task_3()
{
    TF1* energy_minimum = new TF1("energy_minimum", dot_product, -2, 2, 0);

    TF1* result = new TF1("result", wave_function_template, -10, 10, 1);
    parameter = energy_minimum->GetMinimumX(5, 50);
    result->SetParameter(0, parameter); 

    // result->Draw();    
    // result->SetLineColor(kBlue);
    
    // third task
    TGraph* graph = new TGraph();
    for (int i = 0; i < 10; i++)
    {
        TF1* temp = new TF1("temp", dot_product, -2, 2, 0);
        parameter = temp->GetMinimumX(5, 50);
        
        TF1* average_x_squared_func = new TF1("average_x_squared_func", average_x_squared, -10, 10, 0);

        graph->SetPoint(i, U, average_x_squared_func->Eval(parameter));
        
        //case 1
        //U += 0.1;

        //case 2
        //U -= 0.1;

        delete average_x_squared_func;
        delete temp;
    }

    graph->Draw("AP");
    graph->SetMarkerStyle(4);
    graph->SetMarkerSize(1);
    graph->SetMarkerColor(kBlue);

    return 0;
}