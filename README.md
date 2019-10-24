# CreditPricing
Modeling and Calibration of Multiname Credit Product

This repository contains Matlab code that I wrote, with two colleagues, as part of a project for a Financial Engineering course at Politecnico di Milano

<b> Introduction: </b> <br>
The main goal of our project is to price the tranches of an Mortgage Backed Security which has a reference portfolio composed of I=1000 mortgages, considered homogeneous with an average notional of 2MI€, an average recovery of 65% for each mortgage and a default probability of 7% for the period T. We assume for simplicity that mortgages provide a single payment at the end of the interest period T and we neglect the effect due to discount factors. In order to price tranches we have to calibrate the parameters of two different models, double t-Student and t-Student, given the implied correlations for the cumulated tranches:
<table style="width:100%">
  <tr>
    <th>Subordinator (Ku) </th>
    <th>Correlation (Vasicek)</th>
   </tr>
   <tr>
    <th>3% </th>
    <th>22.8%</th>
   </tr>  
  <tr>
    <th>6% </th>
    <th>26.0%</th>
   </tr>  
  <tr>
    <th>9% </th>
    <th>29.1%</th>
   </tr>  
  <tr>
    <th>12% </th>
    <th>32.1%</th>
   </tr>  
  <tr>
    <th>22% </th>
    <th>37.9%</th>
   </tr>      
</table>
<br>
In this project we first calibrate the parameters for the double t-Student model under the assumption of the Large Homogeneous Portfolio (LHP) and then we do the same for the t-Student model. Using these optimal parameters we price the tranches under LHP and HP (Homogeneous Porfolio) assumptions and KL (Kullback-Leibler) approximation and compare the results. Finally we calibrate the parmeters of our models considering the KL approximation.
