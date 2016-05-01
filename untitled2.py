# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 11:06:08 2016

@author: Tushar
"""
class SolarThermal:
    
   I=800 #w/m2
   eta_DNI=0.75;eta_cap=0.5;eta_conv=0.3
   P=500E6 #W
   Dop=24; Dsun=6 #hours per day
   Yop=8000 #hours per year
   RoI=20.0 #%.
   phi_land=0.25
   Cmirror=20000.0
   Cland=250000.0
   Cplant=60000.0
   def Celec(self):
      I=self.I #w/m2
      eta_DNI=self.eta_DNI;eta_cap=self.eta_cap;eta_conv=self.eta_conv
      P=self.P #W
      Dop=self.Dop; Dsun=self.Dsun #hours per day
      Yop=self.Yop #hours per year
      RoI=self.RoI #%.
      phi_land=self.phi_land
      Cmirror=self.Cmirror
      Cland=self.Cland
      Cplant=self.Cplant
     #unit Conversion and Standardization
      Dop*=3600.0; Dsun*=3600.0;Yop*=3600.0
      Cland/=10000.0
      Cplant/=1000.0
     #Calculations
      esp=I*eta_DNI*eta_cap*eta_conv
      Edaily=P*Dop; Psun=Edaily/Dsun
      Amirror=Psun/esp; Aland=Amirror/phi_land
      CT=Amirror*Cmirror+Aland*Cland+P*Cplant; Eyearly=P*Yop
      Celec=RoI/100.0*CT/Eyearly;
      Celec_kWh=Celec*3600.0*1000.0
      return Celec_kWh