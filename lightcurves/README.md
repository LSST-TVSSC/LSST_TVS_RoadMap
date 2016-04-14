Collection of   light curves of different transients.

Formats: 
- asci files with time in days and magnitudes 

FILE NAMES
light curve in multiple bands:   type_name.dat       (eg snIc_94I.dat)
light curve in one band:         type_name_band.dat  (eg snIc_94I_V.dat)

SINGLE BAND LIGHT CURVE
- first line of the fileshould be commented '#' and report time unit (days) and band 
- please report the following information if available
1) source of the file ('# ref= Smith et al 2010')
2) The t0 epoch       ('# jd_expl=   24550059.3')
3) if the light curve is the observed one or a smooth version of the light curve ('# type= smooth') 
example:
#  days   V 
#  ref= Smith et al 2010
#  jd_expl=   24550059.3
#  type= smooth
2.2     14.4
3.1     14.3
4.2     14.1
5.6     13.2
...     ....
___________________________________________________________________

MULTI BAND LIGHT CURVE
- first line of the file should be commented '#' and report time unit (days) and bands 
- please report the following information if available
1) source of the file ('# ref= Smith et al 2010')
2) The t0 epoch       ('# jd_expl=   24550059.3')
3) if the light curve is the observed one or a smooth version of the light curve ('# type= smooth') 
4) report if color correction has been applied  ('#  ebv= 0.4')
example: 
#  days   B     V     R      I 
#  ref= Smith et al 2010
#  jd_expl=   24550059.3
#  type= smooth
#  ebv=  0.4
2.2     14.4        14.3     13.3     15.3    
3.1     14.3        14.4     14.3     14.4    
4.2     14.1        14.2     14.1     14.9  
5.6     13.2        14.3     14.3     14.8  
...     ....

_____________________________________________________________________
NOTES:
1)  smooth light curves are prefered 
2)  sloan magnitudes  are prefered ugriz
