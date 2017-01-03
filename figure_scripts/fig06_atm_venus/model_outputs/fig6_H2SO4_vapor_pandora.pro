PRO H2SO4_vapor_pandora

;calculates saturation vapor pressure of h2so4 in units of atm 

;read in a Venus t-p profile
;readcol, '/astro/users/giada/SMART/profile/Venus/Nov08/atm/venus2_30.atm', p, t, alt

readcol, '/astro/users/giada/pandora/Venus_temp/10.0barCO2_2016-08-14.atm', p, T, alt, co2, h2o, h2so4
p = p * 9.8e-6 ;pascals -> atm

T0 = 340. ;K
Tc = 905. ;K 
R = 8.314e7 ;erg/K/mol
p0 = -(10156./t0)+16.259 ;ln(atm)
W = 85. ; assuming 85% H2SO4 by weight -- btw, the results make no sense if we assume this is supposed to be written as a decimal 
H = 4.184 * 1.e7 * (23624.8 - (1.14208e8)/(4798.69 + (W - 105.315)^2.)) ;ergs / mol
ppm = 4e-6 ;ppm h2so4 in atmosphere (won't actually be constant through atmosphere, but all I really care about is what it is in the lower atmosphere up to the base of the cloud deck to determine where the cloud base would form.)

;ln partial pressure of h2so4 [atm]
h2so4 = p0 + 10156.*(-1./T + 1./T0 + 0.38/(Tc - T0) * ( 1. + alog(T0/T) - T0/T)) - H/(R*T)

;in mmHg to compare to plot in Kulmala 1900 paper for validation -->
;works when we assume 100% H2SO4 
h2so4_hg = exp(h2so4) * 760.

;for plot boundaries: 
bound = where(exp(h2so4) gt 1e-13 and exp(h2so4) lt 1e-4)

pp = p*ppm
tb = t[bound]
eh = exp(h2so4[bound])

!Y.MARGIN=[4,2]
;plot
;set_plot, 'ps'
;device, file='~/pandora/Venus_temp/Venus_h2so4_vapor.eps', /color, xsize='5', ysize='3.5', /encapsulated, /inches
;plot, t, pp, color=0, /ylog, yrange=[1e-4, 1e-13], xtitle='Temperature [K]', ytitle='saturation vapor pressure [atm]', thick=5 , title='10 bar Venus-like Atmo;sphere'
;oplot, t[bound], exp(h2so4[bound]), color=120, thick=5
;;oplot, tclima, pclima*ppm, color=230, thick=5
;legend, ['T profile', 'H!D2!NSO!D4!N SVP'], color=[0, 120], thick=[5,5], textcolor=[0,0], linestyle=[0,0], box=0, position=[270, 1e-13], charsize=0.8
;device, /close
;set_plot, 'x'

!p.multi=[0,2,1]

set_plot, 'ps'
device, file='~/pandora/Venus_temp/Venus_h2so4_vapor.eps', /color, xsize='8', ysize='3.5', /encapsulated, /inches
plot, t, pp/ppm, color=0, /ylog, yrange=[1e2, 1e-8], xtitle='Temperature [K]', ytitle='pressure [atm]', thick=5 , title='10 bar Venus-like Atmosphere'
oplot, t[bound], exp(h2so4[bound])/ppm, color=120, thick=5
;oplot, tclima, pclima*ppm, color=230, thick=5
legend, ['T profile', 'H!D2!NSO!D4!N SVP'], color=[0, 120], thick=[5,5], textcolor=[0,0], linestyle=[0,0], box=0, position=[280, 1e-7], charsize=0.8
;device, /close
;set_plot, 'x'


;plot, t, pp/ppm, color=0, yrange=[10, 8], xtitle='Temperature [K]', ytitle='pressure [atm]', thick=5 , title='10 bar Venus-like Atmosphere'
;oplot, t[bound], exp(h2so4[bound])/ppm, color=120, thick=5




;-----90 bar--------
readcol, '/astro/users/giada/pandora/Venus_temp/90barCO2.atm', p, T, alt, co2, h2o, h2so4
p = p * 9.8e-6 ;pascals -> atm

T0 = 340. ;K
Tc = 905. ;K 
R = 8.314e7 ;erg/K/mol
p0 = -(10156./t0)+16.259 ;ln(atm)
W = 85. ; assuming 85% H2SO4 by weight -- btw, the results make no sense if we assume this is supposed to be written as a decimal 
H = 4.184 * 1.e7 * (23624.8 - (1.14208e8)/(4798.69 + (W - 105.315)^2.)) ;ergs / mol
ppm = 4e-6 ;ppm h2so4 in atmosphere (won't actually be constant through atmosphere, but all I really care about is what it is in the lower atmosphere up to the base of the cloud deck to determine where the cloud base would form.)

;ln partial pressure of h2so4 [atm]
h2so4 = p0 + 10156.*(-1./T + 1./T0 + 0.38/(Tc - T0) * ( 1. + alog(T0/T) - T0/T)) - H/(R*T)

;in mmHg to compare to plot in Kulmala 1900 paper for validation -->
;works when we assume 100% H2SO4 
h2so4_hg = exp(h2so4) * 760.

;for plot boundaries: 
bound = where(exp(h2so4) gt 1e-13 and exp(h2so4) lt 1e-2)

pp = p*ppm
tb = t[bound]
eh = exp(h2so4[bound])


;plot
;set_plot, 'ps'
;device, file='~/pandora/Venus_temp/Venus_h2so4_vapor90.eps', /color, xsize='5', ysize='3.5', /encapsulated, /inches
;plot, t, pp, color=0, /ylog, yrange=[1e-2, 1e-13], xtitle='Temperature [K]', ytitle='saturation vapor pressure [atm]', thick=5 , title='90 bar Venus-like Atmo;sphere'
;oplot, t[bound], exp(h2so4[bound]), color=120, thick=5
;;oplot, tclima, pclima*ppm, color=230, thick=5
;legend, ['T profile', 'H!D2!NSO!D4!N SVP'], color=[0, 120], thick=[5,5], textcolor=[0,0], linestyle=[0,0], box=0, position=[270, 1e-13], charsize=0.8
;device, /close
;set_plot, 'x'


;set_plot, 'ps'
;device, file='~/pandora/Venus_temp/Venus_h2so4_vapor_real_p90.eps', /color, xsize='5', ysize='3.5', /encapsulated, /inches
plot, t, pp/ppm, color=0, /ylog, yrange=[1e3, 1e-8], xtitle='Temperature [K]', ytitle='pressure [atm]', thick=5 , title='90 bar Venus-like Atmosphere'
oplot, t[bound], exp(h2so4[bound])/ppm, color=120, thick=5
;oplot, tclima, pclima*ppm, color=230, thick=5
legend, ['T profile', 'H!D2!NSO!D4!N SVP'], color=[0, 120], thick=[5,5], textcolor=[0,0], linestyle=[0,0], box=0, position=[300, 1e-7], charsize=0.8
device, /close
set_plot, 'x'


plot, t, pp/ppm, color=0, yrange=[20, 8], xtitle='Temperature [K]', ytitle='pressure [atm]', thick=5 , title='10 bar Venus-like Atmosphere'
oplot, t[bound], exp(h2so4[bound])/ppm, color=120, thick=5


stop

END
