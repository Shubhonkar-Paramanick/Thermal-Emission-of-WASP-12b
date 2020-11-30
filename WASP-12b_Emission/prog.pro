; SDSS Z band

min_wave_zband = (0.7730)*10000
max_wave_zband = (1.1230)*10000

; 2MASS J band
min_wave_Jband = (1.06200)*10000
max_wave_Jband = (1.45000)*10000

; 2MASS H band
min_wave_Hband = (1.28900)*10000
max_wave_Hband = (1.91400)*10000

; 2MASS K band
min_wave_Kband = (1.90000)*10000
max_wave_Kband = (2.39900)*10000

; SPITZER 3.6 band
min_wave_3_6 = (3.08106)*10000
max_wave_3_6 = (4.01038)*10000

;SPITZER 4.5 band
min_wave_4_5 = (3.72249)*10000
max_wave_4_5 = (5.22198)*10000


;SPITZER 5.8 band
min_wave_5_8 = (4.74421)*10000
max_wave_5_8 = (6.62251)*10000

; SPITZER 8.0 band
min_wave_8_0 = (6.15115)*10000
max_wave_8_0 = (10.49680)*10000

readcol,'/Volumes/ARVIND_4TB/hot_jupiters/IRR_model/wasp_19b/cloud/models/phot_wasp19_b.txt',w_phot,f_phot,errors,f='f,f',skipline=1
readcol,'/Volumes/ARVIND_4TB/hot_jupiters/IRR_model/wasp_19b/cloud/models/bean_2013.txt',w_phot_1,f_phot_1,errors_1,f='f,f'
readcol,'/Volumes/ARVIND_4TB/hot_jupiters/IRR_model/wasp_19b/cloud/models/others.txt',w_phot_2,f_phot_2,errors_2,f='f,f' 
                                     
;*****************************************************************
; At different Redisitribution with solar metallicity
;model_file = '/Volumes/ARVIND_4TB/hot_jupiters/IRR_model/wasp_19b/cloud/spectra/with_Hminus/list_planet.txt'   ; list_planet.txt
model_file = '/Volumes/ARVIND_4TB/hot_jupiters/IRR_model/wasp_19b/cloud/wasp_19b_paper/spectrum/with_all/list.txt'   ;
n_model = 3
;n_model = file_lines (model_file)
openr,4, model_file
name_synth = strarr (n_model)
readf,4,name_synth
close, 4
free_lun,4

;*******************************
title_name = strmid(name_synth,57,16)
print, 'title_name =',title_name
;*******************************  
color=['red','blue','green','brown']
linestyle=[0,0,0,0]

;cgPS_Open, '/Volumes/ARVIND_4TB/hot_jupiters/IRR_model/wasp_19b/plots/paper/secondary_wasp19b.eps'

cgplot,[0], [0],color='black',yr=[0,0.010],xr=[0.8,10], $
ytitle='Planet-to-Star-Ratio',xtitle = 'Wavelength ($\mu$m)', $
/xlog,ystyle=1, charsize=1.5,xtickv=1,/window,/NODATA


        for n = 0, n_model-1 do begin
                       
        ;readcol,'/Volumes/ARVIND_4TB/hot_jupiters/IRR_model/wasp_19b/cloud/spectra/with_Hminus/' + name_synth (n),wl_1,fl_1,format='(f,f)'
        readcol,'/Volumes/ARVIND_4TB/hot_jupiters/IRR_model/wasp_19b/cloud/wasp_19b_paper/spectrum/with_all/' + name_synth (n),wl_1,fl_1,format='(f,f)'
        wave_trim_start = max(where(wl_1 eq 5000))
        wave_trim_end = (where(wl_1 ge 160000))[0]
                        ;wl_planet=wl_1[wave_trim_start:wave_trim_end]
                        wl_planet = wl_1

                        ;fla = fl_1[wave_trim_start:wave_trim_end]
                        fla = fl_1
                        fl_planet = 10d^(fla - 8d)

            fl_planet_1 = smooth(fl_planet, 1000)
   
            radius_planet = 9.59e9   ; original value from job script 9.69e9


                        rows1 = n_elements(fl_planet)

;stop
;*********** For Planet *******************
; Integrating spectrum over
; different photometric bands
;******************************************

;*** SDSS Z band

wave_planet_zband = where (wl_planet ge min_wave_zband and wl_planet le max_wave_zband)
Planet_wave_zband = wl_planet (wave_planet_zband)
flux_planet_zband = fl_planet_1 (wave_planet_zband)

total_flux_planet_zband = total(flux_planet_zband)
flux_Zband = min(abs(Planet_wave_zband-9124.0),ind)
wl_Zband_fix = Planet_wave_zband[ind]

;*** 2 MASS J band

wave_planet_Jband = where (wl_planet ge min_wave_Jband and wl_planet le max_wave_Jband)
Planet_wave_Jband = wl_planet (wave_planet_Jband)
flux_planet_Jband = fl_planet_1 (wave_planet_Jband)

total_flux_planet_Jband = total(flux_planet_Jband)
flux_Jband = min(abs(Planet_wave_Jband-11900.0),ind)
wl_Jband_fix = Planet_wave_Jband[ind]

;*** 2 MASS H band

wave_planet_Hband = where (wl_planet ge min_wave_Hband and wl_planet le max_wave_Hband)
Planet_wave_Hband = wl_planet (wave_planet_Hband)
flux_planet_Hband = fl_planet_1 (wave_planet_Hband)

total_flux_planet_Hband = total(flux_planet_Hband)
flux_Hband = min(abs(Planet_wave_Hband-16500.0),ind)
wl_Hband_fix = Planet_wave_hband[ind]

;*** 2 MASS K band

wave_planet_Kband = where (wl_planet ge min_wave_Kband and wl_planet le max_wave_Kband)
Planet_wave_Kband = wl_planet (wave_planet_Kband)
flux_planet_Kband = fl_planet_1 (wave_planet_Kband)

total_flux_planet_Kband = total(flux_planet_Kband)
flux_Kband = min(abs(Planet_wave_Kband-20900.0),ind)
wl_Kband_fix = Planet_wave_Kband[ind]

;*** SPITZER 3.6

wave_planet_3_6 = where (wl_planet ge min_wave_3_6 and wl_planet le max_wave_3_6)
Planet_wave_3_6 = wl_planet (wave_planet_3_6)
flux_planet_3_6 = fl_planet_1 (wave_planet_3_6)

total_flux_planet_3_6 = total(flux_planet_3_6)

flux_3_6 = min(abs(Planet_wave_3_6-36000d),ind)
wl_3_6_fix = Planet_wave_3_6[ind]


;*** SPITZER 4.5

wave_planet_4_5 = where (wl_planet ge min_wave_4_5 and wl_planet le max_wave_4_5)
Planet_wave_4_5 = wl_planet (wave_planet_4_5)
flux_planet_4_5 = fl_planet_1 (wave_planet_4_5)

total_flux_planet_4_5 = total(flux_planet_4_5)
flux_4_5 = min(abs(Planet_wave_4_5-45000d),ind)
wl_4_5_fix = Planet_wave_4_5[ind]

;*** SPITZER 5.8

wave_planet_5_8 = where (wl_planet ge min_wave_5_8 and wl_planet le max_wave_5_8)
Planet_wave_5_8 = wl_planet (wave_planet_5_8)
flux_planet_5_8 = fl_planet_1 (wave_planet_5_8)

total_flux_planet_5_8 = total(flux_planet_5_8)
flux_5_8 = min(abs(Planet_wave_5_8-58000d),ind)
wl_5_8_fix = Planet_wave_5_8[ind]

;*** SPITZER 8.0

wave_planet_8_0 = where (wl_planet ge min_wave_8_0 and wl_planet le max_wave_8_0)
Planet_wave_8_0 = wl_planet (wave_planet_8_0)
flux_planet_8_0 = fl_planet_1 (wave_planet_8_0)

total_flux_planet_8_0 = total(flux_planet_8_0)
flux_8_0 = min(abs(Planet_wave_8_0-80000d),ind)
wl_8_0_fix = Planet_wave_8_0[ind]


readcol,'/Volumes/ARVIND_4TB/hot_jupiters/IRR_model/wasp_19b/cloud/spectra/lte055-4.5-0.0a+0.0.BT-Settl.spec.7', wl_star,fl_star,format='(f,f,f)'
fla_star = 10d^(fl_star - 8d)
wl_order_star = sort(wl_star)
Wavelength_star = wl_star[wl_order_star]
flux_star = fla_star[wl_order_star]
new_flux_star = interpol(flux_star, Wavelength_star, wl_planet)
new_flux_star = smooth(new_flux_star, 1000)
radius_star = 6.887e10

;******************************
; Integrating spectrum over
; different photometric bands
;******************************

; *** SDSS Z band


wave_star_zband = where (wl_planet ge min_wave_zband and wl_planet le max_wave_zband)
star_wave_zband = wl_planet (wave_star_zband)
flux_star_zband = new_flux_star (wave_star_zband)
total_flux_star_zband = total(flux_star_zband)

; *** SDSS J band

wave_star_Jband = where (wl_planet ge min_wave_Jband and wl_planet le max_wave_Jband)
star_wave_Jband = wl_planet (wave_star_Jband)
flux_star_Jband = new_flux_star (wave_star_Jband)
total_flux_star_Jband = total(flux_star_Jband)


; *** SDSS H band

wave_star_Hband = where (wl_planet ge min_wave_Hband and wl_planet le max_wave_Hband)
star_wave_Hband = wl_planet (wave_star_Hband)
flux_star_Hband = new_flux_star (wave_star_Hband)
total_flux_star_Hband = total(flux_star_Hband)

; *** SDSS K band

wave_star_Kband = where (wl_planet ge min_wave_Kband and wl_planet le max_wave_Kband)
star_wave_Kband = wl_planet (wave_star_Kband)
flux_star_Kband = new_flux_star (wave_star_Kband)
total_flux_star_Kband = total(flux_star_Kband)


;*** SPITZER 3.6

wave_star_3_6 = where (wl_planet ge min_wave_3_6 and wl_planet le max_wave_3_6)
star_wave_3_6 = wl_planet (wave_star_3_6)
flux_star_3_6 = new_flux_star (wave_star_3_6)
total_flux_star_3_6 = total(flux_star_3_6)

;*** SPITZER 4.5

wave_star_4_5 = where (wl_planet ge min_wave_4_5 and wl_planet le max_wave_4_5)
star_wave_4_5 = wl_planet (wave_star_4_5)
flux_star_4_5 = new_flux_star (wave_star_4_5)
total_flux_star_4_5 = total(flux_star_4_5)


;*** SPITZER 5.8

wave_star_5_8 = where (wl_planet ge min_wave_5_8 and wl_planet le max_wave_5_8)
star_wave_5_8 = wl_planet (wave_star_5_8)
flux_star_5_8 = new_flux_star (wave_star_5_8)
total_flux_star_5_8 = total(flux_star_5_8)

;*** SPITZER 8.0

wave_star_8_0 = where (wl_planet ge min_wave_8_0 and wl_planet le max_wave_8_0)
star_wave_8_0 = wl_planet (wave_star_8_0)
flux_star_8_0 = new_flux_star (wave_star_8_0)
total_flux_star_8_0 = total(flux_star_8_0)

;*************************************************************
radius_ratio = (radius_planet/radius_star)^2
final_ratio = (fl_planet_1/new_flux_star) * radius_ratio

flux_ratio_Zband =  (total_flux_planet_zband / total_flux_star_zband) * radius_ratio
flux_ratio_Jband =  (total_flux_planet_Jband / total_flux_star_Jband) * radius_ratio
flux_ratio_Hband =  (total_flux_planet_Hband / total_flux_star_Hband) * radius_ratio
flux_ratio_Kband =  (total_flux_planet_Kband / total_flux_star_Kband) * radius_ratio

flux_ratio_3_6 =  (total_flux_planet_3_6 / total_flux_star_3_6) * radius_ratio
flux_ratio_4_5 =  (total_flux_planet_4_5 / total_flux_star_4_5) * radius_ratio
flux_ratio_5_8 =  (total_flux_planet_5_8 / total_flux_star_5_8) * radius_ratio
flux_ratio_8_0 =  (total_flux_planet_8_0 / total_flux_star_8_0) * radius_ratio
;******************** Ploting the spectra *************************
ticks = LogLevels([0.5,10],/fine)
                       
nticks = n_elements(ticks)

cgplot,(wl_planet/10000), (final_ratio),color = color[n],yr=[0,0.010],xr=[0.4,1],linestyle=linestyle[n], $
ytitle='Planet-to-Star-Ratio',xtitle = 'Wavelength ($\mu$m)', thick=1.0,$
/xlog,ystyle=1,title='WASP-19b', charsize=1.5, xtickinterval=0.1,/overplot,/addcmd

;*******Observed Photometry *******************
cgplot, (w_phot/10000), (f_phot/100), ERR_XLow=dblarr(n_elements(w_phot)), ERR_XHigh=dblarr(n_elements(w_phot)), $
ERR_YLow=errors/100.0, ERR_YHigh=errors/100.0, ERR_Color='black', $
SymColor='black',Psym=16,symsize=1.2,/overplot,/addcmd

cgplot, (w_phot_1/10000), (f_phot_1/100), ERR_XLow=dblarr(n_elements(w_phot_1)), ERR_XHigh=dblarr(n_elements(w_phot_1)), $
ERR_YLow=errors_1/100.0, ERR_YHigh=errors_1/100.0, ERR_Color='black', $
SymColor='black',Psym=16,symsize=1.2,/overplot,/addcmd

       
cgplot, (w_phot_2/10000), (f_phot_2/100), ERR_XLow=dblarr(n_elements(w_phot_2)), ERR_XHigh=dblarr(n_elements(w_phot_2)), $
ERR_YLow=errors_2/100.0, ERR_YHigh=errors_2/100.0, ERR_Color='black', $
SymColor='black',Psym=16,symsize=1.2,/overplot,/addcmd
;********Synthetic Photometry ******************
cgoplot, (wl_Zband_fix/10000),flux_ratio_Zband, color=color[n],psym=15,symsize=1.5,/addcmd
cgoplot, (wl_Jband_fix/10000),flux_ratio_Jband, color=color[n],psym=15,symsize=1.5,/addcmd
cgoplot, (wl_Hband_fix/10000),flux_ratio_Hband, color=color[n],psym=15,symsize=1.5,/addcmd
cgoplot, (wl_Kband_fix/10000),flux_ratio_Kband, color=color[n],psym=15,symsize=1.5,/addcmd

cgoplot, (wl_3_6_fix/10000),flux_ratio_3_6, color=color[n],psym=15,symsize=1.5,/addcmd
cgoplot, (wl_5_8_fix/10000),flux_ratio_5_8, color=color[n],psym=15,symsize=1.5,/addcmd
cgoplot, (wl_8_0_fix/10000),flux_ratio_8_0, color=color[n],psym=15,symsize=1.5,/addcmd
cgoplot, (wl_4_5_fix/10000),flux_ratio_4_5, color=color[n],psym=15,symsize=1.5,/addcmd

al_legend,['[M/H]=+2.0','[M/H]=0.0','[M/H]=-2.0'], color=['red','green','blue'],linestyle=linestyle[n],/bottom,/right, box=0,thick=1.0,charsize=1.5,/window
;al_legend,['f = 1','f = 2'], color=['red','green'],linestyle=linestyle[n],/bottom,/right, box=0,thick=1.0,charsize=1.5,/window
;al_legend,['f = 1','f = 2'], color = ['grey','black'],linestyle=linestyle[n],/bottom,/right, box=0,thick=1.0,charsize=1.5,/window
;cgPS_Close
endfor

end 
