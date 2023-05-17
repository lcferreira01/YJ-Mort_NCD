///////////////////////////////////////////////////////////////////////////////
/// 	Project: 		YJ-Mort Study - Non-Communicable disease mortality
///		File type: 		Analysis
///		Date: 			15/05/2023
///		Analyst:		Lucas Calais Ferreira
///////////////////////////////////////////////////////////////////////////////


// --------------------------------------------------+
//   Descriptive statistics
// --------------------------------------------------+


/// Loading dataset
use "$datasets/yjmortNCDsetup_individual", clear


/// Median age at death
sum age if NCD_fail==1, detail
sum agefirst if morbseq==1, detail
sum fuyears if morbseq==1, detail
sum ageend if morbseq==1, detail


/// Table 1

tab agefirstdicho sex if yjhist==1, missing col r chi
tab agefirstdicho sex if yjhist==2, missing col r chi
tab agefirstdicho sex if yjhist==3, missing col r chi
tab agefirstdicho sex, missing col r chi
tab agefirstdicho yjhist, missing col r chi

tab indigstat sex if yjhist==1, missing col r chi
tab indigstat sex if yjhist==2, missing col r chi
tab indigstat sex if yjhist==3, missing col r chi
tab indigstat sex, missing col r chi
tab indigstat yjhist, missing col r chi

tab yjchrgtotcat sex if yjhist==1, missing col r chi
tab yjchrgtotcat sex if yjhist==2, missing col r chi
tab yjchrgtotcat sex if yjhist==3, missing col r chi
tab yjchrgtotcat sex, missing col r chi
tab yjchrgtotcat yjhist, missing col r chi

tab yjordtotcat sex if yjhist==1, missing col r chi
tab yjordtotcat sex if yjhist==2, missing col r chi
tab yjordtotcat sex if yjhist==3, missing col r chi
tab yjordtotcat sex, missing col r chi
tab yjordtotcat yjhist, missing col r chi

tab yjdettotcat sex if yjhist==1, missing col r chi
tab yjdettotcat sex if yjhist==2, missing col r chi
tab yjdettotcat sex if yjhist==3, missing col r chi
tab yjdettotcat sex, missing col r chi
tab yjdettotcat yjhist, missing col r chi

tab qcscustotcat sex if yjhist==1, missing col r chi
tab qcscustotcat sex if yjhist==2, missing col r chi
tab qcscustotcat sex if yjhist==3, missing col r chi
tab qcscustotcat sex, missing col r chi
tab qcscustotcat yjhist, missing col r chi

tab yjdettotcat sex if yjhist==1, missing col r chi
tab yjdettotcat sex if yjhist==2, missing col r chi
tab yjdettotcat sex if yjhist==3, missing col r chi
tab yjdettotcat yjhist, missing col r chi

tab ndidied sex if yjhist==1, missing col r chi
tab ndidied sex if yjhist==2, missing col r chi
tab ndidied sex if yjhist==3, missing col r chi
tab ndidied yjhist, missing col r chi

tab ageendcat sex if yjhist==1 & ndidied==1, missing col r chi
tab ageendcat sex if yjhist==2 & ndidied==1, missing col r chi
tab ageendcat sex if yjhist==3 & ndidied==1, missing col r chi
tab ageendcat sex if ndidied==1, missing col r chi
tab ageendcat yjhist if ndidied==1, missing col r chi


//////////////////////////////////////////////////////
// Table 2
//////////////////////////////////////////////////////

/// Loading dataset
use "$datasets/yjmortNCDsetup_survival", clear


// --------------------------------------------------+
//   Counts and percentages
// --------------------------------------------------+

tab yjhist if NCD_fail==1
tab agefirstdicho if NCD_fail==1
tab sex if NCD_fail==1
tab indigstat if NCD_fail==1
tab indigsex if NCD_fail==1


// --------------------------------------------------+
//   Time-at-risk
// --------------------------------------------------+

*All cause deaths
stset datex, fail(fail) id(yjid) origin(time startdate) scale(365.25)

strate, per(100000)


*NCD deaths - stset
stset datex, fail(NCD_fail) id(yjid) origin(time startdate) scale(365.25)

*length of follow-up in each observation
gen length = _t-_t0


*number of deaths and followup for each variable 
foreach var in yjhist agefirstdicho sex indigstat indigsex {
strate `var', per(100000)
tabstat length, by(`var') s(sum)
}


//Calculations

//Yj hist
display 365784.3/653651.2
display 185261.6/653651.2
display 102605.2/653651.2

//Age at first contact
display 92233.79/653651.2
display 561417.4/653651.2

//Indigenous
display 173354.3/653651.2
display 480296.9/653651.2

//Indigsex
display 46193.56/653651.2
display 127160.7/653651.2
display 104707.3/653651.2
display 375589.5/653651.2


// --------------------------------------------------+
//   CMRs 
// --------------------------------------------------+

*All cause - stset
stset datex, fail(fail) id(yjid) origin(time startdate) scale(365.25) 	

*overall
strate, per(100000)

*NCD cause - stset
stset datex, fail(NCD_fail) id(yjid) origin(time startdate) scale(365.25)

*Overall
strate, per(100000)

*By youth justice history
strate yjhist, per(100000)

*By age at first contact
strate agefirstdicho, per(100000)

*By sex
strate sex, per(100000)

*By Indigenous status
strate indigstat, per(100000)

*By Indigenous status and sex
strate indigsex, per(100000)


// --------------------------------------------------+
//   SMRs - standardised by sex and age group
// --------------------------------------------------+

/// Loading dataset
use "$datasets/yjmortNCD_SMRsetup_v3", clear

*All cause - stset
stset datex, fail(fail) id(yjid) origin(time startdate) scale(365.25) 	

*overall
strate, smr(CMR_allcause) per(100000)

*NCD cause - stset
stset datex, fail(NCD_fail) id(yjid) origin(time startdate) scale(365.25) 	

*overall
strate, smr(CMR_ncd) per(100000)

*By youth justice history
strate yjhist, smr(CMR_ncd) per(100000)

*By age at first contact
strate agefirstdicho, smr(CMR_ncd) per(100000)

*By sex
strate sex, smr(CMR_ncd) per(100000)

*By Indigenous status
strate indigstat, smr(CMR_ncd) per(100000)

*By Indigenous status and sex
strate indigsex, smr(CMR_ncd) per(100000)



//////////////////////////////////////////////////////
// Table 3 - Cause specific CMRs and SMRs by sex and Indigenous identification
//////////////////////////////////////////////////////

/// Loading dataset
use "$datasets/yjmortNCD_SMRsetup_v3", clear


// --------------------------------------------------+
//   Neoplasms
// --------------------------------------------------+

//Neoplasms - stset
stset datex, fail(neoplasms_fail) id(yjid) origin(time startdate) scale(365.25) 	

*CMR overall
strate, per(100000)

*CMR
strate indigstat sex, per(100000)

*SMR overall
strate, smr(CMR_neoplasms) per(100000)

*SMR1 by sex
strate indigstat sex, smr(CMR_neoplasms) per(100000)


// --------------------------------------------------+
//   Cardiovascular
// --------------------------------------------------+

//Cardiovascular - stset
stset datex, fail(cardio_fail) id(yjid) origin(time startdate) scale(365.25) 	

*CMR overall
strate, per(100000)

*CMR
strate indigstat sex, per(100000)

*SMR overall
strate, smr(CMR_cardio) per(100000)

*SMR1 by sex
strate indigstat sex, smr(CMR_cardio) per(100000)


// --------------------------------------------------+
//   Respiratory
// --------------------------------------------------+

//Respiratory - stset
stset datex, fail(resp_fail) id(yjid) origin(time startdate) scale(365.25) 	

*CMR overall
strate, per(100000)

*CMR
strate indigstat sex, per(100000)

*SMR overall
strate, smr(CMR_resp) per(100000)

*SMR1 by sex
strate indigstat sex, smr(CMR_resp) per(100000)


// --------------------------------------------------+
//   Digestive
// --------------------------------------------------+

//Digestive - stset
stset datex, fail(digestive_fail) id(yjid) origin(time startdate) scale(365.25) 	

*CMR overall
strate, per(100000)

*CMR
strate indigstat sex, per(100000)

*SMR overall
strate, smr(CMR_digestive) per(100000)

*SMR1 by sex
strate indigstat sex, smr(CMR_digestive) per(100000)


// --------------------------------------------------+
//   Neuro
// --------------------------------------------------+

//Neuro - stset
stset datex, fail(neuro_fail) id(yjid) origin(time startdate) scale(365.25) 	

*CMR overall
strate, per(100000)

*CMR
strate indigstat sex, per(100000)

*SMR overall
strate, smr(CMR_neuro) per(100000)

*SMR1 by sex
strate indigstat sex, smr(CMR_neuro) per(100000)


// --------------------------------------------------+
//   Substance use
// --------------------------------------------------+

//Substance - stset
stset datex, fail(substance_fail) id(yjid) origin(time startdate) scale(365.25) 	

*CMR overall
strate, per(100000)

*CMR
strate indigstat sex, per(100000)

*SMR overall
strate, smr(CMR_substance) per(100000)

*SMR1 by sex
strate indigstat sex, smr(CMR_substance) per(100000)


// --------------------------------------------------+
//   Other NCDs
// --------------------------------------------------+

//Other NCDs - stset
stset datex, fail(otherncds_fail) id(yjid) origin(time startdate) scale(365.25) 	

*CMR overall
strate, per(100000)

*CMR
strate indigstat sex, per(100000)

*SMR overall
strate, smr(CMR_otherncds) per(100000)

*SMR1 by sex
strate indigstat sex, smr(CMR_otherncds) per(100000)



///////////////////////
/// Table 4
///////////////////////


// --------------------------------------------------+
//   SUBHAZARD RATIOS: Fine and Gray competing risk regression
// --------------------------------------------------+

/// Loading dataset
use "$datasets/yjmortNCDsetup_survival", clear

*NCD causes
stset datex, fail(NCD_fail) id(yjid) origin(time startdate) scale(365.25)

///Recoding rolling count of YJ charges

recode yjchrgcntcat 0=1


*length of follow-up in each observation
gen length = _t-_t0


*number of deaths and followup for each variable 
foreach var in indigstat sex agefirstdicho yjpenetration legal_2 qcscuscntcat yjchrgcntcat yjordcntcat yjdetcntcat qcscustimecat yjdettimecat {
strate `var', per(100000)
tabstat length, by(`var') s(sum)
}


///// Competing risks analysis
replace fail=2 if fail==1 & NCD_group>=1 & NCD_group<=8

*All cause
stset datex, fail(fail==2) id(yjid) origin(time startdate) scale(365.25)


//Univariable Competing risks regression
foreach var in ib2.indigstat i.sex i.agefirstdicho i.yjpenetration i.yjhist i.legal_2 i.qcscuscntcat i.yjchrgcntcat i.yjordcntcat i.yjdetcntcat i.qcscustimecat i.yjdettimecat {
	stcrreg `var', compete(fail==1)
}


/// Adjusted for sex and indigenous status
foreach var in ib2.indigstat i.sex i.agefirstdicho i.yjpenetration i.yjhist i.legal_2 i.qcscuscntcat i.yjchrgcntcat i.yjordcntcat i.yjdetcntcat i.qcscustimecat i.yjdettimecat {
	stcrreg `var' i.sex ib2.indigstat, compete(fail==1)
}

gen age_cat=.
replace age_cat=0 if age<19
replace age_cat=1 if age>=19 & age<.

/// Adjusted for sex, indigenous status and age (cat)
foreach var in ib2.indigstat i.sex i.agefirstdicho i.yjpenetration i.yjhist i.legal_2 i.qcscuscntcat i.yjchrgcntcat i.yjordcntcat i.yjdetcntcat i.qcscustimecat i.yjdettimecat {
	stcrreg `var' i.sex ib2.indigstat i.age_cat, compete(fail==1)
}



//////////////////////////////////////////////////////
// Graphs - Cumulative incidence curves
//////////////////////////////////////////////////////


/// Scheme
net install cleanplots, from("https://tdmize.github.io/data/cleanplots")
set scheme cleanplots, perm
graph set window fontface "Arial Narrow"

net install palettes, replace from("https://raw.githubusercontent.com/benjann/palettes/master/")
net install colrspace, replace from("https://raw.githubusercontent.com/benjann/colrspace/master/")

// Example
colorpalette ///
AliceBlue AntiqueWhite Aqua Aquamarine Azure Beige Bisque Black ///
BlanchedAlmond Blue BlueViolet Brown BurlyWood CadetBlue ///
Chartreuse Chocolate Coral CornflowerBlue Cornsilk Crimson ///
Cyan DarkBlue DarkCyan DarkGoldenRod DarkGray DarkGrey ///
DarkGreen DarkKhaki DarkMagenta DarkOliveGreen DarkOrange ///
DarkOrchid DarkRed DarkSalmon DarkSeaGreen DarkSlateBlue /// 
DarkSlateGray DarkTurquoise DarkViolet DeepPink DeepSkyBlue /// 
DimGray DimGrey DodgerBlue FireBrick FloralWhite ForestGreen ///
Fuchsia Gainsboro GhostWhite Gold GoldenRod Gray Grey Green /// 
GreenYellow HoneyDew HotPink IndianRed Indigo Ivory ///
Khaki Lavender LavenderBlush LawnGreen LemonChiffon LightBlue /// 
LightCoral LightCyan LightGoldenRodYellow LightGray LightGrey /// 
LightGreen LightPink LightSalmon LightSeaGreen LightSkyBlue /// 
LightSlateGray LightSlateGrey LightSteelBlue LightYellow Lime /// 
LimeGreen Linen Magenta Maroon MediumAquaMarine MediumBlue /// 
MediumOrchid MediumPurple MediumSeaGreen MediumSlateBlue ///
MediumSpringGreen MediumTurquoise MediumVioletRed MidnightBlue /// 
MintCream MistyRose Moccasin NavajoWhite Navy OldLace Olive /// 
OliveDrab Orange OrangeRed Orchid PaleGoldenRod PaleGreen /// 
PaleTurquoise PaleVioletRed PapayaWhip PeachPuff Peru ///
Pink Plum PowderBlue Purple RebeccaPurple Red RosyBrown ///
RoyalBlue SaddleBrown Salmon SandyBrown SeaGreen SeaShell Sienna /// 
Silver SkyBlue SlateBlue SlateGray SlateGrey Snow SpringGreen /// 
SteelBlue Tan Teal Thistle Tomato Turquoise Violet ///
Wheat White WhiteSmoke Yellow YellowGreen, ///
title("Colrspace named colors")

// Interpolation
colorpalette forest_green gold, ipolate(10)


///////////////////////////////////////
/////// Cumulative incidence curves
///////////////////////////////////////

/// Loading dataset
use "$datasets/yjmortNCDsetup_survival", clear

*NCD causes
stset datex, fail(NCD_fail) id(yjid) origin(time startdate) scale(365.25)


/// Cumulative incidence curve (by youth justice penetration)

sts graph, cumhaz by(yjpenetration) /// 
 risktable(, rowtitle("Charge only") failevents group(1)) ///
 risktable(, rowtitle("Community-based order") failevents group(2)) ///
 risktable(, rowtitle("Detention") failevents group(3)) ///
 risktable(, title("Risk table", size(vsmall)) size(vsmall)) ///
ytitle(Cumulative incidence (NCD deaths)) xtitle(Study time (years)) yscale(range(0.0(0.005)0.015)) ///
		ylabel(0 0.005 0.010 0.015) /// 
			legend(region(lstyle(none) color(none))  ///
				size(small) cols(1) ring(0) bplacement(nw) symysize(0.7) rowgap(0.5) ///
				order(3 2 1) ///
				label(1 "Charge only") /// 
				label(2 "Community-based order") ///
				label(3 "Detention")) ///
				xlabel(0(10)25)  ///
				title("") ///
		graphregion(color(white)) saving("$graphs\NCD_yjp_no_ci.gph", replace)
		
graph export "$graphs\NCD_yjp_no_ci.tif", as(tif) replace


/// Cumulative incidence curve (by Indigenous status and sex)

sts graph, cumhaz by(indigsex) /// 
 risktable(, rowtitle("Indigenous Females") failevents group(1)) ///
 risktable(, rowtitle("Indigenous Males") failevents group(2)) ///
 risktable(, rowtitle("Non-Indigenous Females") failevents group(3)) ///
 risktable(, rowtitle("Non-Indigenous Males") failevents group(4)) ///
 risktable(, title("Risk table", size(vsmall)) size(vsmall)) ///
ytitle(Cumulative incidence (NCD deaths)) xtitle(Study time (years)) yscale(range(0.0(0.005)0.01)) ///
		ylabel(0 0.005 0.01 0.015) /// 
			legend(region(lstyle(none) color(none))  ///
				size(small) cols(1) ring(0) bplacement(nw) symysize(0.7) rowgap(0.5) ///
				order(2 4 1 3) ///
				label(1 "Indigenous Females") /// 
				label(2 "Indigenous Males") ///
				label(3 "Non-Indigenous Females") ///
				label(4 "Non-Indigenous Males")) ///
				xlabel(0(10)25)  ///
				title("") ///
		graphregion(color(white)) saving("$graphs\NCD_indigsex_no_ci.gph", replace)
		
graph export "$graphs\NCD_indigsex_no_ci.tif", as(tif) replace


//////////////////////////////////////////////////////
// Graphs - NCD causes
//////////////////////////////////////////////////////


//Neoplasms - stset
stset datex, fail(neoplasms_fail) id(yjid) origin(time startdate) scale(365.25) 	
*CMR overall
stptime, per(100000) output("$datasets\neoplasms_CMR_graph.dta", replace)

//Cardiovascular - stset
stset datex, fail(cardio_fail) id(yjid) origin(time startdate) scale(365.25) 	
*CMR overall
stptime, per(100000) output("$datasets\cardio_CMR_graph.dta", replace)

//Respiratory - stset
stset datex, fail(resp_fail) id(yjid) origin(time startdate) scale(365.25) 	
*CMR overall
stptime, per(100000) output("$datasets\resp_CMR_graph.dta", replace)

//Respiratory - stset
stset datex, fail(resp_fail) id(yjid) origin(time startdate) scale(365.25) 	
*CMR overall
stptime, per(100000) output("$datasets\resp_CMR_graph.dta", replace)

//Digestive - stset
stset datex, fail(digestive_fail) id(yjid) origin(time startdate) scale(365.25) 	
*CMR overall
stptime, per(100000) output("$datasets\digestive_CMR_graph.dta", replace)

//Neuro - stset
stset datex, fail(digestive_fail) id(yjid) origin(time startdate) scale(365.25) 	
*CMR overall
stptime, per(100000) output("$datasets\neuro_CMR_graph.dta", replace)

//Substance - stset
stset datex, fail(substance_fail) id(yjid) origin(time startdate) scale(365.25) 	
*CMR overall
stptime, per(100000) output("$datasets\substance_CMR_graph.dta", replace)

//Substance - stset
stset datex, fail(otherncds_fail) id(yjid) origin(time startdate) scale(365.25) 	
*CMR overall
stptime, per(100000) output("$datasets\otherncds_CMR_graph.dta", replace)


/// Loading dataset
use "$datasets\neoplasms_CMR_graph.dta", clear
gen NCD_group= 1
la de NCD_group 1 "Neoplasms" 2 "Cardiovascular" 3 "Respiratory" 4 "Digestive" ///
5 "Neurological" 6 "Substance use" 7 "Other NCDs"
la val NCD_group NCD_group

//append on other conditions 
append using "$datasets\cardio_CMR_graph.dta"
recode NCD_group (. = 2)

//append on other conditions 
append using "$datasets\resp_CMR_graph.dta"
recode NCD_group (. = 3)

//append on other conditions 
append using "$datasets\digestive_CMR_graph.dta"
recode NCD_group (. = 4)

//append on other conditions 
append using "$datasets\neuro_CMR_graph.dta"
recode NCD_group (. = 5)

//append on other conditions 
append using "$datasets\substance_CMR_graph.dta"
recode NCD_group (. = 6)

//append on other conditions 
append using "$datasets\otherncds_CMR_graph.dta"
recode NCD_group (. = 7)


/// Droping those who did not die
drop if _group==0


//graph results

twoway con _Rate NCD_group if NCD_group == 1, ///
	xtitle("NCD group", margin(5) si(small))  ///
	xlabel(none) /// 
	ytitle("Crude mortality rate per 100 000 person-years and 95%CI)" " ", bm(right) si(small)) /// 
	ylabel(0(1)10, nogrid labs(vsmall) format(%4.0f)) ///
	lcolor(navy) msize(medium) mcolor(navy) c(navy) ///
    || rcap _Upper _Lower NCD_group if NCD_group==1, lcolor(navy) mcolor(navy) msize(small) ///
	|| con _Rate NCD_group if NCD_group == 2, lcolor(maroon) lwidth(medium) mcolor(maroon) c(maroon) ///
	|| rcap _Upper _Lower NCD_group if NCD_group == 2, lcolor(maroon) msize(medium) ///
	|| con _Rate NCD_group if NCD_group == 3, lcolor(forest_green) lwidth(medium) mcolor(forest_green) c(forest_green) ///
	|| rcap _Upper _Lower NCD_group if NCD_group == 3, lcolor(forest_green) msize(medium) ///
	|| con _Rate NCD_group if NCD_group == 4, lcolor(dkorange) lwidth(medium) mcolor(dkorange) c(dkorange) ///
	|| rcap _Upper _Lower NCD_group if NCD_group == 4, lcolor(dkorange) msize(medium) ///
	|| con _Rate NCD_group if NCD_group == 5, lcolor(teal) lwidth(medium) mcolor(teal) c(teal)  ///
	|| rcap _Upper _Lower NCD_group if NCD_group == 5, lcolor(teal) msize(medium) ///
	|| con _Rate NCD_group if NCD_group == 6, lcolor(grey) lwidth(medium) mcolor(grey) c(grey) ///
	|| rcap _Upper _Lower NCD_group if NCD_group == 6, lcolor(grey) msize(medium) ///
	|| con _Rate NCD_group if NCD_group == 7, lcolor(purple) lwidth(medium) mcolor(purple) c(purple) ///
	|| rcap _Upper _Lower NCD_group if NCD_group == 7, lcolor(purple) msize(medium) ///
	legend(size(vsmall) position(6) region(lc(white)) order(1 5 9 13 3 7 11) cols(4) label(1 "Neoplasms") ///
	label(3 "Cardiovascular diseases") label(5 "Respiratory diseases") label(7 "Digestive diseases") ///
	label(9 "Neurological disorders") label(11 "Substance use disorders") label(13 "Other NCDs")) ///
		graphregion(color(white)) saving("$graphs\NCD_cause_specific_rates.gph", replace)
		
graph export "$graphs\NCD_cause_specific_rates.tif", as(tif) replace
	

/// SMR graphs

/// Loading dataset
use "$datasets/yjmortNCD_SMRsetup_v3", clear


/// Adjusted for sex, indigenous status and age
foreach var in neoplasms cardio resp digestive neuro substance otherncds {
	stset datex, fail(`var'_fail) id(yjid) origin(time startdate) scale(365.25)
	strate, smr(CMR_`var') per(100000) output("`var'_SMR_graph.dta", replace)
}

	
use "neoplasms_SMR_graph.dta", clear

gen NCD_group= 1
la de NCD_group 1 "Neoplasms" 2 "Cardiovascular" 3 "Respiratory" 4 "Digestive" ///
5 "Neurological" 6 "Substance use" 7 "Other NCDs"
la val NCD_group NCD_group

//append on other conditions 
append using "cardio_SMR_graph.dta"
recode NCD_group (. = 2)

//append on other conditions 
append using "resp_SMR_graph.dta"
recode NCD_group (. = 3)

//append on other conditions 
append using "digestive_SMR_graph.dta"
recode NCD_group (. = 4)

//append on other conditions 
append using "neuro_SMR_graph.dta"
recode NCD_group (. = 5)

//append on other conditions 
append using "substance_SMR_graph.dta"
recode NCD_group (. = 6)

//append on other conditions 
append using "otherncds_SMR_graph.dta"
recode NCD_group (. = 7)



//graph results

twoway con _SMR NCD_group if NCD_group == 1, ///
	xtitle("NCD group", margin(5) si(small))  ///
	xlabel(none) /// 
	ytitle("Standardised Mortality Ratio and 95%CI" " ", bm(right) si(small)) /// 
		ylabel(0(1)10, nogrid labs(vsmall) format(%4.0f)) ///
	lcolor(navy) msize(medium) mcolor(navy) c(navy) ///
    || rcap _Upper _Lower NCD_group if NCD_group==1, lcolor(navy) mcolor(navy) msize(small) ///
	|| con _SMR NCD_group if NCD_group == 2, lcolor(maroon) lwidth(medium) mcolor(maroon) c(maroon) ///
	|| rcap _Upper _Lower NCD_group if NCD_group == 2, lcolor(maroon) msize(medium) ///
	|| con _SMR NCD_group if NCD_group == 3, lcolor(forest_green) lwidth(medium) mcolor(forest_green) c(forest_green) ///
	|| rcap _Upper _Lower NCD_group if NCD_group == 3, lcolor(forest_green) msize(medium) ///
	|| con _SMR NCD_group if NCD_group == 4, lcolor(dkorange) lwidth(medium) mcolor(dkorange) c(dkorange) ///
	|| rcap _Upper _Lower NCD_group if NCD_group == 4, lcolor(dkorange) msize(medium) ///
	|| con _SMR NCD_group if NCD_group == 5, lcolor(teal) lwidth(medium) mcolor(teal) c(teal)  ///
	|| rcap _Upper _Lower NCD_group if NCD_group == 5, lcolor(teal) msize(medium) ///
	|| con _SMR NCD_group if NCD_group == 6, lcolor(grey) lwidth(medium) mcolor(grey) c(grey) ///
	|| rcap _Upper _Lower NCD_group if NCD_group == 6, lcolor(grey) msize(medium) ///
	|| con _SMR NCD_group if NCD_group == 7, lcolor(purple) lwidth(medium) mcolor(purple) c(purple) ///
	|| rcap _Upper _Lower NCD_group if NCD_group == 7, lcolor(purple) msize(medium) ///
	legend(size(vsmall) position(6) region(lc(white)) order(1 5 9 13 3 7 11) cols(4) label(1 "Neoplasms") ///
	label(3 "Cardiovascular diseases") label(5 "Respiratory diseases") label(7 "Digestive diseases") ///
	label(9 "Neurological disorders") label(11 "Substance use disorders") label(13 "Other NCDs")) ///
		graphregion(color(white)) saving("$graphs\NCD_cause_specific_SMRs.gph", replace)
		
graph export "$graphs\NCD_cause_specific_SMRs.tif", as(tif) replace
	

	


//////////////////////////////////////////////////////
// Table S3 - Cause specific SMRs by sex, age and Indigenous identification
//////////////////////////////////////////////////////


// --------------------------------------------------+
//   Neoplasms
// --------------------------------------------------+

/// Loading dataset
use "$datasets/yjmortNCD_SMRsetup_v3", clear

//Neoplasms - stset
stset datex, fail(neoplasms_fail) id(yjid) origin(time startdate) scale(365.25) 	

*CMR overall
strate, per(100000)

*CMR
strate indigstat sex, per(100000)

*SMR overall
strate, smr(CMR_neoplasms) per(100000)

*SMR1 by sex
strate indigstat sex, smr(CMR_neoplasms) per(100000)

*SMR2 by sex abd indigenous status

use "$datasets/yjmortNCD_SMRsetup_indig_v3", clear

stset datex, fail(neoplasms_fail) id(yjid) origin(time startdate) scale(365.25) 

strate indigstat sex, smr(CMR_neoplasms) per(100000)


// --------------------------------------------------+
//   Cardiovascular
// --------------------------------------------------+

//Cardiovascular - stset
stset datex, fail(cardio_fail) id(yjid) origin(time startdate) scale(365.25) 	

*CMR overall
strate, per(100000)

*CMR
strate indigstat sex, per(100000)

*SMR overall
strate, smr(CMR_cardio) per(100000)

*SMR1 by sex
strate indigstat sex, smr(CMR_cardio) per(100000)

*SMR2 by sex abd indigenous status

use "$datasets/yjmortNCD_SMRsetup_indig_v3", clear

stset datex, fail(cardio_fail) id(yjid) origin(time startdate) scale(365.25) 

strate indigstat sex, smr(CMR_cardio) per(100000)


// --------------------------------------------------+
//   Respiratory
// --------------------------------------------------+

//Respiratory - stset
stset datex, fail(resp_fail) id(yjid) origin(time startdate) scale(365.25) 	

*CMR overall
strate, per(100000)

*CMR
strate indigstat sex, per(100000)

*SMR overall
strate, smr(CMR_resp) per(100000)

*SMR1 by sex
strate indigstat sex, smr(CMR_resp) per(100000)

*SMR2 by sex abd indigenous status

use "$datasets/yjmortNCD_SMRsetup_indig_v3", clear

stset datex, fail(resp_fail) id(yjid) origin(time startdate) scale(365.25) 

strate indigstat sex, smr(CMR_resp) per(100000)


// --------------------------------------------------+
//   Digestive
// --------------------------------------------------+

//Digestive - stset
stset datex, fail(digestive_fail) id(yjid) origin(time startdate) scale(365.25) 	

*CMR overall
strate, per(100000)

*CMR
strate indigstat sex, per(100000)

*SMR overall
strate, smr(CMR_digestive) per(100000)

*SMR1 by sex
strate indigstat sex, smr(CMR_digestive) per(100000)

*SMR2 by sex abd indigenous status

use "$datasets/yjmortNCD_SMRsetup_indig_v3", clear

stset datex, fail(digestive_fail) id(yjid) origin(time startdate) scale(365.25) 

strate indigstat sex, smr(CMR_digestive) per(100000)


// --------------------------------------------------+
//   Neuro
// --------------------------------------------------+


//Neuro - stset
stset datex, fail(neuro_fail) id(yjid) origin(time startdate) scale(365.25) 	

*CMR overall
strate, per(100000)

*CMR
strate indigstat sex, per(100000)

*SMR overall
strate, smr(CMR_neuro) per(100000)

*SMR1 by sex
strate indigstat sex, smr(CMR_neuro) per(100000)

*SMR2 by sex abd indigenous status

use "$datasets/yjmortNCD_SMRsetup_indig_v3", clear

stset datex, fail(neuro_fail) id(yjid) origin(time startdate) scale(365.25) 

strate indigstat sex, smr(CMR_neuro) per(100000)


// --------------------------------------------------+
//   Substance use
// --------------------------------------------------+


//Substance - stset
stset datex, fail(substance_fail) id(yjid) origin(time startdate) scale(365.25) 	

*CMR overall
strate, per(100000)

*CMR
strate indigstat sex, per(100000)

*SMR overall
strate, smr(CMR_substance) per(100000)

*SMR1 by sex
strate indigstat sex, smr(CMR_substance) per(100000)

*SMR2 by sex abd indigenous status

use "$datasets/yjmortNCD_SMRsetup_indig_v3", clear

stset datex, fail(substance_fail) id(yjid) origin(time startdate) scale(365.25) 

strate indigstat sex, smr(CMR_substance) per(100000)


// --------------------------------------------------+
//   Other NCDs
// --------------------------------------------------+


//Other NCDs - stset
stset datex, fail(otherncds_fail) id(yjid) origin(time startdate) scale(365.25) 	

*CMR overall
strate, per(100000)

*CMR
strate indigstat sex, per(100000)

*SMR overall
strate, smr(CMR_otherncds) per(100000)

*SMR1 by sex
strate indigstat sex, smr(CMR_otherncds) per(100000)

*SMR2 by sex abd indigenous status

use "$datasets/yjmortNCD_SMRsetup_indig_v3", clear

stset datex, fail(otherncds_fail) id(yjid) origin(time startdate) scale(365.25) 

strate indigstat sex, smr(CMR_otherncds) per(100000)
	

// --------------------------------------------------+
//   Table S4 - SMRs - standardised by sex, age grop and Indigenous identification
// --------------------------------------------------+
	
/// Loading dataset
use "$datasets/yjmortNCD_SMRsetup_indig_v3", clear

*NCD cause - stset
stset datex, fail(NCD_fail) id(yjid) origin(time startdate) scale(365.25) 	

*By Indigenous status
strate indigstat, smr(CMR_ncd) per(100000)

*By sex and Indigenous status
strate sex indigstat, smr(CMR_ncd) per(100000)


// --------------------------------------------------+
//   Additional analysis: Request by reviewer
// --------------------------------------------------+

use "$datasets/yjmortNCDsetup_survival", clear

gen attributable_death=.
replace attributable_death=1 if disease_group>1 & disease_group<4
replace attributable_death=0 if disease_group==9

tab attributable_death sex if fail==1, r col chi
tab attributable_death indigstat if fail==1, r col chi
tab attributable_death ageendcat if fail==1, r col chi

ttest ageend if fail==1, by(attributable_death)


// --------------------------------------------------+
//   Additional analysis: Request by reviewer - Goodness of fit
// --------------------------------------------------+

use "$datasets/yjmortNCDsetup_survival", clear

/// Setting up for competing risks
replace fail=2 if fail==1 & NCD_group>=1 & NCD_group<=8

stset datex, fail(NCD_fail) id(yjid) origin(time startdate) scale(365.25)

/// AIC and BIC
stcrreg i.sex ib2.indigstat i.age_cat, compete(fail==1)
estimates store with_age
estat ic

stcrreg i.sex ib2.indigstat, compete(fail==1)
estimates store without_age
estat ic

/// Likelihood ratio test
lrtest with_age without_age



log close
exit
