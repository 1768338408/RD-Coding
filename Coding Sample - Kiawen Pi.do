*******************************************************
* Project: The effect of housing prices on households fertility (Coding Sample)
* Purpose: Clean data, estimate baseline & IV models,
*          explore mechanisms, heterogeneity, robustness
* Author : (Kaiwen Pi)
* Date   : (2025-08-01)
*******************************************************

version 18
clear all
set more off
set seed 20250804

*******************************************************
* 0) Load main micro data (household level)
*******************************************************
use "XXXX", clear

* Standardize key identifiers
capture rename Fid   fid
capture rename fid10 fid10   // keep because we already constructed the panel data CFPS(6 waves/2010, 2012, 2014, 2016, 2018, 2020)
sort fid year

*******************************************************
* 1) Data cleaning & construction 
*******************************************************

* 1.1 Encode gender/health/hukou(residence) from Chinese-labeled strings to binaries
* Gender: 1=male ("男"), 0=female otherwise
gen byte gender_bin = (gender == "男")
label var gender_bin "Male=1, Female=0"

* Health: 1=healthy ("健康"), 0=unhealthy ("不健康"); treat others/missing as 0 after check
gen byte health_bin = .
replace health_bin = 1 if health == "健康"
replace health_bin = 0 if health == "不健康"
label var health_bin "Healthy=1, Unhealthy=0"

* Hukou: 1=agricultural ("农业户口"), 0=non-agricultural ("非农业户口")
gen byte hukou_agri = .
replace hukou_agri = 1 if hukou == "农业户口"
replace hukou_agri = 0 if hukou == "非农业户口"
label var hukou_agri "Agricultural hukou=1, non-agricultural=0"

* Drop original string vars and standardize names used later
drop gender health hukou
rename (gender_bin health_bin hukou_agri) (gender health hukou)

* 1.2 Ensure numerics
destring _all, replace force

* 1.3 New birth indicator: 1 if a child is born in current year
gen byte new_birth = (birth_yc == year) if !missing(birth_yc, year)
label var new_birth "Child birth in current year (=1), else 0"

* 1.4 Construct child/parent age variables
gen child_age = year - birth_yc
gen birth_age = age - child_age
label var child_age "Age of the newborn child (years)"
label var birth_age "Age of mother/father at birth (years)"   // adapt to the definition
drop if missing(birth_age)

* Age-at-birth sample restriction (standard in fertility literature)
drop if birth_age < 18
drop if birth_age > 40

* 1.5 Property and key covariates cleaning
destring property, replace force

foreach v in property income num_f age gender edu_year health hukou ratio_cd new_birth birth_age {
    drop if missing(`v')
}

drop if property < 0

* Restrict to births from 2010 onward (aligns with original script)
keep if birth_yc >= 2010

* 1.6 Transformations and winsorization
gen lnincome   = ln(1 + income)
gen lnproperty = ln(1 + property)

* Winsorize at 1-99% to reduce influence of extreme values
winsor2 lnproperty lnincome num_f age gender edu_year health hukou ratio_cd new_birth birth_age, cut(1 99) replace

* 1.7 Descriptive statistics (export to Word)
logout, dec(3) save("Descriptive_Stats") word replace: ///
    tabstat lnproperty lnincome num_f age gender edu_year health hukou ratio_cd new_birth birth_age, ///
    stats(N mean sd median min max) c(s) f(%10.3f)

*******************************************************
* 2) City-level resources & renaming Chinese variables
*******************************************************
* Original Chinese variables are explicitly renamed to English for clarity.
* See full dictionary at the bottom of this file.
* Handle missing water (years indexed 2–6 in the original panel; adjust if your year is calendar)
capture confirm variable total_water_city_10k_m3
if !_rc {
    replace total_water_city_10k_m3 = . if inrange(year, 2, 6)
    gen ww = total_water_city_10k_m3
    bys fid10 (year): replace ww = ww[1] if missing(ww) & inrange(year, 2, 6)
    * Normalize water by a constant used in the original script
    replace ww = ww/7635
    label var ww "City water resources (normalized)"
}

*******************************************************
* 3) Switch to an aligned balanced panel (household head)
*******************************************************
* The original workflow switches here and continues analysis.
use "XXXX", clear

* Fix community id (cid) by lead/lag within household (fid)
bys fid: replace cid = cid[_n+1]
bys fid: replace cid = cid[_n-1] if cid==.

* Zero to missing for per-sqm residential value
replace resivaluep = . if resivaluep == 0

*******************************************************
* 4) Construct core price variables and panel structure
*******************************************************
* Community mean price per m^2 by (cid, year)
bysort cid year: egen cvaluep = mean(resivaluep)
label var cvaluep "Community avg. market price per m^2"

* Panel declaration
xtset fid year

* Price differences: Δp_{t-1,t-2}
gen lcvaluep  = L.cvaluep
gen l2cvaluep = L.lcvaluep
gen dvaluep   = lcvaluep - l2cvaluep
gen lndvaluep = ln(dvaluep)
label var dvaluep "Δcvaluep = cvaluep(t-1)-cvaluep(t-2)"

* Urban dummy redefinition (clean 0/1)
gen byte urban1 = (urban==1)
drop urban

* Price scaling and logs
gen lncvaluep = ln(cvaluep)
gen cvaluepp  = cvaluep/10000   // in 10k RMB per m^2

*******************************************************
* 5) Baseline regressions (HDFE)
*******************************************************
*******************************************************
* Controls (dimension in brackets) — brief version
* Outcome: nbb | Key regressor: cvaluepp (10k RMB/m^2)
*******************************************************

*The below are three dimensions of the controlled variables
* Household×Year (i,t): income/assets/family + binaries
*local HH    Infinc fml Intotal_asset Innon_house_asset Inhouse_asset      // [HH×Year]
*local HHbin medsure_dum job urban1                                        // [HH×Year]

* Amenities (community or city level, time-varying)
*local AMEN  school hospital yundong ranqi                                 // [Comm/City×Year]

* City×Year controls (for extended/IV specs)
*local CITY  nop nopp GDPPC hukou_pop_city_10k res_land_area_sqkm_urban    // [City×Year]
* Compact summary (so referees don't squint)
*di as txt "HH×Year: `HH' `HHbin'"
*di as txt "Amenities (Comm/City×Year): `AMEN'"
*di as txt "City×Year: `CITY'"

* Baseline HDFE
reghdfe nbb cvaluepp `HH' `HHbin' `AMEN', absorb(cid city year fid) vce(cluster city)
* Baseline HDFE: reghdfe nbb cvaluepp `HH' `HHbin' `AMEN', absorb(cid city year fid) vce(cluster city);
* Interpret cvaluepp as the change in birth probability (percentage points) per +10k RMB/m^2;
* a negative sign is consistent with higher housing costs depressing fertility
* IV example (price endogenous; instrument = RVpc×wwptla)
ivreghdfe nbb `HHbin' Innon_house_asset `AMEN' ///
    (cvaluepp = c.RVpc##c.wwptla), absorb(year cid) cluster(city) first

* Income & assets
gen Intotal_asset      = ln(total_asset)
gen Infinc             = ln(finc)
gen Inhouse_asset      = ln(house_asset)
gen Innon_house_asset  = ln(non_house_asset)

* Baseline specifications (progressively richer FE and controls)
reghdfe nbb cvaluepp Infinc fml Intotal_asset medsure_dum job urban1, absorb(cid city#year fid) vce(cluster city)

reghdfe nbb cvaluepp Infinc fml Intotal_asset medsure_dum job urban1 ///
        school hospital yundong ranqi, absorb(cid city year fid) vce(cluster city)

reghdfe nbb cvaluepp Infinc fml Intotal_asset medsure_dum job urban1 ///
        school hospital yundong ranqi Innon_house_asset Inhouse_asset, absorb(cid city year) vce(cluster city)

reghdfe nbb cvaluepp Infinc fml Intotal_asset medsure_dum job urban1 ///
        school hospital yundong ranqi Innon_house_asset Inhouse_asset, absorb(cid city#year fid) vce(cluster city)

* Economic interpretation:
*  - cvaluepp measures the level of local housing prices. A negative coefficient suggests
*    higher housing costs reduce fertility (via higher cost of childrearing/housing).
*  - Fixed effects (cid / city#year / fid) absorb time-invariant community/household traits
*    and city-specific shocks, sharpening identification from within-unit variation.

*******************************************************
* 6) City-level controls and instruments
*******************************************************

* Scale city population and budget variables
capture confirm variable hukou_pop_city_10k
if !_rc {
    replace hukou_pop_city_10k = hukou_pop_city_10k/100   // rescale to "million persons"
}

capture confirm variable gen_budget_exp_city_10k
if !_rc & !_rc { // guard
    gen RVpercapita = gen_budget_exp_city_10k / hukou_pop_city_10k * 100
    replace RVpercapita = RVpercapita/100000
    gen RVpc = gen_budget_exp_city_10k / hukou_pop_city_10k
    replace RVpc = RVpc/10000
    label var RVpc "Per-capita general budget expenditure (scaled)"
}

capture confirm variable res_land_area_sqkm_urban
if !_rc & !_rc {
    gen tla     = res_land_area_sqkm_urban * 25        // follows original scaling
    label var tla "Scaled residential land area (×25)"
}

capture confirm variable ww
if !_rc & !_rc {
    gen wwptla  = ww / tla                              // water-to-land intensity
    label var wwptla "Water resources per unit residential land"
}

* Schools per capita
capture confirm variable num_middle_schools_city
if !_rc & !_rc {
    gen nop = num_middle_schools_city / hukou_pop_city_10k
    replace nop = nop/100
    label var nop "Middle schools per capita (scaled)"
}
capture confirm variable num_primary_schools_city
if !_rc & !_rc {
    gen nopp = num_primary_schools_city / (hukou_pop_city_10k*100)
    label var nopp "Primary schools per capita (scaled)"
}

* GDP per capita (scaled)
capture confirm variable gdp_city_curprice_bilrmb
if !_rc & !_rc {
    gen GDPPC = gdp_city_curprice_bilrmb*10000 / hukou_pop_city_10k
    replace GDPPC = GDPPC/100000
    replace GDPPC = GDPPC/100
    label var GDPPC "GDP per capita (scaled)"
}

*******************************************************
* 7) Instrumental-variables regressions (HDFE)
*******************************************************
* Endogenous regressor: cvaluepp (local house price level)
* Instruments: interaction c.RVpc##c.wwptla (fiscal capacity × water-to-land)
*
* IV Assumptions (commentary):
*  - Relevance: City fiscal capacity (RVpc) and water-to-land intensity (wwptla)
*    shift the supply/cost side of housing, thus moving local prices.
*    The interaction strengthens variation across cities/years.
*  - Exclusion: Conditional on fixed effects and controls, these city resource measures
*    affect fertility only through their impact on local housing prices—not directly
*    on household fertility decisions. 

ivreghdfe nbb fml medsure_dum health job Innon_house_asset school yundong ranqi ///
          (cvaluepp = c.RVpc##c.wwptla), absorb(year) cluster(city) first

ivreghdfe nbb fml medsure_dum health job Innon_house_asset school yundong ranqi ///
          (cvaluepp = c.RVpc##c.wwptla), absorb(year cid) cluster(city) first

ivreghdfe nbb fml medsure_dum health job Innon_house_asset school yundong ranqi ///
          nop nopp GDPPC hukou_pop_city_10k res_land_area_sqkm_urban ///
          (cvaluepp = c.RVpc##c.wwptla), absorb(year city cid) cluster(city) first

* Interpretation:
*  - The 2SLS coefficient on cvaluepp captures the causal effect of prices on fertility,
*    under the IV assumptions. "first" shows first-stage stats (e.g., F-statistics).

*******************************************************
* 8) Mechanism checks
*******************************************************

* 8.1 Wealth effect (heterogeneous slope by owner)
gen byte owner = (house_amount > 0)
label var owner "Homeowner (=1)"

* OLS with interaction
reghdfe nbb c.cvaluepp##c.owner Infinc fml Intotal_asset medsure_dum job urban1 ///
        school hospital yundong ranqi, absorb(cid city year fid) vce(cluster city)

* Alternative FE structure
reghdfe nbb c.cvaluepp##c.m Infinc fml Intotal_asset medsure_dum job urban1, ///
        absorb(cid city#year fid) vce(cluster city)

* IV with interaction (partialling out owner)
ivreghdfe nbb c.cvaluepp#c.owner owner fml medsure_dum health job Innon_house_asset ///
          school yundong ranqi (cvaluepp = c.RVpc##c.wwptla), absorb(year) cluster(city) first

ivreghdfe nbb c.cvaluepp#c.m m fml medsure_dum health job Innon_house_asset ///
          school yundong ranqi nop nopp GDPPC hukou_pop_city_10k res_land_area_sqkm_urban ///
          (cvaluepp = c.RVpc##c.wwptla), absorb(year city cid) cluster(city)

* 8.2 Homeowners subsample
keep if owner == 1
reghdfe nbb cvaluepp Infinc fml Intotal_asset medsure_dum job urban1 ///
        school hospital yundong ranqi, absorb(cid city year fid) vce(cluster city)

ivreghdfe nbb cvaluepp fml medsure_dum health job Innon_house_asset ///
          school yundong ranqi (cvaluepp = c.RVpc##c.wwptla), absorb(year) cluster(city) first

* 8.3 Credit/financial condition at baseline (2010 median, if year is indexed '1'==2010)
preserve
keep if year == 1    // If your 'year' is calendar, replace with "year==2010".
gen fp = rmb_loans_city_10k / hukou_pop_city_10k
tabstat fp, stats(median)
gen byte f = (fp > 1279767)
restore

* 8.4 Cost effect: smaller per-capita dwelling space proxy
gen mjpc   = zfmj / fml
bysort city year: egen mmjpc = mean(mjpc)
bysort city year: egen memjpc = median(mjpc)
gen byte m = (mjpc < memjpc)

* Restrict to homeowners for cost channel
keep if owner == 1
reghdfe nbb c.cvaluepp##c.m Infinc fml Intotal_asset medsure_dum job urban1 ///
        school hospital yundong ranqi, absorb(cid city year fid) vce(cluster city)

reghdfe nbb c.cvaluepp##c.m Infinc fml Intotal_asset medsure_dum job urban1, ///
        absorb(cid city#year fid) vce(cluster city)

ivreghdfe nbb c.cvaluepp#c.m m fml medsure_dum health job Innon_house_asset ///
          school yundong ranqi (cvaluepp = c.RVpc##c.wwptla), absorb(year) cluster(city) first

ivreghdfe nbb c.cvaluepp#c.m m fml medsure_dum health job Innon_house_asset ///
          school yundong ranqi nop nopp GDPPC hukou_pop_city_10k res_land_area_sqkm_urban ///
          (cvaluepp = c.RVpc##c.wwptla), absorb(year city cid) cluster(city)

*******************************************************
* 9) Heterogeneity analyses
*******************************************************

* Education
gen byte edudd = (xueli > 3)
reghdfe  nbb c.cvaluepp##c.edudd Infinc fml Intotal_asset medsure_dum job urban1 ///
         school hospital yundong ranqi, absorb(cid city year fid) vce(cluster city)
ivreghdfe nbb c.cvaluepp#c.edudd edudd fml medsure_dum health job Innon_house_asset ///
         school yundong ranqi (cvaluepp = c.RVpc##c.wwptla), absorb(year cid) cluster(city) first

* Number of children
gen byte child = (zinv > 1)
reghdfe  nbb c.cvaluepp##c.child Infinc fml Intotal_asset medsure_dum job urban1 ///
         school hospital yundong ranqi, absorb(cid city year fid) vce(cluster city)
ivreghdfe nbb c.cvaluepp#c.child child fml medsure_dum health job Innon_house_asset ///
         school yundong ranqi (cvaluepp = c.RVpc##c.wwptla), absorb(year cid) cluster(city) first

* Hukou residence (assuming 'res' indicates registration type)
reghdfe  nbb c.cvaluepp##c.res Infinc fml Intotal_asset medsure_dum job urban1 ///
         school hospital yundong ranqi, absorb(cid city year fid) vce(cluster city)
ivreghdfe nbb c.cvaluepp#c.res res fml medsure_dum health job Innon_house_asset ///
         school yundong ranqi (cvaluepp = c.RVpc##c.wwptla), absorb(year cid) cluster(city) first

* Child gender (assuming tb2_a_c1 identifies boy)
reghdfe  nbb c.cvaluepp##c.tb2_a_c1 Infinc fml Intotal_asset medsure_dum job urban1 ///
         school hospital yundong ranqi, absorb(cid city year fid) vce(cluster city)
ivreghdfe nbb c.cvaluepp#c.tb2_a_c1 tb2_a_c1 fml medsure_dum health job Innon_house_asset ///
         school yundong ranqi (cvaluepp = c.RVpc##c.wwptla), absorb(year cid) cluster(city) first

*******************************************************
* 10) Robustness checks
*******************************************************

* 10.1 Sample replacement: drop very low-price outliers and winsorize price
drop if resivaluep < 10000
winsor2 cvaluepp, replace cuts(0 98) by(year)

* 10.2 First-difference control
gen lcvaluepp  = L.cvaluepp
gen l2cvaluepp = L.lcvaluepp
gen dcvaluepp  = lcvaluepp - l2cvaluepp
reghdfe nbb cvaluepp dcvaluepp Infinc fml Intotal_asset medsure_dum job urban1 ///
        school hospital yundong ranqi Innon_house_asset Inhouse_asset, ///
        absorb(cid city year) vce(cluster city)

ivreghdfe nbb dcvaluepp fml medsure_dum health job Innon_house_asset ///
          school yundong ranqi (cvaluepp = c.RVpc##c.wwptla), absorb(year) cluster(city) first

*******************************************************
* End of analysis. 
*******************************************************
Define the mapping (Chinese -> English -> Meaning) ---
local N 9

local zh1 "水资源总量_万立方米_全市"
local en1 "total_water_city_10k_m3"
local ds1 "City total water resources (10k m^3); used to build normalized 'ww' later"

local zh2 "年末户籍人口_万人_全市"
local en2 "hukou_pop_city_10k"
local ds2 "City hukou population at year-end (10k persons); later rescaled in some specs"

local zh3 "地方一般公共预算支出_万元_全市"
local en3 "gen_budget_exp_city_10k"
local ds3 "City general budget expenditure (10k RMB)"

local zh4 "居住用地面积_平方公里_市辖区"
local en4 "res_land_area_sqkm_urban"
local ds4 "Urban residential land area (sq. km)"

local zh5 "普通中学_所_全市"
local en5 "num_middle_schools_city"
local ds5 "Number of middle schools in the city"

local zh6 "普通小学_所_全市"
local en6 "num_primary_schools_city"
local ds6 "Number of primary schools in the city"

local zh7 "地区生产总值_当年价格_亿元_全市"
local en7 "gdp_city_curprice_bilrmb"
local ds7 "City GDP at current prices (billion RMB)"

local zh8 "年末金融机构人民币各项贷款余额_万元_全市"
local en8 "rmb_loans_city_10k"
local ds8 "Year-end RMB loans at financial institutions (10k RMB)"

* For completeness: "implicit string values" recoded earlier to dummies
* (we include them in the dic for readers)
local zh9 `"（implicit string values）"男" / "健康" / "不健康" / "农业户口" / "非农业户口""'
local en9 "gender / health / hukou"
local ds9 "Binary dummies: Male=1, Healthy=1, Agricultural hukou=1 (created from Chinese strings)"