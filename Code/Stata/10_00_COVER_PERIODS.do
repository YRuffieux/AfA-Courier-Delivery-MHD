use "$source/tblCOVERPERIODS", clear

tab scheme_code, mi
codebook coverfrom_date    
codebook coverto_date

* periods with no scheme code have no start date and vice versa; they all are isolated cases, i.e. thsese patients have no other cover period
merge m:1 patient using "$source/tblBAS", keep(match) keepusing(enrol_d birth_d) nogen
merge m:1 patient using "$source/tblLTFU", keep(match) keepusing(drop_d death_d) nogen

sort patient coverfrom_date coverto_date
bysort patient: egen coverto_date_max=max(coverto_date)
format coverto_date_max %tdD_m_CY
gen start=max(enrol_d,$open_d)
gen end=min(drop_d,death_d,coverto_date_max,$close_d)
format start end %tdD_m_CY
gunique patient         // 474,076
gunique patient if enrol_d>drop_d   // 6,165
drop if start>=end 

* identifying last scheme code before the patient leaves the study
keep if coverfrom_date<end
gunique patient         // 455,979
bysort patient (coverto_date): keep if _n==_N
rename scheme_code last_scheme_code
tab last_scheme_code
codebook coverto_date
codebook coverfrom_date

* switching off deaths that occur after patient exits from the scheme
gen death_y=0
replace death_y=1 if death_d==end

keep patient start end last_scheme_code death_y
tab last_scheme_code, mi

save "$clean/FUPwide", replace