from cohortextractor import (codelist, codelist_from_csv, combine_codelists)

#######################
### for JCVI groups ###
#######################

# Asthma Diagnosis code
ast_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-ast.csv",
    system="snomed",
    column="code",
)

# Asthma Admission codes
astadm_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-astadm.csv",
    system="snomed",
    column="code",
)

# Asthma systemic steroid prescription codes
astrx_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-astrx.csv",
    system="snomed",
    column="code",
)

# Chronic Respiratory Disease
resp_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-resp_cov.csv",
    system="snomed",
    column="code",
)

# Chronic heart disease codes
chd_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-chd_cov.csv",
    system="snomed",
    column="code",
)

# Chronic kidney disease diagnostic codes
ckd_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-ckd_cov.csv",
    system="snomed",
    column="code",
)

# Chronic kidney disease codes - all stages
ckd15_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-ckd15.csv",
    system="snomed",
    column="code",
)

# Chronic kidney disease codes-stages 3 - 5
ckd35_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-ckd35.csv",
    system="snomed",
    column="code",
)

# Chronic Liver disease codes
cld_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-cld.csv",
    system="snomed",
    column="code",
)

# Diabetes diagnosis codes
diab_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-diab.csv",
    system="snomed",
    column="code",
)

# Immunosuppression diagnosis codes
immdx_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-immdx_cov.csv",
    system="snomed",
    column="code",
)

# Immunosuppression medication codes
immrx_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-immrx.csv",
    system="snomed",
    column="code",
)

# Chronic Neurological Disease including Significant Learning Disorder
cns_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-cns_cov.csv",
    system="snomed",
    column="code",
)

# Asplenia or Dysfunction of the Spleen codes
spln_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-spln_cov.csv",
    system="snomed",
    column="code",
)

# BMI
bmi_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-bmi.csv",
    system="snomed",
    column="code",
)

# All BMI coded terms
bmi_stage_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-bmi_stage.csv",
    system="snomed",
    column="code",
    category_column="term",
)

# Severe Obesity code recorded
sev_obesity_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-sev_obesity.csv",
    system="snomed",
    column="code",
)

# Diabetes resolved codes
dmres_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-dmres.csv",
    system="snomed",
    column="code",
)

# Severe Mental Illness codes
sev_mental_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-sev_mental.csv",
    system="snomed",
    column="code",
)

# Remission codes relating to Severe Mental Illness
smhres_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-smhres.csv",
    system="snomed",
    column="code",
)

# High Risk from COVID-19 code
shield_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-shield.csv",
    system="snomed",
    column="code",
)

# Lower Risk from COVID-19 codes
nonshield_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-nonshield.csv",
    system="snomed",
    column="code",
)

# to represent household contact of shielding individual
hhld_imdef_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-hhld_imdef.csv",
    system="snomed",
    column="code",
)

# Wider Learning Disability
learndis_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-learndis.csv",
    system="snomed",
    column="code",
)

# Carer codes
carer_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-carer.csv",
    system="snomed",
    column="code",
)

# No longer a carer codes
notcarer_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-notcarer.csv",
    system="snomed",
    column="code",
)

# Employed by Care Home codes
carehome_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-carehome.csv",
    system="snomed",
    column="code",
)

# Employed by nursing home codes
nursehome_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-nursehome.csv",
    system="snomed",
    column="code",
)

# Employed by domiciliary care provider codes
domcare_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-domcare.csv",
    system="snomed",
    column="code",
)

# Patients in long-stay nursing and residential care
longres_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-longres.csv",
    system="snomed",
    column="code",
)

# Patients who are housebound
housebound = codelist_from_csv(
    "codelists/opensafely-housebound.csv", 
    system="snomed", 
    column="code"
)
# No longer housebound
no_longer_housebound = codelist_from_csv(
    "codelists/opensafely-no-longer-housebound.csv", 
    system="snomed", 
    column="code"
)

# Pregnancy codes 
preg_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-preg.csv",
    system="snomed",
    column="code",
)

# Pregnancy or Delivery codes
pregdel_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-pregdel.csv",
    system="snomed",
    column="code",
)

#################################
### for demographic variables ###
#################################
# Ethnicity codes
eth2001_primis = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-eth2001.csv",
    system="snomed",
    column="code",
    category_column="grouping_6_id",
)

#############################
### COVID infection codes ###
#############################

## History of covid
covid_codes = codelist_from_csv(
  "codelists/opensafely-covid-identification.csv",
  system = "icd10",
  column = "icd10_code",
)

# probable COVID in primary care
covid_primary_care_positive_test = codelist_from_csv(
    "codelists/opensafely-covid-identification-in-primary-care-probable-covid-positive-test.csv",
    system="ctv3",
    column="CTV3ID",
)
covid_primary_care_code = codelist_from_csv(
    "codelists/opensafely-covid-identification-in-primary-care-probable-covid-clinical-code.csv",
    system="ctv3",
    column="CTV3ID",
)
covid_primary_care_sequalae = codelist_from_csv(
    "codelists/opensafely-covid-identification-in-primary-care-probable-covid-sequelae.csv",
    system="ctv3",
    column="CTV3ID",
)
covid_primary_care_probable_combined=combine_codelists(
    covid_primary_care_positive_test,
    covid_primary_care_code,
    covid_primary_care_sequalae,
)

# for identifying hospitalisations from emergency data
covid_emergency = codelist(
    ["1240751000000100"],
    system="snomed",
)

discharged_to_hospital = codelist(
    ["306706006", "1066331000000109", "1066391000000105"],
    system="snomed",
)


flu_med_codes = codelist_from_csv(
    "codelists/opensafely-influenza-vaccination.csv",  
    system="snomed",  
    column="snomed_id",
)
flu_clinical_given_codes = codelist_from_csv(
    "codelists/opensafely-influenza-vaccination-clinical-codes-given.csv",  
    system="ctv3", 
    column="CTV3ID",
)
flu_clinical_not_given_codes = codelist_from_csv(
    "codelists/opensafely-influenza-vaccination-clinical-codes-not-given.csv",  
    system="ctv3", 
    column="CTV3ID",
)
eol_codes = codelist_from_csv(
    "codelists/nhsd-primary-care-domain-refsets-palcare_cod.csv",
    system="snomed",
    column="code",
)
midazolam_codes = codelist_from_csv(
    "codelists/opensafely-midazolam-end-of-life.csv",
    system="snomed",
    column="dmd_id",   
)
