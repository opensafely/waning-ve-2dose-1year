from cohortextractor import (
    StudyDefinition, 
    patients
)

# Import codelists.py script
from codelists import *

# import json module
import json
#study_parameters
with open("./analysis/lib/study_parameters.json") as f:
  study_parameters = json.load(f)

# define variables explicitly
K=study_parameters["K"]
start_date=study_parameters["start_date"] # start of phase 1
end_date=study_parameters["end_date"] # latest date of data

# set seed so that dummy data can be reproduced
import numpy as np
np.random.seed(study_parameters["seed"])

# read start and end dates from csv
def type_X_date(type, n):
  def var_signature(name, k):
    return {
      name: patients.with_value_from_file(
        f_path='output/data/data_eligible_e.csv', 
        returning=f"{type}_{k}_date", 
        returning_type='date',
        date_format='YYYY-MM-DD'
        ),
    }
  variables=dict()
  for i in range(1, n+1):
    variables.update(var_signature(name=f"{type}_{i}_date", k=i))
  return variables

# any covid test in each comparison period
def anytest_X_date(n):
  def var_signature(name, k):
    return {
      name: patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        test_result="any",
        between=[f"start_{k}_date + 1 days", f"end_{k}_date"],
        restrict_to_earliest_specimen_date=False,
        find_first_match_in_period=True,
        returning="date",
        date_format = "YYYY-MM-DD",
	    ),
    }
  variables=dict()
  for i in range(1, n+1):
    variables.update(var_signature(name=f"anytest_{i}_date", k=i))
  return variables

###
study=StudyDefinition(

    default_expectations={
        "date": {"earliest": start_date, "latest": end_date},
        "rate": "uniform",
        "incidence": 0.8,
    },  

    population=patients.all(),

    # min elig date within subgroup
    svp_start_date=patients.with_value_from_file(
        f_path='output/data/data_eligible_e.csv', 
        returning='svp_start_date', 
        returning_type='date',
        date_format='YYYY-MM-DD'
        ),
    # min elig date within subgroup
    min_elig_date=patients.with_value_from_file(
        f_path='output/data/data_eligible_e.csv', 
        returning='min_elig_date', 
        returning_type='date',
        date_format='YYYY-MM-DD'
        ),
    # comparison start and end dates
    **type_X_date("start", K),
    **type_X_date("end", K),

    # age on svp start date
    age=patients.age_as_of(
        "svp_start_date",
        return_expectations={
            "rate" : "universal",
            "int" : {"distribution" : "population_ages"}
            }
        ),

    ### covid tests as covariates
    # during unvaccinated time (from when tests widely availabe to elig_date)
    test_hist_n=patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        test_result="any",
        between=["2020-05-18", "min_elig_date - 1 day"], # day before 1st vaccine eligibility date
        restrict_to_earliest_specimen_date=False,
        returning="number_of_matches_in_period",
        return_expectations={"int" : {"distribution": "poisson", "mean": 2}, "incidence" : 0.6}
	    ),

    # most recent numeric BMI
    bmi=patients.categorised_as(
    {
      "Not obese": "DEFAULT",
      "Obese I (30-34.9)": """ bmi_value >= 30 AND bmi_value < 35""",
      "Obese II (35-39.9)": """ bmi_value >= 35 AND bmi_value < 40""",
      "Obese III (40+)": """ bmi_value >= 40 AND bmi_value < 100""",
      # set maximum to avoid any impossibly extreme values being classified as obese
    },
    bmi_value=patients.most_recent_bmi(
        between=["svp_start_date - 5 years", "svp_start_date - 1 day"],
        minimum_age_at_measurement=16
    ),
    return_expectations={
      "rate": "universal",
      "category": {
        "ratios": {
          "Not obese": 0.7,
          "Obese I (30-34.9)": 0.1,
          "Obese II (35-39.9)": 0.1,
          "Obese III (40+)": 0.1,
        }
      },
    },
  ),

    # severe asthma
    asthma=patients.satisfying(
        """
        astadm OR
        (ast AND astrxm1 AND astrxm2 AND astrxm3)
        """,
    # Asthma Admission codes
    astadm=patients.with_these_clinical_events(
      astadm_primis,
      returning="binary_flag",
      on_or_before="svp_start_date - 1 day",
    ),
    # Asthma Diagnosis code
    ast = patients.with_these_clinical_events(
      ast_primis,
      returning="binary_flag",
      on_or_before="svp_start_date - 1 day",
    ),
    # Asthma systemic steroid prescription code in month 1
    astrxm1=patients.with_these_medications(
      astrx_primis,
      returning="binary_flag",
      between=["svp_start_date - 30 days", "svp_start_date - 1 day"],
    ),
    # Asthma systemic steroid prescription code in month 2
    astrxm2=patients.with_these_medications(
      astrx_primis,
      returning="binary_flag",
      between=["svp_start_date - 60 days", "svp_start_date - 31 days"],
    ),
    # Asthma systemic steroid prescription code in month 3
    astrxm3=patients.with_these_medications(
      astrx_primis,
      returning="binary_flag",
      between= ["svp_start_date - 90 days", "svp_start_date - 61 days"],
        ),
    ),

    # chronic respiratory disease
    crd=patients.satisfying(
        "asthma OR resp",
        # Chronic Respiratory Disease other than asthma
        resp=patients.with_these_clinical_events(
            resp_primis,
            returning="binary_flag",
            on_or_before="svp_start_date - 1 day",
            return_expectations={"incidence": 0.02},
            ),
    ),

    # Chronic Neurological Disease including Significant Learning Disorder
    cns=patients.with_these_clinical_events(
        cns_primis,
        returning="binary_flag",
        on_or_before="svp_start_date - 1 day",
        return_expectations={"incidence": 0.02},
    ),

    # Chronic kidney disease diagnostic codes
    ckd=patients.satisfying(
        """
        ckd_any OR
        (ckd15_date AND 
        (ckd35_date >= ckd15_date) OR (ckd35_date AND NOT ckd15_date))
        """,
        # Chronic kidney disease codes - all stages
        ckd15_date=patients.with_these_clinical_events(
            ckd15_primis,
            returning="date",
            find_last_match_in_period=True,
            on_or_before="svp_start_date - 1 day",
            date_format="YYYY-MM-DD",
        ),
        # Chronic kidney disease codes-stages 3 - 5
        ckd35_date=patients.with_these_clinical_events(
            ckd35_primis,
            returning="date",
            find_last_match_in_period=True,
            on_or_before="svp_start_date - 1 day",
            date_format="YYYY-MM-DD",
        ),
        # Chronic kidney disease diagnostic codes
        ckd_any=patients.with_these_clinical_events(
            ckd_primis,
            returning="binary_flag",
            on_or_before="svp_start_date - 1 day",
        ),
        return_expectations={"incidence": 0.01},
    ),

    # Diabetes
    diabetes=patients.with_these_clinical_events(
        diab_primis,
        returning="binary_flag",
        on_or_before="svp_start_date - 1 day",
        return_expectations={"incidence": 0.02},
        ),

    # Severe mental illness
    sev_mental=patients.with_these_clinical_events(
        sev_mental_primis,
        returning="binary_flag",
        on_or_before="svp_start_date - 1 day",
        return_expectations={"incidence": 0.02},
        ),

    # Chronic heart disease codes
    chd=patients.with_these_clinical_events(
        chd_primis,
        returning="binary_flag",
        on_or_before="svp_start_date - 1 day",
        return_expectations={"incidence": 0.02},
    ),

    # Chronic Liver disease codes
    cld=patients.with_these_clinical_events(
        cld_primis,
        returning="binary_flag",
        on_or_before="svp_start_date - 1 day",
        return_expectations={"incidence": 0.02},
    ),

    # Immunosuppression
    immunosuppressed=patients.satisfying(
    "immrx OR immdx",
    # Immunosuppression diagnosis codes
    immdx=patients.with_these_clinical_events(
        immdx_primis,
        returning="binary_flag",
        on_or_before="svp_start_date - 1 day",
        return_expectations={"incidence": 0.02},
        ),
    # Immunosuppression medication codes
    immrx=patients.with_these_medications(
        immrx_primis,
        returning="binary_flag",
        between=["svp_start_date - 180 days", "svp_start_date - 1 day"],
        ),
    ),

    # Asplenia or Dysfunction of the Spleen codes
    asplenia=patients.with_these_clinical_events(
        spln_primis,
        returning="binary_flag",
        on_or_before="svp_start_date - 1 day",
        return_expectations={"incidence": 0.02},
    ),

    # Learning Disability
    learndis=patients.with_these_clinical_events(
        learndis_primis,
        returning="binary_flag",
        on_or_before="svp_start_date - 1 day",
        return_expectations={"incidence": 0.02},
    ),

    # flu vaccine in the 5 years before svp start date
    flu_vaccine=patients.satisfying(
        """
        flu_vaccine_tpp_table>0 OR
        flu_vaccine_med>0 OR
        flu_vaccine_clinical>0
        """,
        
        flu_vaccine_tpp_table=patients.with_tpp_vaccination_record(
            target_disease_matches="INFLUENZA",
            between=["svp_start_date - 5 years", "svp_start_date"], 
            returning="binary_flag",
        ),
        
        flu_vaccine_med=patients.with_these_medications(
            flu_med_codes,
            between=["svp_start_date - 5 years", "svp_start_date"], 
            returning="binary_flag",
        ),
        flu_vaccine_clinical=patients.with_these_clinical_events(
            flu_clinical_given_codes,
            ignore_days_where_these_codes_occur=flu_clinical_not_given_codes,
            between=["svp_start_date - 5 years", "svp_start_date"], 
            returning="binary_flag",
        ),
        return_expectations={"incidence": 0.5, },
    ),

    # pregnancy
    pregnancy=patients.satisfying(
        """
        (preg_36wks_date) AND
        (pregdel_pre_date <= preg_36wks_date OR NOT pregdel_pre_date)
        """,
        # date of last pregnancy code in 36 weeks before ref_cev
        preg_36wks_date=patients.with_these_clinical_events(
            preg_primis,
            returning="date",
            find_last_match_in_period=True,
            between=["svp_start_date - 252 days", "svp_start_date - 1 day"],
            date_format="YYYY-MM-DD",
        ),
        # date of last delivery code recorded in 36 weeks before elig_date
        pregdel_pre_date=patients.with_these_clinical_events(
            pregdel_primis,
            returning="date",
            find_last_match_in_period=True,
            between=["svp_start_date - 252 days", "svp_start_date - 1 day"],
            date_format="YYYY-MM-DD",
        ),
    ),

    # first occurence of any covid test in each comparison period
    **anytest_X_date(K),

)
