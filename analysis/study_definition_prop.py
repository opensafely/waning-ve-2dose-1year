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
start_date=study_parameters["start_date"] # start of phase 1
end_date=study_parameters["end_date"] # latest date of data

# set seed so that dummy data can be reproduced
import numpy as np
np.random.seed(study_parameters["seed"])

####################################################################################################
## function to add days to a string date
from datetime import datetime, timedelta
def days(datestring, days):
  
  try: 
     dt = datetime.strptime(datestring, "%Y-%m-%d").date()
     dt_add = dt + timedelta(days)
     datestring_add = datetime.strftime(dt_add, "%Y-%m-%d")
  except ValueError:
     if days > 0: datestring_add = datestring + f" + {days} days"
     else: datestring_add = datestring + f" - {abs(days)} days"

  return datestring_add

####################################################################################################
## function for extracting recurrent positive tests
def postest_X_date(n, index_date, shift):
  def var_signature(index_k_date, name):
    return {
      name: patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        test_result="positive",
        on_or_after=index_k_date,
        restrict_to_earliest_specimen_date=False,
        find_first_match_in_period=True,
        returning="date",
        date_format = "YYYY-MM-DD",
	    ),
    }
  variables= var_signature(
      index_k_date = days(index_date, shift),
      name = "postest_1_date"
      )
  for i in range(2, n+1):
    variables.update(var_signature(
      index_k_date = f"postest_{i-1}_date + 1 day",
      name = f"postest_{i}_date"
      ))
  return variables

####################################################################################################
## function for extracting recurrent hospitalisations
def hospitalisation_X_date(n, index_date, admission_type = None, covid = False):
  if admission_type == "unplanned": 
      with_admission_method = ["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"]
  if admission_type == "planned":
    with_admission_method = ["11", "12", "13", "81"]
  if covid == False:
    with_these_diagnoses = None
  if covid == True:
    with_these_diagnoses = covid_codes
    admission_type = f"covid{admission_type}"
  ##
  def var_signature(
      k, return_type, 
      on_or_before = None, on_or_after = None, 
      find_last_match_in_period = False, find_first_match_in_period = False
      ):
    # latest admission before index date
    if k==0:
      on_or_before = days(index_date, -1)
      find_last_match_in_period = True
    # first admission after index date
    if k==1: 
      on_or_after = index_date
      find_first_match_in_period = True
    # first admission after the previous admission
    if k>1: 
      on_or_after = f"admitted_{admission_type}_{i-1}_date + 1 day"
      find_first_match_in_period = True
    return {
      f"{return_type}_{admission_type}_{k}_date": patients.admitted_to_hospital(
        returning = f"date_{return_type}",
        with_admission_method = with_admission_method,
        with_these_diagnoses = with_these_diagnoses,
        # the dates specify the admission date (even if returning "date_discharged")
        on_or_before = on_or_before,
        on_or_after = on_or_after,
        date_format = "YYYY-MM-DD",
        find_first_match_in_period = find_first_match_in_period,
        find_last_match_in_period = find_last_match_in_period,
        ),
    }
  ##
  variables = dict()
  for i in range(0, n+1):
    variables.update(
      var_signature(i, return_type = "admitted"),
    )
    variables.update(
      var_signature(i, return_type = "discharged"),
    )
  return variables

####################################################################################################
study=StudyDefinition(

    default_expectations={
        "date": {"earliest": days(start_date,-28), "latest": end_date},
        "rate": "uniform",
        "incidence": 0.8,
    },
  
    population=patients.satisfying(
        "eligible_e", 
        eligible_e = patients.which_exist_in_file(
            f_path=f"output/data/data_eligible_e.csv"
        ),
    ),

    start_1_date = patients.with_value_from_file(
        f_path='output/data/data_eligible_e.csv', 
        returning="start_1_date", 
        returning_type='date',
        date_format='YYYY-MM-DD'
        ),

    # positive tests (recurring)
    # shift = -28 because anyone with evidence of covid >28 days before start_1_date already excluded
    **postest_X_date(12, "start_1_date", shift = -28),

    # unplanned hospitalisations (recurring)
    **hospitalisation_X_date(6, "start_1_date", admission_type = "unplanned"),

    # unplanned covid hospitalisations (recurring)
    **hospitalisation_X_date(6, "start_1_date", admission_type = "unplanned", covid = True),

    # planned hospitalisations (recurring)
    **hospitalisation_X_date(6, "start_1_date", admission_type = "planned"),

    # end of life care initiated (once)
    endoflife_date=patients.with_these_clinical_events(
        eol_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        # -28 days because anyone with eol >28 days before start_1_date already excluded
        on_or_after="start_1_date - 28 days",
        find_first_match_in_period=True,
    ),
  
)