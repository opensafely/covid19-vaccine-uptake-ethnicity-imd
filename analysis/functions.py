from cohortextractor import patients
from codelists import *

####################################################################################################
# function to add days to a string date
from datetime import datetime, timedelta
def days(datestring, days):
  
  dt = datetime.strptime(datestring, "%Y-%m-%d").date()
  dt_add = dt + timedelta(days)
  datestring_add = datetime.strftime(dt_add, "%Y-%m-%d")

  return datestring_add

####################################################################################################
# functio to extract a sequence of vacination dates
def vaccination_date_X(name, index_date, n, product_name_matches=None, target_disease_matches=None):
  # vaccination date, given product_name
  def var_signature(
    name,
    on_or_after,
    product_name_matches,
    target_disease_matches
  ):
    return {
      name: patients.with_tpp_vaccination_record(
        product_name_matches=product_name_matches,
        target_disease_matches=target_disease_matches,
        on_or_after=on_or_after,
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
      ),
    }
  variables = var_signature(f"{name}_1_date", index_date, product_name_matches, target_disease_matches)
  for i in range(2, n+1):
    variables.update(var_signature(
      f"{name}_{i}_date", 
      f"{name}_{i-1}_date + 1 days",
      # pick up subsequent vaccines occurring one day or later -- people with unrealistic dosing intervals are later excluded
      product_name_matches,
      target_disease_matches
    ))
  return variables
  
####################################################################################################
# this will need further quality checks to make sure sex='F' and age < 50
def pregnancy(name, index_date, type):
  if type == "value": between = [days(index_date, -252), days(index_date, -1)]
  if type == "column": between = [f"{index_date} - 252 days", f"{index_date} - 1 day"]

  def tmp_preg(name, between):
    return {
      # date of last pregnancy code in 36 weeks before index_date
      f"{name}_36wks_date": patients.with_these_clinical_events(
        preg_primis,
        returning = "date",
        find_last_match_in_period = True,
        between = between,
        date_format = "YYYY-MM-DD",
      )
    }
  
  def tmp_del(name, between):
    return {
      # date of last delivery code recorded in 36 weeks before index_date
      f"{name}_del_date": patients.with_these_clinical_events(
        pregdel_primis,
        returning = "date",
        find_last_match_in_period = True,
        between = between,
        date_format = "YYYY-MM-DD",
      )
    }

  return {
    name: patients.satisfying(
        f"""
        ({name}_36wks_date AND NOT {name}_del_date) OR
        ({name}_36wks_date AND {name}_del_date AND ({name}_del_date < {name}_36wks_date))
        """,
      **tmp_preg(name, between),
      **tmp_del(name, between),
      return_expectations = {"incidence": 0.05},
    ),
  }
