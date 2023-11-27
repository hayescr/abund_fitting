# abund_fitting
Codes that I use for doing a manual abundance analysis

## Set up

1. move into file TOP FILE / synthesis
2. python initialize_params.py
3. test things by running abund_fit_fe_test.py

## Measuring parameters

1. fit vmacro with fit_convol.py
2. initial round of fe with abund_fit_fe.py (can update fe using the convol fits if needed with update_feh_param.py)
3. investigate how the parameters work by moving the fe measurements up from abund_novel and running evaluate_param.py, potentially update fe
4. fit Mg and C with new fe using abund_fit_mg and abund_fit_c, then update each of these with update_param_abunds
5. iterate a couple times on abundances of Fe, Mg, and C then see if the teff and logg are reasonable with evaluate param.  If not update Teff first and see if that fixes things (continuing to iterate on fe mg c until things converge at a given teff)

## Measuring abundances

1. once the params converge run abund_fit_elem_no_velshift on elements reasonable for the S/N
2. look through elements and update the elements in the main param file with update_param_abunds
3. rerun all of the elements and do a final investigation

## Start making tables

NEED TO ADD STEPS
